/*
 * @Author: Hsuan-Wei Fan
 * @Date: 2020-04-23 13:47:13
 * @LastEditors: Hsuan-Wei Fan
 * @LastEditTime: 2020-04-23 13:50:48
 * @Description: The definition of vae subsurface scattering material.
 */

#include "vaescatter.h"
#include "reflection.h"
#include "vaehandlereigen.h"
#include "sampler.h"
#include "stats.h"
#include "scene.h"

namespace pbrt {


VaeScatter::~VaeScatter() {
    // FIXME: Not sure if delete here is ok?
    // if (po.polyStorage) {
    //     delete po.polyStorage;
    // }
}


Spectrum VaeScatter::S(const SurfaceInteraction &pi, const Vector3f &wi) {
    // FIXME: Direct copy yet.
    ProfilePhase pp(Prof::BSSRDFEvaluation);
    Float Ft = FrDielectric(CosTheta(po.wo), 1, eta);
    return (1 - Ft) * Sw(wi); // No sp here.
}


Spectrum VaeScatter::Sample_S(const Scene &scene, Float u1, const Point2f &u2,
                              MemoryArena &arena, SurfaceInteraction *si,
                              Float *pdf) const {
    // TODO: 采样 S 函数的值，将其储存于 pdf 和 si 中, 並初始化完成 VAE 模式的 BSDF
    ProfilePhase pp(Prof::BSSRDFSampling);
    Vector3f refractedD = -this->po.wo; // FIXME: -
    // Ray zeroScatterRay(this->po.p, refractedD); // - refracted?
    // SurfaceInteraction zeroScatterIts;
    // if (!scene.Intersect(zeroScatterRay, &zeroScatterIts)) {
    //     return Spectrum(0.0f);
    // }
    // TODO: consider ray passes through situation.

    Spectrum Sp = Sample_Sp(scene, refractedD, si, pdf, 3);
    if (!Sp.IsBlack()) { // FIXME: valid intersection.
        si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);
        si->bsdf->Add(ARENA_ALLOC(arena, VaeScatterAdapter)(this));
        si->wo = Vector3f(si->shading.n);
    }

    return Sp;
}

void VaeScatter::Sample_Pi(ScatterSamplingRecord *sRecs, const Scene &scene, SurfaceInteraction *res,
                           const Vector3f &w, int nSamples) const {
    const Point3f &po = this->po.p;
    const Vector3f &wo = this->po.wo;
    // TODO: Compute albedo
    const Spectrum &albedo = mVaeHandler->GetMedium().albedo;
    const Spectrum &sigmaT = mVaeHandler->GetMedium().sigmaT;
    Vector3f inDir = -w;
    if (nSamples == 3) {
        for (int i = 0; i < 3; i++) {
            Normal3f polyNormal = PolyUtils::AdjustRayDirForPolynomialTracing(
                    inDir, this->po, 3, PolyUtils::GetFitScaleFactor(mVaeHandler->GetMedium(), i), i);
            sRecs[i] = mVaeHandler->Sample(po, wo, &scene, polyNormal, sigmaT, albedo,
                                           this->g, this->eta, this->po, true, i, &res[i]);
            Spectrum tmp = sRecs[i].throughout;
            sRecs[i].throughout = Spectrum(0.0f);
            sRecs[i].throughout[i] = tmp[i] * 3.0f;
        }
    } else {
        CHECK_EQ(nSamples, 1);
        // Choose a random color channel to sample
        int channel = int(3 * pbrt::GetRandomFloat());
        Normal3f polyNormal = PolyUtils::AdjustRayDirForPolynomialTracing(
                inDir, this->po, 3, PolyUtils::GetFitScaleFactor(mVaeHandler->GetMedium()), channel);
        sRecs[0] = mVaeHandler->Sample(po, wo, &scene, polyNormal, sigmaT, albedo,
                                       this->g, this->eta, this->po, true, channel, res);
    }
}

Spectrum VaeScatter::Sample_Sp(const Scene &scene, const Vector3f &refractedD,
                               SurfaceInteraction *resIsect,
                               Float *pdf, int nSamples) const {
    // FIXME: Check the calculation of Sp function
    ProfilePhase pp(Prof::BSSRDFEvaluation);
    ScatterSamplingRecord sRecs[nSamples];
    DCHECK_GE(nSamples, 0);
    DCHECK_LE(nSamples, 3);
    SurfaceInteraction sIsect[nSamples];
    Sample_Pi(sRecs, scene, sIsect, refractedD, nSamples);
    Spectrum res;
    int nMissed = 0;
    for (int i = 0; i < nSamples; i++) {
        const ScatterSamplingRecord &s = sRecs[i];
        if (!s.isValid) {
            nMissed++;
            continue;
        }
        res += s.throughout;
        *resIsect = sIsect[i];
    }
    *pdf = 1.0f; // FIXME: Maybe the pdf need to be calculated
    return res / (Float) nSamples;
}


} // namespace pbrt