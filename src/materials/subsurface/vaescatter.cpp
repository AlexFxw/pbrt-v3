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
#include "stats.h"

namespace pbrt {


VaeScatter::VaeScatter() {
    mVaeHandler = nullptr;
}

VaeScatter::~VaeScatter() {
    if (mVaeHandler)
        delete mVaeHandler;
}

void VaeScatter::Prepare(const Scene *scene, const std::vector<std::shared_ptr<Shape>> &shapes) {
    DCHECK_EQ(mTriangles.size(), 0);
    DCHECK_GT(shapes.size(), 0);
    for (int i = 0; i < shapes.size(); i++) {
        mTriangles.push_back(shapes[i]);
    }
    // FIXME: set appropriate values
    Float kernelEpsScale = 0.5f;
    mVaeHandler = new VaeHandlerEigen(kernelEpsScale);
    PolyUtils::PolyFitConfig pfConfig;
    std::string modelName = "";
    std::string absModelName = "";
    std::string angularModelName = "";
    std::string outputDir = "";
    int sssSamples = 300;
    mVaeHandler->Prepare(scene, shapes, mSigmaT, mAlbedo,
                         g, eta, modelName, absModelName,
                         angularModelName, outputDir, sssSamples, pfConfig);
}


Spectrum VaeScatter::S(const SurfaceInteraction &pi, const Vector3f &wi) {
    // TODO
    ProfilePhase pp(Prof::BSSRDFEvaluation);
    return pbrt::Spectrum();
}


Spectrum VaeScatter::Sample_S(const Scene &scene, Float u1, const Point2f &u2,
                              MemoryArena &arena, SurfaceInteraction *si,
                              Float *pdf) const {
    // TODO: 采样 S 函数的值，将其储存于 pdf 和 si 中, 並初始化完成 VAE 模式的 BSDF
    ProfilePhase pp(Prof::BSSRDFSampling);
    Spectrum Sp = Sample_Sp(scene, arena, si, pdf);
    if (!Sp.IsBlack()) {
        si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);
        si->wo = Vector3f(si->shading.n);
    }

    return Sp;
}

void VaeScatter::Sample_Pi(ScatterSamplingRecord *sRecs, const Scene &scene, const SurfaceInteraction &itact,
                           const Vector3f &w, Sampler &sampler, int nSamples) const {
    const Point3f &po = this->po.p;
    const Vector3f &wo = this->po.wo;
    SurfaceInteraction sampledIsect; // TODO: Make use of this, it is projected one.
    // TODO: Compute albedo
    Spectrum albedo = mAlbedo;
    Vector3f inDir = -w;
    if (nSamples == 3) {
        for (int i = 0; i < 3; i++) {
            Normal3f polyNormal = PolyUtils::AdjustRayDirForPolynomialTracing(inDir, itact, 3,
                                                                              PolyUtils::GetFitScaleFactor(mMedium, i),
                                                                              i);
            sRecs[i] = mVaeHandler->Sample(po, wo, &scene, polyNormal, mSigmaT, albedo,
                                           this->g, this->eta, sampler, itact, &sampledIsect, true, i);
            Spectrum tmp = sRecs[i].throughout;
            sRecs[i].throughout = Spectrum(0.0f);
            sRecs[i].throughout[i] = tmp[i] * 3.0f;
        }
    } else {
        CHECK_EQ(nSamples, 1);
        // Choose a random color channel to sample
        int channel = int(3 * sampler.Get1D());
        Normal3f polyNormal = PolyUtils::AdjustRayDirForPolynomialTracing(inDir, itact, 3,
                                                                          PolyUtils::GetFitScaleFactor(mMedium),
                                                                          channel);
        sRecs[0] = mVaeHandler->Sample(po, inDir, &scene, polyNormal, mSigmaT, albedo,
                                       this->g, this->eta, sampler, itact, &sampledIsect, true, channel);
    }
}

Spectrum VaeScatter::Sample_Sp(const Scene &scene, MemoryArena &arena,
                               SurfaceInteraction *pi, Float *pdf) const {
    // TODO: Initialize pdf;
    // TODO: Initialize pi, and pi->shading.
    ProfilePhase pp(Prof::BSSRDFEvaluation);
    // FIXME: Sample subsurface scattering. Use 3 nSamples here
    int nSamples = 3;
    ScatterSamplingRecord sRecRgb[nSamples];
    Vector3f refractedW = -this->po.wo; // FIXME: to the world transformation?
    Sample_Pi(sRecRgb, scene, this->po, refractedW, *mSampler, nSamples);
    Spectrum res;
    Point3f pos;
    Normal3f normal;
    int nMissed = 0;
    for (int i = 0; i < nSamples; i++) {
        if (!sRecRgb[i].isValid) {
            nMissed++;
            continue;
        }
        pos = sRecRgb[i].p;
        normal = sRecRgb[i].n;
        res += sRecRgb[i].throughout;
    }
    pi = ARENA_ALLOC(arena, SurfaceInteraction)();
    pi->p = pos;
    pi->n = normal;
    return res / (Float) nSamples;
}

} // namespace pbrt