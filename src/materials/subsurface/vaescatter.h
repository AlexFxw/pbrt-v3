/*
 * @Author: Hsuan-Wei Fan
 * @Date: 2020-04-23 13:47:05
 * @LastEditors: Hsuan-Wei Fan
 * @LastEditTime: 2020-04-23 22:20:53
 * @Description: The declaration of vae subsurface scattering.
 */

#ifndef PBRT_MATERIALS_SUBSURFACE_VAESCATTER_H
#define PBRT_MATERIALS_SUBSURFACE_VAESCATTER_H

#include "pbrt.h"
#include "bssrdf.h"

namespace pbrt {

class VaeScatter : public BSSRDF {
public:
    VaeScatter(const SurfaceInteraction &po, Float eta,
               const std::shared_ptr<VaeHandler> &vaeHandler) : BSSRDF(po, eta) {
        mVaeHandler = vaeHandler;
    }

    ~VaeScatter();
    Spectrum S(const SurfaceInteraction &pi, const Vector3f &wi) override;
    Spectrum Sample_S(const Scene &scene, Float u1, const Point2f &u2,
                      MemoryArena &arena, SurfaceInteraction *si,
                      Float *pdf) const override;
    void Sample_Pi(ScatterSamplingRecord *sRecs, const Scene &scene, const SurfaceInteraction &itact,
                   const Vector3f &w, Sampler &sampler, int nSamples) const;
    Spectrum Sample_Sp(const Scene &scene, ScatterSamplingRecord *sRecs, Float *pdf, int nSamples) const;

private:
    mutable std::shared_ptr<VaeHandler> mVaeHandler;
    std::vector<std::shared_ptr<Shape>> mTriangles;
    // TODO: Initialize those member variables.
    Float eta, g;
    Spectrum mSigmaT, mAlbedo;
    MediumParameters mMedium;
    Sampler *mSampler;
};

};  // namespace pbrt

#endif  // PBRT_MATERIALS_SUBSURFACE_VAESCATTER_H