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
    VaeScatter();  // TODO: implement
    ~VaeScatter();
    Spectrum S(const SurfaceInteraction &pi, const Vector3f &wi);
    Spectrum Sample_S(const Scene &scene, Float u1, const Point2f &u2,
                      MemoryArena &arena, SurfaceInteraction *si,
                      Float *pdf) const;
    void Prepare(const Scene *scene, const std::vector<std::shared_ptr<Shape>> &shapes);
    void Sample_Pi(ScatterSamplingRecord *sRecs, const Scene &scene, const SurfaceInteraction &itact,
                   const Vector3f &w, Sampler &sampler, int nSamples) const;
    Spectrum Sample_Sp(const Scene &scene, MemoryArena &arena, SurfaceInteraction *pi, Float *pdf) const;
private:
    VaeHandler *mVaeHandler;
    std::vector<std::shared_ptr<Shape>> mTriangles;
    // TODO: Initialize those member variables.
    Float eta, g;
    Spectrum mSigmaT, mAlbedo;
    MediumParameters mMedium;
    Sampler *mSampler;
};

};  // namespace pbrt

#endif  // PBRT_MATERIALS_SUBSURFACE_VAESCATTER_H