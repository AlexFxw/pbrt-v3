//
// Created by 范軒瑋 on 2020/5/5.
//

#ifndef PBRT_V3_VAESUBSURFACE_H
#define PBRT_V3_VAESUBSURFACE_H

#include "pbrt.h"
#include "materials/subsurface.h"

namespace pbrt {

class VAESubsurface : public SubsurfaceMaterial {
public:
    VAESubsurface(Float scale,
                  const std::shared_ptr<Texture<Spectrum>> &Kr,
                  const std::shared_ptr<Texture<Spectrum>> &Kt,
                  const std::shared_ptr<Texture<Spectrum>> &sigma_a,
                  const std::shared_ptr<Texture<Spectrum>> &sigma_s,
                  Float g, Float eta,
                  const std::shared_ptr<Texture<Float>> &uRoughness,
                  const std::shared_ptr<Texture<Float>> &vRoughness,
                  const std::shared_ptr<Texture<Float>> &bumpMap,
                  bool remapRoughness) :
            SubsurfaceMaterial(scale, Kr, Kt, sigma_a, sigma_s,
                               g, eta, uRoughness, vRoughness, bumpMap, remapRoughness) {}


    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode, bool allowMultipleLobes) const;

};

}

#endif //PBRT_V3_VAESUBSURFACE_H
