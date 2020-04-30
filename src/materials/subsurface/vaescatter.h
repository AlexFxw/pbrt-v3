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
#include "material.h"

namespace pbrt {

class VaeScatter : public Material {
  public:
    VaeScatter();  // TODO: implement
    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes);
};

};  // namespace pbrt

#endif  // PBRT_MATERIALS_SUBSURFACE_VAESCATTER_H