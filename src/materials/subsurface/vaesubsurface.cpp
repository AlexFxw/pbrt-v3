//
// Created by 范軒瑋 on 2020/5/5.
//

#include "vaesubsurface.h"
#include "vaescatter.h"
#include "memory.h"

namespace pbrt {

void VAESubsurface::ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena, TransportMode mode,
                                               bool allowMultipleLobes) const {
    // TODO
    // Use the BSDF of subsurface material, but replace the bssrdf with vae scattering function implementation.
    SubsurfaceMaterial::ComputeScatteringFunctions(si, arena, mode, allowMultipleLobes);
    si->bssrdf = ARENA_ALLOC(arena,VaeScatter)(*si, eta);
}

}
