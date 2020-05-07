//
// Created by 范軒瑋 on 2020/5/5.
//

#include "vaesubsurface.h"
#include "vaescatter.h"
#include "texture.h"
#include "paramset.h"
#include "memory.h"

namespace pbrt {

void VAESubsurface::ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena, TransportMode mode,
                                               bool allowMultipleLobes) const {
    // TODO
    // Use the BSDF of subsurface material, but replace the bssrdf with vae scattering function implementation.
    SubsurfaceMaterial::ComputeScatteringFunctions(si, arena, mode, allowMultipleLobes);
    si->bssrdf = ARENA_ALLOC(arena,VaeScatter)(*si, eta);
}

VAESubsurface *CreateVaeSubsurfaceMaterial(const TextureParams &mp)
{
    // TODO: Create a vae subsurface material
    Float sig_a_rgb[3] = {.0011f, .0024f, .014f},
            sig_s_rgb[3] = {2.55f, 3.21f, 3.77f};
    Spectrum sig_a = Spectrum::FromRGB(sig_a_rgb),
            sig_s = Spectrum::FromRGB(sig_s_rgb);
    std::string name = mp.FindString("name");
    bool found = GetMediumScatteringProperties(name, &sig_a, &sig_s);
    Float g = mp.FindFloat("g", 0.0f);
    if (name != "") {
        if (!found)
            Warning("Named material \"%s\" not found.  Using defaults.",
                    name.c_str());
        else
            g = 0; /* Enforce g=0 (the database specifies reduced scattering
                      coefficients) */
    }
    Float scale = mp.FindFloat("scale", 1.f);
    Float eta = mp.FindFloat("eta", 1.33f);

    std::shared_ptr<Texture<Spectrum>> sigma_a, sigma_s;
    sigma_a = mp.GetSpectrumTexture("sigma_a", sig_a);
    sigma_s = mp.GetSpectrumTexture("sigma_s", sig_s);
    std::shared_ptr<Texture<Spectrum>> Kr =
            mp.GetSpectrumTexture("Kr", Spectrum(1.f));
    std::shared_ptr<Texture<Spectrum>> Kt =
            mp.GetSpectrumTexture("Kt", Spectrum(1.f));
    std::shared_ptr<Texture<Float>> roughu =
            mp.GetFloatTexture("uroughness", 0.f);
    std::shared_ptr<Texture<Float>> roughv =
            mp.GetFloatTexture("vroughness", 0.f);
    std::shared_ptr<Texture<Float>> bumpMap =
            mp.GetFloatTextureOrNull("bumpmap");
    bool remapRoughness = mp.FindBool("remaproughness", true);
    return new VAESubsurface(scale, Kr, Kt, sigma_a, sigma_s, g, eta,
                                  roughu, roughv, bumpMap, remapRoughness);
}


}
