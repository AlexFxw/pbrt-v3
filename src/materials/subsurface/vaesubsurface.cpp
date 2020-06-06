#include "vaesubsurface.h"
#include "vaehandler.h"
#include "vaehandlereigen.h"
#include "vaescatter.h"
#include "texture.h"
#include "paramset.h"
#include "memory.h"

namespace pbrt {

void VAESubsurface::ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena, TransportMode mode,
                                               bool allowMultipleLobes) const {
    // Use the BSDF of subsurface material, but replace the bssrdf with vae scattering function implementation.
    SubsurfaceMaterial::ComputeScatteringFunctions(si, arena, mode, allowMultipleLobes);
    si->bssrdf = ARENA_ALLOC(arena, VaeScatter)(*si, this->eta, this->mVaeHandler, mode);
}

VAESubsurface *CreateVaeSubsurfaceMaterial(const TextureParams &mp) {
    // Reuse CreateSubsurfaceMaterial code.
    Float sig_a_rgb[3] = {.0011f, .0024f, .014f},
            sig_s_rgb[3] = {2.55f, 3.21f, 3.77f};
    Spectrum sig_a = Spectrum::FromRGB(sig_a_rgb),
            sig_s = Spectrum::FromRGB(sig_s_rgb);
    std::string name = mp.FindString("name");
    bool found = GetMediumScatteringProperties(name, &sig_a, &sig_s);
    Float g = mp.FindFloat("g", 0.f);
    Float kernelEpsScale = mp.FindFloat("kernelEpsScale", 1.0f);
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
    VAESubsurface *vaeSubsurface = new VAESubsurface(scale, Kr, Kt, sigma_a, sigma_s, g, eta,
                                                     roughu, roughv, bumpMap, remapRoughness);
    // FIXME: Configure sigmaT, albedo, etc better
    Float rgb[3] = {52.02, 78.54, 109.12};
    // Spectrum sigmaT = Spectrum::FromRGB(rgb);
    Spectrum sigmaT = sig_a + sig_s;
    Spectrum albedo(0.99f); // Configure albedo.

    std::string modelName = mp.FindString("model_name", "0487_FinalSharedLs7Mixed3_AbsSharedSimComplexMixed3");
    std::string absModelName = modelName;
    std::string outputDir = "data/";

    vaeSubsurface->mVaeHandler = std::make_shared<VaeHandlerEigen>(
            sigmaT, albedo, g, eta, modelName, absModelName, absModelName, outputDir, 64, kernelEpsScale
    );
    return vaeSubsurface;
}

void VAESubsurface::PrepareMaterial(const std::vector<std::shared_ptr<Shape>> &shapes,
                                    const ParamSet &params) {
    Material::PrepareMaterial(shapes, params);
    ProfilePhase p(Prof::MaterialPreparation);
    PolyUtils::PolyFitConfig pfConfig; // FIXME: I just use default here, remember to make use of params
    mVaeHandler->Prepare(shapes, pfConfig);
    LOG(INFO) << "Finish preparing the material.";
}


}
