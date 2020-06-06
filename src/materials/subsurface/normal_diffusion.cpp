
#include "normal_diffusion.h"
#include "stats.h"
#include "spectrum.h"
#include "rng.h"

namespace pbrt {

Float SmallRFilter(Float r) {
    if(r < 1e-6f) {
        return 1e-6f;
    }
    return r;
}

NormalDiffusion::NormalDiffusion(const Spectrum &R, const Spectrum &d, const SurfaceInteraction &po, Float eta,
                                 const Material *material, TransportMode mode):
                                 SeparableBSSRDF(po, eta, material, mode), d(d), R(R){

}

Float NormalDiffusion::Pdf_Sr(int ch, Float r) const {
    r = SmallRFilter(r);
    return (firstTermRatio * std::exp(-r / d[ch]) / (2 * Pi * d[ch] * r) +
            secondTermRatio * std::exp(-r / (3 * d[ch])) / (6 * Pi * d[ch] * r));
}

Float NormalDiffusion::Sample_Sr(int ch, Float u) const {
    if (u < firstTermRatio) {
        // Sample the first exponential term
        u = std::min<Float>(u * 4, OneMinusEpsilon);
        return d[ch] * std::log(1 / (1 - u));
    } else {
        u = std::min<Float>((u - firstTermRatio) / secondTermRatio, OneMinusEpsilon);
        return 3 * d[ch] * std::log(1 / (1 - u));
    }
}

Spectrum NormalDiffusion::Sr(Float r) const {
    ProfilePhase pp(Prof::BSSRDFEvaluation);
    r = SmallRFilter(r);
    Spectrum spectrumR = -Spectrum(r);
    return R * (Exp(spectrumR / this->d) + Exp(spectrumR / (3 * d))) / (8 * Pi * d * r);
}


}
