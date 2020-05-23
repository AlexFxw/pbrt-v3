
#ifndef PBRT_V3_NORMAL_DIFFUSION_H
#define PBRT_V3_NORMAL_DIFFUSION_H

#include "pbrt.h"
#include "bssrdf.h"

namespace pbrt {

class NormalDiffusion : public SeparableBSSRDF {
public:
    NormalDiffusion(const Spectrum &d, const SurfaceInteraction &po,
                    Float eta, const Material *material, TransportMode mode);
    Spectrum Sr(Float r) const override;
    Float Sample_Sr(int ch, Float u) const override;
    Float Pdf_Sr(int ch, Float r) const override;

private:
    Spectrum d;
    constexpr static Float firstTermRatio = 0.25f, secondTermRatio = 0.75f;
};

} // namespace pbrt
#endif //PBRT_V3_NORMAL_DIFFUSION_H
