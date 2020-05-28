#ifndef PBRT_V3_CLASSICAL_DIPOLE_H
#define PBRT_V3_CLASSICAL_DIPOLE_H

#include "pbrt.h"
#include "bssrdf.h"

namespace pbrt {

class ClassicalBSSRDF : public SeparableBSSRDF {
public:
    ClassicalBSSRDF(const SurfaceInteraction &po, Float eta, const Material *material,
                    TransportMode mode, const Spectrum &sigmaA, const Spectrum &sigmaS, Float g);

    Spectrum Sr(Float d) const override;
    Float Sample_Sr(int ch, Float u) const override;
    Float Pdf_Sr(int ch, Float r) const override;

    static Float FresnelDiffuseReflectance(Float eta);
protected:
    Spectrum mZv, mZr;
    Spectrum mSigmaTr;
    Spectrum mSigmaS, mSigmaA;
    Spectrum mSigmaSPrime, mSigmaTPrime;
    Float g, mFdr;

};

}

#endif //PBRT_V3_CLASSICAL_DIPOLE_H
