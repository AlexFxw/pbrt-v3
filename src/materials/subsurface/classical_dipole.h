#ifndef PBRT_V3_CLASSICAL_DIPOLE_H
#define PBRT_V3_CLASSICAL_DIPOLE_H

#include <core/interaction.h>
#include "pbrt.h"
#include "bssrdf.h"
#include "medium.h"

namespace pbrt {

class ClassicalBSSRDF : public SeparableBSSRDF {
public:
    ClassicalBSSRDF(const SurfaceInteraction &po, Float eta, const Material *material,
                    TransportMode mode, const Spectrum &sigmaA, const Spectrum &sigmaS, Float g);

    Spectrum Sp(const SurfaceInteraction &pi) const override {
        Spectrum diffuseTerm = Sr(Distance(po.p, pi.p));
        Spectrum singleTerm = SingleScatterApproximate(pi.wo);
        return diffuseTerm + singleTerm;
    }
    Spectrum Sr(Float d) const override;
    Float Sample_Sr(int ch, Float u) const override;
    Float Pdf_Sr(int ch, Float r) const override;

    static Float FresnelDiffuseReflectance(Float eta);
    Spectrum SingleScatterApproximate(const Vector3f &wi) const {
        Float Ft = (1 - FresnelDiffuseReflectance(eta));
        Float F = Ft * Ft;
        Spectrum sigmaT = mSigmaA + mSigmaS;
        Spectrum alpha = mSigmaS / sigmaT;
        Spectrum singleTerm = alpha * F * phaseFunc.p(po.wo, wi) / (AbsDot(po.wo, po.shading.n) + AbsDot(wi, po.shading.n));
        return singleTerm;
    }
protected:
    Spectrum mZv, mZr;
    Spectrum mSigmaTr;
    Spectrum mSigmaS, mSigmaA;
    Spectrum mSigmaSPrime, mSigmaTPrime;
    Float g, mFdr;
    HenyeyGreenstein phaseFunc;

};

}

#endif //PBRT_V3_CLASSICAL_DIPOLE_H
