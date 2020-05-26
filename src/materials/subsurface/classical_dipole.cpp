#include "classical_dipole.h"
#include "spectrum.h"

namespace pbrt {

Float FresnelDiffuseReflectance(Float eta) {
    /* Fast mode: the following code approximates the
     * diffuse Frensel reflectance for the eta<1 and
     * eta>1 cases. An evalution of the accuracy led
     * to the following scheme, which cherry-picks
     * fits from two papers where they are best.
     * source: https://github.com/RoaldFre/mitsuba-ALVRL-fwddip
     */
    if (eta < 1) {
        /* Fit by Egan and Hilgeman (1973). Works
           reasonably well for "normal" IOR values (<2).

           Max rel. error in 1.0 - 1.5 : 0.1%
           Max rel. error in 1.5 - 2   : 0.6%
           Max rel. error in 2.0 - 5   : 9.5%
        */
        return -1.4399f * (eta * eta)
               + 0.7099f * eta
               + 0.6681f
               + 0.0636f / eta;
    } else {
        /* Fit by d'Eon and Irving (2011)
         *
         * Maintains a good accuracy even for
         * unrealistic IOR values.
         *
         * Max rel. error in 1.0 - 2.0   : 0.1%
         * Max rel. error in 2.0 - 10.0  : 0.2%
         */
        Float invEta = 1.0f / eta,
                invEta2 = invEta * invEta,
                invEta3 = invEta2 * invEta,
                invEta4 = invEta3 * invEta,
                invEta5 = invEta4 * invEta;

        return 0.919317f - 3.4793f * invEta
               + 6.75335f * invEta2
               - 7.80989f * invEta3
               + 4.98554f * invEta4
               - 1.36881f * invEta5;
    }
}

ClassicalBSSRDF::ClassicalBSSRDF(const SurfaceInteraction &po, Float eta, const Material *material, TransportMode mode,
                                 const Spectrum &sigmaA, const Spectrum &sigmaS, Float g) :
        SeparableBSSRDF(po, eta, material, mode), mSigmaA(sigmaA), mSigmaS(sigmaS), g(g) {
    mSigmaSPrime = mSigmaS * (1.0f - g);
    mSigmaTPrime = mSigmaSPrime + mSigmaA;
    Spectrum meanFreePath = Spectrum(1.0f) / mSigmaTPrime;
    mFdr = FresnelDiffuseReflectance(1.0f / eta);
    Float A = (1.0f + mFdr) / (1.0f - mFdr);

    mSigmaTr = Sqrt(mSigmaA * mSigmaTPrime * 3.0f);
    if(mSigmaTr.IsBlack()) {
        LOG(ERROR) << "Classical BSSRDF needs nonzero  sigma Tr";
    }
    mZr = meanFreePath;
    mZv = meanFreePath * (1.0f + 4.0f/3.0f * A);
}

Float ClassicalBSSRDF::Pdf_Sr(int ch, Float r) const {
    return mSigmaTr[ch] * std::exp(-mSigmaTr[ch] * r);
}

Float ClassicalBSSRDF::Sample_Sr(int ch, Float u) const {
    Float d = 1.0f / mSigmaTr[ch];
    return u * d;
}

Spectrum ClassicalBSSRDF::Sr(Float d) const {
    Spectrum spectrumD = Spectrum(d);
    Spectrum dr = Sqrt(spectrumD + mZr * mZr);
    Spectrum dv = Sqrt(spectrumD + mZv * mZv);

    Spectrum C1 = mZr * (mSigmaTr + Spectrum(1.0f) / dr);
    Spectrum C2 = mZv * (mSigmaTr + Spectrum(1.0f) / dv);

    Spectrum dMo = Inv4Pi *
                   ((C1 * Exp(-mSigmaTr * dr) / (dr * dr))
                    + (C2 * Exp(-mSigmaTr * dv) / (dv * dv)));

    for (int i = 0; i < Spectrum::nSamples; i++) {
        if (!std::isfinite(mZr[i])) {
            dMo[i] = 0.0f;
        }
    }

    return dMo;
}


} // namespace pbrt
