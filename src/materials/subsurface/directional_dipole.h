#ifndef PBRT_V3_DIRECTIONAL_DIPOLE_H
#define PBRT_V3_DIRECTIONAL_DIPOLE_H

#include "pbrt.h"
#include "bssrdf.h"

namespace pbrt {

class DirectionalBSSRDF : public SeparableBSSRDF {
public:
    DirectionalBSSRDF(const SurfaceInteraction &po, Float eta,
                      const Material *material, TransportMode mode,
                      const Spectrum &sigmaA, const Spectrum &sigmaS, Float g);
    Spectrum Sr(Float d) const override ;
    Float Sample_Sr(int ch, Float u) const override ;
    Float Pdf_Sr(int ch, Float r) const override ;
    Spectrum Sp(const SurfaceInteraction &pi) const override;

private:
    Spectrum sigmaT, sigmaSp, sigmaTp, albedoP, sigmaTr;
    Spectrum D, de;
    Float cPNorm, cP, cE, A;

    inline Float C1(const Float n) {
        Float r;
        if (n > 1.0) {
            r = -9.23372 + n * (22.2272 + n * (-20.9292 + n * (10.2291 + n * (-2.54396 + 0.254913 * n))));
        } else {
            r = 0.919317 + n * (-3.4793 + n * (6.75335 + n *  (-7.80989 + n *(4.98554 - 1.36881 * n))));
        }
        return r / 2.0;
    }
    inline Float C2(const Float n) {
        Float r = -1641.1 + n * (1213.67 + n * (-568.556 + n * (164.798 + n * (-27.0181 + 1.91826 * n))));
        r += (((135.926 / n) - 656.175) / n + 1376.53) / n;
        return r / 3.0;
    }

    Float Sp_d(const Vector3f &x, const Vector3f &w, Float r, const Vector3f &n, int channel) const;
    Spectrum Sd(const SurfaceInteraction &po,const SurfaceInteraction &pi) const;
};

}

#endif //PBRT_V3_DIRECTIONAL_DIPOLE_H
