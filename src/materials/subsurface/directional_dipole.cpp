#include <core/interaction.h>
#include "directional_dipole.h"
#include "spectrum.h"

namespace pbrt {


Float DirectionalBSSRDF::Pdf_Sr(int ch, Float r) const {
    // TODO
    return sigmaTr[ch] * std::exp(-sigmaTr[ch] * r);
}

Float DirectionalBSSRDF::Sample_Sr(int ch, Float u) const {
    // TODO
    Float d = 1.0f / sigmaTr[ch];
    return u * d;
}

Spectrum DirectionalBSSRDF::Sr(Float d) const {
    // TODO
    // Doesn't matter here because I rewrite the Sp
    return Spectrum();
}

DirectionalBSSRDF::DirectionalBSSRDF(const SurfaceInteraction &po, Float eta, const Material *material,
                                     TransportMode mode, const Spectrum &sigmaA, const Spectrum &sigmaS,
                                     Float g) :
        SeparableBSSRDF(po, eta, material, mode) {
    sigmaT = sigmaA + sigmaS;
    sigmaSp = sigmaS * (1.0f - g);
    sigmaTp = sigmaSp + sigmaA;
    albedoP = sigmaSp / sigmaTp;
    D = Spectrum(1.0f) / (3.0 * sigmaTp);
    sigmaTr = Sqrt(sigmaA / D);
    de = 2.131 * D / Sqrt(albedoP);
    cPNorm = 1.0 / (1.0 - 2.0 * C1(1.0 / eta));
    cP = (1.0 - 2.0 * C1(eta)) / 4.0;
    cE = (1.0 - 3.0 * C2(eta)) / 2.0;
    A = (1.0 - cE) / (2.0 * cP);
}

Spectrum DirectionalBSSRDF::Sd(const SurfaceInteraction &po, const SurfaceInteraction &pi) const {
    const Point3f &xo = po.p, &xi = pi.p;
    const Vector3f &wo = po.wo, wi = -pi.wo;
    const Vector3f no = Vector3f(po.shading.n), ni = Vector3f(pi.shading.n);
    const Vector3f xoxi = xo - xi;
    const Float r = xoxi.Length(), r2 = r * r;
    const Vector3f niS = Cross(Normalize(xoxi), Normalize(Cross(ni, xoxi)));

    const Float nnt = 1.0 / eta, ddn = -Dot(wi, ni);
    const Vector3f wr = Normalize(wi * -nnt - ni * (ddn * nnt + std::sqrt(1.0 - nnt * nnt * (1.0 - ddn * ddn))));
    const Vector3f wv = wr - niS * (2.0 * Dot(wr, niS));


    Spectrum res(0.0f);
    for(int i = 0; i < Spectrum::nSamples; i ++) {
        // Distance to real sources
        Float dotXoxiWr = Dot(xoxi, wr);
        const Float cosBeta = -std::sqrt((r2 - dotXoxiWr * dotXoxiWr) / (r2 + de[i] * de[i]));
        Float dr;
        const Float mu0 = -Dot(no, wr);
        if(mu0 > 0.0f) {
            dr = std::sqrt((D[i] * mu0) * ((D[i] * mu0) - de[i] * cosBeta * 2.0) + r2);
        } else {
            dr = std::sqrt(1.0 / (3.0 * sigmaT[i] * 3.0 * sigmaT[i]) + r2);
        }

        // Distance to virtual light source.
        const Vector3f xoxv = xo - (xi + niS * (2.0 * A * de[i]));
        const Float dv = xoxv.Length();

        const Float iRes = Sp_d(xoxi, wr, dr, no, i) - Sp_d(xoxv, wv, dv, no, i);
        if(iRes > 0.0f)
            res[i] = iRes;
    }
    return res;
}

Float DirectionalBSSRDF::Sp_d(const Vector3f &x, const Vector3f &w, Float r, const Vector3f &n, int channel) const {
    const Float sTrR = sigmaTr[channel] * r;
    const Float sTrROne = 1.0f + sTrR;
    const Float xDotW = Dot(x, w);
    const Float rSqr = r * r;
    const Float t0 = cPNorm * (1.0 * Inv4Pi * InvPi) * std::exp(-sTrR) / (r * rSqr);
    const Float t1 = rSqr / D[channel] + 3.0 * sTrROne * xDotW;
    const Float t2 = 3.0 * D[channel] * sTrROne * xDotW;
    const Float t3 = (sTrROne + 3.0 * D[channel] * (3.0 * sTrROne + sTrR * sTrR) / rSqr * xDotW) * Dot(x, n);
    return t0 * (cP * t1 - cE * (t2 - t3));
}

Spectrum DirectionalBSSRDF::Sp(const SurfaceInteraction &pi) const {
    // FIXME: Reasonable?
    return Sd(pi, po);
    // return Sd(po, pi);
}


}