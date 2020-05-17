/*
 * @Author: Hsuan-Wei Fan
 * @Date: 2020-04-23 13:47:05
 * @LastEditors: Hsuan-Wei Fan
 * @LastEditTime: 2020-04-23 22:20:53
 * @Description: The declaration of vae subsurface scattering.
 */

#ifndef PBRT_MATERIALS_SUBSURFACE_VAESCATTER_H
#define PBRT_MATERIALS_SUBSURFACE_VAESCATTER_H

#include "pbrt.h"
#include "bssrdf.h"
#include "vaehandler.h"

namespace pbrt {

class VaeScatterAdapter;

class VaeScatter : public BSSRDF {
public:
    friend class VaeScatterAdapter;

    VaeScatter(const SurfaceInteraction &po, Float eta,
               const std::shared_ptr<VaeHandler> &vaeHandler, TransportMode mode) :
            BSSRDF(po, eta), mode(mode) {
        mVaeHandler = vaeHandler;
        // eta = vaeHandler->GetMedium().eta;
        this->eta = eta;
        g = vaeHandler->GetMedium().g;
    }

    ~VaeScatter();
    Spectrum S(const SurfaceInteraction &pi, const Vector3f &wi) override;
    Spectrum Sample_S(const Scene &scene, Float u1, const Point2f &u2, MemoryArena &arena,
                      SurfaceInteraction *si, Float *pdf) const override;

    Spectrum Sw(const Vector3f &w) const {
        // Use standard method for outgoing direction
        Float c = 1 - 2 * FresnelMoment1(1 / eta);
        return (1 - FrDielectric(CosTheta(w), 1, eta)) / (c * Pi);
    }


    void Sample_Pi(ScatterSamplingRecord *sRecs, const Scene &scene, SurfaceInteraction *itact,
                   const Vector3f &w, int nSamples) const;

    Spectrum Sample_Sp(const Scene &scene, const Vector3f &refractedD, SurfaceInteraction *resIsect,
                       Float *pdf, int nSamples) const;

private:
    mutable std::shared_ptr<VaeHandler> mVaeHandler;
    Float eta, g;
    const TransportMode mode;
};

class VaeScatterAdapter : public BxDF {
public:
    VaeScatterAdapter(const VaeScatter *bssrdf)
            : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), bssrdf(bssrdf) {}

    Spectrum f(const Vector3f &wo, const Vector3f &wi) const override {
        Spectrum f = bssrdf->Sw(wi);
        // FIXME: No Sw now.
        if (bssrdf->mode == TransportMode::Radiance)
            f *= bssrdf->eta * bssrdf->eta;
        return f;
    }

    std::string ToString() const override { return "[ VaeScatterBSSRDFAdapter ]"; }

private:
    const VaeScatter *bssrdf;
};

};  // namespace pbrt

#endif  // PBRT_MATERIALS_SUBSURFACE_VAESCATTER_H