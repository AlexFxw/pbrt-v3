//
// Created by 范軒瑋 on 2020/5/23.
//

#ifndef PBRT_V3_TWOPASS_DIPOLE_H
#define PBRT_V3_TWOPASS_DIPOLE_H

#include "pbrt.h"
#include "geometry.h"
#include "octree.h"
#include "classical_dipole.h"
#include <mutex>
#include "medium.h"


namespace pbrt {


class TwoPassHelper {
public:
    TwoPassHelper() : octreeBuilt(false) {}

    int PrepareOctree(const Scene &scene, MemoryArena &arena, const PathIntegrator *integrator);
    int Prepare(const std::vector<std::shared_ptr<Shape>> &shapes);

    Spectrum E(const Point3f &p, const TwoPassBSSRDF *bssrdf);

    bool Prepared() const { return octreeBuilt; }


private:
    std::unique_ptr<IrradianceOctree> octree = nullptr;
    std::vector<Interaction> sampledP;
    std::vector<int> indices;
    int triangleSum;
    double areaSum;
    const int nSamples = 100000;
    const static int irrSamples = 4;
    bool octreeBuilt;
};

class TwoPassBSSRDF : public ClassicalBSSRDF {
public:
    TwoPassBSSRDF(const SurfaceInteraction &po, Float eta, const Material *material,
                  TransportMode mode, const Spectrum &sigmaA, const Spectrum &sigmaS,
                  Float g, std::shared_ptr<TwoPassHelper> twoPassHelper, const Spectrum &R);

    // Spectrum Sp(const SurfaceInteraction &pi) const override;

    // FIXME: It looks like just use Sample_S is OK?
    Spectrum Sample_S(const Scene &scene, Float u1, const Point2f &u2,
                      MemoryArena &arena, SurfaceInteraction *si,
                      Float *pdf) const override;

    bool UseCacheCloud() const override { return true; }

    bool Prepared() const override { return twoPassHelper->Prepared(); }

    void SetHelper(const std::shared_ptr<TwoPassHelper> &helper) override {
        twoPassHelper = std::shared_ptr<TwoPassHelper>(helper);
    }

    std::shared_ptr<TwoPassHelper> twoPassHelper;

    Spectrum Mo(const Point3f &pi, const Spectrum &influxI) const {
        Spectrum rSqr = Spectrum((pi - po.p).LengthSquared());

        // Distance to real source
        Spectrum dr = Sqrt(rSqr + mZr * mZr);
        // Distance to virtual source.
        Spectrum dv = Sqrt(rSqr + mZv * mZv);

        Spectrum C1 = mZr * (mSigmaTr + Spectrum(1.0f) / dr);
        Spectrum C2 = mZv * (mSigmaTr + Spectrum(1.0f) / dv);

        // Do not include the reduced albedo
        Spectrum dMo = Inv4Pi *
                       (C1 * (Exp(-mSigmaTr * dr)) / (dr * dr)
                        + C2 * (Exp(-mSigmaTr * dv)) / (dv * dv));

        Float Fdt = 1.0f - FresnelDiffuseReflectance(eta);
        return Fdt * dMo * influxI;
    }

private:
    Spectrum R;
};

class TwoPassBSSRDFAdapter : public BxDF {
public:
    TwoPassBSSRDFAdapter(const TwoPassBSSRDF *bssrdf) :
            BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), bssrdf(bssrdf) {}

    Spectrum f(const Vector3f &wo, const Vector3f &wi) const override {
        // FIXME:
        return Spectrum(0.0f);
    }

    std::string ToString() const { return "[ TwoPassBSSRDFAdapter ]"; }

private:
    const TwoPassBSSRDF *bssrdf;
};


} // namespace pbrt

#endif //PBRT_V3_TWOPASS_DIPOLE_H
