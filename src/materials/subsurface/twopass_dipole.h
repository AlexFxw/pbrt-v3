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

namespace pbrt {


class TwoPassHelper {
public:
    TwoPassHelper() : octreeBuilt(false) {}

    int PrepareOctree(const Scene &scene, MemoryArena &arena, BSDF *bsdf);
    int Prepare(const std::vector<std::shared_ptr<Shape>> &shapes);

    Spectrum E(const Point3f &p);

    bool Prepared() const { return octreeBuilt; }


private:
    Spectrum IrradianceCache(const Scene &scene, const Point3f &p,
                             Sampler &sampler, MemoryArena &arena) const; // Photon mapping.
    std::unique_ptr<IrradianceOctree> octree = nullptr;
    std::vector<Interaction> sampledP;
    std::vector<int> indices;
    int triangleSum, areaSum;
    const int nSamples = 10000;
    const static int irrSamples = 64;
    bool octreeBuilt;
    std::mutex buildLock;
};

class TwoPassBSSRDF : public ClassicalBSSRDF {
public:
    // TwoPassBSSRDF(const SurfaceInteraction &po, Float eta,
    //               std::shared_ptr<TwoPassHelper> twoPassHelper) :
    //         BSSRDF(po, eta), twoPassHelper(twoPassHelper) {

    // }
    TwoPassBSSRDF(const SurfaceInteraction &po, Float eta, const Material *material,
                  TransportMode mode, const Spectrum &sigmaA, const Spectrum &sigmaS,
                  Float g,std::shared_ptr<TwoPassHelper> twoPassHelper);

    Spectrum Sp(const SurfaceInteraction &pi) const override;

    // Spectrum S(const SurfaceInteraction &pi, const Vector3f &wi) override;
    // Spectrum Sample_S(const Scene &scene, Float u1, const Point2f &u2,
    //                   MemoryArena &arena, SurfaceInteraction *si,
    //                   Float *pdf) const override;

    // Spectrum Sr(Float d) const override;
    bool UseCacheCloud() const override { return true; }

    bool Prepared() const override { return twoPassHelper->Prepared(); }

    std::shared_ptr<TwoPassHelper> twoPassHelper;

private:

    Spectrum Mo(const Point3f &pi, const Spectrum &influxI) const {
        // Though call influx, it should be already multiplied by the area here.
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

        return dMo * influxI;
    }
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
