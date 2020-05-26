//
// Created by 范軒瑋 on 2020/5/23.
//

#ifndef PBRT_V3_TWOPASS_DIPOLE_H
#define PBRT_V3_TWOPASS_DIPOLE_H

#include "pbrt.h"
#include "geometry.h"
#include "octree.h"
#include "bssrdf.h"
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
    int triangleSum;
    const int nSamples = 10000;
    const static int irrSamples = 64;
    bool octreeBuilt;
    std::mutex buildLock;
};

class TwoPassBSSRDF : public BSSRDF {
public:
    TwoPassBSSRDF(const SurfaceInteraction &po, Float eta,
                  std::shared_ptr<TwoPassHelper> twoPassHelper) :
            BSSRDF(po, eta), twoPassHelper(twoPassHelper) {

    }

    Spectrum S(const SurfaceInteraction &pi, const Vector3f &wi) override;
    Spectrum Sample_S(const Scene &scene, Float u1, const Point2f &u2,
                      MemoryArena &arena, SurfaceInteraction *si,
                      Float *pdf) const override;

    bool UseCacheCloud() const override { return true; }

    bool Prepared() const override { return twoPassHelper->Prepared(); }

    std::shared_ptr<TwoPassHelper> twoPassHelper;
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
