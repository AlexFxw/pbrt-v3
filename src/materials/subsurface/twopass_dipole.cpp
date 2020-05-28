//
// Created by 范軒瑋 on 2020/5/23.
//

#include <core/scene.h>
#include "twopass_dipole.h"
#include "utils.h"
#include "shapes/triangle.h"
#include "interaction.h"
#include "sampling.h"
#include "octree.h"
#include "samplers/random.h"
#include "reflection.h"
#include "parallel.h"
#include "integrator.h"

namespace pbrt {

int TwoPassHelper::PrepareOctree(const Scene &scene, MemoryArena &arena, BSDF *bsdf) {
    std::lock_guard<std::mutex> lockGuard(buildLock);
    if (octreeBuilt) {
        return 0;
    }
    std::cout << "Preparing Octree" << std::endl;
    std::vector<IrradianceData> data;
#ifdef VISUALIZE_OCTREE
    const std::string octreeData = "../data/octree.txt";
    std::ofstream file;
    file.open(octreeData);
    CHECK(file.is_open());
#endif

    // std::cout << "calculating: " << i << std::endl;
    Float irrSamplesInv = 1.0f / (Float)irrSamples;
    Float sa = areaSum / (Float)nSamples;
    for (int i = 0; i < nSamples; i++) {
        Float pdf;
        Interaction &sP = sampledP[i];
        int trigIndx = indices[i];
        // TODO: Sample a direction of tmp
        Spectrum E = Spectrum(0.0f);
        // Initialize BSDF?
        SurfaceInteraction tmp;
        tmp.p = sP.p;
        tmp.shading.n = sP.n;
        tmp.n = sP.n;
        tmp.time = sP.time;
        Vector3f ss, ts;
        pbrt::CoordinateSystem(Vector3f(tmp.n), &ss, &ts);
        tmp.shading.dpdu = ss;
        tmp.bsdf = bsdf;

#ifdef VISUALIZE_OCTREE
        const auto &p = sP.p;
        // std::cout << p << std::endl;
        file << p.x << " " << p.y << " " << p.z << " " << 1 << " " << 0 << " " << 0 << std::endl;
        bool visualizedWo = false;
#endif
        // FIXME: Should use a new sampler?
        RandomSampler randomSampler(irrSamples << 3);
        randomSampler.StartPixel(Point2i());

        for (int irr = 0; irr < irrSamples; irr++) {
            // TODO: Irradiance caching.
            // tmp.wo = SampleAzimuthVector(Vector3f(tmp.n), randomSampler.Get2D());

#ifdef VISUALIZE_OCTREE
            if (!visualizedWo) {
                auto woP = tmp.p + tmp.wo;
                file << woP.x << " " << woP.y << " " << woP.z << " " << 0 << " " << 1 << " " << 0 << std::endl;
                visualizedWo = true;
            }
#endif
            const auto &lights = scene.lights;
            const int index = Clamp((int) (randomSampler.Get1D() * lights.size()), 0, lights.size());
            const auto &light = *lights[index];
            Float pdf;
            Spectrum directRadiance = light.SampleDirect(scene, tmp, &pdf, randomSampler);

            {
                // Handle indirect; TODO: Add parameter to control
            }

            if (!directRadiance.IsBlack() && pdf != 0) {
                directRadiance /= pdf;
                std::cout << trigIndx << directRadiance << std::endl;
                E += directRadiance;
            }
        }
        E *= irrSamplesInv;
        data.emplace_back(tmp.p, E, sa);
    }
#ifdef VISUALIZE_OCTREE
    file.close();
#endif
    octree = std::unique_ptr<IrradianceOctree>(new IrradianceOctree());
    octree->Build(data);
    octreeBuilt = true;
    std::cout << "Prepared octree." << std::endl;
    return 0;
}

int TwoPassHelper::Prepare(const std::vector<std::shared_ptr<Shape>> &shapes) {
    std::shared_ptr<TriangleMesh> triMesh(Utils::PreprocessTriangles(shapes));
    DCHECK_EQ(triMesh->nTriangles, shapes.size());
    areaSum = triMesh->area;
    triangleSum = triMesh->nTriangles;

    for (int i = 0; i < nSamples; i++) {
        Float pdf;
        int trigIdx = triMesh->areaDistri->SampleDiscrete(GetRandomFloat());
        sampledP.push_back(shapes[trigIdx]->Sample(RandPoint2f(), &pdf));
        indices.push_back(trigIdx);
    }
    return 0;
}

Spectrum TwoPassHelper::E(const Point3f &p) {
    Spectrum E = octree->Search(p);
    // std::cout << "E: " << E << std::endl;
    return E;
}

Spectrum TwoPassHelper::IrradianceCache(const Scene &scene, const Point3f &p,
                                        Sampler &sampler, MemoryArena &arena) const {
    // TODO: Irradiance Cache
    const auto &lights = scene.lights;
    Spectrum irradiance(0.0f);
    if (!lights.empty()) {
        int index = Clamp((int) (sampler.Get1D() * lights.size()), 0, lights.size() - 1);
        const Light &light = *(lights[index]);
        Point2f u1 = sampler.Get2D(), u2 = sampler.Get2D();
        // FIXME: time not sure
        Ray ray;
        Normal3f nLight;
        Float pdfPos, pdfDir, liPdf;
        Spectrum Le = light.Sample_Le(u1, u2, 0, &ray, &nLight, &pdfPos, &pdfDir);
        Vector3f d = Normalize(p - ray.o), wi;
        Ray probeRay(ray.o, d);
        SurfaceInteraction isect;
        if (!scene.Intersect(probeRay, &isect)) {
            LOG(ERROR) << "Miss in irradiance caching";
        }
        // VisibilityTester visibilityTester;
        // Spectrum Li = light.Sample_Li(isect, sampler.Get2D(), &wi, &liPdf, &visibilityTester);
        // if(!Li.IsBlack() && liPdf != 0) {
        //     irradiance += Li / liPdf;
        // }
        irradiance += EstimateDirect(isect, sampler.Get2D(), light, sampler.Get2D(),
                                     scene, sampler, arena);
    }
    return irradiance;
}


// Spectrum
// TwoPassBSSRDF::Sample_S(const Scene &scene, Float u1, const Point2f &u2, MemoryArena &arena, SurfaceInteraction *si,
//                         Float *pdf) const {
//     // TODO: Not sure if it is right here.
//     // In two pass cache method, there is actually no need to sample a output position,
//     // So the si is the SurfaceInteraction object attached this TwoPassBSSRDF itself.
//     // ProfilePhase pp(Prof::BSSRDFSampling);
//     // *pdf = 1.0f;
//     // if (!E.IsBlack()) {
//     //     si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);
//     //     si->bsdf->Add(ARENA_ALLOC(arena, TwoPassBSSRDFAdapter)(this));
//     //     si->wo = Vector3f(si->shading.n);
//     // }
//     return E;
// }

Spectrum TwoPassBSSRDF::Sp(const SurfaceInteraction &pi) const {
    Spectrum E = twoPassHelper->E(pi.p);
    Float Ft = FrDielectric(CosTheta(po.wo), 1, eta);
    return (1 - Ft) * Mo(pi.p, E);
    // E = Mo(pi.p, E) * InvPi;
    // return E;

}

TwoPassBSSRDF::TwoPassBSSRDF(const SurfaceInteraction &po, Float eta, const Material *material,
                             TransportMode mode, const Spectrum &sigmaA, const Spectrum &sigmaS, Float g,
                             std::shared_ptr<TwoPassHelper> twoPassHelper) :
                             ClassicalBSSRDF(po, eta, material, mode, sigmaA, sigmaS, g),
                             twoPassHelper(twoPassHelper){

}

} // namespace pbrt