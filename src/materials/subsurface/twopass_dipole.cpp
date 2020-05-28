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
#include "integrators/path.h"

namespace pbrt {

int TwoPassHelper::PrepareOctree(const Scene &scene, MemoryArena &arena, const PathIntegrator *intergrater) {
    if (octreeBuilt) {
        return 0;
    }
    std::cout << "Preparing Octree" << std::endl;
    std::vector<IrradianceData> data;
    Float sa = areaSum / (Float) nSamples;
#ifdef VISUALIZE_OCTREE
    const std::string octreeData = "../data/octree.txt";
    std::ofstream file;
    file.open(octreeData);
    CHECK(file.is_open());
    for (int i = 0; i < nSamples; i++) {
        Interaction &sP = sampledP[i];
        const auto &p = sP.p;
        // std::cout << p << std::endl;
        file << p.x << " " << p.y << " " << p.z << " " << 1 << " " << 0 << " " << 0 << std::endl;
        bool visualizedWo = false;
#else
    for (int i = 0; i < nSamples; i++) {
        Interaction &sP = sampledP[i];
#endif
        // FIXME: Should use a new sampler?
        RandomSampler randomSampler(irrSamples << 3);
        randomSampler.StartPixel(Point2i());
        Spectrum irradiance = intergrater->Irradiance(scene, sP.p, Vector3f(sP.n), randomSampler, arena, irrSamples);
        data.emplace_back(sP.p, irradiance, sa);
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
    std::cout << "Area: " << areaSum << " " << triangleSum << std::endl;

    for (int i = 0; i < nSamples; i++) {
        Float pdf;
        int trigIdx = triMesh->areaDistri->SampleDiscrete(GetRandomFloat());
        sampledP.push_back(shapes[trigIdx]->Sample(RandPoint2f(), &pdf));
        indices.push_back(trigIdx);
    }
    return 0;
}

Spectrum TwoPassHelper::E(const Point3f &p, const TwoPassBSSRDF *bssrdf) {
    Spectrum E = octree->Search(p, bssrdf);
    return E;
}

// Spectrum TwoPassBSSRDF::Sp(const SurfaceInteraction &pi) const {
//     Spectrum E = twoPassHelper->E(pi.p, this);
//     return E;
// }

TwoPassBSSRDF::TwoPassBSSRDF(const SurfaceInteraction &po, Float eta, const Material *material,
                             TransportMode mode, const Spectrum &sigmaA, const Spectrum &sigmaS, Float g,
                             std::shared_ptr<TwoPassHelper> twoPassHelper) :
        ClassicalBSSRDF(po, eta, material, mode, sigmaA, sigmaS, g),
        twoPassHelper(twoPassHelper) {

}

Spectrum
TwoPassBSSRDF::Sample_S(const Scene &scene, Float u1, const Point2f &u2, MemoryArena &arena, SurfaceInteraction *si,
                        Float *pdf) const {

    Float Ft = FrDielectric(CosTheta(po.wo), 1, eta);
    Spectrum E = twoPassHelper->E(po.p, this);
    Float Fdr = FresnelDiffuseReflectance(eta);
    *pdf = 1.0f;
    // return ((1 - Ft) / Fdr) * (E / Pi);
    return (1 - Ft) * E * InvPi;
}

} // namespace pbrt