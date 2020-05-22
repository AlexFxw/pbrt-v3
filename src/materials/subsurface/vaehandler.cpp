/*
 * @Author: Hsuan-Wei Fan
 * @Date: 2020-04-23 22:23:19
 * @LastEditors: Hsuan-Wei Fan
 * @LastEditTime: 2020-04-23 22:23:19
 * @Description: 
 */

#include "vaehandler.h"
#include "geometry.h"
#include "transform.h"
#include "sampling.h"
#include "shapes/triangle.h"
#include "network_utils.h"
#include "parallel.h"
#include <random>
#include <fstream>
#include <ctime>


namespace pbrt {

VaeHandler::VaeHandler(const Spectrum &sigmaT, const Spectrum &albedo, Float g, Float eta, const std::string &modelName,
                       const std::string &absModelName, const std::string &angularModelName,
                       const std::string &outputDir, int batchSize) {
    mBatchSize = batchSize;
    MediumParameters mediumParas(albedo, g, eta, sigmaT);
    mAvgMedium = mediumParas;
    mModelName = modelName;
    std::string modelPath = outputDir + "models/" + modelName + "/";
    std::string absModelPath = outputDir + "models_abs/" + absModelName + "/";
    std::string angularModelPath = outputDir + "models_angular/" + angularModelName + "/";
    std::string configFile = modelPath + "training-metadata.json";
    std::string angularConfigFile = angularModelPath + "training-metadata.json";

    mConfig = VaeConfig(configFile, angularModelName != "None" ? angularConfigFile : "", outputDir);
    mPolyOrder = mConfig.polyOrder;
}

int VaeHandler::Prepare(const std::vector<std::shared_ptr<Shape>> &shapes, const PolyUtils::PolyFitConfig &pfConfig) {
    LOG(INFO) << "Preparing the vaehandler";
    if (mModelName == "None") {
        // Dont load any ML models
        mPolyOrder = pfConfig.order;
        PrecomputePolynomials(shapes, mAvgMedium, pfConfig);
        return 0;
    }

    PrecomputePolynomials(shapes, mAvgMedium, pfConfig);
    return 0;
}

void VaeHandler::PrecomputePolynomials(const std::vector<std::shared_ptr<Shape>> &shapes,
                                       const pbrt::MediumParameters &mediumPara,
                                       const pbrt::PolyUtils::PolyFitConfig &pfConfig) {
    LOG(INFO) << "Precompute the polynomials";
    // mKDTrees.push_back(std::vector<ConstraintKDTree>());
    if (mediumPara.isRgb()) {
        // mKDTrees.push_back(std::vector<ConstraintKDTree>());
        // mKDTrees.push_back(std::vector<ConstraintKDTree>());
        PrecomputePolynomialsImpl(shapes, 0, mediumPara, pfConfig);
        PrecomputePolynomialsImpl(shapes, 1, mediumPara, pfConfig);
        PrecomputePolynomialsImpl(shapes, 2, mediumPara, pfConfig);
    } else {
        PrecomputePolynomialsImpl(shapes, 0, mediumPara, pfConfig);
    }

}

void VaeHandler::PrecomputePolynomialsImpl(const std::vector<std::shared_ptr<Shape>> &shapes, int channel,
                                           const pbrt::MediumParameters &mediumPara,
                                           const pbrt::PolyUtils::PolyFitConfig &pfConfig) {
    // 1. Sampling max{1024, 2\sigma_n^2 * SurfaceArea} points around the neighborhood.
    Float kernelEps = PolyUtils::GetKernelEps(mediumPara, channel, pfConfig.kernelEpsScale);
    std::shared_ptr<TriangleMesh> triMesh(PreprocessTriangles(shapes));
    int nSamples = std::max(int(triMesh->Area() * 2.0f / kernelEps), 1024);
    nSamples = std::min(nSamples, 1000000); // FIXME: Adjust the kernel eps instead.
    std::vector<Point3f> sampledP;
    std::vector<Normal3f> sampledN;

    DCHECK_EQ(triMesh->nTriangles, shapes.size());


    for (int i = 0; i < nSamples; i++) {
        Float pdf;
        size_t trigIdx = triMesh->areaDistri->SampleDiscrete(GetRandomFloat());
        Interaction isect = shapes[trigIdx]->Sample(RandPoint2f(), &pdf);
        sampledP.push_back(isect.p);
        sampledN.push_back(isect.n);
    }

#ifdef VISUALIZE_SHAPE_DATA
    {
        const std::string fileName = "../data/pointcloud.txt";
        std::ofstream file;
        file.open(fileName);
        DCHECK(file.is_open());
        for(Point3f &p: sampledP) {
            file << p.x << " " << p.y << " " << p.z << " " << 0 << " " << 0 << " " << 1.0 << std::endl;
        }
        file.close();
    }
#endif

    // 2. Build a constraint KD-tree to accelerate the precomputation. Only computed once before rendering.
    ConstraintKDTree kdTree(sampledP.size());
    kdTree.Build(sampledP, sampledN);

    if (!triMesh->HasPolyCoeffs())
        triMesh->CreatePolyCoeffs();
    PolyStorage *polyCoeffs = triMesh->GetPolyeffs();

    // 3. Fit the polynomials in surrounding by solving the 20 * 20 linear systems.
    bool hasN = triMesh->n != NULL;
    ParallelFor([&](int64_t i) {
        PolyUtils::PolyFitRecord pfRec;
        pfRec.p = triMesh->p[i];
        pfRec.d = hasN ? Vector3f(triMesh->n[i]) : Vector3f(1, 0, 0);
        pfRec.n = hasN ? triMesh->n[i] : Normal3f(1, 0, 0);
        pfRec.kernelEps = kernelEps;
        pfRec.config = pfConfig;
        pfRec.config.useLightspace = false;
        PolyUtils::Polynomial res;
        std::vector<Point3f> pts;
        std::vector<Normal3f> dirs;
        std::tie(res, pts, dirs) = PolyUtils::FitPolynomial(pfRec, &kdTree);
        for (int k = 0; k < res.coeffs.size(); k++) {
            polyCoeffs[i].coeffs[channel][k] = res.coeffs[k];
            polyCoeffs[i].kernelEps[channel] = kernelEps;
            polyCoeffs[i].nPolyCoeffs = res.coeffs.size();
        }
    }, triMesh->nVertices, 1024);
}


std::shared_ptr<TriangleMesh> VaeHandler::PreprocessTriangles(const std::vector<std::shared_ptr<Shape>> &shapes) {
    // "shapes" is the collection of triangles.
    // Initialize triangle mesh's area distribution here.
    Float areaSum = 0.0f;
    size_t shapesNum = shapes.size();
    DCHECK_GT(shapesNum, 0);
    Float areas[shapesNum];
    for (size_t i = 0; i < shapesNum; i++) {
        areas[i] = shapes[i]->Area();
        areaSum += areas[i];
    }
    auto tmp = shapes.back();
    std::shared_ptr<Triangle> triangleIter(std::dynamic_pointer_cast<Triangle>(shapes.back()));
    std::shared_ptr<TriangleMesh> triMesh(triangleIter->mesh);
    triMesh->areaDistri = new Distribution1D(areas, shapesNum);
    triMesh->area = areaSum;
    triMesh->invArea = 1.0f / areaSum;
    return triMesh;
}

void VaeHandler::OnbDuff(const Normal3f &n, Vector3f &b1, Vector3f &b2) {
    float sign = copysignf(1.0f, n.z);
    const Float a = -1.0f / (sign + n.z);
    const Float b = n.x * n.y * a;
    b1 = Vector3f(1.0f + sign * n.x * n.x * a, sign * b, -sign * n.x);
    b2 = Vector3f(b, sign + n.y * n.y * a, -n.y);
}

Transform VaeHandler::AzimuthSpaceTransform(const Vector3f &lightDir, const Normal3f &normal) {
    // TODO: Verify
    Normal3f t1(pbrt::Cross(normal, lightDir));
    Normal3f t2(pbrt::Cross(lightDir, t1));
    Normal3f light(lightDir);
    if (std::abs(pbrt::Dot(normal, lightDir)) > 0.99999f) {
        Vector3f s, t;
        OnbDuff(light, s, t);
        Matrix4x4 lsMatrix(s.x, s.y, s.z, 0,
                           t.x, t.y, t.z, 0,
                           light.x, light.y, light.z, 0,
                           0, 0, 0, 1.0f);
        return Transform(lsMatrix);
    } else {
        Matrix4x4 lsMatrix(t1.x, t1.y, t1.z, 0,
                           t2.x, t2.y, t2.z, 0,
                           light.x, light.y, light.z, 0,
                           0, 0, 0, 1.0f);
        return Transform(lsMatrix);
    }
}


}
