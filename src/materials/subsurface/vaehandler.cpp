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
#include "samplers/halton.h"
#include "shapes/triangle.h"
#include "network_utils.h"

namespace pbrt {

int VaeHandler::Prepare(const pbrt::Scene *scene, const std::vector<std::shared_ptr<Shape>> &shapes,
                        const pbrt::Spectrum &sigmaT,
                        const pbrt::Spectrum &albedo, float g, float eta, const std::string &modelName,
                        const std::string &absModelName, const std::string &angularModelName,
                        const std::string &outputDir, int batchSize, const pbrt::PolyUtils::PolyFitConfig &pfConfig) {
    mBatchSize = batchSize;
    MediumParameters mediumParas(albedo, g, eta, sigmaT);
    mAvgMedium = mediumParas;
    if (modelName == "None") {
        // Dont load any ML models
        mPolyOrder = pfConfig.order;
        PrecomputePolynomials(shapes, mediumParas, pfConfig);
        return true;
    }
    std::string modelPath = outputDir + "models/" + modelName + "/";
    std::string absModelPath = outputDir + "models_abs/" + absModelName + "/";
    std::string angularModelPath = outputDir + "models_angular/" + angularModelName + "/";
    std::string configFile = modelPath + "training-metadata.json";
    std::string angularConfigFile = angularModelPath + "training-metadata.json";

    mConfig = VaeConfig(configFile, angularModelName != "None" ? angularConfigFile : "", outputDir);
    mPolyOrder = mConfig.polyOrder;
    PrecomputePolynomials(shapes, mediumParas, pfConfig);
    return true;
}

void VaeHandler::PrecomputePolynomials(const std::vector<std::shared_ptr<Shape>> &shapes,
                                       const pbrt::MediumParameters &mediumPara,
                                       const pbrt::PolyUtils::PolyFitConfig &pfConfig) {
    mKDTrees.push_back(std::vector<ConstraintKDTree>());
    if (mediumPara.isRgb()) {
        mKDTrees.push_back(std::vector<ConstraintKDTree>());
        mKDTrees.push_back(std::vector<ConstraintKDTree>());
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
    // TODO: Modify the implementation. Use primitives or shapes?
    // 1. Sampling max{1024, 2\sigma_n^2 * SurfaceArea} points around the neighborhood.
    Float kernelEps = PolyUtils::GetKernelEps(mediumPara, channel, pfConfig.kernelEpsScale);
    std::shared_ptr<TriangleMesh> triMesh = PreprocessTriangles(shapes);
    int nSamples = std::max(int(triMesh->Area() * 2.0f / kernelEps), 1024);
    std::vector<Point3f> sampledP;
    std::vector<Normal3f> sampledN;
    for (int i = 0; i < nSamples; i++) {
        Float pdf;
        size_t trigIdx = triMesh->areaDistri->SampleDiscrete(mSampler->Get1D());
        Interaction isect = shapes[trigIdx]->Sample(mSampler->Get2D(), &pdf);
        sampledP.push_back(isect.p);
        sampledN.push_back(isect.n);
    }
    // 2. Build a constraint KD-tree to accelerate the precomputation. Only computed once before rendering.
    mKDTrees[channel].push_back(ConstraintKDTree());
    mKDTrees[channel].back().Build(sampledP, sampledN);
    if (!triMesh->HasPolyCoeffs())
        triMesh->CreatePolyCoeffs();
    PolyStorage *polyCoeffs = triMesh->GetPolyeffs();

    // 3. Fit the polynomials in surrounding by solving the 20 * 20 linear systems.
    for (int i = 0; i < triMesh->nVertices; i++) {
        PolyUtils::PolyFitRecord pfRec;
        pfRec.p = triMesh->p[i];
        pfRec.d = (Vector3f) triMesh->n[i];
        pfRec.n = triMesh->n[i];
        pfRec.kernelEps = kernelEps;
        pfRec.config = pfConfig;
        pfRec.config.useLightspace = false;
        PolyUtils::Polynomial res;
        std::vector<Point3f> pts;
        std::vector<Normal3f> dirs;
        std::tie(res, pts, dirs) = PolyUtils::FitPolynomial(pfRec, &(mKDTrees[channel].back()));
        for (int k = 0; k < res.coeffs.size(); k++) {
            polyCoeffs[i].coeffs[channel][k] = res.coeffs[k];
            polyCoeffs[i].kernelEps[channel] = kernelEps;
            polyCoeffs[i].nPolyCoeffs = res.coeffs.size();
        }
    }
}

template<size_t PolyOrder = 3>
std::pair<Eigen::Matrix<float, PolyUtils::nPolyCoeffs(PolyOrder), 1>, Transform>
VaeHandler::GetPolyCoeffsAs(const Point3f &p, const Vector3f &d,
                            const Normal3f &polyNormal,
                            const Interaction &its, int channel) {
    const float *coeffs = its.GetPolyCoeffs(channel);
    Transform transf = AzimuthSpaceTransform(-d, polyNormal);
    Matrix4x4 &m = transf.GetMatrix();
    Vector3f s(m(0, 0), m(0, 1), m(0,2));
    Vector3f t(m(1, 0), m(1, 1), m(1,2));
    Normal3f n(m(2, 0), m(2, 1), m(2,2));
    Eigen::Matrix<float, PolyUtils::nPolyCoeffs(PolyOrder), 1> shapeCoeffs =
            PolyUtils::RotatePolynomialEigen<PolyOrder>(coeffs, s, t, n);
    return std::make_pair(shapeCoeffs, transf);
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
    std::shared_ptr<Triangle> triangleIter = std::dynamic_pointer_cast<Triangle>(shapes.back());
    std::shared_ptr<TriangleMesh> triMesh = triangleIter->mesh;
    triMesh->areaDistri = new Distribution1D(areas, shapesNum);
    triMesh->area = areaSum;
    triMesh->invArea = 1.0f / areaSum;
    return triMesh;
}

void OnbDuff(const Normal3f &n, Vector3f &b1, Vector3f &b2) {
    float sign = copysignf(1.0f, n.z);
    const Float a = -1.0f / (sign + n.z);
    const Float b = n.x * n.y * a;
    b1 = Vector3f(1.0f + sign * n.x * n.x * a, sign * b, -sign * n.x);
    b2 = Vector3f(b, sign + n.y * n.y * a, -n.y);
}

Transform AzimuthSpaceTransform(const Vector3f &lightDir, const Normal3f &normal) {
    // TODO: Verify
    Normal3f t1(pbrt::Cross(normal, lightDir));
    Normal3f t2(pbrt::Cross(lightDir, t1));
    Normal3f light(lightDir);
    if (std::abs(pbrt::Dot(normal, lightDir)) > 0.99999f) {
        Vector3f s, t;
        OnbDuff(light, s, t);
        Matrix4x4 lsMatrix(s.x, s.y, s.z, 0, t.x, t.y, t.z, 0, lightDir.x, lightDir.y, lightDir.z, 0, 0, 0, 0, 1.0f);
        return Transform(lsMatrix);
    } else {
        Matrix4x4 lsMatrix(t1.x, t1.y, t1.z, 0, t2.x, t2.y, t2.z, 0, lightDir.x, lightDir.y, lightDir.z, 0, 0, 0, 0,
                           1.0f);
        return Transform(lsMatrix);
    }
}


}
