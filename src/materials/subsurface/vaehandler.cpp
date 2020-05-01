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
#include "shapes/triangle.h"
#include "network_utils.h"

namespace pbrt {

int VaeHandler::Prepare(const pbrt::Scene *scene, const std::vector<Shape *> &shapes, const pbrt::Spectrum &sigmaT,
                        const pbrt::Spectrum &albedo, float g, float eta, const std::string &modelName,
                        const std::string &absModelName, const std::string &angularModelName,
                        const std::string &outputDir, int batchSize, const pbrt::PolyUtils::PolyFitConfig &pfConfig) {
    // TODO:
    mBatchSize = batchSize;
    MediumParameters mediumParas(albedo, g, eta, sigmaT);
    mAvgMedium = mediumParas;
    if (modelName == "None") {
        mPolyOrder = pfConfig.order;
        precomputePolynomials(shapes, medium, pfConfig);
        return true; // Dont load any ML models
    }
    std::string modelPath = outputDir + "models/" + modelName + "/";
    std::string absModelPath = outputDir + "models_abs/" + absModelName + "/";
    std::string angularModelPath = outputDir + "models_angular/" + angularModelName + "/";
    // std::string configFile = modelPath + "training-metadata.json";
    // std::string angularConfigFile = angularModelPath + "training-metadata.json";

    mConfig = VaeConfig(outputDir);
    mPolyOrder = mConfig.polyOrder;
    PrecomputePolynomials(shapes, medium, pfConfig);
    return true;
}

void VaeHandler::PrecomputePolynomials(const std::vector<Shape *> &shapes, const pbrt::MediumParameters &mediumPara,
                                       const pbrt::PolyUtils::PolyFitConfig &pfConfig) {
    // TODO: push m_trees
    if (mediumPara.isRgb()) {
        PrecomputePolynomialsImpl(shapes, mediumPara, 0, pfConfig);
        PrecomputePolynomialsImpl(shapes, mediumPara, 1, pfConfig);
        PrecomputePolynomialsImpl(shapes, mediumPara, 2, pfConfig);
    } else {
        PrecomputePolynomialsImpl(shapes, mediumPara, 0, pfConfig);
    }

}

void VaeHandler::PrecomputePolynomialsImpl(const std::vector<Shape *> &shapes, int channel,
                                           const pbrt::MediumParameters &mediumPara,
                                           const pbrt::PolyUtils::PolyFitConfig &pfConfig) {
    // 1. Sampling max{1024, 2\sigma_n^2 * SurfaceArea} points around the neighborhood.
    Float kernelEps = PolyUtils::GetKernelEps(mediumPara, channel, pfConfig.kernelEpsScale);
    for (int shapeIdx = 0; shapeIdx < shapes.size(); shapeIdx++) {
        int nSamples = std::max(int(shapes[shapeIdx]->Area() * 2.0f / kernelEps), 1024);
        std::vector<Point3f> sampledP;
        srd::vector <Vector3f> sampledN;
        for (int i = 0; i < nSamples; i++) {
            Float pdf;
            Interaction isect = shapes[shapeIdx]->Sample(mSampler->Get2D(), &pdf);
            sampledP.push_back(isect.p);
            sampledN.push_back(isect.n);
        }
        // 2. Build a constraint KD-tree to accelerate the precomputation. Only computed once before rendering.
        mKDTrees[channel].push_back(ConstraintKDTree());
        mKDTrees[channel].back().build(sampledP, sampledN);
        TriangleMesh *triMesh = dynamic_cast<TriangleMesh *>(shapes[shapeIdx]);
        if (!triMesh->HasPolyCoeffs())
            triMesh->CreatePolyCoeffs();

        PolyStorage *polyCoeffs = triMesh->GetPolyeffs();

        // 3. Fit the polynomials in surrounding by solving the 20 * 20 linear systems.
        for (int i = 0; i < triMesh->nVertices; i++) {
            PolyUtils::PolyFitRecord pfRec;
            pfRec.p = triMesh->p[i];
            pfRee.d = triMesh->n[i];
            pfRee.n = triMesh->n[i];
            pfRec.kernelEps = kernelEps;
            pfRec.config = pfConfig;
            pfRec.config.useLightspace = false;
            PolyUtils::Polynomial res;
            // TODO: Fit polynomial
            std::vector<Point3f> pts;
            std::vector<Vector3f> dirs;
            std::tie(res, pts, dirs) = PolyUtils::FitPolynomial(pfRec, &(mKDTrees[channel].back()));
            for (int k = 0; k < res.coeffs.size(); k++) {
                polyCoeffs[i].coeffs[channel][k] = res.coeffs[k];
                polyCoeffs[i].kernelEps[channel] = kernelEps;
                polyCoeffs[i].nPolyCoeffs = res.coeffs.size();
            }
        }
    }
}

template <size_t PolyOrder = 3>
std::pair<Eigen::Matrix<float, PolyUtils::nPolyCoeffs(PolyOrder), 1>, Transform>
VaeHandler::GetPolyCoeffsAs(const Point3f &p, const Vector3f &d,
                            const Vector3f &polyNormal,
                            const Interaction *its, int channel) {
    const float *coeffs = its->GetPolyCoeffs();
    Transform transf = AzimuthSpaceTransform(-d, polyNormal);
    Matrix4x4 &m = transf.GetMatrix();
    Vector3f s(m[0][0], m[0][1], m[0][2]);
    Vector3f t(m[1][0], m[1][1], m[1][2]);
    Normal3f n(m[2][0], m[2][1], m[2][2]);
    Eigen::Matrix<float, PolyUtils::nPolyCoeffs(PolyOrder), 1> shapeCoeffs =
            PolyUtils::RotatePolynomialEigen<PolyOrder>(coeffs, s, t, n);
    return std::make_pair(shapeCoeffs, transf);
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
    if (std::abs(pbrt::Dot(normal, lightDir)) > 0.99999f) {
        Vector3f s, t;
        OnbDuff(lightDir, s, t);
        Matrix4x4 lsMatrix(s.x, s.y, s.z, 0, t.x, t.y, t.z, 0, lightDir.x, lightDir.y, lightDir.z, 0, 0, 0, 0, 1.0f);
    }
    else {
        Matrix4x4 lsMatrix(t1.x, t1.y, t1.z, 0, t2.x, t2.y, t2.z, 0, lightDir.x, lightDir.y, lightDir.z, 0, 0, 0,0, 1.0f);
    }
    return Transform(lsTransform);
}



}
