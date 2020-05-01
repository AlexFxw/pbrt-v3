/*
 * @Author: Hsuan-Wei Fan
 * @Date: 2020-04-23 23:17:24
 * @LastEditors: Hsuan-Wei Fan
 * @LastEditTime: 2020-04-23 23:17:40
 * @Description: 
 */

#include "network_utils.h"
#include <Eigen/Core>
#include <Eigen/Dense>

namespace pbrt {

void NetworkUtils::onb(const Eigen::Vector3f &n, Eigen::Vector3f &b1, Eigen::Vector3f &b2) {
    float sign = copysignf(1.0f, n[2]);
    const float a = -1.0f / (sign + n[2]);
    const float b = n[0] * n[1] * a;
    b1 = Eigen::Vector3f(1.0f + sign * n[0] * n[0] * a, sign * b, -sign * n[0]);
    b2 = Eigen::Vector3f(b, sign + n[1] * n[1] * a, -n[1]);
}

Eigen::VectorXf NetworkUtils::LoadVectorDynamic(const std::string &filename) {
    std::ifstream f(filename, std::ios::binary);
    if (!f.is_open())
        std::cout << "FILE NOT FOUND: " << filename << std::endl;
    std::cout << "Loading " << filename << std::endl;
    int32_t nDims;
    f.read(reinterpret_cast<char *>(&nDims), sizeof(nDims));
    int32_t size;
    f.read(reinterpret_cast<char *>(&size), sizeof(size));
    Eigen::VectorXf ret(size);
    for (int i = 0; i < size; ++i) {
        f.read(reinterpret_cast<char *>(&ret[i]), sizeof(ret[i]));
    }
    return ret;
}

Eigen::MatrixXf NetworkUtils::LoadMatrixDynamic(const std::string &filename) {

    std::ifstream f(filename, std::ios::binary);
    if (!f.is_open())
        LOG(ERROR) << "FILE NOT FOUND " << filename;
    LOG(INFO) << "Loading " << filename;

    int32_t nDims;
    f.read(reinterpret_cast<char *>(&nDims), sizeof(nDims));

    int32_t rows, cols;
    f.read(reinterpret_cast<char *>(&rows), sizeof(rows));
    f.read(reinterpret_cast<char *>(&cols), sizeof(cols));
    Eigen::MatrixXf ret(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            f.read(reinterpret_cast<char *>(&ret(i, j)), sizeof(ret(i, j)));
        }
    }
    return ret;
}

Eigen::Vector3f NetworkUtils::LocalToWorld(const Eigen::Vector3f &inPos, const Eigen::Vector3f &inNormal,
                                           const Eigen::Vector3f &outPosLocal, bool predictInTangentSpace) {
    if (predictInTangentSpace) {
        Eigen::Vector3f tangent1, tangent2;
        onb(inNormal, tangent1, tangent2);
        return inPos + outPosLocal[0] * tangent1 + outPosLocal[1] * tangent2 + outPosLocal[2] * inNormal;
    } else {
        return inPos + outPosLocal;
    }
}

template<size_t PolyOrder, bool useSimilarityTheory>
Eigen::Matrix<float, NetworkUtils::nInFeatures(PolyOrder), 1>
NetworkUtils::PreprocessFeatures(const Spectrum &albedo, float g, float ior, const Spectrum &sigmaT,
                                 const Eigen::Matrix<float, nPolyCoeffs(PolyOrder), 1> &shapeFeatures, float albedoMean,
                                 float albedoStdInv, float gMean, float gStdInv,
                                 const Eigen::Matrix<float, nPolyCoeffs(PolyOrder), 1> &shapeFeatMean,
                                 const Eigen::Matrix<float, nPolyCoeffs(PolyOrder), 1> &shapeFeatStdInv) {

    float effectiveAlbedo;
    if (useSimilarityTheory) {
        Spectrum sigmaS = albedo * sigmaT;
        Spectrum sigmaA = sigmaT - sigmaS;
        Spectrum albedoP = (1 - g) * sigmaS / ((1 - g) * sigmaS + sigmaA);
        effectiveAlbedo = pbrt::effectiveAlbedo(albedoP).average();
    } else {
        effectiveAlbedo = pbrt::effectiveAlbedo(albedo).average();
    }
    float albedoNorm = (effectiveAlbedo - albedoMean) * albedoStdInv;
    float gNorm = (g - gMean) * gStdInv;
    float iorNorm = 2.0f * (ior - 1.25f);
    Eigen::Matrix<float, nPolyCoeffs(PolyOrder), 1> shapeFeaturesNorm =
            (shapeFeatures - shapeFeatMean).cwiseProduct(shapeFeatStdInv);
    Eigen::Matrix<float, nInFeatures(PolyOrder), 1> features;
    features.segment(0, nPolyCoeffs(PolyOrder)) = shapeFeaturesNorm;
    features[nPolyCoeffs(PolyOrder)] = albedoNorm;
    features[nPolyCoeffs(PolyOrder) + 1] = gNorm;
    features[nPolyCoeffs(PolyOrder) + 2] = iorNorm;
    return features;
}

}