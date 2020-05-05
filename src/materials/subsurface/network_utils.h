/*
 * @Author: Hsuan-Wei Fan
 * @Date: 2020-04-23 22:41:27
 * @LastEditors: Hsuan-Wei Fan
 * @LastEditTime: 2020-04-23 23:17:11
 * @Description: 
 */

#ifndef PBRT_MATERIALS_SUBSURFACE_NETWORKUTILS_H
#define PBRT_MATERIALS_SUBSURFACE_NETWORKUTILS_H

#include "pbrt.h"
#include "json.h"
#include "vaehandler.h"
#include "transform.h"
#include "spectrum.h"
#include <Eigen/Core>
#include <Eigen/Dense>

namespace pbrt {

class NetworkUtils {
public:
    static void onb(const Eigen::Vector3f &n, Eigen::Vector3f &b1, Eigen::Vector3f &b2);

    static inline constexpr int nChooseK(int n, int k) {
        return (k == 0 || n == k) ? 1 : nChooseK(n - 1, k - 1) + nChooseK(n - 1, k);
    }

    static inline constexpr int nPolyCoeffs(int polyOrder) { return nChooseK(3 + polyOrder, polyOrder); }

    static inline constexpr int nInFeatures(int polyOrder) { return nPolyCoeffs(polyOrder) + 3; }

    static inline float sigmoid(float x) { return 1.0f / (1.0f + std::exp(-x)); }

    static Eigen::VectorXf LoadVectorDynamic(const std::string &filename);

    static Eigen::MatrixXf LoadMatrixDynamic(const std::string &filename);

    static Eigen::Vector3f LocalToWorld(const Eigen::Vector3f &inPos, const Eigen::Vector3f &inNormal,
                                        const Eigen::Vector3f &outPosLocal, bool predictInTangentSpace);

    template<size_t PolyOrder = 3, bool useSimilarityTheory = false>
    static Eigen::Matrix<float, nInFeatures(PolyOrder), 1>
    PreprocessFeatures(const Spectrum &albedo, float g, float ior, const Spectrum &sigmaT,
                       const Eigen::Matrix<float, nPolyCoeffs(PolyOrder), 1> &shapeFeatures, float albedoMean,
                       float albedoStdInv, float gMean, float gStdInv,
                       const Eigen::Matrix<float, nPolyCoeffs(PolyOrder), 1> &shapeFeatMean,
                       const Eigen::Matrix<float, nPolyCoeffs(PolyOrder), 1> &shapeFeatStdInv);
};

template<size_t PolyOrder = 3, size_t LayerWidth = 64>
class AbsorptionModel {
public:
    typedef Eigen::Matrix<float, NetworkUtils::nPolyCoeffs(PolyOrder), 1> ShapeVector;

    AbsorptionModel() {}

    // AbsorptionModel(const std::string &variablePath, const VaeConfig &config) {
    AbsorptionModel(const std::string &variablePath, const nlohmann::json &stats,
                    const std::string &shapeFeaturesName) {
        absorption_mlp_fcn_0_biases =
                NetworkUtils::LoadVectorDynamic(variablePath + "/absorption_mlp_fcn_0_biases.bin");
        absorption_mlp_fcn_1_biases =
                NetworkUtils::LoadVectorDynamic(variablePath + "/absorption_mlp_fcn_1_biases.bin");
        absorption_mlp_fcn_2_biases =
                NetworkUtils::LoadVectorDynamic(variablePath + "/absorption_mlp_fcn_2_biases.bin");
        absorption_dense_bias = NetworkUtils::LoadVectorDynamic(variablePath + "/absorption_dense_bias.bin");
        absorption_mlp_fcn_0_weights =
                NetworkUtils::LoadMatrixDynamic(variablePath + "/absorption_mlp_fcn_0_weights.bin");
        absorption_mlp_fcn_1_weights =
                NetworkUtils::LoadMatrixDynamic(variablePath + "/absorption_mlp_fcn_1_weights.bin");
        absorption_mlp_fcn_2_weights =
                NetworkUtils::LoadMatrixDynamic(variablePath + "/absorption_mlp_fcn_2_weights.bin");
        absorption_dense_kernel = NetworkUtils::LoadMatrixDynamic(variablePath + "/absorption_dense_kernel.bin");

        m_gMean = stats["g_mean"][0];
        m_gStdInv = stats["g_stdinv"][0];
        m_albedoMean = stats["effAlbedo_mean"][0];
        m_albedoStdInv = stats["effAlbedo_stdinv"][0];
        std::string degStr = std::to_string(PolyOrder);
        for (int i = 0; i < NetworkUtils::nPolyCoeffs(PolyOrder); ++i) {
            m_shapeFeatMean[i] = stats[shapeFeaturesName + "_mean"][i];
            m_shapeFeatStdInv[i] = stats[shapeFeaturesName + "_stdinv"][i];
        }
    }

    Float Run(Spectrum albedo, float g, float ior, const Spectrum &sigmaT, const ShapeVector &polyCoeffs) const {
        Eigen::Matrix<float, NetworkUtils::nInFeatures(PolyOrder), 1> input =
                NetworkUtils::PreprocessFeatures<PolyOrder, false>(albedo, g, ior, sigmaT, polyCoeffs, m_albedoMean,
                                                                   m_albedoStdInv,
                                                                   m_gMean, m_gStdInv, m_shapeFeatMean,
                                                                   m_shapeFeatStdInv);
        Eigen::Matrix<float, LayerWidth, 1> x =
                (absorption_mlp_fcn_0_weights * input + absorption_mlp_fcn_0_biases).cwiseMax(0.0f);
        x = (absorption_mlp_fcn_1_weights * x + absorption_mlp_fcn_1_biases).cwiseMax(0.0f);
        x = (absorption_mlp_fcn_2_weights * x + absorption_mlp_fcn_2_biases).cwiseMax(0.0f);
        Eigen::Matrix<float, 1, 1> x2 = absorption_dense_kernel * x + absorption_dense_bias;
        return NetworkUtils::sigmoid(x2[0]);
    }

    Eigen::Matrix<float, LayerWidth, 1> absorption_mlp_fcn_0_biases;
    Eigen::Matrix<float, LayerWidth, 1> absorption_mlp_fcn_1_biases;
    Eigen::Matrix<float, LayerWidth, 1> absorption_mlp_fcn_2_biases;
    Eigen::VectorXf absorption_dense_bias;

    Eigen::Matrix<float, LayerWidth, NetworkUtils::nInFeatures(PolyOrder)> absorption_mlp_fcn_0_weights;
    Eigen::Matrix<float, LayerWidth, LayerWidth> absorption_mlp_fcn_1_weights;
    Eigen::Matrix<float, LayerWidth, LayerWidth> absorption_mlp_fcn_2_weights;
    Eigen::Matrix<float, 1, LayerWidth> absorption_dense_kernel;

    ShapeVector m_shapeFeatMean, m_shapeFeatStdInv;
    float m_albedoMean, m_albedoStdInv, m_gMean, m_gStdInv;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class ScatterModelBase {
public:
    virtual std::pair<Eigen::Vector3f, float>
    Run(const Eigen::Vector3f &inPos, const Eigen::Vector3f &inDir, const Spectrum &albedo, float g, float ior,
        const Spectrum &sigmaT, float polyScaleFactor,
        const Eigen::Matrix<float, NetworkUtils::nPolyCoeffs(3), 1> &polyCoeffs, Sampler &sampler,
        const Transform &toAsTransform) const = 0;

    virtual ~ScatterModelBase() {}
};

template<size_t PolyOrder = 3, size_t NLatent = 8, size_t LayerWidth = 64,
        size_t PreLayerWidth = 32>
class ScatterModel : public ScatterModelBase {
public:
    ScatterModel() {}

    ScatterModel(const std::string &variablePath, const std::string &absVariablePath, const nlohmann::json &stats,
                 const std::string &shapeFeaturesName,
                 const std::string &predictionSpace = "LS") {

        absModel = AbsorptionModel<3>(absVariablePath, stats, shapeFeaturesName);

        scatter_decoder_fcn_fcn_0_biases =
                NetworkUtils::LoadVectorDynamic(variablePath + "/scatter_decoder_fcn_fcn_0_biases.bin");
        scatter_decoder_fcn_fcn_1_biases =
                NetworkUtils::LoadVectorDynamic(variablePath + "/scatter_decoder_fcn_fcn_1_biases.bin");
        scatter_decoder_fcn_fcn_2_biases =
                NetworkUtils::LoadVectorDynamic(variablePath + "/scatter_decoder_fcn_fcn_2_biases.bin");
        scatter_shapemlp_fcn_0_biases =
                NetworkUtils::LoadVectorDynamic(variablePath + "/scatter_shapemlp_fcn_0_biases.bin");
        scatter_shapemlp_fcn_1_biases =
                NetworkUtils::LoadVectorDynamic(variablePath + "/scatter_shapemlp_fcn_1_biases.bin");
        scatter_shapemlp_fcn_2_biases =
                NetworkUtils::LoadVectorDynamic(variablePath + "/scatter_shapemlp_fcn_2_biases.bin");
        scatter_dense_2_bias = NetworkUtils::LoadVectorDynamic(variablePath + "/scatter_dense_2_bias.bin");
        scatter_decoder_fcn_fcn_0_weights =
                NetworkUtils::LoadMatrixDynamic(variablePath + "/scatter_decoder_fcn_fcn_0_weights.bin");
        scatter_decoder_fcn_fcn_1_weights =
                NetworkUtils::LoadMatrixDynamic(variablePath + "/scatter_decoder_fcn_fcn_1_weights.bin");
        scatter_decoder_fcn_fcn_2_weights =
                NetworkUtils::LoadMatrixDynamic(variablePath + "/scatter_decoder_fcn_fcn_2_weights.bin");
        scatter_shapemlp_fcn_0_weights =
                NetworkUtils::LoadMatrixDynamic(variablePath + "/scatter_shapemlp_fcn_0_weights.bin");
        scatter_shapemlp_fcn_1_weights =
                NetworkUtils::LoadMatrixDynamic(variablePath + "/scatter_shapemlp_fcn_1_weights.bin");
        scatter_shapemlp_fcn_2_weights =
                NetworkUtils::LoadMatrixDynamic(variablePath + "/scatter_shapemlp_fcn_2_weights.bin");
        scatter_dense_2_kernel = NetworkUtils::LoadMatrixDynamic(variablePath + "/scatter_dense_2_kernel.bin");

        m_gMean = stats["g_mean"][0];
        m_gStdInv = stats["g_stdinv"][0];
        m_albedoMean = stats["effAlbedo_mean"][0];
        m_albedoStdInv = stats["effAlbedo_stdinv"][0];
        std::string degStr = std::to_string(PolyOrder);
        for (int i = 0; i < NetworkUtils::nPolyCoeffs(PolyOrder); ++i) {
            m_shapeFeatMean[i] = stats[shapeFeaturesName + "_mean"][i];
            m_shapeFeatStdInv[i] = stats[shapeFeaturesName + "_stdinv"][i];
        }
        for (int i = 0; i < 3; ++i) {
            m_outPosMean[i] = stats["outPosRel" + predictionSpace + "_mean"][i];
            m_outPosStd[i] = 1.0f / float(stats["outPosRel" + predictionSpace + "_stdinv"][i]);
        }
    }

    std::pair<Eigen::Vector3f, float>
    Run(const Eigen::Vector3f &inPos, const Eigen::Vector3f &inDir, const Spectrum &albedo, float g, float ior,
        const Spectrum &sigmaT,
        const float polyScaleFactor,
        const Eigen::Matrix<float, NetworkUtils::nPolyCoeffs(PolyOrder), 1> &polyCoeffs,
        Sampler &sampler, const Transform &toAsTransform) const override {
        Eigen::Matrix<float, NetworkUtils::nInFeatures(PolyOrder), 1> x =
                NetworkUtils::PreprocessFeatures<PolyOrder, false>(albedo, g, ior, sigmaT, polyCoeffs, m_albedoMean,
                                                                   m_albedoStdInv,
                                                                   m_gMean, m_gStdInv, m_shapeFeatMean,
                                                                   m_shapeFeatStdInv);

        // Apply the preprocessing network
        Eigen::Matrix<float, PreLayerWidth, 1> features =
                (scatter_shapemlp_fcn_0_weights * x + scatter_shapemlp_fcn_0_biases).cwiseMax(0.0f);
        features = (scatter_shapemlp_fcn_1_weights * features + scatter_shapemlp_fcn_1_biases).cwiseMax(0.0f);
        features = (scatter_shapemlp_fcn_2_weights * features + scatter_shapemlp_fcn_2_biases).cwiseMax(0.0f);
        // Concatenate features with random numbers

        Eigen::Matrix<float, NLatent, 1> latent(NLatent);
        VaeHandler::SampleGaussianVector(latent.data(), sampler, NLatent);

        Eigen::Matrix<float, PreLayerWidth + NLatent, 1> featLatent;
        featLatent << latent, features;
        Eigen::Matrix<float, LayerWidth, 1> y =
                (scatter_decoder_fcn_fcn_0_weights * featLatent + scatter_decoder_fcn_fcn_0_biases).cwiseMax(0.0f);

        Eigen::Matrix<float, LayerWidth + PreLayerWidth, 1> stacked;
        stacked << y, features;
        y = (scatter_decoder_fcn_fcn_1_weights * stacked + scatter_decoder_fcn_fcn_1_biases).cwiseMax(0.0f);
        stacked << y, features;
        y = (scatter_decoder_fcn_fcn_2_weights * stacked + scatter_decoder_fcn_fcn_2_biases).cwiseMax(0.0f);
        stacked << y, features;
        Eigen::Vector3f outPos = scatter_dense_2_kernel * stacked + scatter_dense_2_bias;
        outPos = outPos.cwiseProduct(m_outPosStd) + m_outPosMean;
        outPos = NetworkUtils::LocalToWorld(inPos, -inDir, outPos, true);
        outPos = inPos + (outPos - inPos) / sigmaT.Average();
        float absorption = absModel.Run(albedo, g, ior, sigmaT, polyCoeffs);
        return std::make_pair(outPos, absorption);
    }

    // Public member variables.
    Eigen::Matrix<float, LayerWidth, 1> scatter_decoder_fcn_fcn_0_biases;
    Eigen::Matrix<float, LayerWidth, 1> scatter_decoder_fcn_fcn_1_biases;
    Eigen::Matrix<float, LayerWidth, 1> scatter_decoder_fcn_fcn_2_biases;
    Eigen::Matrix<float, PreLayerWidth, 1> scatter_shapemlp_fcn_0_biases;
    Eigen::Matrix<float, PreLayerWidth, 1> scatter_shapemlp_fcn_1_biases;
    Eigen::Matrix<float, PreLayerWidth, 1> scatter_shapemlp_fcn_2_biases;
    Eigen::Matrix<float, 3, 1> scatter_dense_2_bias;

    Eigen::Matrix<float, LayerWidth, PreLayerWidth + NLatent> scatter_decoder_fcn_fcn_0_weights;
    Eigen::Matrix<float, LayerWidth, LayerWidth + PreLayerWidth> scatter_decoder_fcn_fcn_1_weights;
    Eigen::Matrix<float, LayerWidth, LayerWidth + PreLayerWidth> scatter_decoder_fcn_fcn_2_weights;
    Eigen::Matrix<float, 3, LayerWidth + PreLayerWidth> scatter_dense_2_kernel;

    Eigen::Matrix<float, PreLayerWidth, NetworkUtils::nInFeatures(PolyOrder)> scatter_shapemlp_fcn_0_weights;
    Eigen::Matrix<float, PreLayerWidth, PreLayerWidth> scatter_shapemlp_fcn_1_weights;
    Eigen::Matrix<float, PreLayerWidth, PreLayerWidth> scatter_shapemlp_fcn_2_weights;

    Eigen::Matrix<float, NetworkUtils::nPolyCoeffs(PolyOrder), 1> m_shapeFeatMean, m_shapeFeatStdInv;
    Eigen::Vector3f m_outPosMean, m_outPosStd;
    float m_albedoMean, m_albedoStdInv, m_gMean, m_gStdInv;

    AbsorptionModel<3> absModel;

    // EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};


}  // namespace pbrt

#endif // PBRT_MATERIALS_SUBSURFACE_NETWORKUTILS_H