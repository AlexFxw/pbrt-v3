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
#include "spectrum.h"
#include "geometry.h"
#include "transform.h"
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

    static void SampleGaussianVector(float *data, int nVars) {
        bool odd = nVars % 2;
        int idx = 0;
        for (int i = 0; i < nVars / 2; ++i) {
            Point2f uv = pbrt::SquareToStdNormal(pbrt::RandPoint2f());
            data[idx] = uv.x;
            ++idx;
            data[idx] = uv.y;
            ++idx;
        }
        if (odd)
            data[idx] = pbrt::SquareToStdNormal(pbrt::RandPoint2f()).x;
    }

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
                       const Eigen::Matrix<float, nPolyCoeffs(PolyOrder), 1> &shapeFeatStdInv) {
        float effectiveAlbedo;
        if (useSimilarityTheory) {
            Spectrum sigmaS = albedo * sigmaT;
            Spectrum sigmaA = sigmaT - sigmaS;
            Spectrum albedoP = (1 - g) * sigmaS / ((1 - g) * sigmaS + sigmaA);
            effectiveAlbedo = pbrt::EffectiveAlbedo(albedoP).Average();
        } else {
            effectiveAlbedo = pbrt::EffectiveAlbedo(albedo).Average();
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
};


template<size_t PolyOrder, size_t PreLayerWidth>
class ScatterModelBase {
public:
    virtual std::pair<Eigen::Vector3f, Float> Run(const Eigen::Vector3f &inPos, const Eigen::Vector3f &inDir,
                                                  const Spectrum &albedo, float g, float ior, const Spectrum &sigmaT,
                                                  float fitScaleFactor,
                                                  const Eigen::Matrix<float, NetworkUtils::nPolyCoeffs(
                                                          PolyOrder), 1> &polyCoeffs,
                                                  const Transform &toAsTransform) const = 0;

    virtual ~ScatterModelBase() {}
};


template<size_t PolyOrder, size_t NLatent, size_t LayerWidth, size_t PreLayerWidth>
class ScatterModelSimShared : public ScatterModelBase<PolyOrder> {
public:
    ScatterModelSimShared() {}

    ScatterModelSimShared(const std::string &variablePath, const std::string &absVariablePath,
                          const nlohmann::json &stats, const std::string &shapeFeaturesName,
                          const std::string &predictionSpace = "LS", bool useEpsilonSpace = false) {
        m_useEpsilonSpace = useEpsilonSpace;
        absorption_dense_bias =
                NetworkUtils::LoadVectorDynamic(variablePath + "/absorption_dense_bias.bin");
        absorption_mlp_fcn_0_biases =
                NetworkUtils::LoadVectorDynamic(variablePath + "/absorption_mlp_fcn_0_biases.bin");
        scatter_decoder_fcn_fcn_0_biases =
                NetworkUtils::LoadVectorDynamic(variablePath + "/scatter_decoder_fcn_fcn_0_biases.bin");
        scatter_decoder_fcn_fcn_1_biases =
                NetworkUtils::LoadVectorDynamic(variablePath + "/scatter_decoder_fcn_fcn_1_biases.bin");
        scatter_decoder_fcn_fcn_2_biases =
                NetworkUtils::LoadVectorDynamic(variablePath + "/scatter_decoder_fcn_fcn_2_biases.bin");
        shared_preproc_mlp_2_shapemlp_fcn_0_biases =
                NetworkUtils::LoadVectorDynamic(variablePath + "/shared_preproc_mlp_2_shapemlp_fcn_0_biases.bin");
        shared_preproc_mlp_2_shapemlp_fcn_1_biases =
                NetworkUtils::LoadVectorDynamic(variablePath + "/shared_preproc_mlp_2_shapemlp_fcn_1_biases.bin");
        shared_preproc_mlp_2_shapemlp_fcn_2_biases =
                NetworkUtils::LoadVectorDynamic(variablePath + "/shared_preproc_mlp_2_shapemlp_fcn_2_biases.bin");
        scatter_dense_2_bias = NetworkUtils::LoadVectorDynamic(variablePath + "/scatter_dense_2_bias.bin");
        scatter_decoder_fcn_fcn_0_weights =
                NetworkUtils::LoadMatrixDynamic(variablePath + "/scatter_decoder_fcn_fcn_0_weights.bin");
        scatter_decoder_fcn_fcn_1_weights =
                NetworkUtils::LoadMatrixDynamic(variablePath + "/scatter_decoder_fcn_fcn_1_weights.bin");
        scatter_decoder_fcn_fcn_2_weights =
                NetworkUtils::LoadMatrixDynamic(variablePath + "/scatter_decoder_fcn_fcn_2_weights.bin");
        shared_preproc_mlp_2_shapemlp_fcn_0_weights =
                NetworkUtils::LoadMatrixDynamic(variablePath + "/shared_preproc_mlp_2_shapemlp_fcn_0_weights.bin");
        shared_preproc_mlp_2_shapemlp_fcn_1_weights =
                NetworkUtils::LoadMatrixDynamic(variablePath + "/shared_preproc_mlp_2_shapemlp_fcn_1_weights.bin");
        shared_preproc_mlp_2_shapemlp_fcn_2_weights =
                NetworkUtils::LoadMatrixDynamic(variablePath + "/shared_preproc_mlp_2_shapemlp_fcn_2_weights.bin");
        scatter_dense_2_kernel = NetworkUtils::LoadMatrixDynamic(variablePath + "/scatter_dense_2_kernel.bin");

        absorption_dense_kernel = NetworkUtils::LoadMatrixDynamic(variablePath + "/absorption_dense_kernel.bin");
        absorption_mlp_fcn_0_weights = NetworkUtils::LoadMatrixDynamic(
                variablePath + "/absorption_mlp_fcn_0_weights.bin");

        m_gMean = stats["g_mean"][0];
        m_gStdInv = stats["g_stdinv"][0];
        m_albedoMean = stats["effAlbedo_mean"][0];
        m_albedoStdInv = stats["effAlbedo_stdinv"][0];
        std::string degStr = std::to_string(PolyOrder);
        for (int i = 0; i < NetworkUtils::nPolyCoeffs(PolyOrder); ++i) {
            m_shapeFeatMean[i] = stats[shapeFeaturesName + "_mean"][i];
            m_shapeFeatStdInv[i] = stats[shapeFeaturesName + "_stdinv"][i];
        }
        if (predictionSpace != "AS") {
            for (int i = 0; i < 3; ++i) {
                m_outPosMean[i] = stats["outPosRel" + predictionSpace + "_mean"][i];
                m_outPosStd[i] = 1.0f / float(stats["outPosRel" + predictionSpace + "_stdinv"][i]);
            }
        }

        m_useAsSpace = predictionSpace == "AS";
        std::cout << "predictionSpace: " << predictionSpace << std::endl;
        std::cout << "m_useAsSpace: " << m_useAsSpace << std::endl;

    }

    std::pair<Eigen::Vector3f, Float> Run(const Eigen::Vector3f &inPos, const Eigen::Vector3f &inDir,
                                          const Spectrum &albedo, float g, float ior, const Spectrum &sigmaT,
                                          float fitScaleFactor, const Eigen::Matrix<float, NetworkUtils::nPolyCoeffs(
            PolyOrder), 1> &polyCoeffs, const Transform &toAsTransform) const override {
        Spectrum sigmaTp = GetSigmaTP(albedo, g, sigmaT);
        Eigen::Matrix<float, NetworkUtils::nInFeatures(PolyOrder), 1> x =
                NetworkUtils::PreprocessFeatures<PolyOrder, true>(albedo, g, ior, sigmaT, polyCoeffs, m_albedoMean,
                                                                  m_albedoStdInv,
                                                                  m_gMean, m_gStdInv, m_shapeFeatMean,
                                                                  m_shapeFeatStdInv);

        // Apply the preprocessing network
        Eigen::Matrix<float, PreLayerWidth, 1> features =
                (shared_preproc_mlp_2_shapemlp_fcn_0_weights * x + shared_preproc_mlp_2_shapemlp_fcn_0_biases).cwiseMax(
                        0.0f);
        features = (shared_preproc_mlp_2_shapemlp_fcn_1_weights * features +
                    shared_preproc_mlp_2_shapemlp_fcn_1_biases).cwiseMax(0.0f);
        features = (shared_preproc_mlp_2_shapemlp_fcn_2_weights * features +
                    shared_preproc_mlp_2_shapemlp_fcn_2_biases).cwiseMax(0.0f);

        // Compute absorption
        Eigen::Matrix<float, 32, 1> absTmp = (absorption_mlp_fcn_0_weights * features +
                                              absorption_mlp_fcn_0_biases).cwiseMax(0.0f);
        Eigen::Matrix<float, 1, 1> a = absorption_dense_kernel * absTmp + absorption_dense_bias;
        Float absorption = NetworkUtils::sigmoid(a[0]);

        Float randFloat = GetRandomFloat();
        if (randFloat > absorption) {
            absorption = 0.0f; // nothing gets absorbed instead
        } else {
            return std::make_pair(inPos, 1.0f); // all is absorbed
        }

        // Concatenate features with random numbers
        Eigen::Matrix<float, NLatent, 1> latent(NLatent);
        NetworkUtils::SampleGaussianVector(latent.data(), NLatent);

        Eigen::Matrix<float, PreLayerWidth + NLatent, 1> featLatent;
        featLatent << latent, features;

        Eigen::Matrix<float, 64, 1> y = (scatter_decoder_fcn_fcn_0_weights * featLatent +
                                         scatter_decoder_fcn_fcn_0_biases).cwiseMax(0.0f);
        y = (scatter_decoder_fcn_fcn_1_weights * y + scatter_decoder_fcn_fcn_1_biases).cwiseMax(0.0f);
        y = (scatter_decoder_fcn_fcn_2_weights * y + scatter_decoder_fcn_fcn_2_biases).cwiseMax(0.0f);
        Eigen::Vector3f outPos = scatter_dense_2_kernel * y + scatter_dense_2_bias;
        // if (m_useEpsilonSpace) {
        //     if (m_useAsSpace) {
        //         Vector3f tmp = toAsTransform(Vector3f(outPos[0], outPos[1], outPos[2])) / fitScaleFactor;
        //         outPos = Eigen::Vector3f(tmp.x, tmp.y, tmp.z) + inPos;
        //     } else {
        //         outPos = NetworkUtils::LocalToWorld(inPos, -inDir, outPos, true);
        //         outPos = inPos + (outPos - inPos) / fitScaleFactor;
        //     }
        // } else {
        //     outPos = outPos.cwiseProduct(m_outPosStd) + m_outPosMean;
        //     outPos = NetworkUtils::LocalToWorld(inPos, -inDir, outPos, true);
        //     outPos = inPos + (outPos - inPos) / sigmaTp.Average();
        // }
        outPos = outPos.cwiseProduct(m_outPosStd) + m_outPosMean;
        return std::make_pair(outPos, absorption);
    }

    bool m_useEpsilonSpace, m_useAsSpace;

    Eigen::Matrix<float, 32, PreLayerWidth> absorption_mlp_fcn_0_weights;
    Eigen::Matrix<float, 32, 1> absorption_mlp_fcn_0_biases;
    Eigen::Matrix<float, 1, 32> absorption_dense_kernel;
    Eigen::Matrix<float, 1, 1> absorption_dense_bias;

    Eigen::Matrix<float, LayerWidth, 1> scatter_decoder_fcn_fcn_0_biases;
    Eigen::Matrix<float, LayerWidth, 1> scatter_decoder_fcn_fcn_1_biases;
    Eigen::Matrix<float, LayerWidth, 1> scatter_decoder_fcn_fcn_2_biases;
    Eigen::Matrix<float, PreLayerWidth, 1> shared_preproc_mlp_2_shapemlp_fcn_0_biases;
    Eigen::Matrix<float, PreLayerWidth, 1> shared_preproc_mlp_2_shapemlp_fcn_1_biases;
    Eigen::Matrix<float, PreLayerWidth, 1> shared_preproc_mlp_2_shapemlp_fcn_2_biases;
    Eigen::Matrix<float, 3, 1> scatter_dense_2_bias;

    Eigen::Matrix<float, LayerWidth, PreLayerWidth + NLatent> scatter_decoder_fcn_fcn_0_weights;
    Eigen::Matrix<float, LayerWidth, LayerWidth> scatter_decoder_fcn_fcn_1_weights;
    Eigen::Matrix<float, LayerWidth, LayerWidth> scatter_decoder_fcn_fcn_2_weights;
    Eigen::Matrix<float, 3, LayerWidth> scatter_dense_2_kernel;

    Eigen::Matrix<float, PreLayerWidth, NetworkUtils::nInFeatures(
            PolyOrder)> shared_preproc_mlp_2_shapemlp_fcn_0_weights;
    Eigen::Matrix<float, PreLayerWidth, PreLayerWidth> shared_preproc_mlp_2_shapemlp_fcn_1_weights;
    Eigen::Matrix<float, PreLayerWidth, PreLayerWidth> shared_preproc_mlp_2_shapemlp_fcn_2_weights;

    Eigen::Matrix<float, NetworkUtils::nPolyCoeffs(PolyOrder), 1> m_shapeFeatMean, m_shapeFeatStdInv;
    Eigen::Vector3f m_outPosMean, m_outPosStd;
    float m_albedoMean, m_albedoStdInv, m_gMean, m_gStdInv;
    // EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};


}  // namespace pbrt

#endif // PBRT_MATERIALS_SUBSURFACE_NETWORKUTILS_H