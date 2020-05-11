//
// Created by 范軒瑋 on 2020-04-30.
//

#include "vaehandlereigen.h"
#include "vaehandler.h"
#include "interaction.h"
#include "network_utils.h"
#include "bssrdf.h"

namespace pbrt {

int VaeHandlerEigen::Prepare(const std::vector<std::shared_ptr<Shape>> &shapes,
                             const pbrt::PolyUtils::PolyFitConfig &pfConfig) {
    VaeHandler::Prepare(shapes, pfConfig);
    const std::string shapeFeatureName = "mlsPolyLS3";

    featureModel = std::unique_ptr<FeatureModel<3>>(new FeatureModel<3>(
            mScatterVariableDir, mConfig.stats, shapeFeatureName
    ));

    scatterModel = std::unique_ptr<ScatterModelBase<3>>(new ScatterModel<3>(
            mScatterVariableDir, mConfig.stats, shapeFeatureName
    ));

    absorptionModel = std::unique_ptr<AbsorptionModel<3>>(new AbsorptionModel<3>(
            mAbsVariableDir, mConfig.stats, shapeFeatureName
    ));

    LOG(INFO) << "Finish preparing VAEHandlerEigen.";
    return 0;
}

ScatterSamplingRecord VaeHandlerEigen::Sample(const Point3f &po, const Vector3f &wo,
                                              const Scene *scene, const Normal3f &polyNormal, const Spectrum &sigmaT,
                                              const Spectrum &albedo, float g, float eta, Sampler &sampler,
                                              const Interaction &isect, bool projectSamples, int channel) const {
    FeatureModel<3>::ShapeVector shapeCoeffEigen;
    FeatureModel<3>::ShapeVector shapeCoeffEigenWs;
    Transform asTransform;

    // TODO: prediction scope
    const Float *coeffs = isect.GetPolyCoeffs(channel);
    for (int i = 0; i < shapeCoeffEigenWs.size(); i++) {
        shapeCoeffEigenWs[i] = coeffs[i];
    }
    // Calculate the polynomial coefficients
    std::tie(shapeCoeffEigen, asTransform) = VaeHandler::GetPolyCoeffsAs<3>(po, wo, polyNormal, isect, channel);
    const Eigen::Vector3f inPos(po.x, po.y, po.z);
    const Eigen::Vector3f inDir(wo.x, wo.y, wo.z);

    Spectrum albedoChannel(albedo[channel]);
    Spectrum sigmaTChannel(sigmaT[channel]);
    MediumParameters mediumParas(albedoChannel, g, eta, sigmaTChannel);

    Float kernelEps = PolyUtils::GetKernelEps(mAvgMedium, channel, this->mKernelEpsScale);
    Float fitScaleFactor = PolyUtils::GetFitScaleFactor(kernelEps);

    FeatureModel<3>::PreprocessFeatureShape preprocessFeatures =
            featureModel->GetPreprocessFeature(mediumParas.albedo, mediumParas.g, mediumParas.eta, mediumParas.sigmaT,
                                               fitScaleFactor, shapeCoeffEigen);
    FeatureModel<3>::FeatureShape features = featureModel->Run(preprocessFeatures);

    Float absorption = absorptionModel->Run(features);
    Eigen::Vector3f outPos = scatterModel->Run(features);

    outPos = NetworkUtils::LocalToWorld(inPos, -inDir, outPos, true);
    outPos = inPos + (outPos - inPos) / sigmaT.Average();

    ScatterSamplingRecord sRec;
    sRec.throughout = Spectrum(1.0f - absorption);
    sRec.w = Vector3f(1.0, 0, 0); // FIXME
    Point3f sampledP(outPos[0], outPos[1], outPos[2]);
    sRec.isValid = absorption < 1.0f;
    // Project the sampled points to the surface.
    PolyUtils::ProjectPointsToSurface(scene, po, -wo, sRec, shapeCoeffEigen,
                                      mConfig.polyOrder, false, fitScaleFactor, kernelEps);
    return sRec;
}

VaeHandlerEigen::VaeHandlerEigen(const Spectrum &sigmaT, const Spectrum &albedo, float g, float eta,
                                 const std::string &modelName, const std::string &absModelName,
                                 const std::string &angularModelName, const std::string &outputDir, int batchSize,
                                 Float kernelEpsScale) :
        VaeHandler(sigmaT, albedo, g, eta, modelName, absModelName,
                   angularModelName, outputDir, batchSize),
        mKernelEpsScale(kernelEpsScale) {
    std::string modelPath = outputDir + "models/" + modelName + "/";
    std::string absModelPath = outputDir + "models/" + absModelName + "/";
    std::string angularModelPath = outputDir + "models_angular/" + angularModelName + "/";
    std::string graph_path = modelPath + "frozen_model.pb";
    std::string abs_graph_path = absModelPath + "frozen_model.pb";
    std::string angular_graph_path = angularModelPath + "frozen_model.pb";
    std::string configFile = modelPath + "training-metadata.json";
    std::string angularConfigFile = angularModelPath + "training-metadata.json";

    mAbsVariableDir = absModelPath + "variables";
    mScatterVariableDir = modelPath + "variables";
    mEffectiveAlbedo = pbrt::EffectiveAlbedo(albedo);
}


}