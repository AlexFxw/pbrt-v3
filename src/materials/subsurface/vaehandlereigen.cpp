//
// Created by 范軒瑋 on 2020-04-30.
//

#include <core/interaction.h>
#include <core/bssrdf.h>
#include "vaehandler.h"
#include "vaehandlereigen.h"
#include "network_utils.h"

namespace pbrt {

int VaeHandlerEigen::Prepare(const pbrt::Scene *scene, const std::vector<Shape *> &shapes, const pbrt::Spectrum &sigmaT,
                             const pbrt::Spectrum &albedo, float g, float eta, const std::string &modelName,
                             const std::string &absModelName, const std::string &angularModelName,
                             const std::string &outputDir, int batchSize,
                             const pbrt::PolyUtils::PolyFitConfig &pfConfig) {
    std::string modelPath = outputDir + "models/" + modelName + "/";
    std::string absModelPath = outputDir + "models_abs/" + absModelName + "/";
    std::string angularModelPath = outputDir + "models_angular/" + angularModelName + "/";
    std::string graph_path = modelPath + "frozen_model.pb";
    std::string abs_graph_path = absModelPath + "frozen_model.pb";
    std::string angular_graph_path = angularModelPath + "frozen_model.pb";
    std::string configFile = modelPath + "training-metadata.json";
    std::string angularConfigFile = angularModelPath + "training-metadata.json";

    std::string absVariableDir = absModelPath + "/variables/";
    std::string scatterVariableDir = modelPath + "/variables/";

    scatterModel = std::unique_ptr<ScatterModelBase>(new ScatterModel<3, 8>(
            scatterVariableDir, absVariableDir, mConfig.stats, "mlsPolyLS3"
    ));

    mEffectiveAlbedo = pbrt::EffectiveAlbedo(albedo);
    LOG(INFO) << "Finish loading VAE model.";
    return 0;
}

ScatterSamplingRecord VaeHandlerEigen::Sample(const Point3f &po, const Vector3f &wo,
                                              const Scene *scene, const Normal3f &polyNormal, const Spectrum &sigmaT,
                                              const Spectrum &albedo, float g, float eta, Sampler &sampler,
                                              const Interaction &isect, bool projectSamples, int channel) const {
    AbsorptionModel<3>::ShapeVector shapeCoeffEigen;
    AbsorptionModel<3>::ShapeVector shapeCoeffEigenWs;
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

    Float absorption;
    Eigen::Vector3f outPos;
    std::tie(outPos, absorption) = scatterModel->Run(inPos, inDir, mediumParas.albedo, mediumParas.g, mediumParas.eta,
                                                     mediumParas.sigmaT, fitScaleFactor, shapeCoeffEigen, sampler,
                                                     asTransform);
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

VaeHandlerEigen::VaeHandlerEigen(Float kernelEpsScale) {
    mKernelEpsScale = kernelEpsScale;
}


}