//
// Created by 范軒瑋 on 2020-04-30.
//

#include "vaehandlereigen.h"
#include "vaehandler.h"
#include "interaction.h"
#include "network_utils.h"
#include "transform.h"
#include "bssrdf.h"

namespace pbrt {

int VaeHandlerEigen::Prepare(const std::vector<std::shared_ptr<Shape>> &shapes,
                             const pbrt::PolyUtils::PolyFitConfig &pfConfig) {
    VaeHandler::Prepare(shapes, pfConfig);
    const std::string shapeFeatureName = "mlsPolyLS3";
    scatterModel = std::unique_ptr<ScatterModelBase<3>>(new ScatterModelSimShared<3>(
            mScatterVariableDir, mAbsVariableDir, mConfig.stats, shapeFeatureName
    ));

    LOG(INFO) << "Finish preparing VAEHandlerEigen.";
    return 0;
}

ScatterSamplingRecord VaeHandlerEigen::Sample(const Point3f &po, const Vector3f &wo,
                                              const Scene *scene, const Normal3f &polyNormal, const Spectrum &sigmaT,
                                              const Spectrum &albedo, Float g, Float eta,
                                              const SurfaceInteraction &isect,
                                              bool projectSamples, int channel, SurfaceInteraction *res) const {
    typedef Eigen::Matrix<float, NetworkUtils::nPolyCoeffs(3), 1> ShapeVector;
    ShapeVector shapeCoeffEigen;
    ShapeVector shapeCoeffEigenWs;
    Transform asTransform;

    const Float *coeffs = isect.GetPolyCoeffs(channel);
    for (int i = 0; i < shapeCoeffEigenWs.size(); i++) {
        shapeCoeffEigenWs[i] = coeffs[i];
    }

    shapeCoeffEigen = VaeHandler::GetPolyCoeffsEigen<3>(po, wo, polyNormal, &isect,
                                                        mConfig.predictionSpace == "LS", channel);

    const Eigen::Vector3f inPos(po.x, po.y, po.z);
    const Eigen::Vector3f inDir(wo.x, wo.y, wo.z);

    Spectrum albedoChannel(albedo[channel]);
    Spectrum sigmaTChannel(sigmaT[channel]);
    MediumParameters mediumParas(albedoChannel, g, eta, sigmaTChannel);

    Float kernelEps = PolyUtils::GetKernelEps(mAvgMedium, channel, this->mKernelEpsScale);
    Float fitScaleFactor = PolyUtils::GetFitScaleFactor(kernelEps);

    Float absorption;
    Eigen::Vector3f outPos;
    std::tie(outPos, absorption) = scatterModel->Run(inPos, inDir, mediumParas.albedo, g, eta,
                                                     mediumParas.sigmaT, fitScaleFactor, shapeCoeffEigen, asTransform);
    ScatterSamplingRecord sRec;
    sRec.throughout = Spectrum(1.0f - absorption);


    outPos = inPos + (outPos - inPos) / sigmaT.Average();
    Point3f sampledP(outPos[0], outPos[1], outPos[2]);

    sRec.p = sampledP;
    sRec.isValid = absorption < 1.0f;
    // if(!sRec.isValid) {
    //     return sRec;
    // }
    // Project the sampled points to the surface.
    if (mConfig.predictionSpace == "AS") {
        PolyUtils::ProjectPointsToSurface(scene, po, -wo, sRec, shapeCoeffEigen,
                                          mConfig.polyOrder, false, fitScaleFactor, kernelEps, res);
    } else {
        PolyUtils::ProjectPointsToSurface(scene, po, -wo, sRec, shapeCoeffEigen,
                                          mConfig.polyOrder, mConfig.predictionSpace == "LS", fitScaleFactor, kernelEps,
                                          res);
    }
#ifdef VISUALIZE_SCATTER
    {
        // FIXME: To visualize scatter point
        const std::string fileName = "../data/scatterpos.txt";
        std::ofstream file;
        file.open(fileName, std::ios::app);
        DCHECK(file.is_open());
        file << sampledP.x << " " << sampledP.y << " " << sampledP.z << " " << 0 << " " << 1.0 << " " << 0 << std::endl;
        file << po.x << " " << po.y << " " << po.z << " " << 1.0 << " " << 0 << " " << 0 << std::endl;
        auto &p = res->p;
        file << p.x << " " << p.y << " " << p.z << " " << 0 << " " << 0 << " " << 1.0 << std::endl;
    }
#endif
    return sRec;
}

VaeHandlerEigen::VaeHandlerEigen(const Spectrum &sigmaT, const Spectrum &albedo, Float g, Float eta,
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