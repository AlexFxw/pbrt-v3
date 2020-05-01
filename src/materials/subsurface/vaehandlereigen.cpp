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
    // TODO: Load the scatter model.
}

ScatterSamplingRecord VaeHandlerEigen::Sample(const Point3f &pi, const Vector3f &wi, Point3f *po, Vector3f *wo,
                            const Scene *scene, const Vector3f &polyNormal, const Spectrum &sigmaT,
                            const Spectrum &albedo, float g, float eta, Sampler *sampler, const Interaction *isect,
                            bool projectSamples, int channel) const {
    AbsorptionModel<3>::ShapeVector shapeCoeffEigen;
    AbsorptionModel<3>::ShapeVector shapeCoeffEigenWs;
    Transform asTransform;

    // TODO: prediction scope
    const float *coeffs = isect->GetPolyCoeffs()[channel];
    for (int i = 0; i < shapeCoeffEigenWs.size(); i++) {
        shapeCoeffEigenWs[i] = coeffs[i];
    }
    std::tie(shapeCoeffEigen, asTransform) = VaeHandler::GetPolyCoeffsAs<3>(pi, wi, polyNormal, isect, channel);
    const Eigen::Vector3f inPos(pi.x, pi.y, pi.z);
    const Eigen::Vector3f inDir(wi.x, wi.y, wi.z);

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
    sRec.w = Vector3f(); // FIXME
    Point3f sampledP(outPos[0], outPos[1], outPos[2]);
    sRec.isValid = absorption < 1.0f;
    PolyUtils::ProjectPointsToSurface(scene, pi, -wi, sRec, shapeCoeffEigen, mConfig.polyOrder, false, fitScaleFactor, kernelEps);
    return sRec;
}

}