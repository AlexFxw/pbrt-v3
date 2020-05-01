//
// Created by 范軒瑋 on 2020-04-30.
//

#ifndef PBRT_V3_VAEHANDLEREIGEN_H
#define PBRT_V3_VAEHANDLEREIGEN_H

#include "pbrt.h"
#include "vaehandler.h"
#include "network_utils.h"

namespace pbrt {

class VaeHandlerEigen : VaeHandler {
public:
    virtual int Prepare(const pbrt::Scene *scene, const std::vector<Shape *> &shapes, const pbrt::Spectrum &sigmaT,
                    const pbrt::Spectrum &albedo, float g, float eta, const std::string &modelName,
                    const std::string &absModelName, const std::string &angularModelName,
                    const std::string &outputDir, int batchSize,
                    const pbrt::PolyUtils::PolyFitConfig &pfConfig) override;
    virtual ScatterSamplingRecord Sample(const Point3f &pi, const Vector3f &wi, Point3f *po, Vector3f *wo,
                       const Scene *scene, const Vector3f &polyNormal, const Spectrum &sigmaT,
                       const Spectrum &albedo, float g, float eta, Sampler *sampler, const Interaction *isect,
                       bool projectSamples, int channel) const override;
private:
    AbsorptionModel<3> absModel;
    std::unique_ptr<ScatterModelBase> scatterModel;
    Spectrum mEffectiveAlbedo;
    Float mKernelEpsScale; // TODO: need to be initialize
};

}

#endif //PBRT_V3_VAEHANDLEREIGEN_H
