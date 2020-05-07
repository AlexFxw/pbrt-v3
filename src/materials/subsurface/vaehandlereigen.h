//
// Created by 范軒瑋 on 2020-04-30.
//

#ifndef PBRT_V3_VAEHANDLEREIGEN_H
#define PBRT_V3_VAEHANDLEREIGEN_H

#include "pbrt.h"
#include "vaehandler.h"
#include "network_utils.h"
#include "bssrdf.h"

namespace pbrt {

class VaeHandlerEigen : public VaeHandler {
public:
    VaeHandlerEigen(Float kernelEpsScale);
    int Prepare(const pbrt::Scene *scene, const std::vector<Shape *> &shapes, const pbrt::Spectrum &sigmaT,
                const pbrt::Spectrum &albedo, float g, float eta, const std::string &modelName,
                const std::string &absModelName, const std::string &angularModelName,
                const std::string &outputDir, int batchSize,
                const pbrt::PolyUtils::PolyFitConfig &pfConfig);

    ScatterSamplingRecord Sample(const Point3f &po, const Vector3f &wo,
                                 const Scene *scene, const Normal3f &polyNormal, const Spectrum &sigmaT,
                                 const Spectrum &albedo, float g, float eta, Sampler &sampler,
                                 const Interaction &isect, bool projectSamples, int channel) const override;

private:
    AbsorptionModel<3> absModel; // TODO: initialize
    std::unique_ptr<ScatterModelBase> scatterModel;
    Spectrum mEffectiveAlbedo;
    Float mKernelEpsScale;
};

}

#endif //PBRT_V3_VAEHANDLEREIGEN_H
