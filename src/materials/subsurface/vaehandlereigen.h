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
    VaeHandlerEigen(const Spectrum &sigmaT, const Spectrum &albedo, Float g, Float eta,
                    const std::string &modelName, const std::string &absModelName, const std::string &angularModelName,
                    const std::string &outputDir, int batchSize, Float kernelEpsScale);

    int Prepare(const std::vector<std::shared_ptr<Shape>> &shapes,
                const PolyUtils::PolyFitConfig &pfConfig) override;

    ScatterSamplingRecord Sample(const Point3f &po, const Vector3f &wo,
                                 const Scene *scene, const Normal3f &polyNormal, const Spectrum &sigmaT,
                                 const Spectrum &albedo, Float g, Float eta, const SurfaceInteraction &isect,
                                 bool projectSamples, int channel, SurfaceInteraction *res) const override;

private:
    std::unique_ptr<ScatterModelBase<3>> scatterModel;
    // std::unique_ptr<FeatureModel<3>> featureModel;
    // std::unique_ptr<AbsorptionModel<3>> absorptionModel;
    Spectrum mEffectiveAlbedo;
    Float mKernelEpsScale;
    std::string mAbsVariableDir, mScatterVariableDir;
};

}

#endif //PBRT_V3_VAEHANDLEREIGEN_H
