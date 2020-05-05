/*
 * @Author: Hsuan-Wei Fan
 * @Date: 2020-04-23 22:23:13
 * @LastEditors: Hsuan-Wei Fan
 * @LastEditTime: 2020-04-28 09:23:06
 * @Description:
 */

#ifndef PBRT_MATERIALS_SUBSURFACE_VAEHANDLER_H
#define PBRT_MATERIALS_SUBSURFACE_VAEHANDLER_H

#include "pbrt.h"
#include "polynomials.h"
#include "vaeconfig.h"
#include "sampler.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <shapes/triangle.h>

namespace pbrt {

class VaeHandler {
public:
    VaeHandler() {}

    ~VaeHandler() {}

    static void SampleGaussianVector(float *data, Sampler &sampler, int nVars) {
        bool odd = nVars % 2;
        int idx = 0;
        for (int i = 0; i < nVars / 2; ++i) {
            Point2f uv = pbrt::SquareToStdNormal(sampler.Get2D());
            data[idx] = uv.x;
            ++idx;
            data[idx] = uv.y;
            ++idx;
        }
        if (odd)
            data[idx] = pbrt::SquareToStdNormal(sampler.Get2D()).x;
    }

    static void SampleUniformVector(float *data, Sampler *sampler, int nVars) {
        for (int i = 0; i < nVars; ++i) {
            data[i] = sampler->Get1D();
        }
    }

    virtual ScatterSamplingRecord Sample(const Point3f &po, const Vector3f &wo,
                       const Scene *scene, const Normal3f &polyNormal, const Spectrum &sigmaT,
                       const Spectrum &albedo, float g, float eta, Sampler &sampler, const Interaction &isect,
                       bool projectSamples, int channel) const = 0;

    virtual int
    Prepare(const Scene *scene, const std::vector<std::shared_ptr<Shape>> &shapes, const Spectrum &sigmaT,
            const Spectrum &albedo, float g, float eta, const std::string &modelName,
            const std::string &absModelName, const std::string &angularModelName,
            const std::string &outputDir, int batchSize, const PolyUtils::PolyFitConfig &pfConfig);

    void PrecomputePolynomials(const std::vector<std::shared_ptr<Shape>> &shapes, const MediumParameters &mediumPara,
                               const PolyUtils::PolyFitConfig &pfConfig);
    void PrecomputePolynomialsImpl(const std::vector<std::shared_ptr<Shape>> &shapes, int channel, const MediumParameters &mediumPara,
                                   const PolyUtils::PolyFitConfig &pfConfig);

    template<size_t PolyOrder = 3>
    static std::pair<Eigen::Matrix<float, PolyUtils::nPolyCoeffs(PolyOrder), 1>, Transform>
    GetPolyCoeffsAs(const Point3f &p, const Vector3f &d,
                    const Normal3f &polyNormal,
                    const Interaction &its, int channel = 0);

protected:
    Sampler *mSampler; // TODO: Initialize
    VaeConfig mConfig;
    std::vector<std::vector<ConstraintKDTree>> mKDTrees;
    size_t mBatchSize;
    int mPolyOrder;
    MediumParameters mAvgMedium;

    std::shared_ptr<TriangleMesh> PreprocessTriangles(const std::vector<std::shared_ptr<Shape>> &shapes);

};

Transform AzimuthSpaceTransform(const Vector3f &lightDir, const Normal3f &normal);

}  // namespace pbrt

#endif  // PBRT_MATERIALS_SUBSURFACE_VAEHANDLER_H