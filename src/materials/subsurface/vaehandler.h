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
#include "shapes/triangle.h"
#include <Eigen/Core>
#include <Eigen/Dense>

namespace pbrt {


class VaeHandler {
public:
    VaeHandler(const Spectrum &sigmaT,
               const Spectrum &albedo, Float g, Float eta, const std::string &modelName,
               const std::string &absModelName, const std::string &angularModelName,
               const std::string &outputDir, int batchSize);

    ~VaeHandler() {}

    virtual ScatterSamplingRecord Sample(const Point3f &po, const Vector3f &wo,
                                         const Scene *scene, const Normal3f &polyNormal, const Spectrum &sigmaT,
                                         const Spectrum &albedo, Float g, Float eta, const SurfaceInteraction &isect,
                                         bool projectSamples, int channel, SurfaceInteraction *res) const = 0;

    virtual int Prepare(const std::vector<std::shared_ptr<Shape>> &shapes, const PolyUtils::PolyFitConfig &pfConfig);

    void PrecomputePolynomials(const std::vector<std::shared_ptr<Shape>> &shapes, const MediumParameters &mediumPara,
                               const PolyUtils::PolyFitConfig &pfConfig);
    void PrecomputePolynomialsImpl(const std::vector<std::shared_ptr<Shape>> &shapes, int channel,
                                   const MediumParameters &mediumPara,
                                   const PolyUtils::PolyFitConfig &pfConfig);

    template<size_t PolyOrder = 3>
    static std::pair<Eigen::Matrix<float, PolyUtils::nPolyCoeffs(PolyOrder), 1>, Transform>
    GetPolyCoeffsAs(const Point3f &p, const Vector3f &d,
                    const Normal3f &polyNormal,
                    const SurfaceInteraction &its, int channel = 0) {
        const Float *coeffs = its.GetPolyCoeffs(channel);
        Transform transf = AzimuthSpaceTransform(-d, polyNormal);
        const Matrix4x4 &m = transf.GetMatrix();
        Vector3f s(m(0, 0), m(0, 1), m(0, 2));
        Vector3f t(m(1, 0), m(1, 1), m(1, 2));
        Normal3f n(m(2, 0), m(2, 1), m(2, 2));
        Eigen::Matrix<float, PolyUtils::nPolyCoeffs(PolyOrder), 1> shapeCoeffs =
                PolyUtils::RotatePolynomialEigen<PolyOrder>(coeffs, s, t, n);
        return std::make_pair(shapeCoeffs, transf);
    }

    template<size_t PolyOrder = 3>
    Eigen::Matrix<float, PolyUtils::nPolyCoeffs(PolyOrder), 1> GetPolyCoeffsEigen(
            const Point3f &p, const Vector3f &d, const Normal3f &polyNormal,
            const SurfaceInteraction *its, bool useLightSpace, int channel = 0) const {
        if (its) {
            // const Eigen::Vector3fXf &c = its->polyCoeffs;
            const float *coeffs = its->GetPolyCoeffs(channel);
            if (useLightSpace) {
                Vector3f s, t;
                Normal3f n = Normal3f(-d);
                PolyUtils::OnbDuff(Normal3f(n), s, t);
                Eigen::Matrix<float, PolyUtils::nPolyCoeffs(PolyOrder), 1> shapeCoeffs =
                        PolyUtils::RotatePolynomialEigen<PolyOrder>(coeffs, s, t, n);
                return shapeCoeffs;
            } else {
                Eigen::Matrix<float, PolyUtils::nPolyCoeffs(PolyOrder), 1> shapeCoeffs;
                for (int i = 0; i < its->polyStorage.nPolyCoeffs; ++i)
                    shapeCoeffs[i] = coeffs[i];
                return shapeCoeffs;
            }
        }
        return Eigen::Matrix<float, PolyUtils::nPolyCoeffs(PolyOrder), 1>::Zero();
    }

    const MediumParameters &GetMedium() const { return mAvgMedium; }

protected:
    VaeConfig mConfig;
    std::string mModelName;
    size_t mBatchSize;
    int mPolyOrder;
    MediumParameters mAvgMedium;

private:
    static void OnbDuff(const Normal3f &n, Vector3f &b1, Vector3f &b2);
    static Transform AzimuthSpaceTransform(const Vector3f &lightDir, const Normal3f &normal);
    static std::shared_ptr<TriangleMesh> PreprocessTriangles(const std::vector<std::shared_ptr<Shape>> &shapes);


};


}  // namespace pbrt

#endif  // PBRT_MATERIALS_SUBSURFACE_VAEHANDLER_H