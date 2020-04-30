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

namespace pbrt {

class VaeHandler {
public:
    VaeHandler();

    ~VaeHandler();

    virtual int Sample(const Point3f &pi, const Vector3f &wi, Point3f *po, Vector3f *wo) = 0;

    virtual int
    Prepare(const Scene *scene, const std::vector<Shape *> &shapes, Float g, Float eta, const Spectrum &albedo,
            const PolyUtils::PolyFitConfig &pfConfig);

    void PrecomputePolynomials(const std::vector<Shape *> &shapes, const MediumParameters &mediumPara, const PolyUtils::PolyFitConfig &pfConfig);
    void PrecomputePolynomialsImpl(const std::vector<Shape*> &shapes, int channel, const MediumParameters &mediumPara,const PolyUtils::PolyFitConfig &pfConfig);

protected:
    Sampler *mSampler; // TODO: Initialize
    std::vector<std::vector<ConstraintKDTree>> mKDTrees;

};

}  // namespace pbrt

#endif  // PBRT_MATERIALS_SUBSURFACE_VAEHANDLER_H