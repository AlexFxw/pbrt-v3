/*
 * @Author: Hsuan-Wei Fan
 * @Date: 2020-04-23 22:23:19
 * @LastEditors: Hsuan-Wei Fan
 * @LastEditTime: 2020-04-23 22:23:19
 * @Description: 
 */

#include "vaehandler.h"
#include "shapes/triangle.h"

namespace pbrt {

int VaeHandler::Sample(const Point3f &pi, const Vector3f &wi, Point3f *po, Vector3f *wo) {
    // TODO
    // Build acceleration data structure

}

void VaeHandler::PrecomputePolynomials(const std::vector<Shape *> &shapes, const pbrt::MediumParameters &mediumPara,
                                       const pbrt::PolyUtils::PolyFitConfig &pfConfig) {
    // TODO: push m_trees
    if(mediumPara.isRgb()) {
        PrecomputePolynomialsImpl(shapes, mediumPara, 0, pfConfig);
        PrecomputePolynomialsImpl(shapes, mediumPara, 1, pfConfig);
        PrecomputePolynomialsImpl(shapes, mediumPara, 2, pfConfig);
    }
    else {
        PrecomputePolynomialsImpl(shapes, mediumPara, 0, pfConfig);
    }

}

void VaeHandler::PrecomputePolynomialsImpl(const std::vector<Shape *> &shapes, int channel,
                                           const pbrt::MediumParameters &mediumPara,
                                           const pbrt::PolyUtils::PolyFitConfig &pfConfig) {
    // TODO
    // 1. Sampling max{1024, 2\sigma_n^2 * SurfaceArea} points around the neighborhood.
    Float kernelEps = PolyUtils::GetKernelEps(mediumPara, channel, pfConfig.kernelEpsScale);
    for(int shapeIdx = 0; shapeIdx < shapes.size(); shapeIdx ++)
    {
        int nSamples = std::max(int(shapes[shapeIdx]->Area() * 2.0f / kernelEps), 1024);
        std::vector<Point3f> sampledP;
        srd::vector<Vector3f> sampledN;
        for(int i = 0; i < nSamples; i ++)
        {
            Float pdf;
            Interaction isect = shapes[shapeIdx]->Sample(mSampler->Get2D(), &pdf);
            sampledP.push_back(isect.p);
            sampledN.push_back(isect.n);
        }
        // 2. Build a constraint KD-tree to accelerate the precomputation. Only computed once before rendering.
        mKDTrees[channel].push_back(ConstraintKDTree());
        mKDTrees[channel].back().build(sampledP, sampledN);
        TriangleMesh *triMesh = dynamic_cast<TriangleMesh*>(shapes[shapeIdx]);
        if(!triMesh->HasPolyCoeffs())
            triMesh->CreatePolyCoeffs();

        PolyStorage *polyCoeffs = triMesh->GetPolyeffs();

        for(int i = 0; i < triMesh->nVertices; i ++)
        {
            PolyUtils::PolyFitRecord pfRec;
            pfRec.p = triMesh->p[i];
            pfRee.d = triMesh->n[i];
            pfRee.n = triMesh->n[i];
            pfRec.kernelEps = kernelEps;
            pfRec.config = pfConfig;
            pfRec.config.useLightspace = false;
            PolyUtils::Polynomial res;
            // TODO: Fit polynomial
        }
    }
    // 3. Fit the polynomials in surrounding by solving the 20 * 20 linear systems.
}

}
