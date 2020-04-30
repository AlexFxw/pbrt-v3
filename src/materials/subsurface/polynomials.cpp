//
// Created by 范軒瑋 on 2020-04-28.
//

#include "polynomials.h"
#include "geometry.h"

namespace pbrt {


void ConstraintKDTree::Build(const std::vector<Point3f> &sampledP, const std::vector<Noraml3f> &sampledN) {
    int pointsNum = sampledP.size();
    mTree = CTree(pointsNum);

    for (size_t i = 0; i < pointsNum; i++) {
        mTree[i].SetPosition(sampledP[i]);
        ExtraData d;
        d.p = sampledP[i];
        d.n = sampledN[i];
        d.sampleCount = 1;
        d.avgN = Normal3f(-10, -10, -10);
        d.avgP = Point3f(-100, -100, -100);
        m_tree[i].SetData(d);
    }

    mTree.Build(true);
    CalcAvgValues(mTree[0], 0);

    for (size_t i = 0; i < std::min(size_t(32), sampledP.size()); i++) {
        mPoints.push_back(sampledP[i]);
        mNormals.push_back(sampledN[i]);
    }
}

std::tuple<Point3f, Normal3f, size_t> ConstraintKDTree::CalcAvgValues(TreeNode &node, TreeNode::IndexType index) {
    Point3f rp, lp;
    Normal3f rn, ln;
    size_t rSamples = 0, lSamples = 0;
    if (mTree.HasRightChild(index)) {
        auto rightIdx = node.GetRightIndex(index);
        std::tie(rp, rn, rSamples) = CalcAvgValues(mTree[rightIdx], rightIdx);
    }

    if (!node.IsLeaf()) {
        auto leftIdx = node.GetLeftIndex(index);
        std::tie(lp, ln, lSamples) = CalcAvgValues(mTree[leftIdx], leftIdx);
    }

    size_t allSamples = lSamples + rSamples + 1;
    Point3f avgP = (node.GetData().p + ((Float) rSamples) * rp + ((Float) lSamples) * lp) / ((Float) allSamples);
    Normal3f avgN = (node.GetData().n + ((Float) rSamples) * rn + ((Float) lSamples) * ln) / ((Float) allSamples);

    node.GetData().avgP = avgP;
    node.GetData().avgN = avgN;
    node.GetData().sampledNum = allSamples;
    return std::make_tuple(avgP, avgN, allSamples);
}

void ConstraintKDTree::GetConstraints(const Point3f &p, SimpleKDNode<Point3f, pbrt::ConstraintKDTree::ExtraData> node,
                                      TreeNode::IndexType index, const Bounds3f &aabb, std::vector<Point3f> &points,
                                      std::vector<Normal3f> &normals, std::vector<Float> &sampleWeights,
                                      pbrt::Float kernelEps,
                                      const std::function<Float(Float, Float)> &kernelFunc) const {
    Float dist2Threshold = 9.0 * kernelEps;
    if (pbrt::DistanceSquared(p, aabb) > dist2Threshold) {
        return;
    }
    if (pbrt::DistanceSquared(node.GetData().p, p) < dist2Threshold) {
        points.push_back(node.GetData().p)
        normals.push_back(node.GetData().n)
        sampleWeights.push_back(1.0f);
    }

    if (node.IsLeaf())
        return;

    // Compute bounding box of child nodes
    uint16_t ax = node.GetAxis();
    Point3f leftMin = aabb.pMin;
    Point3f leftMax = aabb.pMax;
    leftMax[ax] = node.GetPos()[ax];
    Point3f rightMin = aabb.pMin;
    Point3f rightMax = aabb.pMax;
    rightMin[ax] = node.GetPos()[ax];
    Bounds3f leftBounds(leftMin, leftMax);
    Bounds3f rightBounds(rightMin, rightMax);
    if (mTree.hasRightChild(index)) {
        auto rightIdx = node.GetRightIndex(index);
        GetConstraints(p, mTree[rightIdx], rightIdx, rightBounds, points, normals, sampleWeights, kernelEps,
                       kernelFunc);
    }
    // Node always has left child by construction
    auto leftIdx = node.GetLeftIndex(index);
    GetConstraints(p, mTree[leftIdx], leftIdx, leftBounds, points, normals, sampleWeights, kernelEps, kernelFunc);
}

Float PolyUtils::GetKernelEps(const MediumParameters &mediumParas, int channel,
                              pbrt::Float kernalMultiplier) {
    // The kernel epsilon size is related to the medium's properties
    Float sigmaT = mediumParas.sigmaT[channel];
    Float albedo = medium.albedo[channel];
    Float g = medium.g;
    Float sigmaS = albedo * sigmaT;
    Float sigmaa = sigmaT - sigmaS;
    Float sigmaSp = (1 - g) * sigmaS;
    Float sigmaTp = sigmaSp + sigmaa;
    Float alphaP = sigmaSp / sigmaTp;
    Float effAlphaP = pbrt::effectiveAlbedo(alphaP);
    Float val = 0.25f * g + 0.25 * alphaP + 1.0 * effAlphaP;
    return kernel_multiplier * 4.0f * val * val / (sigmaTp * sigmaTp);
}

std::tuple<Polynomial, std::vector<Point3f>, std::vector<Vector3f>>
PolyUtils::FitPolynomial(const PolyUtils::PolyFitRecord &polyFitRecord, const ConstraintKDTree *kdTree) {
    if (polyFitRecord.config.hardSurfaceConstraint) {
        if (polyFitRecord.config.order == 2) {
            return FitPolynomialImpl<2, true>(polyFitRecord, kdTree);
        } else {
            return FitPolynomialImpl<3, true>(polyFitRecord, kdTree);
        }
    } else {
        LOG(ERROR) << "Polynomial without hard surface constraint is unsupported.";
    }
}

template<size_t polyOrder, bool hardSurfaceConstraint = true>
static std::tuple<Polynomial, std::vector<Point3f>, std::vector<Vector3f>>
PolyUtils::FitPolynomialImpl(const PolyUtils::PolyFitRecord &pfRec, const ConstraintKDTree *kdTree) {
    // TODO: Fit polynomial
    Float kernelEps = pfRec.kernelEps;
    std::function<Float(Float, Float)> kernelFunc = PolyUtils::GaussianKernel;
    std::vector<Point3f> posConstraints;
    std::vector<Normal3f> normConstraints;
    std::vector<Float> sampleWeights;
}

}