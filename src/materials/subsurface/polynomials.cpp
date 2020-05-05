//
// Created by 范軒瑋 on 2020-04-28.
//

#include "polynomials.h"
#include "geometry.h"
#include "scene.h"
#include "bssrdf.h"
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <core/interaction.h>

namespace pbrt {

size_t NumPolynomialCoefficients(size_t deg) {
    return (deg + 1) * (deg + 2) * (deg + 3) / 6;
}

template<size_t polyOrder, bool hardSurfaceConstraint = true>
void basisFunFitBuildA(const Point3f &pos, const Vector3f &inDir, const std::vector<Point3f> &evalP,
                       const Eigen::VectorXi &permX, const Eigen::VectorXi &permY, const Eigen::VectorXi &permZ,
                       const Eigen::VectorXf *weights,
                       Eigen::Matrix<float, Eigen::Dynamic, PolyUtils::nPolyCoeffs(polyOrder,
                                                                                   hardSurfaceConstraint)> &A,
                       const Eigen::VectorXf &weightedB, Float scaleFactor) {

    size_t n = evalP.size();
    Eigen::Matrix<float, Eigen::Dynamic, 3 * (polyOrder + 1)> relPosPow(n, 3 * (polyOrder + 1));
    for (size_t i = 0; i < n; ++i) {
        Vector rel = (evalP[i] - pos) * scaleFactor;
        for (size_t d = 0; d <= polyOrder; ++d) {
            relPosPow(i, d * 3 + 0) = std::pow(rel.x, d);
            relPosPow(i, d * 3 + 1) = std::pow(rel.y, d);
            relPosPow(i, d * 3 + 2) = std::pow(rel.z, d);
        }
    }

    PBRT_CONSTEXPR int nCoeffs = PolyUtils::nPolyCoeffs(polyOrder, false);
    Eigen::Matrix<float, Eigen::Dynamic, nCoeffs> fullA =
            Eigen::Matrix<float, Eigen::Dynamic, nCoeffs>::Zero(n * 4, nCoeffs);
    size_t termIdx = 0;
    for (size_t d = 0; d <= polyOrder; ++d) {
        for (size_t i = 0; i <= d; ++i) {
            size_t dx = d - i;
            for (size_t j = 0; j <= i; ++j) {
                size_t dy = d - dx - j;
                size_t dz = d - dx - dy;
                Eigen::VectorXf col = relPosPow.col(0 + 3 * dx).array() * relPosPow.col(1 + 3 * dy).array() *
                                      relPosPow.col(2 + 3 * dz).array();
                if (weights) {
                    col = col.array() * (*weights).array();
                }
                fullA.block(0, termIdx, n, 1) = col;
                const int pX = permX[termIdx];
                const int pY = permY[termIdx];
                const int pZ = permZ[termIdx];
                if (pX > 0) {
                    fullA.block(n, pX, n, 1) = (dx + 1) * col;
                }
                if (pY > 0) {
                    fullA.block(2 * n, pY, n, 1) = (dy + 1) * col;
                }
                if (pZ > 0) {
                    fullA.block(3 * n, pZ, n, 1) = (dz + 1) * col;
                }
                ++termIdx;
            }
        }
    }
    A = fullA.block(0, 1, fullA.rows(), fullA.cols() - 1);
}

Eigen::VectorXi DerivPermutationEigen(size_t degree, size_t axis) {
    auto numCoeffs = NumPolynomialCoefficients(degree);
    Eigen::VectorXi permutation = Eigen::VectorXi::Constant(numCoeffs, -1);
    for (size_t d = 0; d <= degree; ++d) {
        for (size_t i = 0; i <= d; ++i) {
            size_t dx = d - i;
            for (size_t j = 0; j <= i; ++j) {
                size_t dy = d - dx - j;
                size_t dz = d - dx - dy;
                Vector3i deg(dx, dy, dz);
                Vector3i derivDeg = deg;
                derivDeg[axis] -= 1;
                if (derivDeg[0] < 0 || derivDeg[1] < 0 || derivDeg[2] < 0) {
                    continue;
                }
                // For a valid derivative: add entry to matrix
                permutation[PolyUtils::PowerToIndex(derivDeg[0], derivDeg[1], derivDeg[2])] = PolyUtils::PowerToIndex(
                        dx, dy, dz);
            }
        }
    }
    return permutation;
}


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
    if (mTree.HasRightChild(index)) {
        auto rightIdx = node.GetRightIndex(index);
        GetConstraints(p, mTree[rightIdx], rightIdx, rightBounds, points, normals, sampleWeights, kernelEps,
                       kernelFunc);
    }
    // Node always has left child by construction
    auto leftIdx = node.GetLeftIndex(index);
    GetConstraints(p, mTree[leftIdx], leftIdx, leftBounds, points, normals, sampleWeights, kernelEps, kernelFunc);
}

Float PolyUtils::GetKernelEps(const MediumParameters &mediumParas, int channel,
                              pbrt::Float kernelMultiplier) {
    // The kernel epsilon size is related to the medium's properties
    Float sigmaT = mediumParas.sigmaT[channel];
    Float albedo = mediumParas.albedo[channel];
    Float g = mediumParas.g;
    Float sigmaS = albedo * sigmaT;
    Float sigmaa = sigmaT - sigmaS;
    Float sigmaSp = (1 - g) * sigmaS;
    Float sigmaTp = sigmaSp + sigmaa;
    Float alphaP = sigmaSp / sigmaTp;
    Float effAlphaP = pbrt::effectiveAlbedo(alphaP);
    Float val = 0.25f * g + 0.25 * alphaP + 1.0 * effAlphaP;
    return kernelMultiplier * 4.0f * val * val / (sigmaTp * sigmaTp);
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

    std::tie(posConstraints, normConstraints, sampleWeights) = kdTree->GetConstraints(pfRec.p, kernelEps, kernelFunc);

    size_t n = posConstraints.size();
    float invSqrtN = 1.0f / std::sqrt(n);

    Eigen::VectorXf weights(n);
    Vector3f s, t;
    OnbDuff(pfRec.d, s, t);

    PolyUtils::Frame local(s, t, pfRec.d);
    bool &&useLightspace = pfRec.config.useLightspace;

    for (size_t i = 0; i < n; i++) {
        Float d2 = DistanceSquared(pfRec.p, posConstraints[i]), w;
        if (sampleWeights[i] < 0) {
            w = pfRec.config.globalConstraintWeight * std::sqrt(1.0f / 32.0f);
        } else {
            w = std::sqrt(kernelFunc(d2, kernelEps) * sampleWeights[i]) * invSqrtN;
        }

        weights[i] = w;
        if (useLightspace) {
            auto localPos = local.toLocal(posConstraints[i] - pfRec.p);
            posConstraints[i] = Point3f(localPos.x, localPos.y, localPos.z);
        }
    }

    Eigen::VectorXf weightedB(4 * n);
    for (size_t i = 0; i < n; i++) {
        Normal3f normal = normConstraints[i];
        if (useLightspace) {
            normal = local.toLocal(normal);
        }
        weightedB[i] = 0.0f;
        weightedB[i + 1 * n] = normal.x * weights[i];
        weightedB[i + 2 * n] = normal.y * weights[i];
        weightedB[i + 3 * n] = normal.z * weights[i];
    }

    Eigen::VectorXi pX = DerivPermutationEigen(polyOrder, 0);
    Eigen::VectorXi pY = DerivPermutationEigen(polyOrder, 1);
    Eigen::VectorXi pZ = DerivPermutationEigen(polyOrder, 2);

    PBRT_CONSTEXPR size_t nCoeffs = nPolyCoeffs(polyOrder, hardSurfaceConstraint);

    Eigen::Matrix<float, Eigen::Dynamic, nCoeffs> A(4 * n, nCoeffs);
    Eigen::Matrix<float, nCoeffs, nCoeffs> AtA(nCoeffs, nCoeffs);
    // This scale factor seems to lead to a well behaved fit in many different settings
    float fitScaleFactor = PolyUtils::GetFitScaleFactor(kernelEps);
    Vector3f usedRefDir = useLightSpace ? local.toLocal(pfRec.n) : pfRec.n;

    if (useLightSpace) {
        basisFunFitBuildA<polyOrder, hardSurfaceConstraint>(Point3f(0.0f), usedRefDir, positionConstraints,
                                                            pX, pY, pZ, &weights, A, weightedB, fitScaleFactor);
    } else {
        basisFunFitBuildA<polyOrder, hardSurfaceConstraint>(pfRec.p, usedRefDir, positionConstraints, pX, pY, pZ,
                                                            &weights, A, weightedB, fitScaleFactor);
    }
    Eigen::Matrix<float, nCoeffs, 1> Atb = A.transpose() * weightedB;


}

template<int order, typename T>
static Eigen::Matrix<float, nPolyCoeffs(order), 1> PolyUtils::RotatePolynomialEigen(
        const T &c,
        const Vector3f &s,
        const Vector3f &t,
        const Normal3f &n) {
    Eigen::Matrix<float, nPolyCoeffs(order), 1> c2;
    c2[0] = c[0];
    c2[1] = c[1] * s[0] + c[2] * s[1] + c[3] * s[2];
    c2[2] = c[1] * t[0] + c[2] * t[1] + c[3] * t[2];
    c2[3] = c[1] * n[0] + c[2] * n[1] + c[3] * n[2];
    c2[4] = c[4] * powi(s[0], 2) + c[5] * s[0] * s[1] + c[6] * s[0] * s[2] + c[7] * powi(s[1], 2) + c[8] * s[1] * s[2] +
            c[9] * powi(s[2], 2);
    c2[5] = 2 * c[4] * s[0] * t[0] + c[5] * (s[0] * t[1] + s[1] * t[0]) + c[6] * (s[0] * t[2] + s[2] * t[0]) +
            2 * c[7] * s[1] * t[1] + c[8] * (s[1] * t[2] + s[2] * t[1]) + 2 * c[9] * s[2] * t[2];
    c2[6] = 2 * c[4] * n[0] * s[0] + c[5] * (n[0] * s[1] + n[1] * s[0]) + c[6] * (n[0] * s[2] + n[2] * s[0]) +
            2 * c[7] * n[1] * s[1] + c[8] * (n[1] * s[2] + n[2] * s[1]) + 2 * c[9] * n[2] * s[2];
    c2[7] = c[4] * powi(t[0], 2) + c[5] * t[0] * t[1] + c[6] * t[0] * t[2] + c[7] * powi(t[1], 2) + c[8] * t[1] * t[2] +
            c[9] * powi(t[2], 2);
    c2[8] = 2 * c[4] * n[0] * t[0] + c[5] * (n[0] * t[1] + n[1] * t[0]) + c[6] * (n[0] * t[2] + n[2] * t[0]) +
            2 * c[7] * n[1] * t[1] + c[8] * (n[1] * t[2] + n[2] * t[1]) + 2 * c[9] * n[2] * t[2];
    c2[9] = c[4] * powi(n[0], 2) + c[5] * n[0] * n[1] + c[6] * n[0] * n[2] + c[7] * powi(n[1], 2) + c[8] * n[1] * n[2] +
            c[9] * powi(n[2], 2);
    if (order > 2) {
        c2[10] = c[10] * powi(s[0], 3) + c[11] * powi(s[0], 2) * s[1] + c[12] * powi(s[0], 2) * s[2] +
                 c[13] * s[0] * powi(s[1], 2) + c[14] * s[0] * s[1] * s[2] + c[15] * s[0] * powi(s[2], 2) +
                 c[16] * powi(s[1], 3) + c[17] * powi(s[1], 2) * s[2] + c[18] * s[1] * powi(s[2], 2) +
                 c[19] * powi(s[2], 3);
        c2[11] = 3 * c[10] * powi(s[0], 2) * t[0] + c[11] * (powi(s[0], 2) * t[1] + 2 * s[0] * s[1] * t[0]) +
                 c[12] * (powi(s[0], 2) * t[2] + 2 * s[0] * s[2] * t[0]) +
                 c[13] * (2 * s[0] * s[1] * t[1] + powi(s[1], 2) * t[0]) +
                 c[14] * (s[0] * s[1] * t[2] + s[0] * s[2] * t[1] + s[1] * s[2] * t[0]) +
                 c[15] * (2 * s[0] * s[2] * t[2] + powi(s[2], 2) * t[0]) + 3 * c[16] * powi(s[1], 2) * t[1] +
                 c[17] * (powi(s[1], 2) * t[2] + 2 * s[1] * s[2] * t[1]) +
                 c[18] * (2 * s[1] * s[2] * t[2] + powi(s[2], 2) * t[1]) + 3 * c[19] * powi(s[2], 2) * t[2];
        c2[12] = 3 * c[10] * n[0] * powi(s[0], 2) + c[11] * (2 * n[0] * s[0] * s[1] + n[1] * powi(s[0], 2)) +
                 c[12] * (2 * n[0] * s[0] * s[2] + n[2] * powi(s[0], 2)) +
                 c[13] * (n[0] * powi(s[1], 2) + 2 * n[1] * s[0] * s[1]) +
                 c[14] * (n[0] * s[1] * s[2] + n[1] * s[0] * s[2] + n[2] * s[0] * s[1]) +
                 c[15] * (n[0] * powi(s[2], 2) + 2 * n[2] * s[0] * s[2]) + 3 * c[16] * n[1] * powi(s[1], 2) +
                 c[17] * (2 * n[1] * s[1] * s[2] + n[2] * powi(s[1], 2)) +
                 c[18] * (n[1] * powi(s[2], 2) + 2 * n[2] * s[1] * s[2]) + 3 * c[19] * n[2] * powi(s[2], 2);
        c2[13] = 3 * c[10] * s[0] * powi(t[0], 2) + c[11] * (2 * s[0] * t[0] * t[1] + s[1] * powi(t[0], 2)) +
                 c[12] * (2 * s[0] * t[0] * t[2] + s[2] * powi(t[0], 2)) +
                 c[13] * (s[0] * powi(t[1], 2) + 2 * s[1] * t[0] * t[1]) +
                 c[14] * (s[0] * t[1] * t[2] + s[1] * t[0] * t[2] + s[2] * t[0] * t[1]) +
                 c[15] * (s[0] * powi(t[2], 2) + 2 * s[2] * t[0] * t[2]) + 3 * c[16] * s[1] * powi(t[1], 2) +
                 c[17] * (2 * s[1] * t[1] * t[2] + s[2] * powi(t[1], 2)) +
                 c[18] * (s[1] * powi(t[2], 2) + 2 * s[2] * t[1] * t[2]) + 3 * c[19] * s[2] * powi(t[2], 2);
        c2[14] = 6 * c[10] * n[0] * s[0] * t[0] +
                 c[11] * (2 * n[0] * s[0] * t[1] + 2 * n[0] * s[1] * t[0] + 2 * n[1] * s[0] * t[0]) +
                 c[12] * (2 * n[0] * s[0] * t[2] + 2 * n[0] * s[2] * t[0] + 2 * n[2] * s[0] * t[0]) +
                 c[13] * (2 * n[0] * s[1] * t[1] + 2 * n[1] * s[0] * t[1] + 2 * n[1] * s[1] * t[0]) + c[14] *
                                                                                                      (n[0] * s[1] *
                                                                                                       t[2] +
                                                                                                       n[0] * s[2] *
                                                                                                       t[1] +
                                                                                                       n[1] * s[0] *
                                                                                                       t[2] +
                                                                                                       n[1] * s[2] *
                                                                                                       t[0] +
                                                                                                       n[2] * s[0] *
                                                                                                       t[1] +
                                                                                                       n[2] * s[1] *
                                                                                                       t[0]) +
                 c[15] * (2 * n[0] * s[2] * t[2] + 2 * n[2] * s[0] * t[2] + 2 * n[2] * s[2] * t[0]) +
                 6 * c[16] * n[1] * s[1] * t[1] +
                 c[17] * (2 * n[1] * s[1] * t[2] + 2 * n[1] * s[2] * t[1] + 2 * n[2] * s[1] * t[1]) +
                 c[18] * (2 * n[1] * s[2] * t[2] + 2 * n[2] * s[1] * t[2] + 2 * n[2] * s[2] * t[1]) +
                 6 * c[19] * n[2] * s[2] * t[2];
        c2[15] = 3 * c[10] * powi(n[0], 2) * s[0] + c[11] * (powi(n[0], 2) * s[1] + 2 * n[0] * n[1] * s[0]) +
                 c[12] * (powi(n[0], 2) * s[2] + 2 * n[0] * n[2] * s[0]) +
                 c[13] * (2 * n[0] * n[1] * s[1] + powi(n[1], 2) * s[0]) +
                 c[14] * (n[0] * n[1] * s[2] + n[0] * n[2] * s[1] + n[1] * n[2] * s[0]) +
                 c[15] * (2 * n[0] * n[2] * s[2] + powi(n[2], 2) * s[0]) + 3 * c[16] * powi(n[1], 2) * s[1] +
                 c[17] * (powi(n[1], 2) * s[2] + 2 * n[1] * n[2] * s[1]) +
                 c[18] * (2 * n[1] * n[2] * s[2] + powi(n[2], 2) * s[1]) + 3 * c[19] * powi(n[2], 2) * s[2];
        c2[16] = c[10] * powi(t[0], 3) + c[11] * powi(t[0], 2) * t[1] + c[12] * powi(t[0], 2) * t[2] +
                 c[13] * t[0] * powi(t[1], 2) + c[14] * t[0] * t[1] * t[2] + c[15] * t[0] * powi(t[2], 2) +
                 c[16] * powi(t[1], 3) + c[17] * powi(t[1], 2) * t[2] + c[18] * t[1] * powi(t[2], 2) +
                 c[19] * powi(t[2], 3);
        c2[17] = 3 * c[10] * n[0] * powi(t[0], 2) + c[11] * (2 * n[0] * t[0] * t[1] + n[1] * powi(t[0], 2)) +
                 c[12] * (2 * n[0] * t[0] * t[2] + n[2] * powi(t[0], 2)) +
                 c[13] * (n[0] * powi(t[1], 2) + 2 * n[1] * t[0] * t[1]) +
                 c[14] * (n[0] * t[1] * t[2] + n[1] * t[0] * t[2] + n[2] * t[0] * t[1]) +
                 c[15] * (n[0] * powi(t[2], 2) + 2 * n[2] * t[0] * t[2]) + 3 * c[16] * n[1] * powi(t[1], 2) +
                 c[17] * (2 * n[1] * t[1] * t[2] + n[2] * powi(t[1], 2)) +
                 c[18] * (n[1] * powi(t[2], 2) + 2 * n[2] * t[1] * t[2]) + 3 * c[19] * n[2] * powi(t[2], 2);
        c2[18] = 3 * c[10] * powi(n[0], 2) * t[0] + c[11] * (powi(n[0], 2) * t[1] + 2 * n[0] * n[1] * t[0]) +
                 c[12] * (powi(n[0], 2) * t[2] + 2 * n[0] * n[2] * t[0]) +
                 c[13] * (2 * n[0] * n[1] * t[1] + powi(n[1], 2) * t[0]) +
                 c[14] * (n[0] * n[1] * t[2] + n[0] * n[2] * t[1] + n[1] * n[2] * t[0]) +
                 c[15] * (2 * n[0] * n[2] * t[2] + powi(n[2], 2) * t[0]) + 3 * c[16] * powi(n[1], 2) * t[1] +
                 c[17] * (powi(n[1], 2) * t[2] + 2 * n[1] * n[2] * t[1]) +
                 c[18] * (2 * n[1] * n[2] * t[2] + powi(n[2], 2) * t[1]) + 3 * c[19] * powi(n[2], 2) * t[2];
        c2[19] = c[10] * powi(n[0], 3) + c[11] * powi(n[0], 2) * n[1] + c[12] * powi(n[0], 2) * n[2] +
                 c[13] * n[0] * powi(n[1], 2) + c[14] * n[0] * n[1] * n[2] + c[15] * n[0] * powi(n[2], 2) +
                 c[16] * powi(n[1], 3) + c[17] * powi(n[1], 2) * n[2] + c[18] * n[1] * powi(n[2], 2) +
                 c[19] * powi(n[2], 3);
    }
    return c2;
}

const static int PERMX2[10] = {1, 4, 5, 6, -1, -1, -1, -1, -1, -1};
const static int PERMY2[10] = {2, 5, 7, 8, -1, -1, -1, -1, -1, -1};
const static int PERMZ2[10] = {3, 6, 8, 9, -1, -1, -1, -1, -1, -1};
const static int PERMX3[20] = {1, 4, 5, 6, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
const static int PERMY3[20] = {2, 5, 7, 8, 11, 13, 14, 16, 17, 18, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
const static int PERMZ3[20] = {3, 6, 8, 9, 12, 14, 15, 17, 18, 19, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

const int *DerivPermutation(size_t degree, size_t axis) {
    switch (axis) {
        case 0:
            return degree == 2 ? PERMX2 : PERMX3;
        case 1:
            return degree == 2 ? PERMY2 : PERMY3;
        case 2:
            return degree == 2 ? PERMZ2 : PERMZ3;
        default:
            return nullptr;
    }
}

template<typename T>
std::pair<float, Vector3f> EvalPolyGrad(const Point3f &pos,
                                        const Point3f &evalP, size_t degree,
                                        const int *permX, const int *permY, const int *permZ,
                                        Float scaleFactor, bool useLocalDir, const Vector3f &refDir,
                                        const T &coeffs) {
    Vector3f relPos;
    if (useLocalDir) {
        Vector3f s, t;
        PolyUtils::OnbDuff(Normal3f(refDir), s, t);
        PolyUtils::Frame local(s, t, Normal3f(refDir));
        relPos = local.toLocal(evalP - pos) * scaleFactor;
    } else {
        relPos = (evalP - pos) * scaleFactor;
    }

    size_t termIdx = 0;
    Vector3f deriv(0.0f, 0.0f, 0.0f);
    float value = 0.0f;

    float xPowers[4] = {1.0, relPos.x, relPos.x * relPos.x, relPos.x * relPos.x * relPos.x};
    float yPowers[4] = {1.0, relPos.y, relPos.y * relPos.y, relPos.y * relPos.y * relPos.y};
    float zPowers[4] = {1.0, relPos.z, relPos.z * relPos.z, relPos.z * relPos.z * relPos.z};
    for (size_t d = 0; d <= degree; ++d) {
        for (size_t i = 0; i <= d; ++i) {
            size_t dx = d - i;
            for (size_t j = 0; j <= i; ++j) {
                size_t dy = d - dx - j;
                size_t dz = d - dx - dy;
                float t = xPowers[dx] * yPowers[dy] * zPowers[dz];
                value += coeffs[termIdx] * t;

                int pX = permX[termIdx];
                int pY = permY[termIdx];
                int pZ = permZ[termIdx];
                if (pX > 0)
                    deriv.x += (dx + 1) * t * coeffs[pX];
                if (pY > 0)
                    deriv.y += (dy + 1) * t * coeffs[pY];
                if (pZ > 0)
                    deriv.z += (dz + 1) * t * coeffs[pZ];
                ++termIdx;
            }
        }
    }
    return std::make_pair(value, deriv);
}

void PolyUtils::ProjectPointsToSurface(const Scene *scene, const Point3f &refPoint, const Vector3f &refDir,
                                       ScatterSamplingRecord &rec, const Eigen::VectorXf &polyCoefficients,
                                       size_t polyOrder, bool useLocalDir, Float scaleFactor, Float kernelEps) {
    if (!rec.isValid)
        return;

    Vector3f dir = EvaluateGradient(refPoint, polyCoefficients, rec, polyOrder, scaleFactor, useLocalDir, refDir);
    Float dists[2] = {2 * kernelEps, std::numeric_limits<Float>::infinity()};

    for (int i = 0; i < 2; i++) {
        Float maxDist = dists[i];
        Ray r1(rec.p, dir);
        SurfaceInteraction its;
        Point3f projectedP;
        Normal3f normal;
        Float pointTime = -1.0f;
        bool itsFound = false;
        if (scene->Intersect(r1, &its)) {
            projectedP = its.p;
            normal = its.shading.n; // FIXME: Note the difference between shading.n and n
            pointTime = its.time;
            itsFound = true;
        }
        Float maxT = itsFound ? its.time : maxDist;
        Ray r2(rec.p, -dir);
        if (scene->Intersect(r2, &its)) {
            if (pointTime < 0 || pointTime > its.time) {
                projectedP = its.p;
                normal = its.shading.n;
            }
            itsFound = true;
        }
        rec.isValid = itsFound;
        if (itsFound) {
            rec.p = projectedP;
            rec.n = normal;
            return;
        }
    }
}

Vector3f
PolyUtils::EvaluateGradient(const Point3f &pos, const Eigen::VectorXf &coeffs, const ScatterSamplingRecord &rec,
                            size_t degree, Float scaleFactor,
                            bool useLocalDir, const Vector3f &refDir) {
    const int *permX = DerivPermutation(degree, 0);
    const int *permY = DerivPermutation(degree, 1);
    const int *permZ = DerivPermutation(degree, 2);
    Float polyVal;
    Vector3f gradient;
    std::tie(polyVal, gradient) = EvalPolyGrad(pos, rec.p, degree, permX, permY, permZ, scaleFactor, useLocalDir,
                                               refDir, coeffs);
    if (useLocalDir) {
        Vector3f s, t;
        PolyUtils::OnbDuff(Normal3f(refDir), s, t);
        PolyUtils::Frame local(s, t, Normal3f(refDir));
        gradient = local.toWorld(gradient);
    }
    return pbrt::Vector3f();
}

Normal3f PolyUtils::AdjustRayDirForPolynomialTracing(Vector3f &inDir, const SurfaceInteraction &isect, int polyOrder,
                                                     Float polyScaleFactor, int channel) {
    const int *pX = DerivPermutation(polyOrder, 0);
    const int *pY = DerivPermutation(polyOrder, 1);
    const int *pZ = DerivPermutation(polyOrder, 2);

    float polyValue;
    Vector3f polyNormal;
    std::tie(polyValue, polyNormal) = EvalPolyGrad(isect.p, isect.p, polyOrder, pX, pY, pZ,
                                                   polyScaleFactor, false,
                                                   inDir, isect.polyStorage->coeffs[channel]);
    Vector3f rotationAxis = pbrt::Cross(isect.shading.n, polyNormal);
    if (rotationAxis.Length() < 1e-8f) {
        return (Normal3f) polyNormal;
    }
    Normal3f normalizedTarget = isect.shading.n;
    float angle = acos(std::max(std::min(Dot(polyNormal, normalizedTarget), 1.0f), -1.0f));
    Transform transf = pbrt::Rotate(Degrees(angle), rotationAxis);
    inDir = transf(inDir);
    return (Normal3f) polyNormal;
}

} // namespace pbrt.