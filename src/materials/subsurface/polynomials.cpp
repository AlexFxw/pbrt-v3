//
// Created by 范軒瑋 on 2020-04-28.
//

#include "polynomials.h"
#include "geometry.h"
#include "scene.h"
#include "bssrdf.h"
#include "shapes/triangle.h"
#include "interaction.h"
#include "parallel.h"
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>

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
        Vector3f rel = (evalP[i] - pos) * scaleFactor;
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


void ConstraintKDTree::Build(const std::vector<Point3f> &sampledP, const std::vector<Normal3f> &sampledN) {
    int pointsNum = sampledP.size();
    mTree = CTree(pointsNum);
    LOG(INFO) << "Build constraint kd tree: " << pointsNum;

#ifdef VISUALIZE_SHAPE_DATA
    {
        // FIXME: visualize KD tree.
        const std::string &fileName = "../data/kdtree.txt";
    }
#endif
    ParallelFor([&](int64_t i) {
        mTree[i].SetPosition(sampledP[i]);
        ExtraData d;
        d.p = sampledP[i];
        d.n = sampledN[i];
        d.sampledNum = 1;
        d.avgN = Normal3f(-10, -10, -10);
        d.avgP = Point3f(-100, -100, -100);
        mTree[i].SetData(d);
    }, pointsNum, 1024);
    mTree.Build(true);
    CalcAvgValues(mTree[0], 0);
    LOG(INFO) << "Calculated the average value";

    for (size_t i = 0; i < std::min(size_t(32), sampledP.size()); i++) {
        mPoints.push_back(sampledP[i]);
        mNormals.push_back(sampledN[i]);
    }
}

std::tuple<Point3f, Normal3f, size_t> ConstraintKDTree::CalcAvgValues(TreeNode &node, TreeNode::IndexType index) {
    LOG(INFO) << "Calculating average value of index: " << index;
    Point3f rp, lp;
    Normal3f rn, ln;
    size_t rSamples = 0, lSamples = 0;
    if (mTree.HasRightChild(index)) {
        auto rightIdx = node.GetRightIndex(index);
        if (index != rightIdx)
            std::tie(rp, rn, rSamples) = CalcAvgValues(mTree[rightIdx], rightIdx);
    }

    if (!node.IsLeaf()) {
        auto leftIdx = node.GetLeftIndex(index);
        if (index != leftIdx)
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

std::tuple<std::vector<Point3f>, std::vector<Normal3f>, std::vector<Float>>
ConstraintKDTree::GetConstraints(const Point3f &p, Float kernelEps,
                                 const std::function<Float(Float, Float)> &kernel) const {
    // Extract constraints from KD Tree by traversing from the top
    std::vector<Point3f> positions;
    std::vector<Normal3f> normals;
    std::vector<Float> sampleWeights;

    GetConstraints(p, mTree[0], 0, mTree.GetAABB(), positions, normals, sampleWeights, kernelEps, kernel);
    std::cout << "Constraints num: " << positions.size() << "/" << mPoints.size() << std::endl;

    // // Add the super constraints
    for (size_t i = 0; i < mPoints.size(); ++i) {
        positions.push_back(mPoints[i]);
        normals.push_back(mNormals[i]);
        sampleWeights.push_back(-1.0f);
    }
    return std::make_tuple(positions, normals, sampleWeights);
}

void ConstraintKDTree::GetConstraints(const Point3f &p, SimpleKDNode<Point3f,
        pbrt::ConstraintKDTree::ExtraData> node,
                                      TreeNode::IndexType index, const Bounds3f &aabb,
                                      std::vector<Point3f> &points,
                                      std::vector<Normal3f> &normals, std::vector<Float> &sampleWeights,
                                      pbrt::Float kernelEps,
                                      const std::function<Float(Float, Float)> &kernelFunc) const {
    Float dist2Threshold = 9.0f * kernelEps; // FIXME: Adjust a better threshold.
    Float distSq = pbrt::DistanceSquared(p, aabb);
    if (distSq > dist2Threshold) {
        return;
    }
    if (distSq < dist2Threshold) {
        points.push_back(node.GetData().p);
        normals.push_back(node.GetData().n);
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
        // if (rightIdx != index)
        GetConstraints(p, mTree[rightIdx], rightIdx, rightBounds,
                       points, normals, sampleWeights, kernelEps, kernelFunc);
    }
    // Node always has left child by construction
    auto leftIdx = node.GetLeftIndex(index);
    GetConstraints(p, mTree[leftIdx], leftIdx, leftBounds, points, normals, sampleWeights, kernelEps, kernelFunc);

#ifdef VISUALIZE_SHAPE_DATA
    {
        // FIXME: Make visualize a parameter
        std::ofstream file;
        const std::string fileName = "../data/constraintpct.txt";
        file.open(fileName, std::ios::app);
        DCHECK(file.is_open());
        auto r = GetRandomFloat(), g = GetRandomFloat(), b = GetRandomFloat();
        for (auto &p: points) {
            file << p.x << " " << p.y << " " << p.z << " " << r << " " << g << " " << b << std::endl;
        }
        file.close();
    }
#endif
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
    Float effAlphaP = pbrt::EffectiveAlbedo(alphaP);
    Float val = 0.25f * g + 0.25 * alphaP + 1.0 * effAlphaP;
    return kernelMultiplier * 4.0f * val * val / (sigmaTp * sigmaTp);
}

// FIXME: Use kdtree
std::tuple<PolyUtils::Polynomial, std::vector<Point3f>, std::vector<Normal3f>>
// PolyUtils::FitPolynomial(const PolyFitRecord &polyFitRecord, const ConstraintKDTree *kdTree) {
PolyUtils::FitPolynomial(const PolyFitRecord &polyFitRecord,
                         const std::vector<Point3f> &points, const std::vector<Normal3f> &normals) {
    if (polyFitRecord.config.hardSurfaceConstraint) {
        if (polyFitRecord.config.order == 2) {
            // return FitPolynomialImpl<2, true>(polyFitRecord, kdTree);
            return FitPolynomialImpl<2, true>(polyFitRecord, points, normals);
        } else {
            // return FitPolynomialImpl<3, true>(polyFitRecord, kdTree);
            return FitPolynomialImpl<3, true>(polyFitRecord, points, normals);
        }
    } else {
        LOG(ERROR) << "Polynomial without hard surface constraint is unsupported.";
        exit(1);
    }
}

template<size_t polyOrder, bool hardSurfaceConstraint>
std::tuple<PolyUtils::Polynomial, std::vector<Point3f>, std::vector<Normal3f>>
// PolyUtils::FitPolynomialImpl(const PolyFitRecord &pfRec, const ConstraintKDTree *kdTree) {
PolyUtils::FitPolynomialImpl(const PolyFitRecord &pfRec,
                             const std::vector<Point3f> &points, const std::vector<Normal3f> &normals) {
    Float kernelEps = pfRec.kernelEps;
    std::function<Float(Float, Float)> kernelFunc = PolyUtils::GaussianKernel;
    std::vector<Point3f> posConstraints;
    std::vector<Normal3f> normConstraints;
    std::vector<Float> sampleWeights;

    // Problems may come from here.
    // FIXME: Use kdtree instead
    // std::tie(posConstraints, normConstraints, sampleWeights) = kdTree->GetConstraints(pfRec.p, kernelEps, kernelFunc);
    std::tie(posConstraints, normConstraints, sampleWeights) = PolyUtils::GetConstraints(pfRec.p, kernelEps, points, normals);

#ifdef VISUALIZE_SHAPE_DATA
    {
        // FIXME: Visualize constraints
        std::cout << "collecting constraints" << std::endl;
        const std::string &fileName = "../data/constraints.txt";
        std::ofstream file;
        file.open(fileName, std::ios::app);
        file << pfRec.p.x << " " << pfRec.p.y << " " << pfRec.p.z << " " << 1.0 << " " << 0 << " " << 0 << std::endl;
        for(auto &p: posConstraints) {
            file << p.x << " " << p.y << " " << p.z << " " << 0 << " " << 1.0 << " " << 0 << std::endl;
        }
        file.close();
    }
#endif

    size_t n = posConstraints.size();
    Float invSqrtN = 1.0f / std::sqrt(n);

    Eigen::VectorXf weights(n);
    Vector3f s, t;
    Normal3f pfd = Normal3f(pfRec.d);
    OnbDuff(pfd, s, t);
    PolyUtils::Frame local(s, t, pfd);

    bool useLightspace = pfRec.config.useLightspace;

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
            normal = Normal3f(local.toLocal(Vector3f(normal)));
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
    Vector3f usedRefDir = useLightspace ? local.toLocal((Vector3f) pfRec.n) : (Vector3f) pfRec.n;

    if (useLightspace) {
        basisFunFitBuildA<polyOrder, hardSurfaceConstraint>(Point3f(), usedRefDir, posConstraints,
                                                            pX, pY, pZ, &weights, A, weightedB, fitScaleFactor);
    } else {
        basisFunFitBuildA<polyOrder, hardSurfaceConstraint>(pfRec.p, usedRefDir, posConstraints, pX, pY, pZ,
                                                            &weights, A, weightedB, fitScaleFactor);
    }
    Eigen::Matrix<float, nCoeffs, 1> Atb = A.transpose() * weightedB;

    Eigen::VectorXf coeffs;
    if (pfRec.config.useSvd) {
        Eigen::MatrixXf ADyn = A;
        Eigen::BDCSVD<Eigen::MatrixXf> svd = ADyn.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
        const Eigen::VectorXf &sVal = svd.singularValues();
        float eps = 0.01f;
        coeffs = svd.matrixV() * ((sVal.array() > eps).select(sVal.array().inverse(), 0)).matrix().asDiagonal() *
                 svd.matrixU().transpose() * weightedB;
    } else {
        Eigen::MatrixXf reg = Eigen::MatrixXf::Identity(A.cols(), A.cols()) * pfRec.config.regularization;
        reg(0, 0) = 0.0f;
        reg(1, 1) = 0.0f;
        reg(2, 2) = 0.0f;
        AtA = A.transpose() * A + reg;
        coeffs = AtA.ldlt().solve(Atb);
    }

    std::vector<Float> coeffsVec(NumPolynomialCoefficients(polyOrder));
    size_t cfSize = coeffsVec.size();
    if (hardSurfaceConstraint) {
        coeffsVec[0] = 0.0f;
        for (size_t i = 1; i < coeffs.size(); ++i)
            coeffsVec[i] = coeffs[i - 1];
    } else {
        for (size_t i = 0; i < coeffs.size(); ++i)
            coeffsVec[i] = coeffs[i];
    }

    PolyUtils::Polynomial poly;
    poly.coeffs = coeffsVec;
    poly.refPos = pfRec.p;
    poly.refDir = pfRec.d;
    poly.useLocalDir = pfRec.config.useLightspace;
    poly.scaleFactor = PolyUtils::GetFitScaleFactor(kernelEps);
    poly.order = polyOrder;
    return std::make_tuple(poly, posConstraints, normConstraints);
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
        relPos = local.toLocal(Vector3f(evalP - pos)) * scaleFactor;
    } else {
        relPos = Vector3f(evalP - pos) * scaleFactor;
    }

    size_t termIdx = 0;
    Vector3f deriv(0.0f, 0.0f, 0.0f);
    float value = 0.0f;

    Float xPowers[4] = {1.0, relPos.x, relPos.x * relPos.x, relPos.x * relPos.x * relPos.x};
    Float yPowers[4] = {1.0, relPos.y, relPos.y * relPos.y, relPos.y * relPos.y * relPos.y};
    Float zPowers[4] = {1.0, relPos.z, relPos.z * relPos.z, relPos.z * relPos.z * relPos.z};
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
                    deriv.x += ((float) dx + 1.0f) * t * coeffs[pX];
                if (pY > 0)
                    deriv.y += ((float) dy + 1.0f) * t * coeffs[pY];
                if (pZ > 0)
                    deriv.z += ((float) dz + 1.0f) * t * coeffs[pZ];
                ++termIdx;
            }
        }
    }
    return std::make_pair(value, deriv);
}

void PolyUtils::ProjectPointsToSurface(const Scene *scene, const Point3f &refPoint, const Vector3f &refDir,
                                       ScatterSamplingRecord &rec, const Eigen::VectorXf &polyCoefficients,
                                       size_t polyOrder, bool useLocalDir, Float scaleFactor, Float kernelEps,
                                       SurfaceInteraction *res) {
    if (!rec.isValid)
        return;

    Vector3f dir = EvaluateGradient(refPoint, polyCoefficients, rec, polyOrder, scaleFactor, useLocalDir, refDir);
    Float dists[2] = {2 * kernelEps, std::numeric_limits<Float>::infinity()};

    // FIXME: Adjust the dir.
    if (dir.Length() < 1e-8) {
        rec.isValid = false;
        return;
    }
    dir = Normalize(dir);

    for (int i = 0; i < 2; i++) {
        Float maxDist = dists[i];
        Ray r1(rec.p, dir); // FIXME: Add maxDist?
        SurfaceInteraction its, its2;
        Point3f projectedP;
        Normal3f normal;
        Float pointTime = -1.0f;
        bool itsFound = false;
        if (scene->Intersect(r1, &its)) {
            projectedP = its.p;
            normal = its.shading.n; // FIXME: Note the difference between shading.n and n
            pointTime = its.time;
            itsFound = true;
            *res = its;
        }
        Float maxT = itsFound ? its.time : maxDist;
        Ray r2(rec.p, -dir);
        if (scene->Intersect(r2, &its2)) {
            if (pointTime < 0 || pointTime > its2.time) {
                projectedP = its2.p;
                normal = its2.shading.n;
            }
            itsFound = true;
            *res = its2;
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
    // FIXME: May be bug here.
    return gradient;
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
                                                   inDir, isect.polyStorage.coeffs[channel]);
    Vector3f rotationAxis = pbrt::Cross(isect.shading.n, polyNormal);
    if (rotationAxis.Length() < 1e-8f) {
        return (Normal3f) polyNormal;
    }
    Normal3f normalizedTarget = isect.shading.n;
    float angle = acos(std::max(std::min(Dot(polyNormal, normalizedTarget), 1.0f), -1.0f));
    Transform transf = pbrt::Rotate(Degrees(angle), rotationAxis);
    inDir = transf(inDir);
    return Normal3f(polyNormal);
}

} // namespace pbrt.