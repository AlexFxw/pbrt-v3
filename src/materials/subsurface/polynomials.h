//
// Created by 范軒瑋 on 2020-04-28.
//

#ifndef PBRT_V3_POLYNOMIALS_H
#define PBRT_V3_POLYNOMIALS_H

#include "pbrt.h"
#include "geometry.h"
#include "core/medium.h"
#include <Eigen/Core>

namespace pbrt {

template<typename PointType, typename DataRecord>
struct SimpleKDNode {
    typedef uint32_t IndexType;
    enum {
        ELeafFlag = 0x10,
        EAxisMask = 0x0F
    };

    PointType pos;
    IndexType right;
    DataRecord data;
    uint8_t flags;

    const PointType &GetPos() const { return pos; }

    SimpleKDNode() : pos(PointType()), right(0), data(), flags(0) {}

    DataRecord &GetData() { return data; }

    SimpleKDNode(const DataRecord &data) : pos(PointType()), right(0), data(data), flags(0) {}

    IndexType GetRightIndex(IndexType right) const { return right; }

    IndexType GetLeftIndex(IndexType self) const { return self + 1; }

    inline bool IsLeaf() const { return flags & (uint8_t) ELeafFlag; }

    inline const PointType &GetPosition() const { return pos; }

    inline void SetPosition(const Point3f &pos) { this->pos = pos; }

    inline void SetData(const DataRecord &data) { this->data = data; }

    void SetLeaf(bool value) {
        if (value)
            flags |= (uint8_t) ELeafFlag;
        else
            flags &= (uint8_t)
                    ~ELeafFlag;
    }


    void SetRightIndex(IndexType self, IndexType value) {
        CHECK_NE(self, value);
        right = value;
    }

    void SetLeftIndex(IndexType self, IndexType value) {
        // Do nothing
    }

    // Get the split axis of this node.
    uint16_t GetAxis() const { return flags & (uint8_t) EAxisMask; }

    void SetAxis(uint8_t axis) {
        flags = (flags & (uint8_t)
                ~EAxisMask) | axis;
    }
};

template<typename NodeType>
class PointKDTree {
public:
    typedef typename NodeType::IndexType IndexType;

    PointKDTree(size_t nodes = 0) : mNodes(nodes), mDepth(0) {
        mAABB = Bounds3f(Point3f(), Point3f());
    }

    void Clear() {
        mNodes.clear();
        mAABB = Bounds3f(Point3f(), Point3f());
    }

    inline void Resize(size_t size) { mNodes.resize(size); }

    inline NodeType &operator[](size_t idx) { return mNodes[idx]; }

    inline const NodeType &operator[](size_t idx) const { return mNodes[idx]; }

    inline bool HasRightChild(IndexType index) const {
        // return mNodes[index].GetRightIndex(index) != 0;
        auto rightIndex = mNodes[index].GetRightIndex(index);
        return rightIndex != 0 && rightIndex != index;
    }

    void PushBack(const NodeType &node) {
        mNodes.push_back(node);
        mAABB = pbrt::Union(mAABB, node.pos);
    }

    inline size_t Size() const { return mNodes.size(); }

    inline const Bounds3f &GetAABB() const { return mAABB; }

    void Build(bool recomputeAABB = false) {
        LOG(INFO) << "Building the kd tree.";
        if (recomputeAABB) {
            mAABB = Bounds3f(Point3f(), Point3f());
            for (auto iter = mNodes.begin(); iter != mNodes.end(); iter++) {
                mAABB = pbrt::Union(mAABB, iter->pos);
            }
        }
        int mNodeSize = mNodes.size();
        std::vector<IndexType> indirection(mNodeSize);
        for (size_t i = 0; i < mNodes.size(); i++) {
            indirection[i] = (IndexType) i;
        }
        mDepth = 0;
        BuildImpl(1, indirection.begin(), indirection.begin(), indirection.end());
        permute_inplace(&mNodes[0], indirection);
    }

    void BuildImpl(size_t depth,
                   typename std::vector<IndexType>::iterator base,
                   typename std::vector<IndexType>::iterator rangeStart,
                   typename std::vector<IndexType>::iterator rangeEnd) {
        mDepth = std::max(depth, mDepth);
        IndexType count = (IndexType) (rangeEnd - rangeStart);
        if (count == 1) {
            // Create a leaf
            mNodes[*rangeStart].SetLeaf(true);
            return;
        }
        typename std::vector<IndexType>::iterator split;
        // Balanced strategy
        split = rangeStart + count / 2;
        int axis = mAABB.LargestAxis();
        std::nth_element(rangeStart, split, rangeEnd, CoordinateOrdering(mNodes, axis));
        NodeType &splitNode = mNodes[*split];
        splitNode.SetAxis(axis);
        splitNode.SetLeaf(false);
        if (split + 1 != rangeEnd && split != rangeStart) {
            splitNode.SetRightIndex((IndexType) (rangeStart - base), (IndexType) (split + 1 - base));
        } else {
            splitNode.SetRightIndex((IndexType) (rangeStart - base), 0);
        }
        splitNode.SetLeftIndex((IndexType) (rangeStart - base), (IndexType) (rangeStart + 1 - base));
        std::iter_swap(rangeStart, split);

        Float tmp = mAABB.pMax[axis], splitPos = splitNode.pos[axis];
        mAABB.pMax[axis] = splitPos;
        BuildImpl(depth + 1, base, rangeStart + 1, split + 1);
        mAABB.pMax[axis] = tmp;
        if (split + 1 != rangeEnd) {
            tmp = mAABB.pMin[axis];
            mAABB.pMin[axis] = splitPos;
            BuildImpl(depth + 1, base, split + 1, rangeEnd);
            mAABB.pMin[axis] = tmp;
        }
    }

protected:
    struct CoordinateOrdering : public std::binary_function<IndexType, IndexType, bool> {
    public:
        CoordinateOrdering(const std::vector<NodeType> &nodes, int axis)
                : m_nodes(nodes), m_axis(axis) {}

        bool operator()(const IndexType &i1, const IndexType &i2) const {
            return m_nodes[i1].GetPosition()[m_axis] < m_nodes[i2].GetPosition()[m_axis];
        }

    private:
        const std::vector<NodeType> &m_nodes;
        int m_axis;
    };

    std::vector<NodeType> mNodes;
    size_t mDepth;
    Bounds3f mAABB;
};


class ConstraintKDTree {
public:
    struct ExtraData {
        Point3f p, avgP;
        Normal3f n, avgN;
        int sampledNum;
    };

    typedef SimpleKDNode<Point3f, ExtraData> TreeNode;

    typedef PointKDTree<SimpleKDNode<Point3f, ExtraData>> CTree;

    void Build(const std::vector<Point3f> &sampledP, const std::vector<Normal3f> &sampledN);
    std::tuple<std::vector<Point3f>, std::vector<Normal3f>, std::vector<Float>>
    GetConstraints(const Point3f &p, Float kernelEps, const std::function<Float(Float, Float)> &kernel) const;


    std::tuple<Point3f, Normal3f, size_t> CalcAvgValues(TreeNode &node, TreeNode::IndexType index);

private:
    void GetConstraints(const Point3f &p, TreeNode node, TreeNode::IndexType index, const Bounds3f &aabb,
                        std::vector<Point3f> &points, std::vector<Normal3f> &normals, std::vector<Float> &sampleWeights,
                        Float kernelEps, const std::function<Float(Float, Float)> &kernelFunc) const;
    std::vector<Point3f> mPoints;
    std::vector<Normal3f> mNormals;
    CTree mTree;

};


class PolyUtils {
public:
    struct PolyFitConfig {
        float regularization = 0.0001f;
        bool useSvd = false;
        bool useLightspace = true;
        bool visualize = false;
        int order = 3;
        bool hardSurfaceConstraint = true;
        float globalConstraintWeight = 0.01f;
        float kdTreeThreshold = 0.0f;
        bool extractNormalHistogram = false;
        bool useSimilarityKernel = true;
        float kernelEpsScale = 1.0f;
    };

    struct PolyFitRecord {
        PolyFitConfig config;
        Point3f p;
        Vector3f d;
        Normal3f n;
        MediumParameters medium;
        float kernelEps;
    };

    struct Frame {
        Vector3f s, t;
        Normal3f n;

        Frame() {}

        Frame(const Vector3f &s, const Vector3f &t, const Normal3f &n) : s(s), t(t), n(n) {}

        inline Vector3f toLocal(const Vector3f &v) const {
            return Vector3f(
                    pbrt::Dot(v, s),
                    pbrt::Dot(v, t),
                    pbrt::Dot(v, n)
            );
        }

        /// Convert from local coordinates to world coordinates
        inline Vector3f toWorld(const Vector3f &v) const {
            return s * v.x + t * v.y + Vector3f(n) * v.z;
        }
    };

    struct Polynomial {
        std::vector<float> coeffs;
        Point3f refPos;
        Vector3f refDir;
        bool useLocalDir;
        float scaleFactor;
        int order;
        std::vector<float> normalHistogram;
    };

    static const size_t PowerToIndex(size_t dx, size_t dy, size_t dz) {
        // Converts a polynomial degree to a linear coefficient index
        auto d = dx + dy + dz;
        auto i = d - dx;
        auto j = d - dx - dy;
        return i * (i + 1) / 2 + j + d * (d + 1) * (d + 2) / 6;
    }

    static Float GetKernelEps(const MediumParameters &mediumParas, int channel, Float kernelMultiplier);

    static inline float GetFitScaleFactor(float kernelEps) { return 1.0f / std::sqrt(kernelEps); }

    static inline float
    GetFitScaleFactor(const MediumParameters &medium, int channel = 0, Float kernel_multiplier = 1.0f) {
        return GetFitScaleFactor(GetKernelEps(medium, channel, kernel_multiplier));
    }


    static PBRT_CONSTEXPR int nChooseK(int n, int k) {
        return (k == 0 || n == k) ? 1 : nChooseK(n - 1, k - 1) + nChooseK(n - 1, k);
    }

    static PBRT_CONSTEXPR int nPolyCoeffs(int polyOrder, bool hardSurfaceConstraint = false) {
        return hardSurfaceConstraint ? nChooseK(3 + polyOrder, polyOrder) - 1 : nChooseK(3 + polyOrder, polyOrder);
    }

    // FIXME: use kd tree
    static std::tuple<Polynomial, std::vector<Point3f>, std::vector<Normal3f>>
    // FitPolynomial(const PolyFitRecord &polyFitRecord, const ConstraintKDTree *kdTree);
    FitPolynomial(const PolyFitRecord &polyFitRecord, const std::vector<Point3f> &points, const std::vector<Normal3f> &normals);

    // FIXME: use kd tree
    template<size_t polyOrder, bool hardSurfaceConstraint = true>
    static std::tuple<Polynomial, std::vector<Point3f>, std::vector<Normal3f>>
    // FitPolynomialImpl(const PolyFitRecord &pfRec, const ConstraintKDTree *kdTree);
    FitPolynomialImpl(const PolyFitRecord &pfRec, const std::vector<Point3f> &points, const std::vector<Normal3f> &normals);

    static void OnbDuff(const Normal3f &n, Vector3f &b1, Vector3f &b2) {
        float sign = copysignf(1.0f, n.z);
        const float a = -1.0f / (sign + n.z);
        const float b = n.x * n.y * a;
        b1 = Vector3f(1.0f + sign * n.x * n.x * a, sign * b, -sign * n.x);
        b2 = Vector3f(b, sign + n.y * n.y * a, -n.y);
    }

    static inline float powi(float f, int n) {
        float ret = 1.0f;
        for (int i = 0; i < n; ++i) {
            ret *= f;
        }
        return ret;
    }

    template<int order, typename T>
    static Eigen::Matrix<float, nPolyCoeffs(order, false), 1> RotatePolynomialEigen(const T &c,
                                                                                    const Vector3f &s,
                                                                                    const Vector3f &t,
                                                                                    const Normal3f &n) {
        Eigen::Matrix<float, nPolyCoeffs(order), 1> c2;
        c2[0] = c[0];
        c2[1] = c[1] * s[0] + c[2] * s[1] + c[3] * s[2];
        c2[2] = c[1] * t[0] + c[2] * t[1] + c[3] * t[2];
        c2[3] = c[1] * n[0] + c[2] * n[1] + c[3] * n[2];
        c2[4] = c[4] * powi(s[0], 2) + c[5] * s[0] * s[1] + c[6] * s[0] * s[2] + c[7] * powi(s[1], 2) +
                c[8] * s[1] * s[2] +
                c[9] * powi(s[2], 2);
        c2[5] = 2 * c[4] * s[0] * t[0] + c[5] * (s[0] * t[1] + s[1] * t[0]) + c[6] * (s[0] * t[2] + s[2] * t[0]) +
                2 * c[7] * s[1] * t[1] + c[8] * (s[1] * t[2] + s[2] * t[1]) + 2 * c[9] * s[2] * t[2];
        c2[6] = 2 * c[4] * n[0] * s[0] + c[5] * (n[0] * s[1] + n[1] * s[0]) + c[6] * (n[0] * s[2] + n[2] * s[0]) +
                2 * c[7] * n[1] * s[1] + c[8] * (n[1] * s[2] + n[2] * s[1]) + 2 * c[9] * n[2] * s[2];
        c2[7] = c[4] * powi(t[0], 2) + c[5] * t[0] * t[1] + c[6] * t[0] * t[2] + c[7] * powi(t[1], 2) +
                c[8] * t[1] * t[2] +
                c[9] * powi(t[2], 2);
        c2[8] = 2 * c[4] * n[0] * t[0] + c[5] * (n[0] * t[1] + n[1] * t[0]) + c[6] * (n[0] * t[2] + n[2] * t[0]) +
                2 * c[7] * n[1] * t[1] + c[8] * (n[1] * t[2] + n[2] * t[1]) + 2 * c[9] * n[2] * t[2];
        c2[9] = c[4] * powi(n[0], 2) + c[5] * n[0] * n[1] + c[6] * n[0] * n[2] + c[7] * powi(n[1], 2) +
                c[8] * n[1] * n[2] +
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
                     c[13] * (2 * n[0] * s[1] * t[1] + 2 * n[1] * s[0] * t[1] + 2 * n[1] * s[1] * t[0]) +
                     c[14] * (n[0] * s[1] * t[2] + n[0] * s[2] * t[1] + n[1] * s[0] * t[2] + n[1] * s[2] * t[0] +
                              n[2] * s[0] * t[1] + n[2] * s[1] * t[0]) +
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

    static void ProjectPointsToSurface(const Scene *scene, const Point3f &refPoint, const Vector3f &refDir,
                                       ScatterSamplingRecord &rec, const Eigen::VectorXf &polyCoefficients,
                                       size_t polyOrder, bool useLocalDir, Float scaleFactor, Float kernelEps,
                                       SurfaceInteraction *res);

    static Normal3f AdjustRayDirForPolynomialTracing(Vector3f &inDir, const SurfaceInteraction &isect,
                                                     int polyOrder, Float polyScaleFactor, int channel);


    //FIXME: Just for debug, brute force method
    static std::tuple<std::vector<Point3f>, std::vector<Normal3f>, std::vector<Float>>
    GetConstraints(const Point3f &p, Float kernelEps,
                   const std::vector<Point3f> &points, const std::vector<Normal3f> normals) {
        std::vector<Point3f> pointConstraints;
        std::vector<Normal3f> normalConstraints;
        std::vector<Float> weightConstraints;
        Float dist2Threshold = 9.0f * kernelEps; // FIXME: Adjust a better threshold.
        for (size_t i = 0; i < points.size(); i++) {
            Float distSq = pbrt::DistanceSquared(p, points[i]);
            if (distSq < dist2Threshold) {
                pointConstraints.push_back(points[i]);
                normalConstraints.push_back(normals[i]);
                weightConstraints.push_back(1.0f);
            }
        }

        for (size_t i = 0; i < 32; i++) {
            pointConstraints.push_back(points[i]);
            normalConstraints.push_back(normals[i]);
            weightConstraints.push_back(-1.0f);
        }

        return std::make_tuple(pointConstraints, normalConstraints, weightConstraints);
    }

protected:
    static inline Float GaussianKernel(Float dist2, Float sigma2) {
        return std::exp(-dist2 / (2 * sigma2));
    }

    static Vector3f
    EvaluateGradient(const Point3f &pos, const Eigen::VectorXf &coeffs, const ScatterSamplingRecord &rec,
                     size_t degree, Float scaleFactor, bool useLocalDir, const Vector3f &refDir);


};

}

#endif //PBRT_V3_POLYNOMIALS_H
