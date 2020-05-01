//
// Created by 范軒瑋 on 2020-04-28.
//

#ifndef PBRT_V3_POLYNOMIALS_H
#define PBRT_V3_POLYNOMIALS_H

#include "pbrt.h"
#include "geometry.h"
#include <Eigen/Core>

namespace pbrt {

template<typename PointType, typename DataRecord>
struct SimpleKDNode {
    typedef uint32_t IndexType;
    enum {
        ELeafFlag = 0x10;
        EAxisMask = 0x0F;
    };

    PointType pos;
    IndexType right;
    DataRecord data;
    uint8_t flags;

    PointType &GetPos() const { return pos; }

    SimpleKDNode() : pos(PointType()), right(0), data(), flags(0) {}

    DataRecord &GetData() const { return data; }

    SimpleKDNode(const DataRecord &data) : pos(PointType()), right(0), data(data), flags(0) {}

    IndexType GetRightIndex(IndexType right) const { return right; }

    IndexType GetLeftIndex(IndexType self) const { return self + 1; }

    inline bool IsLeaf() const { return flags & (uint8_t) ELeafFlag; }

    inline void SetPosition(const Point3f &pos) { this->pos = pos; }

    inline void SetData(const DataRecord &data) { this->data = data; }

    void SetLeaf(bool value) {
        if (value)
            flags |= (uint8_t) ELeafFlag;
        else
            flags &= (uint8_t)
                    ~ELeafFlag;
    }


    boid SetRightIndex(IndexType self, IndexType value) { right = value; }

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
    typedef NodeType::IndexType IndexType;

    struct SearchResult {
        Float distSquared;
        IndexType index;

        SearchResult() {}

        SearchResult(Float distSquared, IndexType index) : distSquared(distSquared), index(index) {}

        bool operator==(const SearchResult &r) const {
            return distSquared == r.distSquared && index == r.index;
        }
    };

    PointKDTree(size_t nodes = 0) : mNodes(nodes), mDepth(0) {}

    void Clear() {
        mNodes.clear();
        mAABB = Bounds3f();
    }

    inline void Resize(size_t size) { mNodes.resize(size); }

    inline bool HasRightChild(IndexType index) {
        return mNodes[index].GetRightIndex(index) != 0;
    }

    void PushBack(const NodeType &node) {
        mNodes.push_back(node);
        mAABB = pbrt::Union(mAABB, node.pos);
    }

    inline size_t Size() const { return mNodes.size(); }

    inline NodeType &operator[](size_t idx) { return m_nodes[idx]; }

    void Build(bool recomputeAABB = false) {
        if (recomputeAABB) {
            mAABB = Bounds3f();
            for (auto iter = mNodes.begin(); iter != mNodes.end(); iter++) {
                mAABB = pbrt::Union(mAABB, iter->pos);
            }
        }
        std::vector<IndexType> indirection(mNodes.size());
        for (size_t = 0; i < mNodes.size(); i++) {
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
        if (split + 1 != rangeEnd) {
            splitNode.SetRightIndex((IndexType) (rangeStart - base), (IndexType) (split + 1 - base));
        } else {
            splitNode.SetRightIndex((IndexType) (rangeStart - base), 0);
        }
        // TODO: Set left

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
            return m_nodes[i1].getPosition()[m_axis] < m_nodes[i2].getPosition()[m_axis];
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
    typedef PointKDTree<TreeNode> CTree;
    void Build(const std::vector<Point3f> &sampledP, const std::vector<Noraml3f> &sampledN);
    void GetConstraints(const Point3f &p, TreeNode node, TreeNode::IndexType index, const Bounds3f &aabb,
                        std::vector<Point3f> &points, std::vector<Normal3f> &normals, std::vector<Float> &sampleWeights,
                        Float kernelEps, const std::function<Float(Float, Float)> &kernelFunc) const;
    std::tuple<Point3f, Normal3f, size_t> CalcAvgValues(TreeNode &node, TreeNode::IndexType index);

private:
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
            return s * v.x + t * v.y + n * v.z;
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

    static PBRT_CONSTEXPR size_t powerToIndex(size_t dx, size_t dy, size_t dz) {
        // Converts a polynomial degree to a linear coefficient index
        auto d = dx + dy + dz;
        auto i = d - dx;
        auto j = d - dx - dy;
        return i * (i + 1) / 2 + j + d * (d + 1) * (d + 2) / 6;
    }

    static Float GetKernelEps(const MediumParameters &mediumParas, int channel, Float kernelMultiplier);

    static inline float GetFitScaleFactor(float kernelEps) { return 1.0f / std::sqrt(kernelEps); }

    static PBRT_CONSTEXPR int nChooseK(int n, int k) {
        return (k == 0 || n == k) ? 1 : nChooseK(n - 1, k - 1) + nChooseK(n - 1, k);
    }

    static PBRT_CONSTEXPR int nPolyCoeffs(int polyOrder, bool hardSurfaceConstraint = false) {
        return hardSurfaceConstraint ? nChooseK(3 + polyOrder, polyOrder) - 1 : nChooseK(3 + polyOrder, polyOrder);
    }

    static std::tuple<Polynomial, std::vector<Point3f>, std::vector<Vector3f>>
    FitPolynomial(const PolyFitRecord &polyFitRecord, const ConstraintKDTree *kdTree);

    template<size_t polyOrder, bool hardSurfaceConstraint = true>
    static std::tuple<Polynomial, std::vector<Point3f>, std::vector<Vector3f>>
    FitPolynomialImpl(const PolyFitRecord &pfRec, const ConstraintKDTree *kdTree);

    static void OnbDuff(const Normal3f &n, Vector3f &b1, Vector3f &b2) {
        float sign = copysignf(1.0f, n.z);
        const float a = -1.0f / (sign + n.z);
        const float b = n.x * n.y * a;
        b1 = Vector3f(1.0f + sign * n.x * n.x * a, sign * b, -sign * n.x);
        b2 = Vector3f(b, sign + n.y * n.y * a, -n.y);
    }

    template<int order, typename T>
    static Eigen::Matrix<float, nPolyCoeffs(order), 1> RotatePolynomialEigen(const T &c,
                                                                             const Vector3f &s, const Vector3f &t,
                                                                             const Normal3f &n);
    static void ProjectPointsToSurface(const Scene *scene,
                                       const Point3f &refPoint, const Vector3f &refDir,
                                       ScatterSamplingRecord &rec,
                                       const Eigen::VectorXf &polyCoefficients,
                                       size_t polyOrder, bool useLocalDir, Float scaleFactor, Float kernelEps);

protected:
    static inline Float GaussianKernel(Float dist2, Float sigma2) {
        return std::exp(-dist2 / (2 * sigma2));
    }

    static Vector3f EvaluateGradient(const Point3f &pos, const Eigen::VectorXf &coeffs, const ScatterSamplingRecord &rec,
                                     size_t degree, Float scaleFactor, bool useLocalDir, const Vector3f &refDir);


};

}

#endif //PBRT_V3_POLYNOMIALS_H
