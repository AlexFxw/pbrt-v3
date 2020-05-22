/*
 * @Author: Hsuan-Wei Fan
 * @Date: 2020-05-22 11:25:18
 * @LastEditors: Hsuan-Wei Fan
 * @LastEditTime: 2020-05-22 17:01:19
 * @Description: The kdtree for computing shape descriptor
 */

#ifndef PBRT_SUBSURFACE_KDTREE_H
#define PBRT_SUBSURFACE_KDTREE_H

#include "geometry.h"
#include "pbrt.h"

namespace pbrt {


struct KDNode {
    enum class Type {
        Leaf, Axis, Invalid
    };
    Point3f pos, avgP;
    Normal3f norm, avgN;
    int rightIndex, self;
    int sampleCount, axis;
    Type type;

    KDNode() : type(Type::Invalid), axis(0) {}

    KDNode(const Point3f &pos, const Point3f &avgP, const Normal3f &norm, const Normal3f &avgN) :
            pos(pos), avgP(avgP), norm(norm), avgN(avgN), rightIndex(0), sampleCount(1), axis(0), type(Type::Invalid) {}

    bool IsLeaf() const { return type == Type::Leaf; }

    int LeftIndex() const { return self + 1; }

};

inline std::ostream &operator<<(std::ostream &output, const KDNode &node) {
    output << "[pos: " << node.pos;
    output << ", right: " << node.rightIndex;
    std::string type = "invalid";
    if (node.type == KDNode::Type::Leaf) {
        type = "Leaf";
    } else if (node.type == KDNode::Type::Axis) {
        type = "Axis";
    }
    output << ",type: " << type << " ]";
    return output;
}

class ConstraintKDTree {
public:
    explicit ConstraintKDTree(size_t pointNum = 0) : mNodes(pointNum), maxDepth(0) {
        aabb = Bounds3f(Point3f(), Point3f());
    }

    void Build(const std::vector<Point3f> &positions, const std::vector<Normal3f> &normals) {
        DCHECK_EQ(positions.size(), normals.size());
        aabb = Bounds3f(Point3f(), Point3f());
        for (int i = 0; i < positions.size(); i++) {
            aabb = Union(aabb, positions[i]);
            mNodes[i] = KDNode(positions[i], Point3f(), normals[i], Normal3f());
            mNodes[i].self = i;
        }

        Build(1, mNodes.begin(), mNodes.begin(), mNodes.end());

        for (int i = 0; i < positions.size(); i++) {
            mNodes[i].self = i;
        }

        for (size_t i = 0; i < std::min(size_t(32), positions.size()); i++) {
            globalPoints.push_back(positions[i]);
            globalNormals.push_back(normals[i]);
        }
    }

    std::tuple<std::vector<Point3f>, std::vector<Normal3f>, std::vector<Float>>
    GetConstraints(const Point3f &p, Float kernelEps) const {
        std::vector<Point3f> positions;
        std::vector<Normal3f> normals;
        std::vector<Float> weights;

        GetConstraints(p, mNodes[0], aabb, positions, normals, weights, kernelEps);

        auto sP = positions.size();
        for (size_t i = 0; i < globalPoints.size(); i++) {
            positions.push_back(globalPoints[i]);
            normals.push_back(globalNormals[i]);
            weights.push_back(-1.0f);
        }

        return std::make_tuple(positions, normals, weights);
    }

    void Traverse(int curIndex) const {
        const KDNode &node = mNodes[curIndex];
        std::cout << node.self << " " << curIndex << " " << node << std::endl;
        if (node.IsLeaf())
            return;
        if (node.rightIndex != 0) {
            std::cout << "--Right--" << std::endl;
            Traverse(node.rightIndex);
        }
        std::cout << "--Left--" << std::endl;
        Traverse(node.LeftIndex());
    }

private:
    std::vector<KDNode> mNodes;
    int maxDepth;
    std::vector<Point3f> globalPoints;
    std::vector<Normal3f> globalNormals;

    Bounds3f aabb;

    void Build(int depth, std::vector<KDNode>::iterator base,
               std::vector<KDNode>::iterator start, std::vector<KDNode>::iterator end) {
        maxDepth = std::max(maxDepth, depth);
        int count = (int) (end - start);
        CHECK_GT(count, int(0));
        if (count == 1) {
            // Is leaf node.
            start->type = KDNode::Type::Leaf;
            return;
        }

        int axis = 0;
        std::vector<KDNode>::iterator splitNode;
        splitNode = start + count / 2;
        axis = aabb.LargestAxis();
        std::nth_element(start, splitNode, end, [&axis](KDNode &n1, KDNode &n2) {
            return n1.pos[axis] < n2.pos[axis];
        });

        splitNode->type = KDNode::Type::Axis;
        splitNode->axis = axis;

        if (splitNode + 1 != end) {
            // Not the last one of the nodes.
            splitNode->rightIndex = (int) (splitNode + 1 - base);
        } else {
            splitNode->rightIndex = 0;
        }

        // Put the split point (node) to the start of the array.
        int rightIndex = splitNode->rightIndex;
        std::iter_swap(start, splitNode);

        Float tmp = aabb.pMax[axis], splitPos = splitNode->pos[axis];
        aabb.pMax[axis] = splitPos;
        Build(depth + 1, base, start + 1, splitNode + 1); // Left subtree
        aabb.pMax[axis] = tmp;

        if (rightIndex != 0) {
            // Has right subtree.
            tmp = aabb.pMin[axis];
            aabb.pMin[axis] = splitPos;
            Build(depth + 1, base, splitNode + 1, end);
            aabb.pMin[axis] = tmp;
        }
    }

    void GetConstraints(const Point3f &p, KDNode node, const Bounds3f &aabb,
                        std::vector<Point3f> &points, std::vector<Normal3f> &normals,
                        std::vector<Float> &weights, Float kernelEps) const {
        Float dist2Threshold = 9.0 * kernelEps;

        auto d = DistanceSquared(p, aabb);

        if (DistanceSquared(p, aabb) > dist2Threshold) {
            return;
        }

        if (DistanceSquared(p, node.pos) < dist2Threshold) {
            points.push_back(node.pos);
            normals.push_back(node.norm);
            weights.push_back(1.0f);
        }

        if (node.IsLeaf()) {
            return;
        }

        int axis = node.axis;
        Point3f leftMin = aabb.pMin;
        Point3f leftMax = aabb.pMax;
        leftMax[axis] = node.pos[axis];
        Point3f rightMin = aabb.pMin;
        Point3f rightMax = aabb.pMax;
        rightMin[axis] = node.pos[axis];

        Bounds3f leftBounds(leftMin, leftMax);
        Bounds3f rightBounds(rightMin, rightMax);

        int rightIndex = node.rightIndex;
        if (rightIndex != 0) {
            GetConstraints(p, mNodes[rightIndex], rightBounds, points, normals, weights, kernelEps);
        }

        int leftIndex = node.self + 1;
        GetConstraints(p, mNodes[leftIndex], leftBounds, points, normals, weights, kernelEps);
    }
};

} // namespace pbrt
#endif  // PBRT_SUBSURFACE_KDTREE_H