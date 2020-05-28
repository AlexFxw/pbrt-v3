/*
 * The acceleration structure used in two-pass dipole rendering.
 */

#ifndef PBRT_V3_OCTREE_H
#define PBRT_V3_OCTREE_H

#include "core/spectrum.h"
#include "pbrt.h"
#include "geometry.h"
#include "sampling.h"

namespace pbrt {

// Used for
struct IrradianceData {
    Point3f pos;
    Spectrum E;
    Float area;

    IrradianceData() : pos(Point3f()), E(Spectrum(0.0f)), area(0.0f) {}

    IrradianceData(const Point3f &pos, const Spectrum &E, Float area) : pos(pos), E(E), area(area) {}

    const Point3f &GetPos() const { return pos; }
};

struct OctreeNode {
    const static int threshold = 8;

    enum class Type {
        Axis, Leaf, Invalid
    };

    OctreeNode() : type(Type::Invalid) {}

    Type type;
    IrradianceData avgData;
    std::shared_ptr<OctreeNode> children[8];
    std::vector<IrradianceData> relatedNodes;

    void Build(int depth, std::vector<IrradianceData> &points, const Bounds3f &bound) {
        if (points.size() <= threshold) {
            type = Type::Leaf;
            relatedNodes = std::move(points);
            return;
        }

        type = Type::Axis;

        std::vector<IrradianceData> subspacePoints[8];

        // FIXME: Rewrite in IrradianceData format.
        for (IrradianceData &data: points) {
            const Point3f &p = data.pos;
            for (int i = 0; i < 8; i++) {
                Bounds3f subBound = GetSubspace(bound, i);
                Float dist = pbrt::DistanceSquared(p, subBound);
                if (dist == 0.0f) {
                    // In the subspace.
                    subspacePoints[i].push_back(std::move(data));
                    break;
                }
            }
        }

        for (int i = 0; i < 8; i++) {
            if (!subspacePoints[i].empty()) {
                Bounds3f subBound = GetSubspace(bound, i);
                children[i] = std::make_shared<OctreeNode>();
                children[i]->Build(depth + 1, subspacePoints[i], subBound);
            } else {
                children[i] = nullptr;
            }
        }
    }

    static Bounds3f GetSubspace(const Bounds3f &bound, int index) {
        Point3f middle = bound.pMin + (bound.Diagonal() / 2.0f);
        Point3f pMin = bound.pMin, pMax = bound.pMax;
        if (index == 0) {
            pMax = Point3f(middle.x, middle.y, middle.z);
        } else if (index == 1) {
            pMin.x = middle.x;
            pMax.y = middle.y;
            pMax.z = middle.z;
        } else if (index == 2) {
            pMin.y = middle.y;
            pMax.x = middle.x;
            pMax.z = middle.z;
        } else if (index == 3) {
            pMax.z = middle.z;
            pMin.y = middle.y;
            pMin.x = middle.x;
        } else if (index == 4) {
            pMin.y = middle.y;
            pMax.x = middle.x;
            pMax.z = middle.z;
        } else if (index == 5) {
            pMax.y = middle.y;
            pMin.x = middle.x;
            pMin.z = middle.z;
        } else if (index == 6) {
            pMin.y = middle.y;
            pMin.z = middle.z;
            pMax.x = middle.x;
        } else if (index == 7) {
            pMin = Point3f(middle.x, middle.y, middle.z);
        }
        Bounds3f res;
        res.pMin = pMin;
        res.pMax = pMax;
        return res;
    }


    bool IsLeaf() {
        return type == Type::Leaf;
    }
};

class IrradianceOctree {
public:

    IrradianceOctree() : nSamples(0) {}

    void Build(std::vector<IrradianceData> &samplePoints) {
        aabb = Bounds3f((Point3f()), Point3f());
        for (const IrradianceData &data: samplePoints) {
            const Point3f &p = data.GetPos();
            aabb = Union(aabb, p);
        }
        root = std::make_shared<OctreeNode>();
        root->Build(0, samplePoints, aabb);
        Propagate(root);
    }

    Spectrum Search(const Point3f &p) {
        return Search(p, root, aabb);
    }

private:
    Spectrum Search(const Point3f &p, const std::shared_ptr<OctreeNode> &node, const Bounds3f &aabb) {
        // FIXME: Bugs must be here!
        bool contained = Distance(p, aabb) == 0.0f;
        Float approxSolidAngle = node->avgData.area / (p - node->avgData.pos).LengthSquared();
        if (!contained && approxSolidAngle < solidAngleThreshold) {
            return node->avgData.E * node->avgData.area;
        } else {
            Spectrum E(0.0f);
            if (node->IsLeaf()) {
                for(auto &rNode: node->relatedNodes) {
                    E += rNode.E * rNode.area;
                }
            } else {
                for (int i = 0; i < 8; i++) {
                    if (node->children[i] == nullptr) {
                        continue;
                    }
                    E += Search(p, node->children[i], OctreeNode::GetSubspace(aabb, i));
                }
            }
            return E;
        }
    }

    void Propagate(const std::shared_ptr<OctreeNode> &node) {
        IrradianceData &clusterData = node->avgData;
        clusterData = IrradianceData();
        Float weightSum = 0.0f;
        if (node->IsLeaf()) {
            // Leaf node
            for (IrradianceData &sample: node->relatedNodes) {
                clusterData.E += sample.E * sample.area;
                clusterData.area += sample.area;
                // FIXME: Use luminance?
                Float weight = sample.E.y() * sample.area;
                clusterData.pos += sample.pos * weight;
                weightSum += weight;
            }
        } else {
            // Inner node.
            for (int i = 0; i < 8; i++) {
                if (node->children[i] != nullptr) {
                    Propagate(node->children[i]);
                    IrradianceData &childAvg = node->children[i]->avgData;
                    clusterData.E += childAvg.E * childAvg.area;
                    clusterData.area += childAvg.area;
                    Float weight = childAvg.E.y() * childAvg.area;
                    clusterData.pos += childAvg.pos * weight;
                    weightSum += weight;
                }
            }
        }

        if (clusterData.area != 0) {
            clusterData.E /= clusterData.area;
        }

        if (weightSum != 0) {
            clusterData.pos /= weightSum;
        }
    }

    bool Criterion(const Point3f &p, const std::shared_ptr<OctreeNode> &node) {
        Float distSq = (p - node->avgData.pos).LengthSquared();
        if (distSq == 0.0f) {
            return false;
        }
        Float approxSolidAngle = node->avgData.area / distSq;
        return approxSolidAngle < solidAngleThreshold;
    }

    std::shared_ptr<OctreeNode> root;
    int nSamples;
    Bounds3f aabb;
    static constexpr Float solidAngleThreshold = 0.01f;
};


} // namespace pbrt

#endif //PBRT_V3_OCTREE_H
