/*
 * The acceleration structure used in two-pass dipole rendering.
 */

#ifndef PBRT_V3_OCTREE_H
#define PBRT_V3_OCTREE_H

#include "pbrt.h"
#include "geometry.h"
#include "utils.h"
#include "interaction.h"
#include "sampling.h"

namespace pbrt {

// Used for

template<typename ExtraData>
struct OctreeNode {
    const static int threshold = 8;

    enum class Type {
        Axis, Leaf, Invalid
    };

    OctreeNode() : type(Type::Invalid) {}

    Type type;
    ExtraData avgData;
    Bounds3f mBound;
    int mDepth;
    std::unique_ptr<OctreeNode> children[8];
    std::vector<Point3f> relatedNodes;

    void Build(int depth, std::vector<Point3f> &points, const Bounds3f &bound) {
        mDepth = depth;
        mBound = bound;
        if (points.size() <= threshold) {
            type = Type::Leaf;
            relatedNodes = std::move(points);
            printf("Depth and size: %d/%d\n", depth, relatedNodes.size());
            for(auto &p: relatedNodes) {
                std::cout << "  " << p << std::endl;
            }
            return;
        }

        type = Type::Axis;

        std::vector<Point3f> subspacePoints[8];

        for (Point3f &p: points) {
            for (int i = 0; i < 8; i++) {
                Bounds3f subBound = GetSubspace(bound, i);
                Float dist = pbrt::DistanceSquared(p, subBound);
                if (dist == 0.0f) {
                    // In the subspace.
                    subspacePoints[i].push_back(std::move(p));
                    break;
                }
            }
        }

        for (int i = 0; i < 8; i++) {
            if(subspacePoints[i].size() > 0) {
                Bounds3f subBound = GetSubspace(bound, i);
                children[i] = std::unique_ptr<OctreeNode>(new OctreeNode());
                children[i]->Build(depth + 1, subspacePoints[i], subBound);
            }
        }
    }

    Bounds3f GetSubspace(const Bounds3f &bound, int index) {
        Point3f middle = bound.pMin + (bound.Diagonal() / 2.0f);
        Point3f pMin = bound.pMin, pMax = bound.pMax;
        if(index == 0) {
            pMax = Point3f(middle.x, middle.y, middle.z);
        } else if(index == 1) {
            pMin.x = middle.x;
            pMax.y = middle.y;
            pMax.z = middle.z;
        } else if(index == 2) {
            pMin.y = middle.y;
            pMax.x = middle.x;
            pMax.z = middle.z;
        } else if(index == 3) {
            pMax.z = middle.z;
            pMin.y = middle.y;
            pMin.x = middle.x;
        } else if(index == 4) {
            pMin.y = middle.y;
            pMax.x = middle.x;
            pMax.z = middle.z;
        } else if(index == 5) {
            pMax.y = middle.y;
            pMin.x = middle.x;
            pMin.z = middle.z;
        } else if(index == 6) {
            pMin.y = middle.y;
            pMin.z = middle.z;
            pMax.x = middle.x;
        } else if(index == 7) {
            pMin = Point3f(middle.x, middle.y, middle.z);
        }
        Bounds3f res;
        res.pMin = pMin;
        res.pMax = pMax;
        return res;
    }
};

template<typename ExtraData>
class Octree {
public:
    explicit Octree(const std::vector<std::shared_ptr<Shape>> &shapes, int nSamples) : nSamples(nSamples) {
        // TODO:
        std::shared_ptr<TriangleMesh> triMesh(Utils::PreprocessTriangles(shapes));
        DCHECK_EQ(triMesh->nTriangles, shapes.size());

        for (int i = 0; i < nSamples; i++) {
            Float pdf;
            size_t trigIdx = triMesh->areaDistri->SampleDiscrete(GetRandomFloat());
            Interaction tmp = shapes[trigIdx]->Sample(RandPoint2f(), &pdf);
            sampledP.push_back(tmp.p);
        }
    }

    Octree(const std::vector<Point3f> &sampleP) : sampledP(sampleP), nSamples(sampleP.size()) {

    }

    void Build() {
        Bounds3f bound((Point3f()), Point3f());
        for (const Point3f &p: sampledP) {
            bound = Union(bound, p);
        }
        root = std::unique_ptr<OctreeNode<ExtraData>>(
                new OctreeNode<ExtraData>()
        );

        root->Build(0, sampledP, bound);
    }

private:
    std::unique_ptr<OctreeNode<ExtraData>> root;
    std::vector<Point3f> sampledP;
    int nSamples;
};


} // namespace pbrt

#endif //PBRT_V3_OCTREE_H
