#include "octree.h"
#include "twopass_dipole.h"

namespace pbrt {

Spectrum IrradianceOctree::Search(const Point3f &p, const std::shared_ptr<OctreeNode> &node, const Bounds3f &aabb,
                                  const TwoPassBSSRDF *bssrdf) {
    // FIXME: Bugs must be here!
    bool contained = Distance(p, aabb) == Float(0);
    Float approxSolidAngle = node->avgData.area / (p - node->avgData.pos).LengthSquared();
    if (!contained && approxSolidAngle < solidAngleThreshold) {
        // return bssrdf->Mo(p, node->avgData.E) * node->avgData.area;
        // Abandom the region too far away.
        return Spectrum(0.0f);
    } else {
        Spectrum E(0.0f);
        if (node->IsLeaf()) {
            for (auto &rNode: node->relatedNodes) {
                E += bssrdf->Mo(p, rNode.E) * rNode.area;
            }
        } else {
            // Inner node.
            for (int i = 0; i < 8; i++) {
                if (node->children[i] == nullptr) {
                    continue;
                }
                E += Search(p, node->children[i], OctreeNode::GetSubspace(aabb, i), bssrdf);
            }
        }
        return E;
    }
}

void IrradianceOctree::Propagate(const std::shared_ptr<OctreeNode> &node) {
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

}
