#include "utils.h"
#include "shape.h"
#include "shapes/triangle.h"
#include "sampling.h"

namespace pbrt {

std::shared_ptr<TriangleMesh> Utils::PreprocessTriangles(const std::vector<std::shared_ptr<Shape>> &shapes) {
    // "shapes" is the collection of triangles.
    // Initialize triangle mesh's area distribution here.
    Float areaSum = 0.0f;
    size_t shapesNum = shapes.size();
    DCHECK_GT(shapesNum, 0);
    Float areas[shapesNum];
    for (size_t i = 0; i < shapesNum; i++) {
        areas[i] = shapes[i]->Area();
        areaSum += areas[i];
    }
    auto tmp = shapes.back();
    std::shared_ptr<Triangle> triangleIter(std::dynamic_pointer_cast<Triangle>(shapes.back()));
    std::shared_ptr<TriangleMesh> triMesh(triangleIter->mesh);
    triMesh->areaDistri = new Distribution1D(areas, shapesNum);
    triMesh->area = areaSum;
    triMesh->invArea = 1.0f / areaSum;
    return triMesh;
}

}