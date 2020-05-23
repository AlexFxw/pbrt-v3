#ifndef PBRT_V3_UTILS_H
#define PBRT_V3_UTILS_H

#include "pbrt.h"
#include "shapes/triangle.h"

namespace pbrt {

class Utils {
public:
    static std::shared_ptr<TriangleMesh>
    PreprocessTriangles(const std::vector<std::shared_ptr<Shape>> &shapes);
};

} // namespace pbrt
#endif //PBRT_V3_UTILS_H
