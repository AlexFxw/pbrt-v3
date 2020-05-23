//
// Created by 范軒瑋 on 2020/5/23.
//

#ifndef PBRT_V3_TWOPASS_DIPOLE_H
#define PBRT_V3_TWOPASS_DIPOLE_H

#include "pbrt.h"

namespace pbrt
{

class TwoPassHelper {
public:
    int Prepare(const std::vector<std::shared_ptr<Shape>> &shapes);

};

} // namespace pbrt

#endif //PBRT_V3_TWOPASS_DIPOLE_H
