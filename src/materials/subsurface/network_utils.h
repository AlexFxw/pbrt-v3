/*
 * @Author: Hsuan-Wei Fan
 * @Date: 2020-04-23 22:41:27
 * @LastEditors: Hsuan-Wei Fan
 * @LastEditTime: 2020-04-23 23:17:11
 * @Description: 
 */

#ifndef PBRT_MATERIALS_SUBSURFACE_NETWORKUTILS_H
#define PBRT_MATERIALS_SUBSURFACE_NETWORKUTILS_H

#include "pbrt.h"
#include <Eigen/Core>
#include <Eigen/Dense>

namespace pbrt {

class NetworkUtils
{
  public:
    static Eigen::VectorXf LoadVectorDynamic(const std::string &filename);
};

}  // namespace pbrt

#endif // PBRT_MATERIALS_SUBSURFACE_NETWORKUTILS_H