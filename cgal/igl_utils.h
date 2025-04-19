#pragma once

#include <Eigen/Eigen>
#include "Config.h"
#include <tuple>


std::tuple<Eigen::MatrixX3<msf::Real>, Eigen::MatrixX3i> readMesh(const std::string& filename);
