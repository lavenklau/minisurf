#pragma once

#include "Config.h"
#include <Eigen/Eigen>
#include <tuple>

namespace msf {
	void reorient_mesh(const Eigen::MatrixX3<Real>& V, Eigen::MatrixX3i& F);
};
