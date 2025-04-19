#pragma once

#include "Config.h"
#include <Eigen/Eigen>



std::pair<Eigen::MatrixX3<msf::Real>, Eigen::MatrixX3i> isosurf_mesh(
	int reso, const std::pair<Eigen::Vector3<msf::Real>, Eigen::Vector3<msf::Real>>& bbox,
	std::function<msf::Real(msf::Real, msf::Real, msf::Real)> sdf, const double isolevel);
// |sample| = (reso+1)^3
std::pair<Eigen::MatrixX3<msf::Real>, Eigen::MatrixX3i> isosurf_mesh(
	int reso, const Eigen::MatrixX3<msf::Real>& sampleLocation,
	const Eigen::VectorX<msf::Real>& samplevalue);
// |sample| = (reso+1)^3
