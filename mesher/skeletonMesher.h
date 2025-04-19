#pragma once

#include "geometry/surface/mesh.h"
#include "Eigen/Eigen"
#include <tuple>

using namespace msf;

struct SkeletonMesher
	: public MeshSurface
{
	std::vector<Eigen::Vector3<Real>> sk_vertex;
	std::vector<Eigen::Vector2i> sk_line;
	std::vector<Eigen::Vector3i> sk_tri;
	Eigen::VectorX<Real> grid_dist;
	Eigen::MatrixX<Real> grid_pos;
public:
	SkeletonMesher(const std::string& filename, int reso = 128, Real period = 2);
	bool read_skeleton(const std::string& filename);
	void buildSDF(int reso, Real period = 0);
	std::tuple<Eigen::MatrixX3<Real>, Eigen::MatrixX3i> extract_level_set(Real s, Real edge_len, bool closed);
	std::tuple<Eigen::MatrixX3<Real>, Eigen::MatrixX3i> extract_isovolume(Real s);
	void exportLelveSet(const std::vector<Real>& slist, const std::string& fileprefix);
};

