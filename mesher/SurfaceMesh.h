#pragma once
#include "geometry/surface/mesh.h"
#include <string>
#include "Flag.h"

BEGIN_MINSURF_NAMESPACE

class SurfaceMesh
	: public MeshSurface
{
	std::vector<VertexFlag> _vflags;

public:
	auto& getVFlags(void) { return _vflags; }
	VertexFlag& getVFlag(OM::SmartVertexHandle vh) { return _vflags[vh.idx()]; }

	int vDofid(OM::VertexHandle vh) const { return vh.idx(); }

	void faceVdof(OM::FaceHandle fh, int gid[3]) const {
		int counter = 0;
		for (auto fv : fv_range(fh)) { gid[counter++] = vDofid(fv); }
	}
	int n_vdof(void) const { return n_vertices(); }

	void saveVdofVector(std::string filename, const Eigen::VectorXd& u, Eigen::Vector3i beg, Eigen::Vector3i stride);
};

END_MINSURF_NAMESPACE
