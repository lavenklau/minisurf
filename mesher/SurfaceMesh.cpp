#include "SurfaceMesh.h"
#include <fstream>

void msf::SurfaceMesh::saveVdofVector(std::string filename, const Eigen::VectorXd& u, Eigen::Vector3i beg, Eigen::Vector3i stride)
{
	std::ofstream ofs(filename, std::ios::binary);
	std::vector<Eigen::Vector3d> vdata(n_vertices());
	for (auto vh : vertices()) {
		int vdofid = vDofid(vh);
		Eigen::Vector3d u_v;
		for (int axis = 0; axis < 3; axis++) {
			u_v[axis] = u[beg[axis] + stride[axis] * vdofid];
		}
		vdata[vh.idx()] = u_v;
	}
	ofs.write((const char*)vdata.data(), sizeof(Eigen::Vector3d) * n_vertices());
	ofs.close();
}

