#include "PeriodicMesher.h"
#include <fstream>

std::string getPath(std::string str);

void msf::PeriodSurfaceMesh::meanCurvatureFlow(
	Eigen::MatrixX3<Real>& vlist, const std::vector<OM::SmartVertexHandle>& vcut, const std::vector<OM::SmartHalfedgeHandle>& hetrans, const std::vector<HETag>& hetags)
{
	for (int iter = 0; iter < 200; iter++) {
		Eigen::MatrixX3<Real> Hv(vlist.rows(), 3);
		Hv.setZero();
		for (auto vh : vertices()) {
			auto [o, ring] = find1ring(vlist, vh, hetags);
			Eigen::Vector3<Real> hn;
			meanH(o, ring, hn);
			Hv.row(vh.idx()) = hn.transpose();
		}
		vlist += Hv;
		{std::ofstream ofs(getPath("Hv")); ofs << Hv; }
		savePeriodicMesh(getPath("meshiter.obj"), vlist, vcut, hetags);
	}
}

