#include "PeriodicMesher.h"


using namespace msf;

OM::Vec3d close_point(const OM::Vec3d& p0, const OM::Vec3d& p1, const OM::Vec3d& p) {
	auto v = p1 - p0;
	OM::Vec3d tL = p0 + v * (p - p0).dot(v) / v.sqrnorm();
	return tL;
}

auto cross_point(const OM::Vec3d& p0, const OM::Vec3d& p1, const OM::Vec3d& pL, const OM::Vec3d& pR) {
	auto cL = (p0 + p1 + pL) / 3;
	auto cR = (p0 + p1 + pR) / 3;
	auto v = p1 - p0;
	OM::Vec3d tL = p0 + v * (cL - p0).dot(v) / v.sqrnorm();
	OM::Vec3d tR = p0 + v * (cR - p0).dot(v) / v.sqrnorm();
	double tc = (tL - pL).norm() / (tR - pR).norm() + 1;
	OM::Vec3d itc = tR * tc + tL * (1 - tc);
	return std::make_tuple(itc, tc);
}

Eigen::MatrixX3d msf::PeriodSurfaceMesh::second_fundamental_form(int type)
{
	Eigen::MatrixX3d bform(3 * n_faces(), 3);
	auto vnlist = getVertexNormal(2, 0.7);
	auto fflist = getFaceFrame(2, 0.7);
	if (type == 0) {
		Eigen::SparseMatrix<Real> A(n_faces() * 3, n_faces() * 3);
		std::vector<Eigen::Triplet<Real>> Awise;
		for (auto eh : edges()) {
			auto h0 = eh.halfedge(0);
			auto d0 = h0.from();
			auto d1 = h0.to();
			auto lh = h0.next().to();
			auto rh = h0.opp().next().to();
			auto fL = h0.face();
			auto fR = h0.opp().face();
			auto p0 = point(d0), p1 = point(d1);
			auto pL = point(lh), pR = point(rh);
			p1 = sync_period(p1, p0); pL = sync_period(pL, p0); pR = sync_period(pR, p0);
			auto [itc, t] = cross_point(p0, p1, pL, pR);
			Eigen::Vector3d nL = fflist.row(fL.idx() * 3 + 2), nR = fflist.row(fR.idx() * 3 + 2);
			Eigen::Matrix<Real, 3, 2> gL, gR;
			gL << fflist.row(fL.idx() * 3).transpose(), fflist.row(fL.idx() * 3 + 1).transpose();
			gR << fflist.row(fR.idx() * 3).transpose(), fflist.row(fR.idx() * 3 + 1).transpose();
			Eigen::Vector2<Real> vL = gL.transpose() * toEigen(itc - pL);
			Eigen::Vector2<Real> vR = gL.transpose() * toEigen(pR - itc);
			//Awise.emplace_back(eh.idx())
		}
	}
	else if (type == 1) {
	
	}
	else {
	
	}
	return bform;
}

