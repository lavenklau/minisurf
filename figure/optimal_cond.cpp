#include "mesher/asymptotic_analysis.h"
#include "Eigen/PardisoSupport"
#include "material/materail.h"
#include "fmt/core.h"
#include "boost/algorithm/string.hpp"

extern int max_iter;
extern double converge_tol;
extern double converge_kc;
extern double converge_sc;
extern double step_tol;
extern double precondition_strength;
extern bool asym_no_remesh;
extern bool asym_no_surgery;
extern double weight_area;
extern int delaunay_remesh_outer_iter;
extern int delaunay_remesh_inner_iter;

using namespace msf;


auto face_load(PeriodSurfaceMesh& m, const std::vector<Eigen::Matrix3d>& frlist, const std::vector<Compile1ring>& vrings, double lam0, double mu, const Eigen::Matrix3d& eps) {
	std::vector<Eigen::Vector3d> fload(m.n_faces());
	for (auto fh : m.faces()) {
		auto tri = m.getFacePeriodVertex(fh, 1);
		auto vring = retrieve_element_vertex_data(m, fh, vrings);
		auto fr = frlist[fh.idx()];
		auto bform = second_fundamental_form(tri, fr, vring.data());
		Eigen::Vector3d n = fr.col(2).normalized();
		Eigen::Matrix3d P = Eigen::Matrix3d::Identity() - n * n.transpose();
		double H = bform.trace() / 2;
		Eigen::Vector3d Ft = (2 * (lam0 + mu) * bform + 4 * mu * H * P) * eps * n;
		Eigen::Vector3d F3 = 2 * ((lam0 * H * P + mu * bform) * eps).trace() * n;
		fload[fh.idx()] = Ft + F3;
	}
	return fload;
}

auto distribute_vertex_load(PeriodSurfaceMesh& m, const std::vector<Eigen::Vector3d>& fload, const std::vector<Eigen::Matrix3d>& frlist) {
	std::vector<Eigen::Vector3d> vload(m.n_vertices());
	for (auto& v : vload) v.setZero();
	for (auto fh : m.faces()) {
		auto A = frlist[fh.idx()].col(2).norm();
		for (auto vh : fh.vertices()) {
			vload[vh.idx()] += A / 3 * fload[fh.idx()];
		}
	}
	return vload;
}

auto get_mass_matrix(const std::vector<Compile1ring>& vrings) {
	Eigen::VectorXd Av(vrings.size());
	for (int i = 0; i < vrings.size(); i++) {
		Av[i] = vrings[i].As;
	}
	return Av;
}

auto homogeneous_membrane_energy(PeriodSurfaceMesh& m, const std::vector<Eigen::Matrix3d>& frlist, const Eigen::Matrix3d& eps, double lam0, double mu) {
	using ADScalar = Eigen::AutoDiffScalar<Eigen::Vector<double, 9>>;
	double Em = 0;
	Eigen::VectorXd dEm(m.n_vertices() * 3, 1); dEm.setZero();
	for (auto fh : m.faces()) {
		auto tri = m.getFacePeriodVertex(fh, 1);
		auto fvh = m.getFaceVertexHandle(fh);
		Eigen::Matrix3<ADScalar> tri_ad;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				tri_ad(i, j).value() = tri(i, j);
				tri_ad(i, j).derivatives().setZero();
				tri_ad(i, j).derivatives()[i + j * 3] = 1;
			}
		}
		Eigen::Matrix<ADScalar, 2, 3> fr;
		auto e1 = (tri_ad.col(1) - tri_ad.col(0)).normalized();
		auto e2 = (tri_ad.col(2) - tri_ad.col(0)).eval();
		ADScalar A_ad = (tri_ad.col(2) - tri_ad.col(0)).cross(tri_ad.col(1) - tri_ad.col(0)).norm() / 2;
		e2 = (e2 - e2.dot(e1) * e1).normalized();
		fr.row(0) = e1.transpose();
		fr.row(1) = e2.transpose();
		Eigen::Matrix2<ADScalar> eps_ad = fr * eps * fr.transpose();
		
		ADScalar E = (lam0 * Eigen::pow(eps_ad.trace(), 2) + 2 * mu * (eps_ad * eps_ad.transpose()).trace()) * A_ad;
		Em += E.value();
		for (int i = 0; i < fvh.size(); i++) {
			dEm.block<3, 1>(fvh[i].idx() * 3, 0) += E.derivatives().block<3, 1>(i * 3, 0);
		}
	}
	return std::make_tuple(Em, dEm);
}

double eval_ads(
	PeriodSurfaceMesh& m,
	const std::vector<Eigen::Matrix3d>& frlist, const std::vector<Compile1ring>& vrings,
	const Eigen::Matrix3d& eps, double lam0, double mu
) {
	auto [fm, um, Enmem, As] = asym_stif_fem(m, frlist, vrings, lam0, mu);
	auto CA = eval_asym_elastic_tensor(fm, um, Enmem, As);
	auto eps_v = voigt(eps);
	return eps_v.dot(CA * eps_v);
}

void optimal_cond(std::string meshfile, std::string obj_type) {
	PeriodSurfaceMesh m;
	m.readMergePeriodBoundary(meshfile, false);

	auto [lam, mu] = mtl::lameCoeff(1, 0.3);
	double lam0 = 2 * lam * mu / (lam + 2 * mu);

	Eigen::Matrix3d eps; eps.setZero(); 
	eps << 1, 0, 1, 0, -2, 0, 1, 0, 1;

	for (int iter = 0; iter < 1000; iter++) {
		m.saveUnitCell(getPath("iter.obj"));

		if (iter % 4 == 0) { if (!asym_no_surgery) m.surgery(); }

		m.delaunayRemesh(delaunay_remesh_inner_iter, 0.2, 0.02, delaunay_remesh_outer_iter);

		m.saveUnitCell(getPath("iter.obj"));

		auto [frlist, vrings] = generate_face_frame_and_vring(m);

		if (iter % 5 == 0) {
			double ads = eval_ads(m, frlist, vrings, eps, lam0, mu);
			std::cout << "ads = " << ads << std::endl;
		}

		auto Av = get_mass_matrix(vrings);

		auto flist = face_load(m, frlist, vrings, lam0, mu, eps);

		auto Fvlist = distribute_vertex_load(m, flist, frlist);

		auto [Em, dEm] = homogeneous_membrane_energy(m, frlist, eps, lam0, mu);

		auto L = m.getPeriodicLaplacian(2, 1);

		// preconditioner
		Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> so(L);

		double mem_weight = 0.5;
		// non-homogeneous part
		Eigen::MatrixX3d F = /*Av.asDiagonal() **/ Eigen::Matrix3Xd::Map((const double*)Fvlist.data(), 3, Fvlist.size()).transpose();
		// homogeneous part	
		F += dEm.reshaped(3, F.rows()).transpose() * mem_weight;

		auto Fmean = F.colwise().mean().eval(); F.rowwise() -= Fmean;
		std::cout << "Fmean = " << F.colwise().mean() << std::endl;
		Eigen::MatrixX3d Fnew = -so.solve(F);
		Fmean = Fnew.colwise().mean().eval(); Fnew.rowwise() -= Fmean;
		for (int i = 0; i < Fvlist.size(); i++) { Fvlist[i] = Fnew.row(i).transpose(); }
		{
			std::ofstream ofs(getPath("fvlist"), std::ios::binary);
			for (auto vh : m.vertices()) {
				auto p = m.point(vh); ofs.write((const char*)p.data(), sizeof(p));
				auto fv = Fvlist[vh.idx()]; ofs.write((const char*)fv.data(), sizeof(fv));
			}
			ofs.close(); ofs.open(getPath("Em"), std::ios::binary);
			for (auto vh : m.vertices()) {
				auto p = m.point(vh); ofs.write((const char*)p.data(), sizeof(p));
				Eigen::Vector3d fv = dEm.block<3, 1>(vh.idx() * 3, 0); 
				ofs.write((const char*)fv.data(), sizeof(fv));
			}
		}

		double step = 0.01;
		for (auto vh : m.vertices()) {
			auto p = m.point(vh); auto fv = Fvlist[vh.idx()];
			auto pnew = p + step * toOM(fv); m.set_point(vh, pnew);
		}
	}
}
