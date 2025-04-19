//#include "solver/ShellFEM.h"
#include "fundamental_forms.h"
#include "PeriodicMesher.h"
#include <Eigen/PardisoSupport>
#include "matlab/matlab_utils.h"
#include "asymptotic_analysis.h"
#include <Eigen/Eigenvalues>
#include "material/materail.h"
#include <boost/algorithm/string.hpp>

#define REMOVE_RIGID_MOTION

using namespace msf;

void opti_asym_stif(PeriodSurfaceMesh& m);

std::string getPath(std::string s);

extern std::vector<Real> aux_number;

std::tuple<double, double> lamu(double Y, double nu) {
	double lam = Y * nu / (1 + nu) / (1 - 2 * nu);
	double mu = Y / 2 / (1 + nu);
	return { lam, mu };
}

void asymptotic_stiffnss(PeriodSurfaceMesh& m, const Eigen::Matrix3d& Em, double lam, double mu);

void test_asymptotic_stiffness(std::string meshfile) {
#if 0
	ShellFEMSolver fem;
	fem.loadMesh(meshfile, false, true, true);
	fem.mesher->mergePeriodBoundary();

	Eigen::Matrix3<Real> eps_M;
	eps_M <<0.3333333, 0, 0,
			0, 0.3333333, 0,
			0, 0, 0.3333333;
	Real E_A = fem.asymptoticStiffness(1, 0.3, eps_M);
	std::cout << "eps_M = \n" << eps_M << std::endl;
	std::cout << "E_A = " << E_A << std::endl;
#elif 0
	PeriodSurfaceMesh m;
	m.readMergePeriodBoundary(meshfile);
	Eigen::Matrix3d em;
	em <<
		0.33333333, 0, 0,
		0, 0.33333333, 0,
		0, 0, 0.33333333;
	em = (em + em.transpose()).eval() / 2;
	std::cout << "epsM = \n" << em << std::endl;
	auto [lam, mu] = lamu(1, 0.3);
	std::cout << "lam = " << lam << "; mu = " << mu << std::endl;
	asymptotic_stiffnss(m, em, lam, mu);
#else
	PeriodSurfaceMesh m;
	m.readMergePeriodBoundary(meshfile);
	opti_asym_stif(m);
#endif
}



void asymptotic_stiffnss(PeriodSurfaceMesh& m, const Eigen::Matrix3d& Em, double lam, double mu) {
	double lam0 = 2 * lam * mu / (lam + 2 * mu);
	
	std::vector<Compile1ring> vrings;
	double Aring = 0;
	for (auto vh : m.vertices()) {
		auto [o, ring] = m.find1ring(vh, 1);
		vrings.emplace_back(o, ring);
		Aring += vrings.rbegin()->As;
	}
	std::cout << "Aring = " << Aring << std::endl;
	{
		std::ofstream ofs("c1ring");
		for (auto vh : m.vertices()) {
			auto c = vrings[vh.idx()];
			ofs << c.H << " " << c.K << " " << c.nv.transpose() << " " << m.point(vh) << std::endl;
		}
	}

	std::vector<Eigen::Matrix<double, 3, 4>> edge_fr(m.n_edges() * 2);

	for (auto fh : m.faces()) {
		auto tri = m.getFacePeriodVertex(fh, 1);
		Eigen::Vector3d n = (tri.col(1) - tri.col(0)).cross(tri.col(2) - tri.col(0)).normalized();
		auto he = fh.halfedge();
		for (int i = 0; i < 3; i++) {
			auto ev = toEigen(m.calc_edge_vector(he));
			Eigen::Vector3d e1 = make_period(ev, 2, 1).normalized();
			Eigen::Vector3d e2 = n.cross(e1).normalized();
			if (he.edge().halfedge(0) == he) {
				edge_fr[he.edge().idx()].block<3, 1>(0, 0) = e1;
				edge_fr[he.edge().idx()].block<3, 1>(0, 1) = e2;
			}
			else {
				edge_fr[he.edge().idx()].block<3, 1>(0, 2) = -e1;
				edge_fr[he.edge().idx()].block<3, 1>(0, 3) = -e2;
			}
			he = he.next();
		}
	}

#if 0
	size_t ndof = m.n_vertices() + 2 * m.n_edges();
	Eigen::SparseMatrix<double> K(ndof, ndof);
	Eigen::VectorXd f(ndof); f.setZero();
	Eigen::VectorXd me(m.n_faces());
	double As = 0, Es = 0;
	std::vector<Eigen::Triplet<double>> triplist;
	std::vector<Eigen::Matrix<double, 3, 2>> face_fr(m.n_faces());
	std::vector<Eigen::Matrix<double, 2, 6>> translist(m.n_faces());
	std::vector<Eigen::Matrix<double, 9, 9>> Kelist(m.n_faces());
	std::vector<Eigen::Vector<double, 9>> Felist(m.n_faces());
	for (auto fh : m.faces()) {
		auto tri = m.getFacePeriodVertex(fh, 1);

		Compile1ring v[3];
		auto fvh = m.getFaceVertexHandle(fh);
		for (int i = 0; i < 3; i++) {
			v[i] = vrings[fvh[i].idx()];
		}
		auto [Ke, fe, e1, e2] = membrane_element_stif_matrix_vector(lam0, mu, tri, v, Em);
		face_fr[fh.idx()] << e1, e2;
		Kelist[fh.idx()] = Ke;
		Felist[fh.idx()] = fe;
		auto [Ee, Ae] = membrane_strain_energy(lam0, mu, tri, Em);
		As += Ae; Es += Ee;
		me[fh.idx()] = Ee;
		Eigen::Matrix<double, 3, 2> elocal;
		elocal << e1, e2;

		auto he = m.findHe(fvh[0], fvh[1]);
		Eigen::Matrix2d trans[3];
		int eid[3];
		for (int i = 0; i < 3; i++) {
			Eigen::Matrix<double, 3, 2> eglobal;
			if (he.edge().halfedge(0) == he) {
				eglobal = edge_fr[he.edge().idx()].leftCols(2);
			} else {
				eglobal = edge_fr[he.edge().idx()].rightCols(2);
			}
			Eigen::Matrix2d tran = elocal.transpose() * eglobal;
			trans[i] = tran;
			eid[i] = he.edge().idx();
			he = he.next();
		}
		translist[fh.idx()] << trans[0], trans[1], trans[2];
		for (int i = 3; i < 9; i += 2) {
			auto& tran_i = trans[(i - 3) / 2];
			Ke.block<9, 2>(0, i) *= tran_i;
			Ke.block<2, 9>(i, 0).applyOnTheLeft(tran_i.transpose());
			fe.block<2, 1>(i, 0).applyOnTheLeft(tran_i.transpose());
			//for (int j = 3; j < 9; j += 2) {
			//	auto& tran_j = trans[(j - 3) / 2];
			//	K.block<2, 2>(i, j) = (tran_i.transpose() * K.block<2, 2>(i, j) * tran_j).eval();
			//}
		}
		for (int i = 0; i < 9; i++) {
			int row = i < 3 ? fvh[i].idx() : eid[(i - 3) / 2] * 2 + (i - 3) % 2 + m.n_vertices();
			for (int j = 0; j < 9; j++) {
				int col = j < 3 ? fvh[j].idx() : eid[(j - 3) / 2] * 2 + (j - 3) % 2 + m.n_vertices();
				triplist.emplace_back(row, col, Ke(i, j));
			}
			f[row] += fe[i];
		}
		
	}

	K.setFromTriplets(triplist.begin(), triplist.end());
#else
	size_t ndof = m.n_vertices() * 3;
	Eigen::SparseMatrix<double> K(ndof, ndof);
	Eigen::MatrixXd fm(ndof, 6); fm.setZero();
	Eigen::VectorXd f(ndof); f.setZero();
	Eigen::VectorXd me(m.n_faces());
	double As = 0, Es = 0;
	Eigen::Matrix<double, 6, 6> Enmem; Enmem.setZero();
	std::vector<Eigen::Matrix<double, 9, 9>> Kelist(m.n_faces());
	std::vector<Eigen::Vector<double, 9>> Felist(m.n_faces());
	std::vector<Eigen::Triplet<double>> triplist;
	auto Emvgt = voigt(Em);
	for (auto fh : m.faces()) {
		auto tri = m.getFacePeriodVertex(fh, 1);
		Compile1ring v[3];
		auto fvh = m.getFaceVertexHandle(fh);
		for (int i = 0; i < 3; i++) {
			v[i] = vrings[fvh[i].idx()];
		}
		auto face_fr = face_frame(tri);
		auto [Ke, feem] = membrane_element_vertex_stif_matrix_vector(lam0, mu, tri, face_fr, v);
		Eigen::Vector<double, 9> fe = feem * Emvgt;
		auto [Eem, Ae] = membrane_strain_energy(lam0, mu, tri, face_fr, v);
		Enmem += Eem;
		double Ee = Emvgt.dot(Eem * Emvgt);
		As += Ae; Es += Ee;
		me[fh.idx()] = Ee;

		for (int i = 0; i < 9; i++) {
			int vi = fvh[i / 3].idx();
			for (int j = 0; j < 9; j++) {
				int vj = fvh[j / 3].idx();
				triplist.emplace_back(vi * 3 + i % 3, vj * 3 + j % 3, Ke(i, j));
			}
			f[vi * 3 + i % 3] += fe[i];
			fm.row(vi * 3 + i % 3) += feem.row(i);
		}
		Kelist[fh.idx()] = Ke;
		Felist[fh.idx()] = fe;
	}
	K.setFromTriplets(triplist.begin(), triplist.end());
#endif

	
	eigen2ConnectedMatlab("K", K);
	eigen2ConnectedMatlab("f", f);

	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> so(K.transpose() * K);
	Eigen::VectorXd kf = K.transpose() * f;
	Eigen::VectorXd u = so.solve(kf);

	{
		std::ofstream ofs("usol", std::ios::binary);
#if 1
		ofs.write((const char*)u.data(), m.n_vertices() * sizeof(double) * 3);
#else
		ofs.write((const char*)u.data(), m.n_vertices() * sizeof(double));
		ofs.close();

		ofs.open("esol", std::ios::binary);
		for (auto fh : m.faces()) {
			auto ke = Kelist[fh.idx()];
			auto fe = Felist[fh.idx()];
			Eigen::Matrix2d trans[3];
			for (int i = 0; i < 3; i++) {
				trans[i] = translist[fh.idx()].block<2, 2>(0, 2 * i);
			}
			Eigen::Vector<double, 9> ue;
			auto fvh = m.getFaceVertexHandle(fh);
			auto tri = m.getFacePeriodVertex(fh, 1);
			auto he = m.findHe(fvh[0], fvh[1]);
			auto fr = face_fr[fh.idx()];
			for (int i = 0; i < 3; i++) {
				ue[i] = u[fvh[i].idx()];
				int eid = he.edge().idx();
				ue.block<2, 1>(3 + 2 * i, 0) = trans[i] * u.block<2, 1>(m.n_vertices() + 2 * eid, 0);
				Eigen::Vector3d uv = fr * ue.block<2, 1>(3 + 2 * i, 0);
				Eigen::Vector3d p = (tri.col((i + 1) % 3) + tri.col(i)) / 2;
				ofs.write((const char*)p.data(), sizeof(p));
				ofs.write((const char*)uv.data(), sizeof(uv));
				he = he.next();
			}
		}
#endif
	}

	eigen2ConnectedMatlab("u", u);

	double em = (Es - u.dot(f)) / As;

	std::cout << "As = " << As << std::endl;
	std::cout << "Es = " << Es << std::endl;
	std::cout << "ea = " << em << ";  em = " << Es / As << std::endl;

	std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;

	Eigen::VectorXd enlist(m.n_faces()); enlist.setZero();
	for (auto fh : m.faces()) {
		auto ke = Kelist[fh.idx()];
		auto fe = Felist[fh.idx()];
		Eigen::Vector<double, 9> ue;
		auto fvh = m.getFaceVertexHandle(fh);
#if 0
		Eigen::Matrix2d trans[3];
		for (int i = 0; i < 3; i++) { trans[i] = translist[fh.idx()].block<2, 2>(0, 2 * i); }
		auto he = m.findHe(fvh[0], fvh[1]);
		for (int i = 0; i < 3; i++) {
			ue[i] = u[fvh[i].idx()];
			int eid = he.edge().idx();
			ue.block<2, 1>(3 + 2 * i, 0) = trans[i] * u.block<2, 1>(m.n_vertices() + 2 * eid, 0);
			he = he.next();
		}

#else
		for (int i = 0; i < 3; i++) {
			ue.block<3, 1>(i * 3, 0) = u.block<3, 1>(fvh[i].idx() * 3, 0);
		}
#endif
		double e = ue.dot(ke * ue) - 2 * ue.dot(fe) + me[fh.idx()];
		enlist[fh.idx()] = e;
	}
	eigen2ConnectedMatlab("enlist", enlist);

	std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
	// solve for asymptotic matrix
	Eigen::MatrixXd kfm = K.transpose() * fm;
	Eigen::MatrixXd um = so.solve(kfm);
	eigen2ConnectedMatlab("fm", fm);
	eigen2ConnectedMatlab("um", um);

	std::cout << "Pbar = \n" << Enmem / As << std::endl;

	Eigen::Matrix<double, 6, 6> CA;
	for (int i = 0; i < 6; i++) {
		for (int j = i; j < 6; j++) {
			CA(i, j) = (Enmem(i, j) - um.col(i).dot(fm.col(j))) / As;
			CA(j, i) = CA(i, j);
		}
	}
	std::cout << "CA = \n" << CA << std::endl;
}

std::tuple<
	Eigen::Matrix<double, -1, 6>,
	Eigen::Matrix<double, -1, 6>,
	Eigen::Matrix<double, 6, 6>,
	double
> asym_stif_fem(PeriodSurfaceMesh& m, const std::vector<Eigen::Matrix3d>& frlist, const std::vector<Compile1ring>& vrings, double lam0, double mu) {

	size_t ndof = m.n_vertices() * 3;
	Eigen::SparseMatrix<double> K(ndof, ndof);
	Eigen::Matrix<double, -1, 6> fm(ndof, 6); fm.setZero();
	Eigen::VectorXd me(m.n_faces());
	double As = 0, Es = 0;
	Eigen::Matrix<double, 6, 6> Enmem; Enmem.setZero();
	std::vector<Eigen::Matrix<double, 9, 9>> Kelist(m.n_faces());
	std::vector<Eigen::Matrix<double, 9, 6>> Felist(m.n_faces());
	std::vector<Eigen::Triplet<double>> triplist;
	//frlist.resize(m.n_faces());
	std::vector<OM::SmartFaceHandle> fhlist;
	for (auto fh : m.faces()) fhlist.emplace_back(fh);
#pragma omp parallel for
	for (int fid = 0; fid < m.n_faces(); fid++) {
		auto fh = fhlist[fid];
		auto tri = m.getFacePeriodVertex(fh, 1);
		auto& fr = frlist[fh.idx()];
		//fr = face_frame(tri);
		Compile1ring v[3];
		auto fvh = m.getFaceVertexHandle(fh);
		for (int i = 0; i < 3; i++) { v[i] = vrings[fvh[i].idx()]; }
		auto [Ke, feem] = membrane_element_vertex_stif_matrix_vector(lam0, mu, tri, fr, v);
		auto [Eem, Ae] = membrane_strain_energy(lam0, mu, tri, fr, v);

#pragma omp critical
		{
			Enmem += Eem; As += Ae;
			for (int i = 0; i < 9; i++) {
				int vi = fvh[i / 3].idx();
				for (int j = 0; j < 9; j++) {
					int vj = fvh[j / 3].idx();
					triplist.emplace_back(vi * 3 + i % 3, vj * 3 + j % 3, Ke(i, j));
				}
				fm.row(vi * 3 + i % 3) += feem.row(i);
			}
		}
		Kelist[fh.idx()] = Ke;
		Felist[fh.idx()] = feem;
	}

	K.setFromTriplets(triplist.begin(), triplist.end());
	
	//eigen2ConnectedMatlab("K", K);
	//eigen2ConnectedMatlab("fm", fm);

#ifdef REMOVE_RIGID_MOTION
	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> so(K);
	//fm = (fm.rowwise() - fm.colwise().mean().eval()).eval();
	Eigen::Matrix<double, -1, 6> um = so.solve(fm);
#else
	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> so(K.transpose() * K);
	Eigen::Matrix<double, -1, 6> kf = K.transpose() * fm;
	Eigen::Matrix<double, -1, 6> um = so.solve(kf);
#endif

	//eigen2ConnectedMatlab("um", um);

	return std::make_tuple(fm, um, Enmem, As);
}

std::tuple<
	Eigen::SparseMatrix<double>,
	Eigen::Matrix<double, -1, 6>,
	Eigen::Matrix<double, 6, 6>,
	double
> assembly_stif_fem(msf::PeriodSurfaceMesh& m, const std::vector<Eigen::Matrix3d>& frlist, const std::vector<Compile1ring>& vrings, double lam0, double mu)
{
	size_t ndof = m.n_vertices() * 3;
	Eigen::SparseMatrix<double> K(ndof, ndof);
	Eigen::Matrix<double, -1, 6> fm(ndof, 6); fm.setZero();
	Eigen::VectorXd me(m.n_faces());
	double As = 0, Es = 0;
	Eigen::Matrix<double, 6, 6> Enmem; Enmem.setZero();
	std::vector<Eigen::Matrix<double, 9, 9>> Kelist(m.n_faces());
	std::vector<Eigen::Matrix<double, 9, 6>> Felist(m.n_faces());
	std::vector<Eigen::Triplet<double>> triplist;
	//frlist.resize(m.n_faces());
	std::vector<OM::SmartFaceHandle> fhlist;
	for (auto fh : m.faces()) fhlist.emplace_back(fh);
#pragma omp parallel for
	for (int fid = 0; fid < m.n_faces(); fid++) {
		auto fh = fhlist[fid];
		auto tri = m.getFacePeriodVertex(fh, 1);
		auto& fr = frlist[fh.idx()];
		//fr = face_frame(tri);
		Compile1ring v[3];
		auto fvh = m.getFaceVertexHandle(fh);
		for (int i = 0; i < 3; i++) { v[i] = vrings[fvh[i].idx()]; }
		auto [Ke, feem] = membrane_element_vertex_stif_matrix_vector(lam0, mu, tri, fr, v);
		auto [Eem, Ae] = membrane_strain_energy(lam0, mu, tri, fr, v);

#pragma omp critical
		{
			Enmem += Eem; As += Ae;
			for (int i = 0; i < 9; i++) {
				int vi = fvh[i / 3].idx();
				for (int j = 0; j < 9; j++) {
					int vj = fvh[j / 3].idx();
					triplist.emplace_back(vi * 3 + i % 3, vj * 3 + j % 3, Ke(i, j));
				}
				fm.row(vi * 3 + i % 3) += feem.row(i);
			}
		}
		Kelist[fh.idx()] = Ke;
		Felist[fh.idx()] = feem;
	}

	K.setFromTriplets(triplist.begin(), triplist.end());

	return std::make_tuple(K, fm, Enmem, As);
}


std::tuple<
	Eigen::Matrix<double, -1, 6>,
	Eigen::Matrix<double, -1, 6>,
	Eigen::Matrix<double, 6, 6>,
	double
> asym_stif_fem(PeriodSurfaceMesh& m, const std::vector<Eigen::Matrix3d>& frlist, const std::vector<Compile1ring>& vrings, double lam0, double mu, const std::vector<Eigen::Matrix3d>& bvlist) {
	size_t ndof = m.n_vertices() * 3;
	Eigen::SparseMatrix<double> K(ndof, ndof);
	Eigen::Matrix<double, -1, 6> fm(ndof, 6); fm.setZero();
	Eigen::VectorXd me(m.n_faces());
	double As = 0, Es = 0;
	Eigen::Matrix<double, 6, 6> Enmem; Enmem.setZero();
	std::vector<Eigen::Triplet<double>> triplist;
	//frlist.resize(m.n_faces());
	for (auto fh : m.faces()) {
		auto tri = m.getFacePeriodVertex(fh, 1);
		auto& fr = frlist[fh.idx()];
		//fr = face_frame(tri);
		Compile1ring v[3];
		auto fvh = m.getFaceVertexHandle(fh);
		for (int i = 0; i < 3; i++) { v[i] = vrings[fvh[i].idx()]; }
		auto bv_e = retrieve_element_vertex_data(m, fh, bvlist);
		auto [Ke, feem] = membrane_element_vertex_stif_matrix_vector(lam0, mu, tri, fr, v, bv_e.data());
		auto [Eem, Ae] = membrane_strain_energy(lam0, mu, tri, fr, v);
		Enmem += Eem;
		As += Ae;

		for (int i = 0; i < 9; i++) {
			int vi = fvh[i / 3].idx();
			for (int j = 0; j < 9; j++) {
				int vj = fvh[j / 3].idx();
				triplist.emplace_back(vi * 3 + i % 3, vj * 3 + j % 3, Ke(i, j));
			}
			fm.row(vi * 3 + i % 3) += feem.row(i);
		}
	}

	K.setFromTriplets(triplist.begin(), triplist.end());

	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> so(K.transpose() * K);
	Eigen::Matrix<double, -1, 6> kf = K.transpose() * fm;
	Eigen::Matrix<double, -1, 6> um = so.solve(kf);

	return std::make_tuple(fm, um, Enmem, As);
}

Eigen::SparseMatrix<double> asym_stif_fem_matrix(msf::PeriodSurfaceMesh& m, const std::vector<Eigen::Matrix3d>& frlist, const std::vector<Compile1ring>& vrings, double lam0, double mu, const std::vector<Eigen::Matrix3d>& bvlist)
{
	size_t ndof = m.n_vertices() * 3;
	Eigen::SparseMatrix<double> K(ndof, ndof);
	std::vector<Eigen::Triplet<double>> triplist;
	for (auto fh : m.faces()) {
		auto tri = m.getFacePeriodVertex(fh, 1);
		auto fvh = m.getFaceVertexHandle(fh);
		auto& fr = frlist[fh.idx()];
		auto v = retrieve_element_vertex_data(m, fh, vrings);
		auto bv_e = retrieve_element_vertex_data(m, fh, bvlist);
		auto [Ke, feem] = membrane_element_vertex_stif_matrix_vector(lam0, mu, tri, fr, v.data(), bv_e.data());
		for (int i = 0; i < 9; i++) {
			int vi = fvh[i / 3].idx();
			for (int j = 0; j < 9; j++) {
				int vj = fvh[j / 3].idx();
				triplist.emplace_back(vi * 3 + i % 3, vj * 3 + j % 3, Ke(i, j));
			}
		}
	}
	K.setFromTriplets(triplist.begin(), triplist.end());
	return K;
}

Eigen::SparseMatrix<double> asym_stif_fem_matrix(msf::PeriodSurfaceMesh& m, const std::vector<Eigen::Matrix3d>& frlist, const std::vector<Compile1ring>& vrings, double lam0, double mu)
{
	size_t ndof = m.n_vertices() * 3;
	Eigen::SparseMatrix<double> K(ndof, ndof);
	std::vector<Eigen::Triplet<double>> triplist;
	for (auto fh : m.faces()) {
		auto tri = m.getFacePeriodVertex(fh, 1);
		auto fvh = m.getFaceVertexHandle(fh);
		auto& fr = frlist[fh.idx()];
		auto v = retrieve_element_vertex_data(m, fh, vrings);
		auto [Ke, feem] = membrane_element_vertex_stif_matrix_vector(lam0, mu, tri, fr, v.data());
		for (int i = 0; i < 9; i++) {
			int vi = fvh[i / 3].idx();
			for (int j = 0; j < 9; j++) {
				int vj = fvh[j / 3].idx();
				triplist.emplace_back(vi * 3 + i % 3, vj * 3 + j % 3, Ke(i, j));
			}
		}
	}
	K.setFromTriplets(triplist.begin(), triplist.end());
	return K;
}



std::tuple<Eigen::Matrix<double, -1, 21>, Eigen::VectorXd> asym_stif_shape_sensitivity(
	PeriodSurfaceMesh& m, const std::vector<Eigen::Matrix3d>& fr, const std::vector<Compile1ring>& vrings,
	const Eigen::Matrix<double, -1, 6>& um, const Eigen::Matrix<double, -1, 6>& fm,
	double lam0, double mu
) {
	Eigen::Matrix<double, -1, 21> Vn_en(m.n_vertices(), 21); Vn_en.setZero();
	Eigen::VectorXd Vn_A(m.n_vertices()); Vn_A.setZero();
	double As = 0;
	std::vector<OM::SmartFaceHandle> fhlist;
	for (auto fh : m.faces()) fhlist.emplace_back(fh);
#pragma omp parallel for
	for (int fid = 0; fid < fhlist.size(); fid++) {
		auto fh = fhlist[fid];
		auto tri = m.getFacePeriodVertex(fh, 1);
		auto fvh = m.getFaceVertexHandle(fh);
		Compile1ring v[3];
		for (int i = 0; i < 3; i++) v[i] = vrings[fvh[i].idx()];

		Eigen::Matrix<double, 9, 6> ue;
		for (int k = 0; k < 3; k++) {
			ue.middleRows(k * 3, 3) = um.middleRows(fvh[k].idx() * 3, 3);
		}

		// integral part (nominator)
		auto Vne = asym_stif_element_shape_derivative(tri, fr[fh.idx()], v, lam0, mu, ue, nullptr);

		// area part (denominator)
		//As += fr[fh.idx()].col(2).norm();
		auto dAdv = area_shape_derivative(tri, fr[fh.idx()], v);

#pragma omp critical
		{
			for (int k = 0; k < 3; k++) { Vn_en.row(fvh[k].idx()) += Vne.row(k); }
			for (int i = 0; i < 3; i++) { Vn_A[fvh[i].idx()] += dAdv[i]; }
		}
	}


	return std::make_tuple(Vn_en, Vn_A);
}


void show_asym_stif(std::string mfile) {
	PeriodSurfaceMesh m;
	m.readMergePeriodBoundary(mfile);
	auto [lam, mu] = lamu(1, 0.3);
	double lam0 = 2 * lam * mu / (lam + 2 * mu);

	auto [frlist, vrings] = generate_face_frame_and_vring(m);

	auto [fm, um, Enmem, As] = asym_stif_fem(m, frlist, vrings, lam0, mu);

	Eigen::Matrix<double, 6, 6> CA;

	for (int i = 0; i < 6; i++) {
		for (int j = i; j < 6; j++) {
			CA(i, j) = (Enmem(i, j) - um.col(i).dot(fm.col(j))) / As;
			CA(j, i) = CA(i, j);
		}
	}
	
	{
		std::ofstream ofs(getPath("enmem"));
		ofs << Enmem / As;
		ofs.close(); ofs.open(getPath("ads"));
		ofs << CA;
	}

	std::cout << "EM = \n" << Enmem << std::endl;
	std::cout << "EA = \n" << CA << std::endl;
}

void show_asym_stif(std::string mfile);
void output_abaqus_format(std::ostream& os, const Eigen::Matrix<double, 6, 6>& C) {
	std::map<int, int> id = {  };
	id[11] = 0; id[22] = 1; id[33] = 2;
	id[12] = 5; id[23] = 3; id[13] = 4;
	os << "abaqus format CA =\n" <<
		C(id[11], id[11]) << '\n' <<
		C(id[11], id[22]) << '\n' <<
		C(id[22], id[22]) << '\n' <<
		C(id[11], id[33]) << '\n' <<
		C(id[22], id[33]) << '\n' <<
		C(id[33], id[33]) << '\n' <<
		C(id[11], id[12]) << '\n' <<
		C(id[22], id[12]) << '\n' <<
		C(id[33], id[12]) << '\n' <<
		C(id[12], id[12]) << '\n' <<
		C(id[11], id[13]) << '\n' <<
		C(id[22], id[13]) << '\n' <<
		C(id[33], id[13]) << '\n' <<
		C(id[12], id[13]) << '\n' <<
		C(id[13], id[13]) << '\n' <<
		C(id[11], id[23]) << '\n' <<
		C(id[22], id[23]) << '\n' <<
		C(id[33], id[23]) << '\n' <<
		C(id[12], id[23]) << '\n' <<
		C(id[13], id[23]) << '\n' <<
		C(id[23], id[23]) << '\n' << std::endl;
}

void test_gamma_accuracy(void);
extern void test_asym_stif_shape_sensitive(void);

void opti_asym_stif(PeriodSurfaceMesh& m) {
	auto [lam, mu] = lamu(1, 0.3);
	double lam0 = 2 * lam * mu / (lam + 2 * mu);
	
	//test_asym_stif_shape_sensitive(); exit(0);
	//test_gamma_accuracy(); exit(0);

	Eigen::Matrix3d etgt;
	etgt <<
		1. / 3, 0, 0,
		0, 1. / 3, 0,
		0, 0, 1. / 3;
	etgt = (etgt + etgt.transpose()).eval() / 2;
	std::cout << "etgt = \n" << etgt << std::endl;
	auto Evgt = voigt(etgt);
	std::cout << "Evgt = " << Evgt.transpose() << std::endl;

	Eigen::VectorXd dfdC(21); dfdC.setZero();

	int counter = 0;
	for (int i = 0; i < 6; i++) {
		for (int j = i; j < 6; j++) {
			dfdC[counter] = -Evgt[i] * Evgt[j];
			if (i != j) dfdC[counter] *= 2;
			counter++;
		}
	}


	for (int iter = 0; iter < 200; iter++) {
		m.savePeriodicMesh("iter.obj", std::vector<OM::SmartVertexHandle>{}, 1.);
		std::cout << "\033[32m* Iter " << iter << "\033[0m" << std::endl;
		// update vring information
		auto [frlist, vrings] = generate_face_frame_and_vring(m);

		{
			std::ofstream ofs("Vp", std::ios::binary);
			for (auto vh : m.vertices()) {
				auto p = m.point(vh); ofs.write((const char*)p.data(), sizeof(p)); 
				auto n = vrings[vh.idx()].nv; ofs.write((const char*)n.data(), sizeof(n));
			}
		}

		auto [fm, um, Enmem, As] = asym_stif_fem(m, frlist, vrings, lam0, mu);

		{
			std::ofstream ofs("um", std::ios::binary);
			Eigen::VectorXd U = um * voigt(etgt);
			ofs.write((const char*)U.data(), U.size() * sizeof(double)); ofs.close();

			ofs.open("fm", std::ios::binary);
			Eigen::VectorXd F = fm * voigt(etgt);
			ofs.write((const char*)F.data(), F.size() * sizeof(double)); ofs.close();
		}

		Eigen::Matrix<double, 6, 6> CA;

		for (int i = 0; i < 6; i++) {
			for (int j = i; j < 6; j++) {
				CA(i, j) = (Enmem(i, j) - um.col(i).dot(fm.col(j))) / As;
				CA(j, i) = CA(i, j);
			}
		}

		std::cout << "As = " << As << std::endl;
		std::cout << "Em = \n" << Enmem / As << std::endl;
		std::cout << "Ea = \n" << CA << std::endl;
		output_abaqus_format(std::cout, CA);

		auto [Vn_ij, Vn_A] = asym_stif_shape_sensitivity(m, frlist, vrings, um, fm, lam0, mu);
		auto [Av, Nv] = eval_vertex_mass(m);
		Vn_ij.applyOnTheLeft(Av.cwiseInverse().asDiagonal()); Vn_A.applyOnTheLeft(Av.cwiseInverse().asDiagonal());

		double f = -Evgt.dot(CA * Evgt);

		std::cout << "*obj = " << f << std::endl;

		Eigen::VectorXd Vn_E = Vn_ij * dfdC;

		{
			std::ofstream ofs("VnA", std::ios::binary); ofs.write((const char*)Vn_A.data(), Vn_A.size() * sizeof(double)); ofs.close();
			ofs.open("VnE", std::ios::binary); ofs.write((const char*)Vn_E.data(), Vn_E.size() * sizeof(double));
		}

		Eigen::VectorXd Vn(m.n_vertices());

		Vn = Vn_E / As - Vn_A * f / As; // f already divides As

		//eigen2ConnectedMatlab("VnA", Vn_A);
		//eigen2ConnectedMatlab("VnE", Vn_E);
		//eigen2ConnectedMatlab("Vn", Vn);
		//Vn *= -1;

		// precondition
		if (1) {
			Eigen::SparseMatrix<double> L = -m.getPeriodicLaplacian(2, 1);
			//eigen2ConnectedMatlab("L", L);
			//Eigen::SparseMatrix<double> eye(L.rows(), L.rows()); eye.setIdentity();
			Eigen::VectorXd mass(L.rows());
			for (auto vh : m.vertices()) mass[vh.idx()] = vrings[vh.idx()].mass;
			Eigen::VectorXd fv = Vn.cwiseProduct(mass);
			//eigen2ConnectedMatlab("m", mass);
			//eigen2ConnectedMatlab("Vm", Vn.cwiseProduct(mass));
			if (!aux_number.empty()) {
				mass *= aux_number[0];
			}
			L += mass.asDiagonal();
			//eigen2ConnectedMatlab("Lm", L);
			Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> Ls(L);
			Vn = Ls.solve(fv);
			//eigen2ConnectedMatlab("PVn", Vn);
		}

		Eigen::VectorXd VnV(m.n_vertices());
		Eigen::VectorXd VnH(m.n_vertices());
		double step_E = 5;
		double step_H = 0.3; // 0.5
		if (aux_number.size() > 1) step_E = aux_number[1];
		if (aux_number.size() > 2) step_H = aux_number[2];
		for (auto vh : m.vertices()) {
			VertexFlag vflag;
			auto pv = m.point(vh);
			vflag.set_period_boundary(1 - 1e-5, pv[0], pv[1], pv[2]);
			Eigen::Vector3d vn_e = -Vn[vh.idx()] * vrings[vh.idx()].nv;
			Eigen::Vector3d vn_h = vrings[vh.idx()].nv * vrings[vh.idx()].H * vrings[vh.idx()].As;
			Eigen::Vector3d vn = vn_e * step_E + vn_h * step_H;
			VnH[vh.idx()] = vrings[vh.idx()].H * vrings[vh.idx()].As * step_H;
			VnV[vh.idx()] = vn.dot(vrings[vh.idx()].nv);
			for (int k = 0; k < 3; k++) {
				if (!vflag.is_period_boundary(k)) { pv[k] += vn[k]; }
			}
			m.set_point(vh, toOM(pv));
		}
		{
			std::ofstream ofs("VnV", std::ios::binary);
			ofs.write((const char*)VnV.data(), sizeof(double) * VnV.size()); ofs.close();
			ofs.open("VnH", std::ios::binary);
			ofs.write((const char*)VnH.data(), sizeof(double) * VnH.size()); ofs.close();
		}
		std::cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n";
	}
}

void test_second_fundamental_form(std::string mfile) {
	PeriodSurfaceMesh m;
	m.readMergePeriodBoundary(mfile);
	
	std::vector<Compile1ring> vrings;
	std::vector<Eigen::Matrix3d> frlist(m.n_faces());
	for (auto vh : m.vertices()) {
		auto [o, ring] = m.find1ring(vh, 1);
		vrings.emplace_back(o, ring); 
	}
	std::vector<std::tuple<Eigen::Vector3d, Eigen::Matrix<double, 3, 2>>> eigs(m.n_faces());
	for (auto fh : m.faces()) {
		auto tri = m.getFacePeriodVertex(fh, 1);
		auto fvh = m.getFaceVertexHandle(fh);
		Compile1ring v[3];
		for (int k = 0; k < 3; k++) v[k] = vrings[fvh[k].idx()];
		auto fr = face_frame(tri);
		auto Bn = strain_matrix_edge_stretch(tri, fr.col(0), fr.col(1));
		auto be = second_fundamental_form(tri, v[0], v[1], v[2]);
		auto bform = fromvoigt((Bn * be).eval());
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> so(bform);
		auto d = so.eigenvalues();
		auto V = so.eigenvectors();
		Eigen::Matrix<double, 3, 2> V3;
		V3 <<
			d[0] * fr.leftCols(2) * V.col(0),
			d[1] * fr.leftCols(2)* V.col(1);
		Eigen::Vector3d fc = tri.rowwise().mean();
		eigs[fh.idx()] = std::make_tuple(fc, V3);
	}

	std::ofstream ofs("bform", std::ios::binary);
	for (int k = 0; k < eigs.size(); k++) {
		auto [fc, V] = eigs[k];
		ofs.write((const char*)fc.data(), sizeof(fc));
		ofs.write((const char*)V.data(), sizeof(V));
	}
}

void test_energy_min(std::string mfile) {
	PeriodSurfaceMesh m;
	m.readMergePeriodBoundary(mfile);
	auto [lam, mu] = lamu(1, 0.3);
	double lam0 = 2 * lam * mu / (lam + 2 * mu);

	auto [frlist, vrings] = generate_face_frame_and_vring(m);

	auto [fm, um, Enmem, As] = asym_stif_fem(m, frlist, vrings, lam0, mu);

	Eigen::Matrix<double, 6, 6> CA;

	for (int i = 0; i < 6; i++) {
		for (int j = i; j < 6; j++) {
			CA(i, j) = (Enmem(i, j) - um.col(i).dot(fm.col(j))) / As;
			CA(j, i) = CA(i, j);
		}
	}

	std::cout << "EM = \n" << Enmem << std::endl;
	std::cout << "EA = \n" << CA << std::endl;

	Eigen::Matrix3d eps_M; eps_M.setZero();
	eps_M(0, 0) = 1;

	double C11 = 0;
	for (auto fh : m.faces()) {
		auto fvh = m.getFaceVertexHandle(fh);
		auto tri = m.getFacePeriodVertex(fh, 1);
		Compile1ring v[3] = { vrings[fvh[0].idx()], vrings[fvh[1].idx()],vrings[fvh[2].idx()] };
		Eigen::Vector<double, 9> ue;
		for (int k = 0; k < 3; k++) {
			ue.block<3, 1>(k * 3, 0) = um.block<3, 1>(fvh[k].idx() * 3, 0);
		}
		double em = membrane_strain_energy(lam0, mu, tri, frlist[fh.idx()], v, ue, eps_M);
		C11 += em;
	}

	{
		std::ofstream ofs("Vp", std::ios::binary);
		for (auto vh : m.vertices()) {
			auto p = m.point(vh);
			ofs.write((const char*)p.data(), sizeof(p));
		}
		ofs.close(); ofs.open("um", std::ios::binary);
		ofs.write((const char*)um.col(0).data(), sizeof(double) * um.rows());
	}

	auto [Av, Nv] = eval_vertex_mass(m);

	std::cout << "C11 = " << C11 / As << std::endl;
}

template<typename T1, typename T2>
auto sym_inner_prod(const Eigen::Matrix<T1, 6, 6>& A1, const Eigen::Matrix<T2, 6, 6>& A2) {
	Eigen::Matrix<double, 6, 6> M; M.setOnes();
	M.rightCols(3) *= 2; M.bottomRows(3) *= 2;
	return A1.cwiseProduct(M).cwiseProduct(A2).sum();
}
template<typename T1, typename T2>
auto sym_inner_prod(const Eigen::Matrix<T1, 6, 1>& e1, const Eigen::Matrix<T2, 6, 1>& e2) {
	Eigen::Matrix<double, 6, 1> M = { 1, 1, 1, 0.5, 0.5, 0.5 };
	return e1.dot(M.cwiseProduct(e2));
}

using ADScalar = Eigen::AutoDiffScalar<Eigen::Vector<double, 21>>;

Eigen::Matrix<double, 6, 6> elastic_matrix(double Y, double nu) {
	auto [lam, mu] = mtl::lameCoeff(Y, nu);
	Eigen::Matrix<double, 6, 6> D; D.setZero();
	D.block<3, 3>(0, 0).setConstant(lam);
	D.block<3, 3>(0, 0).diagonal() += Eigen::Vector3d(2 * mu, 2 * mu, 2 * mu);
	D.block<3, 3>(3, 3).diagonal().setConstant(mu);
	return D;
}

std::vector<Eigen::Vector3d> parse_point_list(std::string s) {
	std::vector<std::string> coordstr;
	boost::split(coordstr, s, boost::is_any_of("{}, "));
	coordstr.erase(std::remove_if(coordstr.begin(), coordstr.end(), [](const std::string& s) {return s.empty(); }), coordstr.end());
	if (coordstr.size() % 3 != 0) {
		std::cout << "Invalid coordinate string " << s << std::endl;
		throw std::runtime_error("invalid point list");
	}
	std::vector<Eigen::Vector3d> plist;
	for (int i = 0; i < coordstr.size(); i += 3) {
		plist.emplace_back(std::stod(coordstr[i]), std::stod(coordstr[i + 1]), std::stod(coordstr[i + 2]));
	}
	return plist;
}

inline Eigen::Matrix3d r3_e_ij(int i, int j) {
	return Eigen::Vector3d::Unit(i) * Eigen::Vector3d::Unit(j).transpose();
}

ADScalar parse_one_ads_objective(std::string subtype_str, const Eigen::Matrix<ADScalar, 6, 6> CA_ad, const Eigen::Matrix<double, 6, 6>& CA) {
	double weight = 1;
	std::string sub;
	std::string alphabet = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
	int afirst = subtype_str.find_first_of(alphabet);
	if (afirst != 0) weight = std::stod(subtype_str.substr(0, afirst));
	sub = subtype_str.substr(afirst);
	if (sub[0] == 'c') {
		if (std::isdigit(sub[1]) && std::isdigit(sub[2])) {
			int i = sub[1] - '0', j = sub[2] - '0';
			if (i < 6 && j < 6) { return weight * CA_ad(i, j); }
		}
	}
	else if (sub == "bulk") {
		return weight * CA_ad.block<3, 3>(0, 0).reshaped().sum() / 9;
	}
	else if (sub == "norm") {
		return  weight * Eigen::sqrt(sym_inner_prod(CA_ad, CA_ad));
	}
	else if (sub[0] == 'E') {
		auto nustr = sub.find_first_of("n");
		if (nustr != std::string::npos) {
			double E = std::stod(sub.substr(1, nustr - 1));
			double nu = std::stod(sub.substr(nustr + 2));
			auto D = elastic_matrix(E, nu);
			auto dist = (D - CA_ad).eval();
			return 1 - weight * Eigen::sqrt(sym_inner_prod(dist, dist));
		}
	}
	else if (sub.substr(0, 3) == "npr") {
		auto vlist = parse_point_list(sub.substr(3));
		if (vlist.size() == 2) {
			Eigen::Matrix3d R;
			Eigen::Vector3d z = vlist[0].normalized();
			Eigen::Vector3d y = z.cross(vlist[1]).normalized();
			Eigen::Vector3d x = y.cross(z).normalized();
			R << x, y, z;
			Eigen::Matrix3d F11 = R * r3_e_ij(0, 0) * R.transpose(); auto f11 = voigt(F11);
			Eigen::Matrix3d F22 = R * r3_e_ij(1, 1) * R.transpose(); auto f22 = voigt(F22);
			Eigen::Matrix3d F33 = R * r3_e_ij(2, 2) * R.transpose(); auto f33 = voigt(F33);
			Eigen::Matrix2<ADScalar> AtA;
			AtA << f11.dot(CA_ad * f11), f11.dot(CA_ad * f22),
				f22.dot(CA_ad * f11), f22.dot(CA_ad * f22);
			Eigen::Vector2<ADScalar> Atb; Atb << f11.dot(CA_ad * f33), f22.dot(CA_ad * f33);
			Eigen::Vector2<ADScalar> nu2 = AtA.lu().solve(Atb);
			return 1 - nu2.sum() / 2;
		}
	}
	else if (sub == "iso" || sub == "sqrtiso" || sub == "logiso") {
		// Identify isotropic projection
		Eigen::Matrix<double, 6, 6> lamBasis, muBasis;
		lamBasis.setZero(); muBasis.setZero(); 
		lamBasis.block<3, 3>(0, 0).setOnes();
		lamBasis /= std::sqrt(sym_inner_prod(lamBasis, lamBasis));
		muBasis.block<3, 3>(0, 0).diagonal() = Eigen::Vector3d::Ones() * 2;
		muBasis.block<3, 3>(3, 3).diagonal() = Eigen::Vector3d::Ones();
		muBasis /= std::sqrt(sym_inner_prod(muBasis, muBasis));
		Eigen::Matrix2d AtA;
		AtA <<
			sym_inner_prod(lamBasis, lamBasis), sym_inner_prod(lamBasis, muBasis),
			sym_inner_prod(muBasis, lamBasis), sym_inner_prod(muBasis, muBasis);
		Eigen::Vector2d b;
		b[0] = sym_inner_prod(lamBasis, CA); b[1] = sym_inner_prod(muBasis, CA);
		Eigen::Vector2d x = AtA.lu().solve(b);
#if 1
		Eigen::Matrix<double, 6, 6> iso_proj = lamBasis * x[0] + muBasis * x[1];
		//std::cout << "CA = \n" << CA << std::endl;
		//std::cout << "CA_iso = \n" << iso_proj << std::endl;
		Eigen::Matrix<ADScalar, 6, 6> iso_dist = iso_proj - CA_ad;
		if (sub == "iso") {
			return 1. - weight * sym_inner_prod(iso_dist, iso_dist);
		} else if (sub == "sqrtiso") {
			return 1 - weight * Eigen::sqrt(sym_inner_prod(iso_dist, iso_dist));
		} else if (sub == "logiso") {
			return -weight * Eigen::log(5e-2 + sym_inner_prod(iso_dist, iso_dist));
		}
#elif 0
		// compliance projection
		auto Siso = mtl::isotropic_projection_s(CA.inverse().eval(), x[0], x[1]);
		std::cout << "CA = \n" << CA << std::endl;
		std::cout << "Ciso = \n" << Siso.inverse() << std::endl;
		Eigen::Matrix<ADScalar, 6, 6> iso_dist = Siso - CA_ad.inverse();
		if (sub == "iso") {
			return 1. - weight * sym_inner_prod(iso_dist, iso_dist);
		} else if (sub == "sqrtiso") {
			return 1 - weight * Eigen::sqrt(sym_inner_prod(iso_dist, iso_dist));
		} else if (sub == "logiso") {
			return -weight * Eigen::log(5e-2 + sym_inner_prod(iso_dist, iso_dist));
		}
#else
		ADScalar err1 = CA_ad.block<3, 3>(0, 3).squaredNorm() * 2;
		ADScalar err2_a = (CA_ad(0, 1) + CA_ad(1, 2) + CA_ad(0, 2)) / 3;
		ADScalar err2 = 2 * (
			Eigen::pow(CA_ad(0, 1) - err2_a, 2) +
			Eigen::pow(CA_ad(1, 2) - err2_a, 2) +
			Eigen::pow(CA_ad(0, 2) - err2_a, 2));
		err2 = Eigen::pow(Eigen::min(Eigen::min(CA_ad(0, 1), CA_ad(1, 2)), CA_ad(0, 2)) - err2_a, 2);
		ADScalar err3_a = CA_ad.block<3, 3>(0, 0).diagonal().mean();
		ADScalar err3 =
			Eigen::pow(CA_ad(0, 0) - err3_a, 2) +
			Eigen::pow(CA_ad(1, 1) - err3_a, 2) +
			Eigen::pow(CA_ad(2, 2) - err3_a, 2);
		err3 = Eigen::pow(Eigen::min(Eigen::min(CA_ad(0, 0), CA_ad(1, 1)), CA_ad(2, 2)) - err3_a, 2);
		ADScalar err4_a = CA_ad.block<3, 3>(3, 3).diagonal().mean();
		ADScalar err4 =
			Eigen::pow(CA_ad(3, 3) - err4_a, 2) +
			Eigen::pow(CA_ad(4, 4) - err4_a, 2) +
			Eigen::pow(CA_ad(5, 5) - err4_a, 2);
		err4 = Eigen::pow(Eigen::min(Eigen::min(CA_ad(3, 3), CA_ad(4, 4)), CA_ad(5, 5)) - err4_a, 2);
		ADScalar err5 = Eigen::pow((err3_a - err2_a) / 2 - err4_a, 2);
		ADScalar err = Eigen::sqrt(err1 + err2 + err3 + err4 + 0.2 * err5);
		return 1 - weight * err;
#endif
	}
	else if (sub[0] == 'Y') {
		Eigen::Vector3d z; z << sub[1] - '0', sub[2] - '0', sub[3] - '0';
		z.normalize();
		//Eigen::BDCSVD<Eigen::Vector3d> svd(z, Eigen::ComputeFullU);
		//Eigen::Matrix3d zxy = svd.matrixU();
		//zxy.applyOnTheRight(Eigen::PermutationMatrix<3>(Eigen::Vector3d(1, 2, 0)));
		auto sigma_z = voigt_stress((z * z.transpose()).eval());
		auto S_ad = CA_ad.inverse().eval();
		auto eps = (S_ad * sigma_z).eval();
		ADScalar Y_z = 1. / sym_inner_prod(voigt((z * z.transpose()).eval()), eps); // unit z-pressure / eps_z
		return Y_z;
	}
	else if (sub[0] == 'S') {
		auto v3list = parse_point_list(sub.substr(1));
		if (v3list.size() != 2) { std::cout << "invalid objective " << sub << std::endl; throw std::runtime_error("invalid objective"); }
		Eigen::Vector<double, 6> eps_v;
		eps_v << v3list[0], v3list[1];
		return eps_v.dot(CA_ad * eps_v);
	}

	std::cout << "Invalid ADS objective " << subtype_str << std::endl;
	throw std::runtime_error("invalid ads objective");
}

Eigen::AutoDiffScalar<Eigen::Vector<double, 21>> ads_objective(std::string type, const Eigen::Matrix<double, 6, 6>& CA)
{	
	Eigen::Matrix<ADScalar, 6, 6> CA_ad;
	auto [id2ij, ij2id] = elastic_matrix_flatthen_map();
	for (int i = 0; i < 6; i++) {
		for (int j = i; j < 6; j++) {
			CA_ad(i, j).value() = CA(i, j);
			CA_ad(i, j).derivatives() = Eigen::Vector<double, 21>::Unit(ij2id(i, j));
			CA_ad(j, i) = CA_ad(i, j);
		}
	}

	std::vector<std::string> subtypes;
	boost::split(subtypes, type, boost::is_any_of("+"));
	ADScalar f(0); f.derivatives().setZero();
	for (int i = 0; i < subtypes.size(); i++) {
		f += parse_one_ads_objective(subtypes[i], CA_ad, CA);
	}

	return f;
}

Eigen::Vector<double, 21> extract_elastic_matrix_coeff(const Eigen::Matrix<double, 6, 6>& CA)
{
	Eigen::Vector<double, 21> CAcoef;
	int counter = 0;
	for (int i = 0; i < 6; i++) {
		for (int j = i; j < 6; j++) {
			CAcoef[counter] = CA(i, j);
			counter++;
		}
	}
	return CAcoef;
}

std::tuple<std::array<std::pair<int, int>, 21>, Eigen::Matrix<int, 6, 6>> elastic_matrix_flatthen_map(void)
{
	std::array<std::pair<int, int>, 21> id2ij;
	Eigen::Matrix<int, 6, 6> ij2id;

	Eigen::Vector<double, 21> CAcoef;
	int counter = 0;
	for (int i = 0; i < 6; i++) {
		for (int j = i; j < 6; j++) {
			id2ij[counter] = std::make_pair(i, j);
			ij2id(i, j) = counter;
			ij2id(j, i) = counter;
			counter++;
		}
	}
	return std::make_tuple(id2ij, ij2id);
}


std::tuple<std::vector<Eigen::Matrix3d>, std::vector<Compile1ring>> generate_face_frame_and_vring(const PeriodSurfaceMesh& m)
{
	return std::make_tuple(generate_face_frame(m), generate_vring(m));
}

std::vector<Compile1ring> generate_vring(const msf::PeriodSurfaceMesh& m)
{
	std::vector<Compile1ring> vrings(m.n_vertices());
	double Aring = 0;
	
#pragma omp parallel for reduction(+:Aring)
	for (int vid = 0; vid < m.n_vertices(); vid++) {
		auto vh = OM::SmartVertexHandle(vid, &m);
		auto [o, ring] = m.find1ring(vh, 1);
		vrings[vid] = Compile1ring(o, ring); 
		Aring += vrings[vid].As;
	}
	return vrings;
}

std::vector<Eigen::Matrix3d> generate_face_frame(const msf::PeriodSurfaceMesh& m)
{
	std::vector<Eigen::Matrix3d> frlist(m.n_faces());
#pragma omp parallel for
	for (int fid = 0; fid < m.n_faces(); fid++) {
		auto fh = OM::FaceHandle(fid);
		auto tri = m.getFacePeriodVertex(fh, 1);
		frlist[fid] = face_frame(tri);
	}
	return frlist;
}


Eigen::Matrix<double, 6, 6> eval_asym_elastic_tensor(const Eigen::Matrix<double, -1, 6>& fm, const Eigen::Matrix<double, -1, 6>& um, const Eigen::Matrix<double, 6, 6>& Enmem, double As)
{
	Eigen::Matrix<double, 6, 6> CA;

	for (int i = 0; i < 6; i++) {
		for (int j = i; j < 6; j++) {
			CA(i, j) = (Enmem(i, j) - um.col(i).dot(fm.col(j))) / As;
			CA(j, i) = CA(i, j);
		}
	}
	return CA;
}

double search_ads_step(bool enable, PeriodSurfaceMesh& m, const Eigen::VectorXd& step_vector, double obj_last, double gTp, double max_step,
	std::function<Eigen::AutoDiffScalar<Eigen::Vector<double, 21>>(const Eigen::Matrix<double, 6, 6>&)> objfunc) {
	auto [lam, mu] = mtl::lameCoeff(1, 0.3);
	double lam0 = 2 * lam * mu / (lam + 2 * mu);

	double step = maximal_unflip_step(m, step_vector, max_step, M_PI / 3);

	std::vector<OM::Vec3d> v0;
	for (auto vh : m.vertices()) { v0.push_back(m.point(vh)); }

	double c = 0.1;

	for (int iter = 0; iter < 30; iter++) {
		Eigen::VectorXd dk = step_vector * step;
		for (auto vh : m.vertices()) { m.point(vh) = v0[vh.idx()] + toOM(dk.block<3, 1>(vh.idx() * 3, 0)); }
		if (!enable) return step;
		auto [frlist, vrings] = generate_face_frame_and_vring(m);
		auto [fm, um, Enmem, As] = asym_stif_fem(m, frlist, vrings, lam0, mu);
		auto CA = eval_asym_elastic_tensor(fm, um, Enmem, As);
		auto obj = objfunc(CA);
		// we used negative objective, thus it is '>' here 
		if (obj.value() > obj_last - c * step * gTp || iter > 15) {
			//return std::make_tuple(frlist, vrings, fm, um, CA);
			return step;
		}
		step *= 0.7;
	};
	return step;
}



Eigen::VectorXd assemble_mean_curvature_flow(const std::vector<Compile1ring>& vrings) {
	Eigen::VectorXd Hv(vrings.size() * 3);
	for (int i = 0; i < vrings.size(); i++) {
		Hv.block<3, 1>(i * 3, 0) = vrings[i].Lx;
	}
	return Hv;
}


std::ofstream logger(void) {
	return std::ofstream(getPath("log"), std::ios::app); 
}
