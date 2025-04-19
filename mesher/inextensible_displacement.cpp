//#include "solver/ShellFEM.h"
#include "fundamental_forms.h"
#include "PeriodicMesher.h"
#include <Eigen/PardisoSupport>
#include "matlab/matlab_utils.h"
#include "asymptotic_analysis.h"
#include <Eigen/Eigenvalues>
#include "material/materail.h"
//#include "Spectra/SymEigsShiftSolver.h"
//#include "Spectra/GenEigsRealShiftSolver.h"
//#include "Spectra/MatOp/SparseSymMatProd.h"
//#include "Spectra/MatOp/SparseGenMatProd.h"
//#include "Spectra/MatOp/SparseSymShiftSolve.h"
#include <igl/cotmatrix.h>
#include <igl/read_triangle_mesh.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/gaussian_curvature.h>
#include "flow/triangle_quadrature.h"


using namespace msf;

auto translation_basis(const Eigen::VectorXd& Mv) {
	Eigen::MatrixX3d R(Mv.rows() * 3, 3);
	for (int k = 0; k < 3; k++) {
		Eigen::VectorXd Mtri(Mv.rows() * 3, 1); Mtri.setZero();
		for (int i = 0; i < Mv.rows(); i++) { Mtri[i * 3 + k] = std::sqrt(Mv[i]); }
		R.col(k) = Mtri.normalized();
	}
	return R;
}

auto rigid_motion_basis(const Eigen::MatrixX3d& V, const Eigen::VectorXd& Mv) {
	Eigen::Matrix<double, -1, 6> R(V.rows() * 3, 6); R.setZero();
	for (int i = 0; i < V.rows(); i++) {
		Eigen::Vector3d r = V.row(i).transpose();
		Eigen::Matrix3d Rv;
		Rv << 0, -r[2], r[1],
			r[2], 0, -r[0],
			-r[1], r[0], 0;
		Eigen::Matrix3d T; T.setIdentity();
		R.block<3, 3>(i * 3, 0) = T /** Mv[i]*/;
		R.block<3, 3>(i * 3, 3) = Rv /** Mv[i]*/;
	}
	Eigen::VectorXd v_mass = Mv.replicate(1, 3).transpose().reshaped();
	for (int i = 0; i < 6; i++) {
		for (int k = 0; k < i; k++) {
			R.col(i) -= R.col(i).dot(v_mass.cwiseProduct(R.col(k))) * R.col(k);
		}
		R.col(i) /= std::sqrt(R.col(i).dot(v_mass.cwiseProduct(R.col(i))));
	}
	return R;
}

std::tuple<Eigen::VectorXd, double, double> inextensible_displacement(
	PeriodSurfaceMesh& m,
	double lam0, double mu,
	const std::vector<Eigen::Matrix3d>& frlist,
	const std::vector<Compile1ring>& vrings
) {
	auto K = asym_stif_fem_matrix(m, frlist, vrings, lam0, mu);

	//eigen2ConnectedMatlab("K", K);

	//srand((unsigned int) time(0));

	auto [Av, Nv] = eval_vertex_mass(m);
	Eigen::VectorXd Mhalf = Av.cwiseSqrt();
	//eigen2ConnectedMatlab("Mhafl", Mhalf);
	Eigen::MatrixX3d R = translation_basis(Av);
	Eigen::VectorXd u(m.n_vertices() * 3); u.setRandom();
	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> so(K);
	//Mhalf = Mhalf.cwiseInverse().replicate(1, 3).transpose().reshaped().eval();
	Eigen::VectorXd M = Av.replicate(1, 3).transpose().reshaped();
	auto L = m.getPeriodicLaplacian(2, 1);
	auto ulast = u;
	double err;
	for (int iter = 0; iter < 200; iter++) {
		Eigen::VectorXd Pu = u - R * (R.transpose() * u);
		// For H1 space eigenvector
		//Eigen::VectorXd y = M.asDiagonal() * Pu - (L * (Pu.reshaped(3, Pu.size() / 3).transpose())).transpose().reshaped();
		// For L2 space eigenvector
		Eigen::VectorXd y = M.asDiagonal() * Pu;
		//eigen2ConnectedMatlab("y", y);
		u = so.solve(y);
		u /= std::sqrt(u.dot(M.cwiseProduct(u)));
		err = (ulast - u).norm();
		ulast = u;
		//std::cout << "ch = " << err << std::endl;
	}
	{
		//std::ofstream ofs("iextu", std::ios::binary);
		//for (auto vh : m.vertices()) {
		//	auto p = m.point(vh);
		//	ofs.write((const char*)p.data(), sizeof(p));
		//	Eigen::Vector3d uv = u.middleRows(vh.idx() * 3, 3);
		//	ofs.write((const char*)uv.data(), sizeof(uv));
		//}
	}

	double lam_eig = u.dot(K * u);

	return std::make_tuple(u, lam_eig, err);
}

Eigen::Matrix<double, 4, 9> tgnt_displacement_grad_matrix(
	const Eigen::Matrix3d& tri, const Compile1ring v[3],
	const Eigen::Vector3d& n, const Eigen::Matrix<double, 3, 2>& P, bool vertex_proj = true
);

inline double dot2(const Eigen::Matrix2d& m1, const Eigen::Matrix2d& m2) { return m1.cwiseProduct(m2).sum(); }

Eigen::Vector3d abap_sensitivity_element(
	const Eigen::Matrix3d& tri, const Eigen::Matrix3d& fr,
	const Compile1ring v[3],
	double lam0, double mu,
	const Eigen::Vector<double, 9>& ue,
	double lam_min
) {
	Eigen::Vector3d n = fr.col(2);
	double A = n.norm(); n /= A;

	Eigen::Vector3d e1 = fr.col(0), e2 = fr.col(1);

	Eigen::Matrix<double, 3, 2> face_fr; face_fr << e1, e2;

	Eigen::Matrix<double, 3, 9> B0 = tgnt_membrane_strain_displacement_matrix(tri, v, n, face_fr);
	auto Bn = strain_matrix_edge_stretch(tri, e1, e2);
	auto be = second_fundamental_form(tri, v[0], v[1], v[2]);

	Eigen::Vector2d dn[3] = {};
	Eigen::Vector2d dV[3] = {};
	Eigen::Vector<double, 6> utm;
	Eigen::Vector3d unm;
	Eigen::Vector2d dun;
	Eigen::Matrix2d dut;

	utm = (ue.reshaped(3, 3).transpose() * face_fr).transpose().reshaped();
	for (int j = 0; j < 3; j++) { unm[j] = ue.block<3, 1>(j * 3, 0).dot(v[j].nv); }
	dun = dV[0] * unm[0] + dV[1] * unm[1] + dV[2] * unm[2];
	dut = (tgnt_displacement_grad_matrix(tri, v, n, face_fr) * ue).reshaped(2, 2);

	for (int i = 0; i < 3; i++) {
		Eigen::Matrix<double, 3, 2> dVp;
		dVp << tri.col(i) - tri.col((i + 1) % 3), tri.col(i) - tri.col((i + 2) % 3);
		dV[i] = (dVp.transpose() * face_fr).lu().solve(Eigen::Vector2d::Ones());
		dn[i] = dV[i];
	}

	static auto qlist = qd::tri::symq::get(4);

	Eigen::Vector3d Vn; Vn.setZero();

	for (auto qp : qlist) {
		Eigen::Vector3d q(qp.c[0], qp.c[1], qp.c[2]);
		double w = qp.w;
		auto B = B0;
		Eigen::Vector<double, 9> qn; for (int j = 0; j < 3; j++) { qn.block<3, 1>(j * 3, 0) = q[j] * v[j].nv; }
		Eigen::Matrix<double, 3, 9> Bb = Bn * be * qn.transpose();
		// negative?
		Eigen::Matrix2d bform = fromvoigt((Bn * be).eval());
		Eigen::Matrix2d cform = bform * bform;

		double H = bform.trace() / 2;

		B -= Bb;

		Eigen::Vector3d vn = q;

		int counter = 0;
		double u_3 = unm.dot(q);
		Eigen::Vector2d u_t = utm.reshaped(2, 3) * q;
		Eigen::Matrix2d gam = fromvoigt((B * ue).eval());
		// B : gam(u) : gam(u)
		double Bw = (lam0 * (gam.trace() * dot2(gam, bform) + gam.trace() * dot2(gam, bform))
			+ 4 * mu * dot2((gam * gam).eval(), bform));
		// A : gam(u) : gam(u)
		double Aw = H * (lam0 * gam.trace() * gam.trace() + 2 * mu * dot2(gam, gam));
		for (int k = 0; k < 3; k++) {
			Eigen::Matrix2d dgam;
			dgam = bform * u_t * dV[k].transpose() + vn[k] * bform * dut + dun * dV[k].transpose() + u_3 * vn[k] * cform;
			dgam = (dgam + dgam.transpose()).eval() / 2; 
			// A : gam(u) : zeta(u)
			double C = lam0 * gam.trace() * dgam.trace() + 2 * mu * dot2(gam, dgam);
			// u . (vn b uw - u3 Dvn) + u3 (uw . Dvn)
			double D = u_t.dot(vn[k] * bform * u_t - u_3 * dV[k]) + u_3 * u_t.dot(dV[k]);
			// vn H u^2
			double E = vn[k] * H * (u_t.squaredNorm() + u_3 * u_3);
			Vn[k] += 2 * vn[k] * (Bw - Aw) * w;
			Vn[k] += 2 * C * w;
			Vn[k] += 2 * lam_min * D;
			Vn[k] += 2 * lam_min * E;
		}
		counter++;
	}
	return Vn;
}

Eigen::VectorXd abap_sensitivity(
	PeriodSurfaceMesh& m, double lam0, double mu,
	const std::vector<Eigen::Matrix3d>& frlist,
	const std::vector<Compile1ring>& vrings,
	const Eigen::VectorXd& usol,
	double lam_min
) {
	Eigen::VectorXd sens(m.n_vertices()); sens.setZero();

	std::vector<OM::SmartFaceHandle> fhlist;

	for (auto fh : m.faces()) fhlist.emplace_back(fh);

#pragma omp parallel for
	for (int fid = 0; fid < fhlist.size(); fid++) {
		auto fh = fhlist[fid];
		auto tri = m.getFacePeriodVertex(fh, 1);
		auto fvh = m.getFaceVertexHandle(fh);
		auto v = retrieve_element_vertex_data(m, fh, vrings);
		Eigen::Vector<double, 9> ue;
		for (int k = 0; k < 3; k++) ue.block<3, 1>(k * 3, 0) = usol.block<3, 1>(fvh[k].idx() * 3, 0);
		auto fv_sens = abap_sensitivity_element(tri, frlist[fid], v.data(), lam0, mu, ue, lam_min);
#pragma omp critical
		{
			for (int k = 0; k < 3; k++) { sens[fvh[k].idx()] += fv_sens[k]; }
		}
	}

	return sens;
}


void inextensible_displacement_general(std::string mfile) {
	Eigen::MatrixX3d V;
	Eigen::MatrixX3i F;
	igl::read_triangle_mesh(mfile, V, F);
	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);
	Eigen::SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
	Eigen::MatrixX3d Nv;
	igl::per_vertex_normals(V, F, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE, Nv);
	Eigen::MatrixX3d Hv = L * V;
	Eigen::VectorXd Kv;
	igl::gaussian_curvature(V, F, Kv);
	
	std::vector<Eigen::Matrix3d> frlist;
	std::vector<Eigen::Matrix3d> trilist;
	std::vector<Eigen::Matrix2d> blist;
	for (int f = 0; f < F.rows(); f++) {
		Eigen::Matrix3d tri;
		tri << V.row(F(f, 0)), V.row(F(f, 1)), V.row(F(f, 2)); tri.transposeInPlace();
		trilist.emplace_back(tri);
		auto fr = face_frame(tri);
		frlist.emplace_back(fr);
		
		Eigen::Matrix3d nv;
		nv << Nv.row(F(f, 0)), Nv.row(F(f, 1)), Nv.row(F(f, 2)); nv.transposeInPlace();
		auto b = second_fundamental_form_2d(tri, fr, nv);
		blist.emplace_back(b);
	}

	auto [lam, mu] = mtl::lameCoeff(1, 0.3);
	double lam0 = 2 * lam * mu / (lam + 2 * mu);

	Eigen::SparseMatrix<double> K(V.size(), V.size());
	std::vector<Eigen::Triplet<double>> triplist;
	for (auto fid = 0; fid < F.rows(); fid++) {
		Compile1ring v[3] = {};
		Eigen::Vector3i fv = F.row(fid).transpose();
		for (int i = 0; i < 3; i++) {
			v[i].As = M.coeffRef(fv[i], fv[i]);
			v[i].nv = Nv.row(fv[i]).transpose();
		}
		auto tri = trilist[fid];
		auto fr = frlist[fid];
		auto [Ke, fe] = membrane_element_vertex_stif_matrix_vector(lam0, mu, tri, fr, v);
		for (int i = 0; i < 9; i++) {
			int vi = fv[i / 3];
			for (int j = 0; j < 9; j++) {
				int vj = fv[j / 3];
				triplist.emplace_back(vi * 3 + i % 3, vj * 3 + j % 3, Ke(i, j));
			}
		}
	}
	K.setFromTriplets(triplist.begin(), triplist.end());

	//eigen2ConnectedMatlab("K", K);

	Eigen::VectorXd u(V.rows() * 3); u.setRandom();

	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> so(K);

	auto ulast = u;

	Eigen::VectorXd mass = M.diagonal().eval();
	Eigen::VectorXd v_mass = mass.replicate(1, 3).transpose().reshaped();
	auto R = rigid_motion_basis(V, mass);

	//eigen2ConnectedMatlab("mass", mass);
	//eigen2ConnectedMatlab("vmass", v_mass);
	//eigen2ConnectedMatlab("R", R);

	Eigen::VectorXd Pu;

	for (auto iter = 0; iter < 200; iter++) {
		Pu = u - R * (R.transpose() * v_mass.asDiagonal() * u);
		// For H1 space eigenvector
		//Eigen::VectorXd y = M.asDiagonal() * Pu - (L * (Pu.reshaped(3, Pu.size() / 3).transpose())).transpose().reshaped();
		// For L2 space eigenvector
		Eigen::VectorXd y = v_mass.asDiagonal() * Pu;
		//eigen2ConnectedMatlab("y", y);
		u = so.solve(y);
		u.normalize();
		double err = (ulast - u).norm(); ulast = u;
		std::cout << "ch = " << err << std::endl;
	}

	//eigen2ConnectedMatlab("Pu", Pu);
	//eigen2ConnectedMatlab("u", u);

	{
		//std::ofstream ofs("iextu", std::ios::binary);
		//for (int k = 0; k < V.rows(); k++) {
		//	Eigen::Vector3d p = V.row(k).transpose();
		//	ofs.write((const char*)p.data(), sizeof(p));
		//	//Eigen::Vector3d uv = u.middleRows(k * 3, 3);
		//	Eigen::Vector3d uv = Pu.middleRows(k * 3, 3);
		//	ofs.write((const char*)uv.data(), sizeof(uv));
		//}
	}
}

void test_inextensible_displacement(std::string meshfile) {
#if 0
	//std::string mfile = "D:/projects/minisurf/image/siggraph/refine-mesh/psub04.stl";
	//std::string mfile = "D:/projects/minisurf/image/siggraph/thicknesslimit/rfn_0.0400.stl";
	//std::string mfile = "D:/projects/minisurf/image/siggraph/thicknesslimit/rfn_0.0200.stl";
	auto mfile = meshfile;
	PeriodSurfaceMesh m;
	m.readMergePeriodBoundary(mfile);
	inextensible_displacement(m);
#else
	inextensible_displacement_general(meshfile);
#endif
}


