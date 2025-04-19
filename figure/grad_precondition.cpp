//#include "solver/ShellFEM.h"
#include "mesher/fundamental_forms.h"
#include "mesher/PeriodicMesher.h"
#include <Eigen/PardisoSupport>
#include <unsupported/Eigen/AutoDiff>
#include "matlab/matlab_utils.h"
#include "material/materail.h"
#include "function_pools.h"
#include "mesher/asymptotic_analysis.h"
#include "fmt/core.h"

using namespace msf;

extern std::vector<Real> aux_number;


void test_grad_precondition(void) {

	double Y = 1, nu = 0.3;

	auto [lam, mu] = mtl::lameCoeff(Y, nu);
	double lam0 = 2 * lam * mu / (lam + 2 * mu);

	std::cout << "lam = " << lam << ", lam0 = " << lam0 << ", mu = " << mu << std::endl;

	PeriodSurfaceMesh m;
	std::string mfile = "D:/projects/minisurf/image/siggraph/precondition/IWP0.5-23.stl";
	m.readMergePeriodBoundary(mfile);
	int ndof = m.n_vertices() * 3;

	auto [frlist, vrings] = generate_face_frame_and_vring(m);

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

	m.savePeriodicMesh("temp.obj", std::vector<OM::SmartVertexHandle>{}, 1.);

	{
		std::ofstream ofs("Vp", std::ios::binary);
		for (auto vh : m.vertices()) {
			auto p = m.point(vh); ofs.write((const char*)p.data(), sizeof(p));
			auto n = vrings[vh.idx()].nv; ofs.write((const char*)n.data(), sizeof(n));
		}
	}

	auto [fm, um, Enmem, As] = asym_stif_fem(m, frlist, vrings, lam0, mu);

	{
		//std::ofstream ofs("um", std::ios::binary);
		//Eigen::VectorXd U = um * voigt(etgt);
		//ofs.write((const char*)U.data(), U.size() * sizeof(double)); ofs.close();

		//ofs.open("fm", std::ios::binary);
		//Eigen::VectorXd F = fm * voigt(etgt);
		//ofs.write((const char*)F.data(), F.size() * sizeof(double)); ofs.close();
	}

	Eigen::Matrix<double, 6, 6> CA;

	for (int i = 0; i < 6; i++) {
		for (int j = i; j < 6; j++) {
			CA(i, j) = (Enmem(i, j) - um.col(i).dot(fm.col(j))) / As;
			CA(j, i) = CA(i, j);
		}
	}

	std::cout << "Em = \n" << Enmem / As << std::endl;
	std::cout << "Ea = \n" << CA << std::endl;

	auto [Vn_ij, Vn_A] = asym_stif_shape_sensitivity(m, frlist, vrings, um, fm, lam0, mu);
	auto [Av, Nv] = eval_vertex_mass(m);
	Vn_ij.applyOnTheLeft(Av.cwiseInverse().asDiagonal()); Vn_A.applyOnTheLeft(Av.cwiseInverse().asDiagonal());

	double f = -Evgt.dot(CA * Evgt);

	std::cout << "*obj = " << f << std::endl;

	Eigen::VectorXd Vn_E = Vn_ij * dfdC;


	Eigen::VectorXd Vn(m.n_vertices());

	Vn = Vn_E / As - Vn_A * f / As; // f already divides As

	{
		std::ofstream ofs("Vnraw", std::ios::binary);
		ofs.write((const char*)Vn.data(), Vn.size() * sizeof(double));
	}

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

	{
		std::ofstream ofs("Vnpre", std::ios::binary);
		ofs.write((const char*)Vn.data(), Vn.size() * sizeof(double));
	}

	//Eigen::VectorXd VnV(m.n_vertices());
	//Eigen::VectorXd VnH(m.n_vertices());
	//double step_E = 5;
	//double step_H = 0.3; // 0.5
	//if (aux_number.size() > 1) step_E = aux_number[1];
	//if (aux_number.size() > 2) step_H = aux_number[2];
	//for (auto vh : m.vertices()) {
	//	VertexFlag vflag;
	//	auto pv = m.point(vh);
	//	vflag.set_period_boundary(1 - 1e-5, pv[0], pv[1], pv[2]);
	//	Eigen::Vector3d vn_e = -Vn[vh.idx()] * vrings[vh.idx()].nv;
	//	Eigen::Vector3d vn_h = vrings[vh.idx()].nv * vrings[vh.idx()].H * vrings[vh.idx()].As;
	//	Eigen::Vector3d vn = vn_e * step_E + vn_h * step_H;
	//	VnH[vh.idx()] = vrings[vh.idx()].H * vrings[vh.idx()].As * step_H;
	//	VnV[vh.idx()] = vn.dot(vrings[vh.idx()].nv);
	//	for (int k = 0; k < 3; k++) {
	//		if (!vflag.is_period_boundary(k)) { pv[k] += vn[k]; }
	//	}
	//	m.set_point(vh, toOM(pv));
	//}
}



extern int max_iter;
extern double converge_tol;
extern double step_tol;
extern double precondition_strength;
extern bool asym_no_remesh;
extern bool asym_no_surgery;
extern double weight_area;
extern int delaunay_remesh_outer_iter;
extern int delaunay_remesh_inner_iter;

void test_stiffness_matrix_precondition() {
	std::vector<double> obj_rec, time_rec;
	PeriodSurfaceMesh m;
	std::string obj_type = "c33";
	
	m.readMergePeriodBoundary("D:/projects/minisurf/image/thermal/SMdata/non-TPMS/samples/D0.1-2.stl");

	auto [lam, mu] = mtl::lameCoeff(1, 0.3);
	double lam0 = 2 * lam * mu / (lam + 2 * mu);
		

	CheckConverge conv(converge_tol, step_tol * max_time_step, precondition_strength, 0);

	double time = 0;

	if (obj_type.empty()) { std::cout << "empty objective\n"; return; }

	Eigen::VectorXd step_vector(m.n_vertices() * 3);

	for (int iter = 0; iter < max_iter; iter++) {
		if (iter % 4 == 0) { if (!asym_no_surgery) m.surgery(); }

		if (!asym_no_remesh) m.delaunayRemesh(delaunay_remesh_inner_iter, 0.2, 0.02, delaunay_remesh_outer_iter);
	
		if (iter % 50 == 0) m.saveUnitCell(getPath("iter.obj"));
		if (log_detail > 2) m.saveUnitCell(getPath(fmt::format("iter{:04d}.obj", iter)));

		auto [frlist, vrings] = generate_face_frame_and_vring(m);

		//auto [fm, um, Enmem, As] = asym_stif_fem(m, frlist, vrings, lam0, mu);
		auto [Ks, fm, Enmem, As] = assembly_stif_fem(m, frlist, vrings, lam0, mu);
		Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> so(Ks);
		auto um = so.solve(fm).eval();

		auto [Av, Nv] = eval_vertex_mass(m); Eigen::MatrixXd nmat = Eigen::MatrixXd::Map((const double*)Nv.data(), 3, m.n_vertices());
		
		auto CA = eval_asym_elastic_tensor(fm, um, Enmem, As);

		if (iter % 50 == 0) { std::cout << "CA = \n" << CA << std::endl; }

		if (log_detail > 0) { std::ofstream ofs(getPath(fmt::format("CA")), std::ios::app); ofs << "iter = " << iter << ", CA = " << CA.reshaped().transpose() << std::endl; }

		// Compute objective
		auto obj = ads_objective(obj_type, CA); auto CAwise = extract_elastic_matrix_coeff(CA);

		// Solve L2 gradient
		auto [v_sens, A_sens] = asym_stif_shape_sensitivity(m, frlist, vrings, um, fm, lam0, mu);
		v_sens.applyOnTheLeft(Av.cwiseInverse().asDiagonal()); A_sens.applyOnTheLeft(Av.cwiseInverse().asDiagonal());

		// Solve objective gradient
		Eigen::VectorXd dfdvn = (v_sens / As - A_sens * CAwise.transpose() / As) * -obj.derivatives();

		// Precondition
		double c = precondition_strength;
		double fweight = c + 1;
#if 0
		auto L = m.getPeriodicLaplacian(2, 1);
		Eigen::SparseMatrix<double> G = -c * L; G += Av.asDiagonal();
		Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> so(G);
		Eigen::VectorXd dn = -so.solve(Av.cwiseProduct(fweight * dfdvn).eval());
		Eigen::VectorXd dv = (nmat * dn.asDiagonal()).reshaped();
#else
		Eigen::VectorXd dv = fweight * (nmat * dfdvn.cwiseProduct(Av).asDiagonal()).reshaped();
		{ std::ofstream ofs(getPath("Vnbef")); ofs << dv; }
		eigen2ConnectedMatlab("Ks", Ks);
		eigen2ConnectedMatlab("dv", dv);
		Ks *= c;
		for (int k = 0; k < Ks.rows(); k += 3) {
			for (int i = 0; i < 3; i++) {
				Ks.coeffRef(k + i, k + i) += Av[k / 3];
			}
		}
		so.compute(Ks);
		dv = so.solve(dv).eval();
		{ std::ofstream ofs(getPath("Vp")); for (auto vh : m.vertices()) { ofs << m.point(vh) << std::endl; } }
		{ std::ofstream ofs(getPath("Vnaft")); ofs << dv; }
#endif

		// regularization weight
		double w_obj = 1, w_willmore = weight_willmore, w_area = weight_area;

		// Willmore energy
		Eigen::VectorXd willdn(m.n_vertices());
		Eigen::MatrixX3d g_willmore(m.n_vertices(), 3); g_willmore.setZero();
		if (w_willmore != 0) { g_willmore = willmore_H2_gradient(m); g_willmore *= -1; }

		// area energy
		auto g_area = assemble_mean_curvature_flow(vrings);

		if (log_detail > 3) log_sensitivity(getPath(fmt::format("sens{:04d}", iter)), m, dfdvn, g_willmore);

		// backtracking linesearch
		step_vector = w_obj * dv + w_willmore * g_willmore.transpose().reshaped() + w_area * g_area;
		auto gTp = expect_descent(Nv, step_vector, dfdvn, Av);
		double tbar = conv.estimate_next_step(max_time_step);
		double step = search_ads_step(!disable_line_search, m, step_vector, obj.value(), gTp, tbar, [=](const Eigen::Matrix<double, 6, 6>& CA) {return ads_objective(obj_type, CA); });
		std::cout << "obj = " << obj.value() << ", iter = " << iter << ", step = " << step << std::endl;

		// check convergence
		if (conv(obj.value(), step)) break;

		obj_rec.push_back(obj.value());
		time_rec.push_back(w_obj * step);
	}

	{ std::ofstream ofs(getPath("objrec")); for (auto ob : obj_rec) { ofs << ob << std::endl; } }
	{ std::ofstream ofs(getPath("timerec")); for (auto ti : time_rec) { ofs << ti << std::endl; } }
	m.saveUnitCell(getPath("final.obj"));

}

void plot_precondition_gradient(std::string mfile, std::string obj_type) {
	PeriodSurfaceMesh m;
	m.readMergePeriodBoundary(mfile, false);

	auto [lam, mu] = mtl::lameCoeff(1, 0.3);
	double lam0 = 2 * lam * mu / (lam + 2 * mu);

	CheckConverge conv(converge_tol, step_tol * max_time_step, precondition_strength, 0);

	double time = 0;

	if (obj_type.empty()) { std::cout << "empty objective\n"; return; }

	m.delaunayRemesh(delaunay_remesh_inner_iter, 0.2, 0.02, delaunay_remesh_outer_iter);

	m.saveUnitCell(getPath("iter.obj"));

	auto [frlist, vrings] = generate_face_frame_and_vring(m);

	auto [fm, um, Enmem, As] = asym_stif_fem(m, frlist, vrings, lam0, mu);
	//auto [Ks, fm, Enmem, As] = assembly_stif_fem(m, frlist, vrings, lam0, mu);
	//Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> so(Ks);
	//auto um = so.solve(fm).eval();

	auto [Av, Nv] = eval_vertex_mass(m); Eigen::MatrixXd nmat = Eigen::MatrixXd::Map((const double*)Nv.data(), 3, m.n_vertices());

	auto CA = eval_asym_elastic_tensor(fm, um, Enmem, As);

	// Compute objective
	auto obj = ads_objective(obj_type, CA); auto CAwise = extract_elastic_matrix_coeff(CA);

	// Solve L2 gradient
	auto [v_sens, A_sens] = asym_stif_shape_sensitivity(m, frlist, vrings, um, fm, lam0, mu);
	v_sens.applyOnTheLeft(Av.cwiseInverse().asDiagonal()); A_sens.applyOnTheLeft(Av.cwiseInverse().asDiagonal());

	// Solve objective gradient
	Eigen::VectorXd dfdvn = (v_sens / As - A_sens * CAwise.transpose() / As) * -obj.derivatives();

	// Precondition
	double c = precondition_strength;
	double fweight = c + 1;
	c = conv.estimate_precondition(20);
	auto L = m.getPeriodicLaplacian(2, 1);
	Eigen::SparseMatrix<double> G = -c * L; G += Av.asDiagonal();
	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> so(G);
#if 1
	Eigen::VectorXd dn = -so.solve(Av.cwiseProduct(fweight * dfdvn).eval());
	Eigen::VectorXd dv = (nmat * dn.asDiagonal()).reshaped();
#else
	Eigen::VectorXd dv = -fweight * so.solve((Av.cwiseProduct(dfdvn).asDiagonal() * nmat.transpose()).eval()).transpose().reshaped();
#endif
	{
		std::ofstream ofs(getPath("dv"), std::ios::binary);
		for (auto vh : m.vertices()) {
			auto p = m.point(vh);
			ofs.write((const char*)p.data(), sizeof(p));
			Eigen::Vector3d vN = dfdvn[vh.idx()] * nmat.col(vh.idx());
			ofs.write((const char*)vN.data(), sizeof(p));
			ofs.write((const char*)&dv[vh.idx() * 3], sizeof(p));
		}
	}
	
	
}
