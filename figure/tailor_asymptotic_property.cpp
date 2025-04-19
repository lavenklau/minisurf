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

void detectSurgery(PeriodSurfaceMesh& m, int iter) {
	std::cout << "Checking if surgery ...";
	auto mold = m;
	if (m.surgery()) {
		if (log_detail > 1) mold.saveUnitCell(getPath(fmt::format("befsur{:04d}.obj", iter)));
		if (log_detail > 1) m.saveUnitCell(getPath(fmt::format("aftsur{:04d}.obj", iter)));
		std::cout << "Done" << std::endl; 
	}
	else { std::cout << "No neck" << std::endl; }

}


void log_fem_solution(int iter, const Eigen::Matrix<double, -1, 6>& fm, const Eigen::Matrix<double, -1, 6>& um, const Eigen::Matrix<double, 6, 6>& enmem) {
	logger() << "Enmem = \n" << enmem << std::endl;
	std::ofstream ofs(getPath(fmt::format("fm{:04d}", iter)), std::ios::binary);
	ofs.write((const char*)fm.data(), fm.size() * sizeof(double)); ofs.close();
	ofs.open(getPath(fmt::format("um{:04d}", iter)), std::ios::binary);
	ofs.write((const char*)um.data(), um.size() * sizeof(double)); ofs.close();
}

void tailor_ads(PeriodSurfaceMesh& m, std::string obj_type) {
	auto [lam, mu] = mtl::lameCoeff(1, 0.3);
	double lam0 = 2 * lam * mu / (lam + 2 * mu);
		
	std::vector<double> obj_rec;
	std::vector<double> time_rec;

	CheckConverge conv(converge_tol, step_tol * max_time_step, precondition_strength, converge_sc);

	double time = 0;

	if (obj_type.empty()) { std::cout << "empty objective\n"; return; }

	Eigen::VectorXd step_vector(m.n_vertices() * 3);

	for (int iter = 0; iter < max_iter; iter++) {
		if (iter % 4 == 0) { if (!asym_no_surgery) detectSurgery(m, iter); }

		if (!asym_no_remesh) m.delaunayRemesh(delaunay_remesh_inner_iter, 0.2, 0.02, delaunay_remesh_outer_iter);
	
		if (iter % 50 == 0) m.saveUnitCell(getPath("iter.obj"));
		if (log_detail > 2) m.saveUnitCell(getPath(fmt::format("iter{:04d}.obj", iter)));

		auto [frlist, vrings] = generate_face_frame_and_vring(m);

		auto [fm, um, Enmem, As] = asym_stif_fem(m, frlist, vrings, lam0, mu);

		if (log_detail > 3) { log_fem_solution(iter, fm, um, Enmem); }

		auto [Av, Nv] = eval_vertex_mass(m); Eigen::MatrixXd nmat = Eigen::MatrixXd::Map((const double*)Nv.data(), 3, m.n_vertices());
		
		auto CA = eval_asym_elastic_tensor(fm, um, Enmem, As);

		if (iter % 50 == 0) { std::cout << "CA = \n" << CA << std::endl; }

		if (log_detail > 0) { std::ofstream ofs(getPath(fmt::format("CA")), std::ios::app); ofs << "iter = " << iter << ", CA = " << CA.reshaped().transpose() << std::endl; }

		// Compute objective
		auto obj = ads_objective(obj_type, CA); auto CAwise = extract_elastic_matrix_coeff(CA);

		if (std::isnan(obj.value())) break;

		// Solve L2 gradient
		auto [v_sens, A_sens] = asym_stif_shape_sensitivity(m, frlist, vrings, um, fm, lam0, mu);
		v_sens.applyOnTheLeft(Av.cwiseInverse().asDiagonal()); A_sens.applyOnTheLeft(Av.cwiseInverse().asDiagonal());

		// Solve objective gradient
		Eigen::VectorXd dfdvn = (v_sens / As - A_sens * CAwise.transpose() / As) * -obj.derivatives();

		// Precondition
		double c = precondition_strength;
		c = conv.estimate_precondition(20);
		auto L = m.getPeriodicLaplacian(2, 1);
		Eigen::SparseMatrix<double> G = -c * L; G += Av.asDiagonal();
		Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> so(G);
		double fweight = c + 1;
#if 1
		Eigen::VectorXd dn = -so.solve(Av.cwiseProduct(fweight * dfdvn).eval());
		Eigen::VectorXd dv = (nmat * dn.asDiagonal()).reshaped();
#else
		Eigen::VectorXd dv = -fweight * so.solve((Av.cwiseProduct(dfdvn).asDiagonal() * nmat.transpose()).eval()).transpose().reshaped();
#endif
		//{
		//	std::ofstream ofs(getPath("dv"), std::ios::binary);
		//	for (auto vh : m.vertices()) {
		//		auto p = m.point(vh);
		//		ofs.write((const char*)p.data(), sizeof(p));
		//		Eigen::Vector3d vN = dfdvn[vh.idx()] * nmat.col(vh.idx());
		//		ofs.write((const char*)vN.data(), sizeof(p));
		//		ofs.write((const char*)&dv[vh.idx() * 3], sizeof(p));
		//	}
		//}

		// regularization weight
		double w_obj = 1, w_willmore = weight_willmore, w_area = weight_area;

		// Willmore energy
		Eigen::VectorXd willdn(m.n_vertices());
		Eigen::MatrixX3d g_willmore(m.n_vertices(), 3); g_willmore.setZero();
		if (w_willmore != 0) { g_willmore = willmore_H2_gradient(m); g_willmore *= -1; }

		// area energy
		Eigen::VectorXd g_area(dv.rows()); g_area.setZero();
		if (w_area != 0) g_area = assemble_mean_curvature_flow(vrings);

		if (log_detail > 3) log_sensitivity(getPath(fmt::format("sens{:04d}", iter)), m, dfdvn, g_willmore);

		// backtracking linesearch
		step_vector = w_obj * dv + w_willmore * g_willmore.transpose().reshaped() + w_area * g_area;
		auto gTp = expect_descent(Nv, step_vector, dfdvn, Av);
		double tbar = conv.estimate_next_step(max_time_step);
		double step = search_ads_step(!disable_line_search, m, step_vector, obj.value(), gTp, tbar, [=](const Eigen::Matrix<double, 6, 6>& CA) {return ads_objective(obj_type, CA); });
		std::cout << "obj = " << obj.value() << ", iter = " << iter <<", c = " << c << ", tbar = " << tbar << ", step = " << step << std::endl;
		logger() << fmt::format("iter = {}, obj = {}, c = {}, tbar = {}, step = {}\n", iter, obj.value(), c, tbar, step);

		// check convergence
		if (conv(obj.value(), step)) break;

		obj_rec.push_back(obj.value());
		time_rec.push_back(w_obj * step);
	}

	{ std::ofstream ofs(getPath("objrec")); for (auto ob : obj_rec) { ofs << ob << std::endl; } }
	{ std::ofstream ofs(getPath("timerec")); for (auto ti : time_rec) { ofs << ti << std::endl; } }
	m.saveUnitCell(getPath("final.obj"));
}

void tailor_ads(std::string meshfile, std::string obj_type) {
	PeriodSurfaceMesh m;
	m.readMergePeriodBoundary(meshfile);
	tailor_ads(m, obj_type);
}


