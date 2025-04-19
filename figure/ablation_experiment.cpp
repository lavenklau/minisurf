#include "mesher/asymptotic_analysis.h"
#include "Eigen/PardisoSupport"
#include "material/materail.h"

using namespace msf;

void abalation_ads_willmore_remesh(void) {
	std::string mfile = "D:/projects/minisurf/image/siggraph/abalation_study/adc_willm_rem/P0.5-31.stl";
	PeriodSurfaceMesh m;
	m.readMergePeriodBoundary(mfile);

	auto [lam, mu] = mtl::lameCoeff(1, 0.3);
	double lam0 = 2 * lam * mu / (lam + 2 * mu);

	bool use_remesh = true;
	bool use_willmore_energy = false;
		
	std::vector<double> obj_rec;
	std::vector<double> time_rec;

	double time = 0;

	for (int iter = 0; iter < 300; iter++) {
		if (use_remesh) { m.delaunayRemesh(5, 0.1, 0.02, 0.6); }
	
		if (iter % 3 == 0) m.saveUnitCell("iter.obj");

		auto [frlist, vrings] = generate_face_frame_and_vring(m);

		auto [fm, um, Enmem, As] = asym_stif_fem(m, frlist, vrings, lam0, mu);

		auto [Av, Nv] = eval_vertex_mass(m);
		
		Eigen::Matrix<double, 6, 6> CA;

		for (int i = 0; i < 6; i++) {
			for (int j = i; j < 6; j++) {
				CA(i, j) = (Enmem(i, j) - um.col(i).dot(fm.col(j))) / As;
				CA(j, i) = CA(i, j);
			}
		}

		// Solve L2 gradient
		auto [v_sens, A_sens] = asym_stif_shape_sensitivity(m, frlist, vrings, um, fm, lam0, mu);
		v_sens.applyOnTheLeft(Av.cwiseInverse().asDiagonal()); A_sens.applyOnTheLeft(Av.cwiseInverse().asDiagonal());

		// Solve objective gradient
		auto obj = ads_objective("11", CA); auto CAwise = extract_elastic_matrix_coeff(CA);
		Eigen::VectorXd dfdvn = (v_sens / As - A_sens * CAwise.transpose() / As) * -obj.derivatives();

		// Precondition
		double c = 1;
		auto L = m.getPeriodicLaplacian(2, 1);
		Eigen::SparseMatrix<double> G = -c * L; G += Av.asDiagonal();
		Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> so(G);
		Eigen::VectorXd dn = -so.solve(Av.cwiseProduct(dfdvn).eval());

		// regularization
		double w_willmore = 0.2;
		Eigen::VectorXd willdn(m.n_vertices());
		auto g_willmore = willmore_H2_gradient(m);
		//for (int k = 0; k < m.n_vertices(); k++) willdn[k] = Nv[k].dot(g_willmore.row(k).transpose());
		g_willmore *= -1;

		// time step
		double dt = 3;

		for (auto vh : m.vertices()) {
			auto p = m.point(vh);
			m.point(vh) += dt * toOM((dn[vh.idx()] * Nv[vh.idx()]).eval());
			if (use_willmore_energy) m.point(vh) += dt * w_willmore * toOM(g_willmore.row(vh.idx()));
		}

		std::cout << "obj = " << obj.value() << ", iter = " << iter << std::endl;

		obj_rec.push_back(obj.value());
		time_rec.push_back(time);
		time += dt;
	}
}

