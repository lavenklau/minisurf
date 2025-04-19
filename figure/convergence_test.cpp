#include "mesher/PeriodicMesher.h"
#include "mesher/isosurface_generation.h"
#include "matlab/matlab_utils.h"
#include "cgal/cgal_utils.h"
#include "function_pools.h"
#include "mesher/asymptotic_analysis.h"
#include "material/materail.h"
#include "Eigen/PardisoSupport"

using namespace msf;

void clamp_vertices(Eigen::MatrixX3<Real>& V, Eigen::MatrixX3i& F);


Eigen::Matrix<double, 6, 6> eval_ads_matrix(PeriodSurfaceMesh& mesh, double lam0, double mu) {

	auto [frlist, vrings] = generate_face_frame_and_vring(mesh);

	auto [fm, um, Enmem, As] = asym_stif_fem(mesh, frlist, vrings, lam0, mu);

	auto CA = eval_asym_elastic_tensor(fm, um, Enmem, As);

	return CA;
}


template<typename Vector3>
double circumRadius(const Vector3& p1, const Vector3& p2, const  Vector3& p3) {
	double a2 = (p1 - p2).squaredNorm();
	double b2 = (p1 - p3).squaredNorm();
	double c2 = (p3 - p2).squaredNorm();
	double ab = (p1 - p2).dot(p1 - p3);
	return std::sqrt(a2 * b2 * c2 / 4 / (a2 * b2 - ab * ab));
}

double eval_diameter(PeriodSurfaceMesh& m) {
	double rsum = 0;
	for (auto fh : m.faces()) {
		auto tri = m.getFacePeriodVertex(fh, 1);
		rsum += circumRadius(tri.col(0), tri.col(1), tri.col(2));
	}
	rsum /= m.n_faces();
	return rsum * 2;
}

template<typename F,typename DF>
void project_vertex_to_zero(PeriodSurfaceMesh& mesh, F flev, DF df) {
#pragma omp parallel for
	for (int vid = 0; vid < mesh.n_vertices(); vid++) {
		auto p = toEigen(mesh.points()[vid]);
		auto pz = search_for_zero(flev, p, df(p), 7e-2);
		auto vh = OM::VertexHandle(vid);
		mesh.set_point(vh, toOM(pz));
	}
	return;
}

auto select_objective(std::string obj) {
	return sfn::period::get(obj);

	throw std::runtime_error("unknown objective");
}

void test_asym_convergence(std::string obj) {

	auto [flev, df] = select_objective(obj);

	Eigen::Vector3d pl, pr;
	pl.setConstant(-1); pr.setConstant(1);
	auto [V, F] = isosurf_mesh(128, std::make_pair(pl, pr), [=](double x, double y, double z) {return flev({ x,y,z }); }, 0);
	std::tie(V, F) = cgal_remesh(V, F);
	clamp_vertices(V, F);

	PeriodSurfaceMesh mesh;
	mesh.read(V, F, false, false, false);

	auto vcut = mesh.mergePeriodBoundary();

	mesh.periodic_remesh(20, vcut, 0.02, 0);

	project_vertex_to_zero(mesh, flev, df);

	mesh.period_shift();

	mesh.clamp_period_boundary(2e-4);

	mesh.savePeriodicMesh(getPath("input.obj"), std::set<OM::SmartVertexHandle>{}, 1);

	double h0 = 0.1;

	auto [lam, mu] = mtl::lameCoeff(1, 0.3);
	double lam0 = 2 * lam * mu / (lam + 2 * mu);


	for (int iter = 0; iter < 15; iter++) {
		std::cout << "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = \n";
		double h = h0 * std::pow(0.8, iter);
		mesh.delaunayRemesh(5, h, h * 0.8, 1e-10, 5);
		//mesh.saveUnitCell(getPath("rem.obj"));

		project_vertex_to_zero(mesh, flev, df);

		//mesh.saveUnitCell(getPath("proj.obj"));

		auto CA = eval_ads_matrix(mesh, lam0, mu);

		double hm = eval_diameter(mesh);

		std::cout << "h = " << hm << std::endl;
		std::cout << "CA = \n" << CA << std::endl;

		std::ofstream ofs(getPath("conv.log"), std::ios::app);
		ofs << "iter " << iter << std::endl;
		ofs << "h = " << hm << std::endl;
		ofs << "CA = \n" << CA << std::endl;
	}
}


