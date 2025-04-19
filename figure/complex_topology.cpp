#include "mesher/asymptotic_analysis.h"
#include "Eigen/PardisoSupport"
#include "material/materail.h"
#include "fmt/core.h"
#include "function_pools.h"
#include "boost/algorithm/string.hpp"
#include "mesher/isosurface_generation.h"
#include "cgal/cgal_utils.h"

using namespace msf;

auto dup_period_surface(std::string fname, int dup) {
	auto [f, df] = msf::sfn::period::get(fname);

	auto fk = [=](const Eigen::Vector3d& p) {
		return f(dup * p);
		};

	auto dfk = [=](const Eigen::Vector3d& p) {
		return dup * df(dup * p);
		};
	return std::make_tuple(fk, dfk);
}


void clamp_vertices(Eigen::MatrixX3d& V, Eigen::MatrixX3i& F);

template<typename F,typename DF>
inline void project_vertex_to_zero(PeriodSurfaceMesh& mesh, F flev, DF df) {
#pragma omp parallel for
	for (int vid = 0; vid < mesh.n_vertices(); vid++) {
		auto p = toEigen(mesh.points()[vid]);
		auto pz = search_for_zero(flev, p, df(p), 2e-2);
		auto vh = OM::VertexHandle(vid);
		mesh.set_point(vh, toOM(pz));
	}
	return;
}



void tailor_ads(PeriodSurfaceMesh& m, std::string obj_type);

void test_complex_topology(std::string obj_type) {
	std::vector<std::string> segs;
	boost::split(segs, obj_type, boost::is_any_of("@"));
	auto fname = segs.at(0);
	auto cors = segs.at(1);
	auto obj = segs.at(2);

	int dup = 1;
	if (std::isdigit(fname[0])) {
		dup = fname[0] - '0'; fname = fname.substr(1);
	}
	auto [f, df] = dup_period_surface(fname, dup);

	Eigen::Vector3d pl, pr;
	pl.setConstant(-1); pr.setConstant(1);
	auto [V, F] = isosurf_mesh(128, std::make_pair(pl, pr), [=](double x, double y, double z) {return f({ x,y,z }); }, 0);
	std::tie(V, F) = cgal_remesh(V, F);
	clamp_vertices(V, F);

	PeriodSurfaceMesh mesh;
	mesh.read(V, F, false, false, false);

	auto vcut = mesh.mergePeriodBoundary();

	mesh.periodic_remesh(20, vcut, 0.02, 0);

	project_vertex_to_zero(mesh, f, df);

	mesh.period_shift();

	mesh.clamp_period_boundary(2e-4);

	mesh.savePeriodicMesh(getPath("input.obj"), std::set<OM::SmartVertexHandle>{}, 1);

	if (cors == "ads") {
		tailor_ads(mesh, obj);
	} else {
		throw std::runtime_error("unknown cors type");
	}
	
	
}
