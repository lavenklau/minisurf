//#define CGAL_USE_BASIC_VIEWER
//#include "cgal_utils.h"
#include <iostream>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Linear_cell_complex_incremental_builder_3.h>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/draw_face_graph_with_paths.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/Random.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_3_to_lcc.h>
#include <CGAL/Triangulation_2_to_lcc.h>
typedef CGAL::Simple_cartesian<double>               Kernel;
typedef Kernel::Point_3                              Point_3;
typedef CGAL::Polyhedron_3<Kernel>                   Polyhedron;
typedef Polyhedron::Facet_iterator                   Facet_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;

using LCC_3 = CGAL::Linear_cell_complex_for_combinatorial_map<2, 3>;
using IB = CGAL::Linear_cell_complex_incremental_builder_3<LCC_3>;
using Path_on_surface = CGAL::Surface_mesh_topology::Path_on_surface<LCC_3>;
double cycle_length(const LCC_3& lcc, const ::Path_on_surface& cycle)
{ // Compute the length of the given cycle.
	double res = 0;
	for (std::size_t i = 0; i < cycle.length(); ++i)
	{
		res += std::sqrt
		(CGAL::squared_distance(lcc.point(cycle[i]),
			lcc.point(lcc.other_extremity(cycle[i]))));
	}
	return res;
}
void display_cycle_info(const LCC_3& lcc, const ::Path_on_surface& cycle)
{ // Display information about the given cycle.
	if (cycle.is_empty()) { std::cout << "Empty." << std::endl; return; }
	std::cout << "Root: " << lcc.point(cycle[0]) << "; "
		<< "Number of edges: " << cycle.length() << "; "
		<< "Length: " << cycle_length(lcc, cycle) << std::endl;
}
int cgal_main(int argc, char* argv[])
{
	std::string filename(argc == 1 ? CGAL::data_file_path("meshes/3torus.off") : argv[1]);
	bool draw = (argc < 3 ? false : std::string(argv[2]) == "-draw");
	LCC_3 lcc;
	if (!CGAL::load_off(lcc, filename.c_str())) // Load the off file.
	{
		std::cout << "Cannot read file '" << filename << "'. Exiting program" << std::endl;
		return EXIT_FAILURE;
	}
	std::cout << "File '" << filename << "' loaded. Finding shortest non contractible cycle..." << std::endl;
	CGAL::Surface_mesh_topology::Curves_on_surface_topology<LCC_3> cst(lcc);
	LCC_3::Dart_const_descriptor root = lcc.dart_descriptor
	(CGAL::get_default_random().get_int(0, static_cast<int>(lcc.number_of_darts()))); // One dart of the mesh
	::Path_on_surface cycle1 =
		cst.compute_shortest_non_contractible_cycle_with_base_point(root);
	CGAL::Surface_mesh_topology::Euclidean_length_weight_functor<LCC_3> wf(lcc);
	::Path_on_surface cycle2 =
		cst.compute_shortest_non_contractible_cycle_with_base_point(root, wf);
	std::cout << "Cycle 1 (pink): "; display_cycle_info(lcc, cycle1);
	std::cout << "Cycle 2 (green): "; display_cycle_info(lcc, cycle2);
	if (draw) { CGAL::draw(lcc, { cycle1, cycle2 }); }
	return EXIT_SUCCESS;
}


std::vector<std::vector<size_t>> extract_shortest_homotopy_path(const Eigen::Matrix<double, -1, 3>& V, const Eigen::Matrix<int, -1, 3>& F, const std::vector<int>& baselist)
{
	LCC_3 lcc;
	IB ib(lcc);
	std::map<LCC_3::Point, int> p2id;
	for (int i = 0; i < V.rows(); i++) {
		ib.add_vertex(LCC_3::Point(V(i, 0), V(i, 1), V(i, 2)));
		p2id[LCC_3::Point(V(i, 0), V(i, 1), V(i, 2))] = i;
	}
	ib.begin_surface();
	for (int i = 0; i < F.rows(); i++) {
		ib.begin_facet();
		ib.add_vertex_to_facet(F(i, 0));
		ib.add_vertex_to_facet(F(i, 1));
		ib.add_vertex_to_facet(F(i, 2));
		ib.end_facet();
	}
	ib.end_surface();
	
	CGAL::Surface_mesh_topology::Curves_on_surface_topology<LCC_3> cst(lcc);
	std::vector<std::vector<size_t>> cyclelist(baselist.size());
	for (int i = 0; i < baselist.size(); i++) {
		// find dart for basepont
		LCC_3::Dart_const_descriptor root = lcc.dart_descriptor(baselist[i]); // One dart of the mesh
		//std::cout << "Initial dart = \n";
		//lcc.display_dart(root);
		for (size_t k = 0; k < lcc.number_of_darts(); k++) {
			auto dar_descr = lcc.dart_descriptor(k);
			if (lcc.is_free<2>(dar_descr))
				continue;
			if (lcc.point(dar_descr) == LCC_3::Point(V(baselist[i], 0), V(baselist[i], 1), V(baselist[i], 2))) {
				root = dar_descr;
				// lcc.display_dart(dar_descr);
				//std::cout << "Found a dart\n";
				break;
			}
		}
		CGAL::Surface_mesh_topology::Euclidean_length_weight_functor<LCC_3> wf(lcc);
		::Path_on_surface cycle2 =
			cst.compute_shortest_non_contractible_cycle_with_base_point(root, wf);
		display_cycle_info(lcc, cycle2);
		//CGAL::draw(lcc, cycle2);
		// std::vector<LCC_3::Point> path_vlist;
		std::vector<LCC_3::Point> clist;
		for (int k = 0; k < cycle2.length(); k++) {
			auto p1 = lcc.point(cycle2[k]);
			auto p2 = lcc.point(lcc.other_extremity(cycle2[k]));
			if (k == 0) {
				clist.emplace_back(p1);
				clist.emplace_back(p2);
			} else if (k == 1) {
				if (p1 == clist[0] || p1 == clist[1]) { 
					clist.emplace_back(p2);
				} else {
					clist.emplace_back(p1);
				}
			} else if (k == cycle2.length() - 1) {
				if (p1 != clist[0] && p2 != clist[0]) {
					std::swap(clist[0], clist[1]);
				}
			} else {
				if (p1 == *clist.rbegin()) {
					clist.emplace_back(p2);
				} else {
					clist.emplace_back(p1);
				}
			}
		}
		for (int k = 0; k < clist.size(); k++) {
			cyclelist[i].push_back(p2id[clist[k]]);
		}
	}
	return cyclelist;
}

using namespace CGAL;
using namespace CGAL::Surface_mesh_topology;
// Euclidean distance weight functor: each edge has as weight its Euclidean length.
template<typename Mesh>
struct FlatTorus_length_weight_functor
{
	using Weight_t = double;
	using Dart_const_descriptor = typename Get_map<Mesh, Mesh>::type::Dart_const_descriptor;

	FlatTorus_length_weight_functor(const Mesh& m, double period, double border) : m_mesh(m), m_map(m), _period(period), _border(border)
	{}

	Weight_t operator() (Dart_const_descriptor dh) const
	{
		using Point_t = decltype(Get_traits<Mesh>::get_point(m_mesh, dh));
		Point_t p1 = Get_traits<Mesh>::get_point(m_mesh, dh);
		Point_t p2 = Get_traits<Mesh>::get_point(m_mesh, m_map.other_extremity(dh));
		double vec[3] = { p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2] };
		for (int k = 0; k < 3; k++) {
			if (vec[k] < -_border) {
				vec[k] += _period;
			}
			if (vec[k] > _border) {
				vec[k] -= _period;
			}
		}
		//return CGAL::sqrt(CGAL::squared_distance(p1, p2));
		return std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
	}

protected:
	double _period;
	double _border;
	const Mesh& m_mesh;
	const typename Get_map<Mesh, Mesh>::storage_type m_map;
};


std::vector<size_t> extract_shortest_homotopy_path(const Eigen::Matrix<double, -1, 3>& V, const Eigen::Matrix<int, -1, 3>& F, double period, double deduce)
{
	LCC_3 lcc;
	IB ib(lcc);
	std::map<LCC_3::Point, int> p2id;
	for (int i = 0; i < V.rows(); i++) {
		ib.add_vertex(LCC_3::Point(V(i, 0), V(i, 1), V(i, 2)));
		p2id[LCC_3::Point(V(i, 0), V(i, 1), V(i, 2))] = i;
	}
	ib.begin_surface();
	for (int i = 0; i < F.rows(); i++) {
		ib.begin_facet();
		ib.add_vertex_to_facet(F(i, 0));
		ib.add_vertex_to_facet(F(i, 1));
		ib.add_vertex_to_facet(F(i, 2));
		ib.end_facet();
	}
	ib.end_surface();
	
	CGAL::Surface_mesh_topology::Curves_on_surface_topology<LCC_3> cst(lcc);
	LCC_3::Dart_const_descriptor root = lcc.dart_descriptor
	(CGAL::get_default_random().get_int(0, static_cast<int>(lcc.number_of_darts()))); // One dart of the mesh
	::Path_on_surface cycle1 =
		cst.compute_shortest_non_contractible_cycle_with_base_point(root);
	FlatTorus_length_weight_functor<LCC_3> wf(lcc, period, deduce);
	::Path_on_surface cycle = cst.compute_shortest_non_contractible_cycle(wf);
	display_cycle_info(lcc, cycle);
	//CGAL::draw(lcc, cycle);
	// std::vector<LCC_3::Point> path_vlist;
	std::vector<size_t> vloop;
	std::vector<LCC_3::Point> clist;
	for (int k = 0; k < cycle.length(); k++) {
		auto p1 = lcc.point(cycle[k]);
		auto p2 = lcc.point(lcc.other_extremity(cycle[k]));
		if (k == 0) {
			clist.emplace_back(p1);
			clist.emplace_back(p2);
		} else if (k == 1) {
			if (p1 == clist[0] || p1 == clist[1]) {
				clist.emplace_back(p2);
			} else {
				clist.emplace_back(p1);
			}
		} else if (k == cycle.length() - 1) {
			if (p1 != clist[0] && p2 != clist[0]) {
				std::swap(clist[0], clist[1]);
			}
		} else {
			if (p1 == *clist.rbegin()) { clist.emplace_back(p2); }
			else { clist.emplace_back(p1); }
		}
	}
	for (int k = 0; k < clist.size(); k++) {
		vloop.push_back(p2id[clist[k]]);
	}
	return vloop;
}

