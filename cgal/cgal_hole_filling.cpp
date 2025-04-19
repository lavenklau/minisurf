#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <Eigen/Eigen>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef boost::graph_traits<Mesh>::halfedge_descriptor   halfedge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor       face_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor     vertex_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;
bool is_small_hole(halfedge_descriptor h, Mesh& mesh,
	double max_hole_diam, int max_num_hole_edges)
{
	int num_hole_edges = 0;
	CGAL::Bbox_3 hole_bbox;
	for (halfedge_descriptor hc : CGAL::halfedges_around_face(h, mesh))
	{
		const Point& p = mesh.point(target(hc, mesh));
		hole_bbox += p.bbox();
		++num_hole_edges;
		// Exit early, to avoid unnecessary traversal of large holes
		if (num_hole_edges > max_num_hole_edges) return false;
		if (hole_bbox.xmax() - hole_bbox.xmin() > max_hole_diam) return false;
		if (hole_bbox.ymax() - hole_bbox.ymin() > max_hole_diam) return false;
		if (hole_bbox.zmax() - hole_bbox.zmin() > max_hole_diam) return false;
	}
	return true;
}

bool is_small_hole(halfedge_descriptor h, Mesh& mesh, const std::set<int> holeedge_set)
{
	int num_hole_edges = 0;
	for (halfedge_descriptor hc : CGAL::halfedges_around_face(h, mesh)) {
		++num_hole_edges;
	}
	return holeedge_set.count(num_hole_edges);
}

// Incrementally fill the holes that are no larger than given diameter
// and with no more than a given number of edges (if specified).
int cgal_hole_filling_example(int argc, char* argv[])
{
	const char* filename = (argc > 1) ? argv[1] : "data/mech-holes-shark.off";
	// Both of these must be positive in order to be considered
	double max_hole_diam = (argc > 2) ? boost::lexical_cast<double>(argv[2]) : -1.0;
	int max_num_hole_edges = (argc > 3) ? boost::lexical_cast<int>(argv[3]) : -1;
	std::ifstream input(filename);
	Mesh mesh;
	if (!input || !(input >> mesh)) {
		std::cerr << "Not a valid off file." << std::endl;
		return 1;
	}
	unsigned int nb_holes = 0;
	std::vector<halfedge_descriptor> border_cycles;
	// collect one halfedge per boundary cycle
	CGAL::Polygon_mesh_processing::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));
	for (halfedge_descriptor h : border_cycles)
	{
		if (max_hole_diam > 0 && max_num_hole_edges > 0 &&
			!is_small_hole(h, mesh, max_hole_diam, max_num_hole_edges))
			continue;
		std::vector<face_descriptor>  patch_facets;
		std::vector<vertex_descriptor> patch_vertices;
		bool success = std::get<0>(PMP::triangulate_refine_and_fair_hole(mesh,
			h,
			CGAL::parameters::face_output_iterator(std::back_inserter(patch_facets))
			.vertex_output_iterator(std::back_inserter(patch_vertices))));
		std::cout << "* Number of facets in constructed patch: " << patch_facets.size() << std::endl;
		std::cout << "  Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
		std::cout << "  Is fairing successful: " << success << std::endl;
		++nb_holes;
	}
	std::cout << std::endl;
	std::cout << nb_holes << " holes have been filled" << std::endl;
	std::string outfile = "filled_SM.off";
	std::ofstream out(outfile.c_str());
	std::cout << "Mesh written to: " << outfile << std::endl;
	out.precision(17);
	out << mesh << std::endl;
	return 0;
}

std::tuple<Eigen::MatrixX3d, Eigen::MatrixX3i> patch_holes(const Eigen::MatrixX3d& V, const Eigen::MatrixX3i& F, const std::set<int>& cand_hole_edges) {
	Mesh mesh;
	using VH = Mesh::Vertex_index;
	for (int i = 0; i < V.rows(); i++) { mesh.add_vertex(Point(V(i, 0), V(i, 1), V(i, 2))); }
	for (int i = 0; i < F.rows(); i++) { mesh.add_face(VH(F(i, 0)), VH(F(i, 1)), VH(F(i, 2))); }

	std::map<Mesh::Vertex_index, int> vh2id;
	std::vector<Eigen::Vector3i> Flist;

	unsigned int nb_holes = 0;
	std::vector<halfedge_descriptor> border_cycles;
	// collect one halfedge per boundary cycle
	CGAL::Polygon_mesh_processing::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));
	for (halfedge_descriptor h : border_cycles)
	{
		if (!is_small_hole(h, mesh, cand_hole_edges)) continue;
		std::vector<face_descriptor>  patch_facets;
		std::vector<vertex_descriptor> patch_vertices;
		bool success = std::get<0>(PMP::triangulate_refine_and_fair_hole(mesh,
			h,
			CGAL::parameters::face_output_iterator(std::back_inserter(patch_facets))
			.vertex_output_iterator(std::back_inserter(patch_vertices))));
		//std::cout << "* Number of facets in constructed patch: " << patch_facets.size() << std::endl;
		//std::cout << "  Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
		//std::cout << "  Is fairing successful: " << success << std::endl;
		//for (auto vh : patch_vertices) {
		//	if (!vh2id.count(vh)) { vh2id[vh] = vh2id.size(); }
		//}
		for (auto fh : patch_facets) {
			Eigen::Vector3i fv(-1, -1, -1);
			int counter = 0;
			for (auto vh : mesh.vertices_around_face(mesh.halfedge(fh))) {
				int vid;
				if (vh2id.count(vh)) vid = vh2id[vh];
				else { vh2id[vh] = vid = vh2id.size(); }
				fv[counter++] = vid;
			}
			Flist.emplace_back(fv);
		}
		++nb_holes;
	}

	Eigen::MatrixX3d Vpatch(vh2id.size(), 3);
	Eigen::MatrixX3i Fpatch(Flist.size(), 3);

	for (auto [vh, id] : vh2id) {
		auto p = mesh.point(vh);
		for (int k = 0; k < 3; k++) { Vpatch(id, k) = p[k]; }
	}

	Fpatch = Eigen::Matrix3Xi::Map((const int*)Flist.data(), 3, Flist.size()).transpose();

	return std::make_tuple(Vpatch, Fpatch);
}