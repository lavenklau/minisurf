#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <boost/iterator/function_output_iterator.hpp>
#include <iostream>
#include <string>
#include <vector>
#include "Eigen/Eigen"
#include "igl/remove_unreferenced.h"
typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
using Point = K::Point_3;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;
typedef boost::graph_traits<Mesh>::halfedge_descriptor        halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor            edge_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;
struct halfedge2edge
{
	halfedge2edge(const Mesh& m, std::vector<edge_descriptor>& edges)
		: m_mesh(m), m_edges(edges)
	{}
	void operator()(const halfedge_descriptor& h) const
	{
		m_edges.push_back(edge(h, m_mesh));
	}
	const Mesh& m_mesh;
	std::vector<edge_descriptor>& m_edges;
};

std::string getPath(std::string);

int cgal_remesh_main(void)
{
	const std::string filename = "D:/projects/minisurf/image/temp/octet/octet_0.1.obj";
	Mesh mesh;
	if (!PMP::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh))
	{
		std::cerr << "Invalid input." << std::endl;
		return 1;
	}
	double target_edge_length = 0.03;
	unsigned int nb_iter = 10;
	std::cout << "Split border...";
	std::vector<edge_descriptor> border;
	PMP::border_halfedges(faces(mesh), mesh, boost::make_function_output_iterator(halfedge2edge(mesh, border)));
	PMP::split_long_edges(border, target_edge_length, mesh);
	std::cout << "done." << std::endl;
	std::cout << "Start remeshing of " << filename
		<< " (" << num_faces(mesh) << " faces)..." << std::endl;
	PMP::isotropic_remeshing(faces(mesh), target_edge_length, mesh,
		CGAL::parameters::number_of_iterations(nb_iter)
		.protect_constraints(true)); //i.e. protect border, here
	CGAL::IO::write_polygon_mesh(getPath("out.off"), mesh, CGAL::parameters::stream_precision(17));
	std::cout << "Remeshing done." << std::endl;
	return 0;
}


double aspect_ratio(double a, double b, double c) {
	double s = (a + b + c) / 2.0;
	double AR = (a * b * c) / (8.0 * (s - a) * (s - b) * (s - c) + 1e-30);
	return AR;
}

template<typename Vec>
auto aspect_ratio(const Vec& v0, const Vec& v1, const Vec& v2) {
	return aspect_ratio((double)(v0 - v1).norm(), (double)(v1 - v2).norm(), (double)(v2 - v0).norm());
}


template<typename T>
std::tuple<Eigen::MatrixX3<T>, Eigen::MatrixX3i> cgal_remesh(const Eigen::MatrixX3<T>& V, const Eigen::MatrixX3i& F, T len = 0.03) {
	Mesh mesh;
	for (int i = 0; i < V.rows(); i++) {
		Point p(V(i, 0), V(i, 1), V(i, 2));
		mesh.add_vertex(p);
	}
	for (int i = 0; i < F.rows(); i++) {
		Mesh::Vertex_index fv[3];
		for (int j = 0; j < 3; j++) { fv[j] = Mesh::Vertex_index(F(i, j)); }
		mesh.add_face(fv[0], fv[1], fv[2]);
	}

	double target_edge_length = len;
	unsigned int nb_iter = 10;
	//std::cout << "Split border...";
	std::vector<edge_descriptor> border;
	PMP::border_halfedges(faces(mesh), mesh, boost::make_function_output_iterator(halfedge2edge(mesh, border)));
	PMP::split_long_edges(border, target_edge_length, mesh);
	PMP::isotropic_remeshing(faces(mesh), target_edge_length, mesh,
		CGAL::parameters::number_of_iterations(nb_iter)
		.protect_constraints(true)); //i.e. protect border, here
	//CGAL::IO::write_polygon_mesh(getPath("out.off"), mesh, CGAL::parameters::stream_precision(17));
	if (mesh.has_garbage()) { mesh.collect_garbage(); }
	Eigen::MatrixX3<T> Vnew(mesh.number_of_vertices(), 3);
	Eigen::MatrixX3i Fnew(mesh.number_of_faces(), 3);
	int max_vid = 0;
	for (auto v : mesh.vertices()) {
		//if (max_vid < v.idx()) max_vid = v.idx();
		//std::cout << "v = " << v.idx() << std::endl;
		auto p = mesh.point(v);
		Vnew.row(v.idx()) = Eigen::Vector3<T>(p[0], p[1], p[2]).transpose();
	}
	//std::cout << "vnew finish\n";
	for (auto f : mesh.faces()) {
		auto he = mesh.halfedge(f);
		int counter = 0;
		for (auto fv : mesh.vertices_around_face(he)) {
			Fnew(f.idx(), counter++) = fv.idx();
		}
	}

	// remove degenerate face
	//int counter = 0;
	//for (int i = 0; i < Fnew.rows(); i++) {
	//	if (aspect_ratio(Vnew.row(Fnew(i, 0)), Vnew.row(Fnew(i, 1)), Vnew.row(Fnew(i, 2))) > 10000) {
	//		continue;
	//	}
	//	Fnew.row(counter++) = Fnew.row(i);
	//}
	//Fnew.conservativeResize(counter, 3);

	// remove unreferenced points
	Eigen::MatrixX3<T> NV;
	Eigen::MatrixX3i NF;
	Eigen::VectorXi I;
	igl::remove_unreferenced(Vnew, Fnew, NV, NF, I);

	return { NV, NF };
}

template std::tuple<Eigen::MatrixX3<double>, Eigen::MatrixX3i> cgal_remesh<double>(const Eigen::MatrixX3<double>& V, const Eigen::MatrixX3i& F, double len/* = 0.03*/);
template std::tuple<Eigen::MatrixX3<float>, Eigen::MatrixX3i> cgal_remesh<float>(const Eigen::MatrixX3<float>& V, const Eigen::MatrixX3i& F, float len/* = 0.03*/);
