#include "mesh_intersection.h"
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <igl/connected_components.h>
#include <igl/adjacency_matrix.h>

extern std::string getPath(std::string);

typedef CGAL::Epick K;
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor   face_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;


std::map<msf::UDEdge, char> msf::extractNonManifoldUDEdges(const Eigen::MatrixX3i& Flist)
{
	std::map<UDEdge, char> edge_refs;
	for (int i = 0; i < Flist.rows(); i++) {
		for (int k = 0; k < 3; k++) {
			edge_refs[UDEdge(Flist(i, (k + 1) % 3), Flist(i, k % 3))] += 1;
		}
	}
	for (auto it = edge_refs.begin(); it != edge_refs.end();) {
		if (it->second < 3) { edge_refs.erase(it++); }
		else { it++; }
	}
	std::cout << "Found " << edge_refs.size() << " non-manifold edges" << std::endl;
	return edge_refs;
}

auto genCGALMesh(const Eigen::MatrixX3<msf::Real>& Vlist, const Eigen::MatrixX3i& Flist) {
	Mesh m1;
	for (int i = 0; i < Vlist.rows(); i++) {
		m1.add_vertex(Point(Vlist(i, 0), Vlist(i, 1), Vlist(i, 2)));
	}
	using VH = Mesh::vertex_index;
	for (int i = 0; i < Flist.rows(); i++) {
		m1.add_face((VH)Flist(i, 0), (VH)Flist(i, 1), (VH)Flist(i, 2));
	}
	return m1;
}

std::vector<std::vector<int>> msf::connected_components(const Eigen::MatrixX3<Real>& Vlist, const Eigen::MatrixX3i& Flist)
{
#if 0
	auto m = genCGALMesh(Vlist, Flist);
	std::map<face_descriptor, size_t> face2component;
	Mesh::Property_map<face_descriptor, std::size_t> fccmap =
		m.add_property_map<face_descriptor, std::size_t>("f:CC").first;
	std::size_t num = PMP::connected_components(m, fccmap);
	std::vector<std::vector<int>> components(num);
	for (auto f : m.faces()) {
		int fid = f.idx();
		size_t cid = fccmap[f];
		components[cid].push_back(fid);
	}
	return components;
#else
	Eigen::SparseMatrix<Real> A;
	igl::adjacency_matrix(Flist, A);
	Eigen::VectorXi Cid, Ks;
	igl::connected_components(A, Cid, Ks); // vertex connected components
	std::vector<std::set<int>> connected_vset(Ks.size());
	for (int k = 0; k < Cid.size(); k++) {
		connected_vset[Cid[k]].insert(k);
	}
	std::vector<std::vector<int>> components(connected_vset.size());
	for (int i = 0; i < Flist.rows(); i++) {
		for (int k = 0; k < connected_vset.size(); k++) {
			if (connected_vset[k].count(Flist(i, 0))) {
				components[k].push_back(i);
			}
		}
	}
	return components;
#endif
}

void msf::mesh_split(Eigen::MatrixX3<Real>& V1, Eigen::MatrixX3i& F1, const  Eigen::MatrixX3<Real>& V2, const Eigen::MatrixX3i& F2)
{
	auto vlist1 = V1; auto flist1 = F1;
	auto vlist2 = V2; auto flist2 = F2;
	Mesh m1, m2;
	for (int i = 0; i < vlist1.rows(); i++) {
		m1.add_vertex(Point(vlist1(i, 0), vlist1(i, 1), vlist1(i, 2)));
	}
	using VH = Mesh::vertex_index;
	for (int i = 0; i < flist1.rows(); i++) {
		m1.add_face((VH)flist1(i, 0), (VH)flist1(i, 1), (VH)flist1(i, 2));
	}
	for (int i = 0; i < vlist2.rows(); i++) {
		m2.add_vertex(Point(vlist2(i, 0), vlist2(i, 1), vlist2(i, 2)));
	}
	for (int i = 0; i < flist2.rows(); i++) {
		m2.add_face((VH)flist2(i, 0), (VH)flist2(i, 1), (VH)flist2(i, 2));
	}
	bool m1selit = PMP::does_self_intersect(m1, CGAL::parameters::vertex_point_map(get(CGAL::vertex_point, m1)));
	bool m2selit = PMP::does_self_intersect(m2, CGAL::parameters::vertex_point_map(get(CGAL::vertex_point, m2)));
	if (m1selit || m2selit) {
		std::cout << "splitting self intersected mesh" << std::endl;
		std::vector<std::pair<face_descriptor, face_descriptor> > intersected_tris;
		PMP::self_intersections(faces(m1), m1, std::back_inserter(intersected_tris));
		std::cout << intersected_tris.size() << " pairs of triangles intersect on m1." << std::endl;
		intersected_tris.clear();
		PMP::self_intersections(faces(m2), m2, std::back_inserter(intersected_tris));
		std::cout << intersected_tris.size() << " pairs of triangles intersect on m2." << std::endl;
		throw std::runtime_error("mesh has self intersection");
	}
	PMP::split(m1, m2);
	V1.resize(m1.number_of_vertices(), 3);
	int counter = 0;
	for (auto v : m1.vertices()) {
		auto p = m1.point(v);
		for (int i = 0; i < 3; i++) {
			V1(counter, i) = p[i];
		}
		counter++;
	}
	counter = 0;
	F1.resize(m1.number_of_faces(), 3);
	for (auto f : m1.faces()) {
		int k = 0;
		for (auto v : m1.vertices_around_face(m1.halfedge(f))) {
			F1(counter, k++) = v.idx();
		}
		counter++;
	}
}


void msf::mesh_split_non_manifold_vertices(Eigen::MatrixX3<Real>& V, Eigen::MatrixX3i& F)
{
#if 1
	Mesh m;
	for (int i = 0; i < V.rows(); i++) {
		m.add_vertex(Point(V(i, 0), V(i, 1), V(i, 2)));
	}
	using VH = Mesh::vertex_index;
	for (int i = 0; i < F.rows(); i++) {
		m.add_face((VH)F(i, 0), (VH)F(i, 1), (VH)F(i, 2));
	}
	PMP::duplicate_non_manifold_vertices(m);
	V.resize(m.number_of_vertices(), 3);
	F.resize(m.number_of_faces(), 3);
	int counter = 0;
	for (auto v : m.vertices()) {
		auto p = m.point(v);
		for (int k = 0; k < 3; k++) { V(counter, k) = p[k]; }
		counter++;
	}
	counter = 0;
	for (auto f : m.faces()) {
		auto fv = m.vertices_around_face(m.halfedge(f));
		int k = 0;
		for (auto v : fv) { F(counter, k) = v.idx(); k++; }
		counter++;
	}
#else
	auto edgeRefs = extractNonManifoldUDEdges(F);
	for (int i = 0; i < F.rows(); i++) {
		for (int k = 0; k < 3; k++) {
			auto e = UDEdge(F(i, (k + 1) % 3), F(i, k));
			if (edgeRefs.count(e)) {

			}
		}
	}
#endif
}
