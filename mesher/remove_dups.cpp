#include "PeriodicMesher.h"
#include <igl/remove_duplicate_vertices.h>
#include <igl/remove_unreferenced.h>


void msf::PeriodSurfaceMesh::remove_dup_vertices(Real eps)
{
	auto [v, f] = getVFlist();
	// clear old geometry
	OM_Mesh::clear(); OM_Mesh::clean(); OM_Mesh::reset_status();
	// remove duplicated vertices
	auto [vlist, flist] = remove_dup_vertices(v, f, eps);
	// rebuild meshes
	using V = OM::VertexHandle;
	for (int i = 0; i < vlist.rows(); i++) {
		OM::Vec3d p(vlist(i, 0), vlist(i, 1), vlist(i, 2));
		add_vertex(p);
	}
	for (int i = 0; i < flist.rows(); i++) {
		add_face(V(flist(i, 0)), V(flist(i, 1)), V(flist(i, 2)));
	}
	request_normal_status();
}

template<typename Scalar>
inline int countUnique(Scalar a, Scalar b, Scalar c) {
	if (a == b && b == c) {
		return 1; // All three numbers are the same
	} else if (a == b || b == c || a == c) {
		return 2; // Two numbers are the same
	} else {
		return 3; // All numbers are unique
	}
}

std::tuple<Eigen::MatrixX3<msf::Real>, Eigen::MatrixX3i> removeDupVertices(const Eigen::MatrixX3<msf::Real>& v, const Eigen::MatrixX3i& f, msf::Real eps /*= 1e-5*/) {
	// remove unreferencedd vertices and faces
	Eigen::MatrixX3<msf::Real> vlist_refed;
	Eigen::MatrixX3i flist_refed;
	Eigen::VectorXi referid;
	igl::remove_unreferenced(v, f, vlist_refed, flist_refed, referid);

	// remove duplicated vertices
	Eigen::MatrixX3<msf::Real> vlist;
	Eigen::MatrixX3i flist;
	Eigen::VectorXi svi, svj;
	igl::remove_duplicate_vertices(vlist_refed, flist_refed, eps, vlist, svi, svj, flist);

	// remove degenerate faces
	int counter = 0;
	for (int i = 0; i < flist.rows(); i++) {
		if (countUnique(flist(i, 0), flist(i, 1), flist(i, 2)) < 3) { continue; }
		flist.row(counter) = flist.row(i);
		counter++;
	}
	flist.conservativeResize(counter, 3);
	return { vlist,flist };

}

std::tuple<Eigen::MatrixX3<msf::Real>, Eigen::MatrixX3i> msf::PeriodSurfaceMesh::remove_dup_vertices(const Eigen::MatrixX3<Real>& v, const Eigen::MatrixX3i& f, Real eps /*= 1e-5*/)
{
	return removeDupVertices(v, f, eps);
}


