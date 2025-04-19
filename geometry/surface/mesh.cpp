#include "Config.h"
#include "mesh.h"
#include "grid/PeriodicGrid.h"
#include "facet/BezierTri.h"
#include "fmt/core.h"

using namespace msf;

std::pair<Eigen::Matrix<msf::Real, -1, 3>, Eigen::Matrix<int, -1, 3>> MeshSurface::getVFlist(void)
{
	Eigen::Matrix<Real, -1, 3> vlist(MeshSurface::n_vertices(), 3);
	Eigen::Matrix<int, -1, 3> flist(n_faces(), 3);
	int counter = 0;
	for (auto v : MeshSurface::vertices()) {
		auto pos = point(v);
		for (int k = 0; k < 3; k++) { vlist(counter, k) = pos[k]; }
		counter++;
	}
	counter = 0;
	for (auto f : MeshSurface::faces()) {
		if (status(f).deleted()) {
			printf("\033[31mWarning ! Face has been deleted\033[0m\n");
		}
		//if (!f.is_valid()) continue;
		auto fvlist = getFaceVertexHandle(f);
		if (fvlist[0].is_valid() && fvlist[1].is_valid() && fvlist[2].is_valid()) {
			for (int k = 0; k < 3; k++) flist(counter, k) = fvlist[k].idx();
			counter++;
		}
	}
	flist.conservativeResize(counter, 3);
	return { vlist,flist };
}

msf::OM::SmartHalfedgeHandle MeshSurface::findHe(OM::VertexHandle from, OM::VertexHandle to)
{
	for (auto voh : voh_range(from)) {
		if (voh.to() == to) {
			return voh;
		}
	}
	OM::SmartHalfedgeHandle invalid_he;
	return invalid_he;
}

int MeshSurface::heid(OM::HalfedgeHandle h) const
{
	int counter = 0;
	auto fh = face_handle(h);
	for (auto he : fh_ccw_range(fh)) {
		if (he == h) break;
		counter++;
	}
	return counter;
}

int BezierMesh::heid(OM::HalfedgeHandle h) const {
	auto fh = face_handle(h);
	auto fv = getFaceVertexHandle(fh);
	for (int i = 0; i < 3; i++) {
		if (fv[i] == from_vertex_handle(h)) return i;
	}
	throw std::runtime_error("cannot find halfedge in face");
}

bool MeshSurface::read(std::string filepath, bool normalize_mesh)
{
	OM::IO::read_mesh(*this, filepath);
	if (normalize_mesh) normalize();
	request_normal_status();
	return true;
}

void MeshSurface::write(std::string filename)
{
	if (!OM::IO::write_mesh(*this, filename)) {
		std::cout << "\033[31m" << "Cannot write mesh to " << filename << "\033[0m \n";
		throw std::runtime_error("write mesh failed");
	}
}

void MeshSurface::request_normal_status(void)
{
	MeshSurface::request_vertex_normals();
	MeshSurface::request_face_normals();
	MeshSurface::update_normals();
	MeshSurface::request_vertex_status();
	MeshSurface::request_halfedge_status();
	MeshSurface::request_face_status();
	MeshSurface::request_edge_status();
}



Real MeshSurface::cothe(OM::HalfedgeHandle he) const
{
	if (is_boundary(he)) {
		throw std::runtime_error("evaluating cot weight on boundary edge");
	}
	Real ta = (std::tan)(calc_sector_angle(next_halfedge_handle(he)));
	Real cota = (std::abs)(ta) > 1e30 ? 0 : 1 / ta;
	return cota;
}

Real MeshSurface::cote(OM::EdgeHandle eh) const
{
	if (is_boundary(eh)) {
		throw std::runtime_error("evaluating cot weight on boundary edge");
	}
	Real ta = (std::tan)(calc_sector_angle(next_halfedge_handle(halfedge_handle(eh,0))));
	Real cota = (std::abs)(ta) > 1e30 ? 0 : 1 / ta;
	Real tb = (std::tan)(calc_sector_angle(next_halfedge_handle(halfedge_handle(eh, 1))));
	Real cotb = (std::abs)(tb) > 1e30 ? 0 : 1 / tb;
	return (cota + cotb) / 2;
}

bool MeshSurface::read(const Eigen::MatrixX3<Real>& V, const Eigen::MatrixX3i& F)
{
	for (int i = 0; i < V.rows(); i++) {
		OM::Vec3d p(V(i, 0), V(i, 1), V(i, 2));
		add_vertex(p);
	}
	using VH = OM::VertexHandle;
	for (int i = 0; i < F.rows(); i++) {
		add_face((VH)F(i, 0), (VH)F(i, 1), (VH)F(i, 2));
	}
	MeshSurface::request_vertex_normals();
	MeshSurface::request_face_normals();
	MeshSurface::update_normals();
	return true;
}

void MeshSurface::normalize(void)
{
	auto [pmin, pmax] = boundingBox();
	Real dmax = (pmax - pmin).maxCoeff();
	Eigen::Vector3<Real> c = (pmin + pmax) / 2;
	Real eps = 1e-5;
	for (auto vh : vertices()) {
		OM::Vec3d p = point(vh);
		for (int k = 0; k < 3; k++) {
			p[k] = 2 / dmax * (p[k] - c[k]);
			if (-1 + eps > p[k] || p[k] > 1 - eps) { p[k] = p[k] > 0 ? 1 : -1; }
		}
		set_point(vh, p);
	}
}

void MeshSurface::faceVhandle(OM::FaceHandle fh, OM::VertexHandle vhlist[3]) const
{	
	int counter = 0;
	for (auto vh : fv_range(fh)) {
		vhlist[counter++] = vh;
	}
	return;
}

std::pair<msf::OM::SmartHalfedgeHandle, msf::OM::SmartHalfedgeHandle> MeshSurface::adjacent_he(OM::FaceHandle f1, OM::FaceHandle f2) const
{
	auto he1 = OM::SmartHalfedgeHandle(halfedge_handle(f1).idx(), this);
	auto he2 = OM::SmartHalfedgeHandle(halfedge_handle(f2).idx(), this);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (he1.opp() == he2) { return { he1,he2 }; }
			he2 = he2.next();
		}
		he1 = he1.next();
	}
	OM::SmartHalfedgeHandle invhe(-1, this);
	return { invhe,invhe };
}

std::vector<std::vector<msf::OM::SmartHalfedgeHandle>> MeshSurface::extract_boundary_loop(void) const
{
	std::set<OM::SmartHalfedgeHandle> heset;
	for (auto he : halfedges()) {
		if (he.is_boundary()) { heset.insert(he); }
	}
	std::vector<std::vector<OM::SmartHalfedgeHandle>> loops;
	std::vector<std::vector<OM::SmartVertexHandle>> vloops;
	while (!heset.empty()) {
		auto he0 = *heset.begin();
		auto he = he0;
		std::vector<OM::SmartHalfedgeHandle> loop;
		std::vector<OM::SmartVertexHandle> vloop;
		do {
			vloop.push_back(he.from());
			heset.erase(he);
			loop.push_back(he);
			he = he.next();
		} while (he != he0);
		loops.emplace_back(loop);
		vloops.emplace_back(vloop);
	}
	return loops;
}

Real MeshSurface::average_edge_length(void)
{
	int counter = 0;
	Real l_sum = 0;
	for (auto eh : edges()) {
		l_sum += OM_Mesh::calc_edge_length(eh);
		counter++;
	}
	l_sum /= counter;
	return l_sum;
}

Real MeshSurface::totalArea(void) const
{
	Real s = 0;
	for (auto f : faces()) {
		s += calc_face_area(f);
	}
	return s;
}

std::vector<msf::OM::FaceHandle> MeshSurface::facelist(void) const
{
	std::vector<msf::OM::FaceHandle> flist;
	for (auto f : faces()) {
		flist.push_back(f);
	}
	return flist;
}

std::vector<msf::OM::VertexHandle> MeshSurface::vertexlist(void) const
{
	std::vector<OM::VertexHandle> vlist;
	for (auto v : vertices()) {
		vlist.push_back(v);
	}
	return vlist;
}

int BezierMesh::_zid(Real u, Real v) const
{
	if (u < 1e-5) {
		return 0;
	}
	else if (v < 1e-5) {
		return 1;
	}
	return 2;
}

int BezierMesh::_zid(int i, int j, int k) const
{
	if (i == 0) { return 0; }
	else if (j == 0) { return 1; }
	return 2;

}

Eigen::Vector2<msf::Real> BezierMesh::transition(OM::HalfedgeHandle hefrom, OM::HalfedgeHandle heto, Real u, Real v) const
{
	auto fh = MeshSurface::face_handle(hefrom);
	int hid = heid(hefrom);
	int hid2 = heid(heto);
	int match = (hid2 - hid + 3) % 3;
	Real uthis[3] = { u, v, 1 - u - v };
	//int z = _zid(u, v);
	int z = (hid + 2) % 3;
	if (uthis[z] != 0) throw std::runtime_error("unexpected transition");
	Real uother[3] = { 0 };
	uother[(z + match + 1) % 3] = uthis[(z + 2) % 3];
	uother[(z + match + 2) % 3] = uthis[(z + 1) % 3];
	return Eigen::Vector2<Real>(uother[0], uother[1]);
}

Eigen::Vector3i BezierMesh::transition(OM::HalfedgeHandle hefrom, OM::HalfedgeHandle heto, int i, int j, int k) const
{
	auto fh = MeshSurface::face_handle(hefrom);
	//std::cout << getFaceVertex(fh).transpose() << std::endl;
	int hid = heid(hefrom);
	//std::cout << "he  = " << point(from_vertex_handle(hefrom)) << " " << point(to_vertex_handle(hefrom)) << ", id = " << hid << std::endl;
	int hid2 = heid(heto);
	//std::cout << "heto = " << point(from_vertex_handle(heto)) << " " << point(to_vertex_handle(heto)) << ", id = " << hid2 << std::endl;
	int match = (hid2 - hid + 3) % 3;
	//std::cout << "match = " << match << std::endl;
	int uthis[3] = { i, j, order - i - j };
	//std::cout << "uthis = " << uthis[0] << " " << uthis[1] << " " << uthis[2] << std::endl;
	//int z = _zid(i, j, order - i - j);
	int z = (hid + 2) % 3;
	if (uthis[z] != 0) {
		fmt::print("uthis = {}, {}, {};  hid = {};  z = {}", uthis[0], uthis[1], uthis[2], hid, z);
		throw std::runtime_error("unexpected transition");
	}
	int uother[3] = { 0 };
	uother[(z + match + 1) % 3] = uthis[(z + 2) % 3];
	uother[(z + match + 2) % 3] = uthis[(z + 1) % 3];
	//std::cout << "uother = " << uother[0] << " " << uother[1] << " " << uother[2] << std::endl;
	return Eigen::Vector3i(uother[0], uother[1], uother[2]);
}

int BezierMesh::start_vid(OM::HalfedgeHandle h, bool global /*= false*/)
{
	int hid0 = heid(h);
	int o_ijk[3] = { 0 };
	o_ijk[hid0] = order;
	int oid = pid(o_ijk[0], o_ijk[1], o_ijk[2]);
	if (global) {
		oid = f_cid(oid, face_handle(h).idx());
	}
	return oid;
}

int BezierMesh::end_vid(OM::HalfedgeHandle h, bool global /*= false*/)
{
	int hid0 = (heid(h) + 1) % 3;
	int o_ijk[3] = { 0 };
	o_ijk[hid0] = order;
	int oid = pid(o_ijk[0], o_ijk[1], o_ijk[2]);
	if (global) {
		oid = f_cid(oid, face_handle(h).idx());
	}
	return oid;
}

void BezierMesh::set_order(int order_)
{
	order = order_;
}

int BezierMesh::n_ctrlpoints(void)
{
	return ctrl_points.rows();
}

void BezierMesh::buildBezierMesh(void)
{
	int np_tri = (order + 2) * (order + 1) / 2;
	auto [pmin, pmax] = boundingBox();
	Eigen::Vector3<Real> d = pmax - pmin;
	d *= 1 + 2e-5;
	pmin -= d * 1e-5;
	
	// barycentric coordinate list
	std::vector<Eigen::Vector3<Real>> fplist;
	std::vector<int> fpid;
	{
		Eigen::Vector3<Real> fp;
		for (int i = 0; i <= order; i++) {
			fp[0] = Real(i) / order;
			for (int j = 0; j <= order - i; j++) {
				int k = order - i - j;
				fp[1] = Real(j) / order;
				fp[2] = Real(k) / order;
				fplist.emplace_back(fp);
				fpid.emplace_back(pid(i, j, k));
			}
		}
	}

	// bucket points
	PeriodicGridIndex indexer(pmin, d);
	for (auto fh : faces()) {
		auto V = getFaceVertex(fh);
		for (int i = 0; i < fplist.size(); i++) {
			Eigen::Vector3<Real> p = V * fplist[i];
			indexer.insert(p);
		}
	}

	// query index
	ctrl_points = indexer.dumpPoints();
	f_cid.resize(fplist.size(), n_faces());
	for (auto fh : faces()) {
		auto V = getFaceVertex(fh);
		for (int i = 0; i < fplist.size(); i++) {
			Eigen::Vector3<Real> p = V * fplist[i];
			int c_id = indexer.query(p);
			int p_id = fpid[i];
			f_cid(p_id, fh.idx()) = c_id;
		}
	}

	
}

msf::BezierTri<msf::Real, -1> BezierMesh::getBezierTri(OM::FaceHandle fh)
{
	std::vector<Eigen::Vector3<Real>> cpts(f_cid.rows());
	for (int i = 0; i < f_cid.rows(); i++) {
		cpts[i] = ctrl_points.row(f_cid(i, fh.idx())).transpose();
	}
	BezierTri<msf::Real, -1> tri(cpts.data(), order);
	return tri;
}



