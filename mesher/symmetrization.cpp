#include "PeriodicMesher.h"
#include "igl/bfs_orient.h"
#include "igl/write_triangle_mesh.h"

using namespace msf;

extern std::string getPath(std::string);

enum SymType{
	Orthogonal = 0,
	Cubic = 1
};

// 1 - 

bool in_unit_cell(const Eigen::Vector3d& p, int type) {
	if (type == Orthogonal) {
		return p.minCoeff() >= 0;
	}
	else if (type == Cubic) {
		return p[0] <= p[1] && p[1] <= p[2] && p.minCoeff() >= 0;
	}
}

std::vector<Eigen::Vector4d> split_plane(int type) {
	if (type == Orthogonal) {
		return {
			{1,0,0,0},{1,0,0,0},{1,0,0,0},
			{-1,0,0,1},{0,-1,0,1},{0,0,-1,1} };
	}
	else if (type == Cubic) {
		return { {-1,1,0,0},{0,-1,1,0},{-1,1,0,0} ,
			{1,0,0,0},{1,0,0,0},{1,0,0,0},
			{-1,0,0,1},{0,-1,0,1},{0,0,-1,1}
		};
	}
}

// Comparator for Eigen::Matrix3d
struct Matrix3dComparator {
	bool operator()(const Eigen::Matrix3d& lhs, const Eigen::Matrix3d& rhs) const {
		// Compare lexicographically by rows
		for (int row = 0; row < lhs.rows(); ++row) {
			for (int col = 0; col < lhs.cols(); ++col) {
				if (lhs(row, col) < rhs(row, col)) {
					return true;
				}
				if (lhs(row, col) > rhs(row, col)) {
					return false;
				}
			}
		}
		return false; // Matrices are equal
	}
};

auto group_recursive(const std::vector<Eigen::Matrix3d>& generators, int max_iter) {
	std::set<Eigen::Matrix3d, Matrix3dComparator> matrixSet{ generators.begin(),generators.end() };
	matrixSet.insert(Eigen::Matrix3d::Identity());
	for (int i = 0; i < max_iter; i++) {
		auto oldset = matrixSet;
		int oldsize = oldset.size();
		for (auto m : generators) {
			for (auto o : oldset) {
				matrixSet.insert(m * o);
				matrixSet.insert(o * m);
			}
		}
		int newsize = matrixSet.size();
		if (newsize == oldsize) break;
	}
	return matrixSet;
}

const auto& generate_group(int type) {
	static std::map<int, std::vector<Eigen::Matrix3d>> groups;
	if (groups.count(type)) return  groups[type];
	if (type == Orthogonal) {
		std::vector<Eigen::Matrix3d> gen;
		for (int k = 0; k < 3; k++) {
			Eigen::Vector3d refl; refl.setConstant(1);
			refl[k] = -1;
			gen.emplace_back(refl.asDiagonal());
		}
		auto matset = group_recursive(gen, 10);
		groups[type] = std::vector<Eigen::Matrix3d>(matset.begin(), matset.end());
	}
	else if (type == Cubic) {
		std::vector<Eigen::Matrix3d> gen;
		for (int k = 0; k < 3; k++) {
			Eigen::Vector3d refl; refl.setConstant(1);
			refl[k] = -1;
			gen.emplace_back(refl.asDiagonal());
		}
		gen.emplace_back(Eigen::PermutationMatrix<3>(Eigen::Vector3i(1, 0, 2)));
		gen.emplace_back(Eigen::PermutationMatrix<3>(Eigen::Vector3i(0, 2, 1)));
		gen.emplace_back(Eigen::PermutationMatrix<3>(Eigen::Vector3i(2, 1, 0)));
		auto matset = group_recursive(gen, 20);
		groups[type] = std::vector<Eigen::Matrix3d>(matset.begin(), matset.end());
	}
	return groups.at(type);
}

std::tuple<Eigen::MatrixX3<msf::Real>, Eigen::MatrixX3i> removeDupVertices(const Eigen::MatrixX3<msf::Real>& v, const Eigen::MatrixX3i& f, msf::Real eps /*= 1e-5*/);

auto orbit_transform(const Eigen::MatrixX3d& V, const Eigen::MatrixX3i& F, int type) {
	Eigen::MatrixX3d Vnew;
	Eigen::MatrixX3i Fnew;
	if (type == Orthogonal) {
		Vnew.resize(V.rows() * 8, 3);
		Fnew = F.replicate(8, 1);
		auto rs = generate_group(type);
		int counter = 0;
		for (auto g : rs) {
			Vnew.block(counter * V.rows(), 0, V.rows(), 3) = Vnew * g.transpose();
			Fnew.block(counter * F.rows(), 0, F.rows(), 3).array() += V.rows();
			counter++;
		}
	}
	else if (type == Cubic) {
		Vnew.resize(V.rows() * 48, 3);
		Fnew = F.replicate(48, 1);
		std::cout << "Vnew size = " << Vnew.rows() << std::endl;
		std::cout << "F.size = " << F.rows() << ", Fnew size = " << Fnew.rows() << std::endl;
		auto rs = generate_group(type);
		int counter = 0;
		for (auto g : rs) {
			Vnew.block(counter * V.rows(), 0, V.rows(), 3) = V * g.transpose();
			Fnew.block(counter * F.rows(), 0, F.rows(), 3).array() += V.rows() * counter;
			counter++;
		}
	}
	igl::write_triangle_mesh(getPath("dup.obj"), Vnew, Fnew);

	std::tie(Vnew, Fnew) = removeDupVertices(Vnew, Fnew, 2e-5);
	
	Eigen::MatrixX3i FF;
	Eigen::VectorXi C;
	igl::bfs_orient(Fnew, FF, C);
	return std::make_tuple(Vnew, FF);
}


extern std::tuple<
	PeriodSurfaceMesh,
	std::map<OM::SmartVertexHandle, OM::SmartVertexHandle>,
	std::vector<OM::SmartVertexHandle>
> extract_submesh(PeriodSurfaceMesh& m, const std::set<OM::SmartFaceHandle>& fhlist);

void symmetrize_period_mesh(PeriodSurfaceMesh& m, int type) {
	auto splist = split_plane(type);
	std::vector<OM::SmartVertexHandle> vhsplit;
	for (auto& pl : splist) {
		for (auto eh : m.edges()) {
			auto p0 = toEigen(m.point(eh.v0()));
			auto p1 = toEigen(m.point(eh.v1()));
			p1 = p0 + make_period((p1 - p0).eval());
			for (int k = 0; k < 3; k++) { if (p1[k] < -1) { p1[k] += 2; p0[k] += 2; } }
			double w0 = p0.homogeneous().dot(pl);
			double w1 = p1.homogeneous().dot(pl);
			if (w0 * w1 > 0) { continue; }
			w0 = std::abs(w0); w1 = std::abs(w1);
			Eigen::Vector3d p = w1 / (w0 + w1) * p0 + w0 / (w0 + w1) * p1;
			vhsplit.emplace_back(m.split(eh, toOM(p)));
		}
	}

	std::set<OM::SmartFaceHandle> fhset;
	for (auto vh : vhsplit) {
		for (auto fh : vh.faces()) {
			auto tri = m.getFacePeriodVertex(fh, vh, 1);
			if (in_unit_cell(tri.rowwise().mean(), type)) { fhset.insert(fh); }
		}
	}

	for (auto fh : m.faces()) {
		auto tri = m.getFacePeriodVertex(fh, 1);
		if (in_unit_cell(tri.rowwise().mean(), type)) { fhset.insert(fh); }
	}

	auto [mnew, vh2sub, sub2vh] = extract_submesh(m, fhset);

	mnew.save(getPath("sub.obj"));
	
	// collapse short edge on boundary
	for (auto he : mnew.halfedges()) {
		if (!he.is_boundary() || !he.is_valid() || he.deleted()) continue;
		auto prev = he.prev();
		auto angle = mnew.calc_sector_angle(prev);
		if (angle < M_PI * 0.9) continue;
		double len = mnew.calc_edge_length(he);
		if (len < 0.005) { if (mnew.is_collapse_ok(he)) mnew.collapse(he); continue; }
		angle = mnew.calc_sector_angle(he.opp().next());
		if (angle < M_PI / 180 * 15) { if (mnew.is_collapse_ok(he)) mnew.collapse(he); continue; }
	}

	// collapse short incident edge on boundary
	for (auto vh : mnew.vertices()) {
		if (!vh.is_boundary() || !vh.is_valid() || vh.deleted()) continue;
		// Maybe dangerous due to delete of he 
		for (auto vih : vh.incoming_halfedges()) {
			if (vih.is_boundary()) continue;
			if (vih.deleted()) continue;
			if (!vih.is_valid()) continue;
			if (vih.from().is_boundary()) continue;
			auto len = mnew.calc_edge_length(vih);
			if (len < 0.005 && mnew.is_collapse_ok(vih)) mnew.collapse(vih);
		}
	}
	
	mnew.garbage_collection();

	auto [V, F] = mnew.getVFlist();
	std::tie(V, F) = orbit_transform(V, F, type);
	igl::write_triangle_mesh(getPath("xform.obj"), V, F);
	
	PeriodSurfaceMesh mm;
	for (int i = 0; i < V.rows(); i++) { mm.add_vertex(toOM(V.row(i))); }
	for (int i = 0; i < F.rows(); i++) {
		OM::VertexHandle fvh[3];
		for (int k = 0; k < 3; k++) fvh[k] = OM::VertexHandle(F(i, k));
		auto fh = mm.add_face(fvh[0], fvh[1], fvh[2]);
		if (!fh.is_valid()) {
			for (int k = 0; k < 3; k++) {
				std::cout << V.row(fvh[k].idx()) << " ";
			}
			std::cout << std::endl;
		}
	}
	//mm.read(V, F, false, false, false);
	mm.save(getPath("VF.obj"));

	m.clear();
	m.read(V, F, false, false, false);
}

void test_symmetrization(void) {
	PeriodSurfaceMesh m;
	m.readMergePeriodBoundary("D:/projects/minisurf/image/thermal/SMdata/non-TPMS/samples/P0.2-8.stl");
	symmetrize_period_mesh(m, Cubic);

	m.save(getPath("sym.obj"));
}
