#include "PeriodicMesher.h"
#include "fundamental_forms.h"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"
#include "igl/min_quad_with_fixed.h"
#include "igl/connected_components.h"
#include "igl/adjacency_matrix.h"
#include "igl/write_triangle_mesh.h"
#include "igl/boundary_facets.h"
#include "matlab/matlab_utils.h"
#include <Eigen/PardisoSupport>
#include <queue>
#include <unordered_set>
#include "cgal/cgal_utils.h"
#include "igl/connected_components.h"
#include "igl/adjacency_matrix.h"
#include "igl/facet_adjacency_matrix.h"
#include "igl/edges.h"
#include "igl/remove_unreferenced.h"
#include "igl/facet_adjacency_matrix.h"
#include "asymptotic_analysis.h"

using namespace msf;

extern std::vector<double> aux_number;

extern double singularity_surgery_tol;
extern int surgery_type;

extern double island_cull;

//#define SINGULARITY_THRESHOLD 30

extern std::string getPath(std::string s);

auto extract_one_connected_component(std::set<OM::SmartFaceHandle>& fpool) {
	std::set<OM::SmartFaceHandle> conn;
	std::queue<OM::SmartFaceHandle> frti;
	auto v0 = *fpool.begin();
	conn.insert(v0);
	frti.push(v0); fpool.erase(v0);
	while (!frti.empty()) {
		auto fh = frti.front(); frti.pop();
		conn.insert(fh);
		for (auto ff : fh.faces()) {
			if (fpool.count(ff)) { frti.push(ff); fpool.erase(ff); }
		}
	}
	return conn;
}

std::vector<std::vector<OM::SmartHalfedgeHandle>> extract_patch_boundary_loop(const std::set<OM::SmartFaceHandle>& patch)
{
	std::set<OM::SmartEdgeHandle> ehset;
	for (auto fh : patch) {
		for (auto eh : fh.edges()) {
			if (ehset.count(eh)) {
				ehset.erase(eh);
			} else {
				ehset.insert(eh);
			}
		}
	}

	std::vector<std::vector<OM::SmartHalfedgeHandle>> heloops;

	while (!ehset.empty()) {
		std::vector<OM::SmartHalfedgeHandle> heloop;
		auto he = ehset.rbegin()->h0();
		if (patch.count(he.face())) he = he.opp();
		ehset.erase(*ehset.rbegin()); heloop.emplace_back(he);
		OM::SmartVertexHandle v1 = he.to(), v0 = he.from();
		while (v1 != v0) {
			bool found_next = false;
			for (auto voh : v1.outgoing_halfedges()) {
				if (ehset.count(voh.edge())) {
					heloop.emplace_back(voh);
					ehset.erase(voh.edge());
					v1 = voh.to();
					found_next = true;
					break;
				}
			}
			if (!found_next) throw std::runtime_error("found no adjecent edges");
		}
		heloops.emplace_back(heloop);
	}

	return heloops;
}

extern std::tuple<Eigen::MatrixX3<msf::Real>, Eigen::MatrixX3i> removeDupVertices(const Eigen::MatrixX3<msf::Real>& v, const Eigen::MatrixX3i& f, msf::Real eps /*= 1e-5*/);

template<typename Flist>
void savePatches(std::string filename, PeriodSurfaceMesh& m, const Flist& fset) {
	Eigen::MatrixX3d V(fset.size() * 3, 3);
	Eigen::MatrixX3i F(fset.size(), 3);
	int fcounter = 0;
	for (auto fh : fset) {
		int counter = 0;
		for (auto vh : fh.vertices()) {
			auto p = m.point(vh);
			V.row(fcounter * 3 + counter) = toEigen(p).transpose();
			counter++;
		}
		F.row(fcounter) = Eigen::Vector3i(fcounter * 3, fcounter * 3 + 1, fcounter * 3 + 2).transpose();
		fcounter++;
	}
	
	std::tie(V, F) = removeDupVertices(V, F, 1e-5);
	igl::write_triangle_mesh(filename, V, F);
}

template<typename Elist>
void saveEdges(std::string filename, PeriodSurfaceMesh& m, const Elist& eset) {
	std::ofstream ofs(filename, std::ios::binary);
	for (auto eh : eset) {
		auto v0 = eh.h0().from();
		auto vec = make_period(m.calc_edge_vector(eh.h0()));
		auto p0 = m.point(v0);
		auto p1 = p0 + vec;
		ofs.write((const char*)p0.data(), sizeof(p0));
		ofs.write((const char*)p1.data(), sizeof(p1));
	}
}

void save_vertex_vector(std::string filename, PeriodSurfaceMesh& m, const Eigen::MatrixX3d& V) {
	std::ofstream ofs(filename, std::ios::binary);
	for (auto vh : m.vertices()) {
		auto p = m.point(vh);
		ofs.write((const char*)p.data(), sizeof(p));
		Eigen::Vector3d v = V.row(vh.idx()).transpose();
		ofs.write((const char*)v.data(), sizeof(v));
	}
}

void save_vertex_vector(std::string filename, PeriodSurfaceMesh& m, const Eigen::VectorXd& V, const Eigen::MatrixX3d& N) {
	std::ofstream ofs(filename, std::ios::binary);
	for (auto vh : m.vertices()) {
		auto p = m.point(vh);
		ofs.write((const char*)p.data(), sizeof(p));
		Eigen::Vector3d v = N.row(vh.idx()).transpose() * V[vh.idx()];
		ofs.write((const char*)v.data(), sizeof(v));
	}
}

void save_vertex_vector(std::string filename, PeriodSurfaceMesh& m, const Eigen::VectorXd& V) {
	std::ofstream ofs(filename, std::ios::binary);
	for (auto vh : m.vertices()) {
		auto p = m.point(vh);
		ofs.write((const char*)p.data(), sizeof(p));
		ofs.write((const char*)&V[vh.idx()], sizeof(double));
	}
}

std::set<OM::SmartFaceHandle> expand_region(const msf::PeriodSurfaceMesh& m,
	const Eigen::PardisoLDLT<Eigen::SparseMatrix<double>>& so, const Eigen::SparseMatrix<double>& L,
	const Eigen::VectorXd& Mv, std::vector<std::set<OM::SmartFaceHandle>>& patches) {
	if (patches.empty()) return {};
#if 0
	// expansion three times
	for (int iter = 0; iter < 2; iter++) {
		std::set<OM::SmartFaceHandle> toAdd;
		for (auto fh : patch) {
			for (auto vh : fh.vertices()) { for (auto ff : vh.faces()) toAdd.insert(ff); }
		}
		for (auto fh : toAdd) patch.insert(fh);
	}

	// erosion three times
	auto heloops = extract_patch_boundary_loop(patch);
	std::set<OM::SmartVertexHandle> todel;
	for (auto heloop : heloops) { for (auto he : heloop) { todel.insert(he.from()); } }
	for (int iter = 0; iter < 2; iter++) {
		std::set<OM::SmartVertexHandle> newdel;
		for (auto vh : todel) {
			for (auto vv : vh.vertices()) { newdel.insert(vv); }
		}
		for (auto vh : newdel) todel.insert(vh);
	}
	// del boundary faces
	for (auto vh : todel) { for (auto vf : vh.faces()) { patch.erase(vf); } }
#else
	std::set<OM::SmartFaceHandle> toDeleteF;

	for (auto patch : patches) {
		Eigen::VectorXd rho(m.n_vertices()); rho.setZero();
		for (auto fh : patch) { for (auto vh : fh.vertices()) { rho[vh.idx()] = 1; } }

		Eigen::VectorXd x = so.solve((L.transpose() * Mv.cwiseProduct(rho)).eval());
		x.array() -= Mv.dot(x) / Mv.sum();

		double xmin = 1e30, xmax = -1e30;
		for (auto fh : patch) { for (auto vh : fh.vertices()) { xmin = std::min(x[vh.idx()], xmin); xmax = std::max(x[vh.idx()], xmax); } }

		std::cout << "xmin = " << xmin << ", xmax = " << xmax << std::endl;

		//eigen2ConnectedMatlab("x", x);
		//{
		//	std::ofstream ofs(getPath("xrho"), std::ios::binary);
		//	for (auto vh : m.vertices()) { auto p = m.point(vh); ofs.write((const char*)p.data(), sizeof(p)); ofs.write((const char*)&x[vh.idx()], sizeof(double)); }
		//}

		std::set<OM::SmartFaceHandle> pext;
		std::queue<OM::SmartFaceHandle> fr;
		std::vector<bool> infr(m.n_faces(), false);
		for (auto fh : patch) { fr.push(fh); infr[fh.idx()] = true; }
		while (!fr.empty()) {
			auto fh = fr.front(); fr.pop();
			pext.insert(fh);
			for (auto he : fh.halfedges()) {
				if (infr[he.opp().face().idx()]) continue;
				auto vop = he.opp().next().to();
				if (x[vop.idx()] >= xmin && x[vop.idx()] <= xmax) {
					fr.push(he.opp().face()); infr[he.opp().face().idx()] = true;
				}
			}
		}
		// remove  tip face
		for (auto fh : pext) {
			int counter = 0;
			for (auto ff : fh.faces()) { counter += pext.count(ff); }
			if (counter == 1) { continue; }
			else { toDeleteF.insert(fh); }
		}
	}
	return toDeleteF;
#endif
}

extern std::tuple<std::vector<double>, std::vector<Eigen::Vector3d>> eval_cot_edges(PeriodSurfaceMesh& m);

Eigen::VectorXd compute_singularity_measure(msf::PeriodSurfaceMesh& m, const std::vector<Compile1ring>& vrings, const Eigen::SparseMatrix<double>& L) {
	Eigen::MatrixX3d N(vrings.size(), 3);
	Eigen::VectorXd Av(vrings.size());
	for (int i = 0; i < vrings.size(); i++) {
		N.row(i) = vrings[i].nv.transpose();
		Av[i] = vrings[i].As;
	}
#if 0
	//Eigen::MatrixX3d Ln = so.solve((L.transpose() * N).eval());
	Eigen::MatrixX3d Ln = L * N;
	return { Ln.rowwise().norm().cwiseQuotient(Av), N };
#else
	//auto [ecot, evec] = eval_cot_edges(m);
	//Eigen::Matrix3Xd D2x(3, m.n_vertices());
	Eigen::VectorXd delx(m.n_vertices());
	if (surgery_type == 1) {
		for (int k = 0; k < vrings.size(); k++) delx[k] = std::abs(vrings[k].H);
	} else if (surgery_type == 2) {
		for (int k = 0; k < vrings.size(); k++) {
			double k1difk2 = std::sqrt(std::abs(std::pow(vrings[k].H, 2) - vrings[k].K));
			delx[k] = std::max(std::abs(vrings[k].H - k1difk2), std::abs(vrings[k].H + k1difk2));
		}
	}
	return delx;
#endif
}

auto maximal_connected_face_region(msf::PeriodSurfaceMesh& m, OM::SmartVertexHandle vhseed, const std::set<OM::SmartEdgeHandle>& banEdges) {
	std::vector<bool> planeFaces(m.n_faces(), false);
	std::queue<OM::SmartFaceHandle> frti;
	for (auto vf : vhseed.faces()) { planeFaces[vf.idx()] = true; frti.push(vf); }
	while (!frti.empty()) {
		auto fh = frti.front(); frti.pop();
		planeFaces[fh.idx()] = true;
		for (auto he : fh.halfedges()) {
			auto ff = he.opp().face();
			if (!planeFaces[ff.idx()]) {
				if (banEdges.count(he.edge())) continue;
				frti.push(ff);  planeFaces[ff.idx()] = true;
			}
		}
	}

	std::set<OM::SmartFaceHandle> isofaces;
	for (auto fh : m.faces()) { if (!planeFaces[fh.idx()])isofaces.insert(fh); }
	return isofaces;
}

auto filter_singular_faces(const std::vector<std::set<OM::SmartFaceHandle>>& patches) {
#if 0
	std::set<OM::SmartFaceHandle> toDeleteF;
	for (int k = 0; k < patches.size(); k++) {
		auto heloops = extract_patch_boundary_loop(patches[k]);
		// If not a disk
		if (heloops.size() > 1) {
			for (auto pt : patches[k]) toDeleteF.insert(pt);
		}
	}
	return toDeleteF;
#else
	return patches;
#endif
}

auto extract_boundary_edges(msf::PeriodSurfaceMesh& m, const std::set<OM::SmartFaceHandle>& toDeleteF) {
	std::set<OM::SmartEdgeHandle> banEdges;
	for (auto fh : toDeleteF) {
		for (auto eh : fh.edges()) {
			if (banEdges.count(eh)) banEdges.erase(eh); else banEdges.insert(eh);
		}
	}
	return banEdges;
}

auto search_singular_strips(msf::PeriodSurfaceMesh& m) {
	std::vector<Compile1ring> vrings;
	Eigen::VectorXd Av(m.n_vertices());
	for (auto vh : m.vertices()) {
		auto [o, ring] = m.find1ring(vh, 1);
		vrings.emplace_back(o, ring);
		Av[vh.idx()] = vrings.rbegin()->As;
	}

	std::set<OM::SmartVertexHandle> vcan;
	std::set<OM::SmartFaceHandle> fcan;

	//// sift candidate vertices
	//for (auto vh : m.vertices()) {
	//	if (vrings[vh.idx()].H > 15 && vrings[vh.idx()].K < 10) {
	//		vcan.insert(vh);
	//		for (auto vvh : vh.vertices()) { vcan.insert(vvh); }
	//		for (auto fh : vh.faces()) { fcan.insert(fh); }
	//	}
	//}

	auto L = m.getPeriodicLaplacian(2, 1);
	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> so(L.transpose() * L);

	Eigen::VectorXd deln = compute_singularity_measure(m, vrings, L);

	//eigen2ConnectedMatlab("deln", deln);

	//save_vertex_vector(getPath("Ln"), m, Ln);
	//save_vertex_vector(getPath("deln"), m, deln);

	std::vector<OM::SmartVertexHandle> toDeleteV;
	std::set<OM::SmartFaceHandle> toDeleteF;
	for (auto vh : m.vertices()) {
		if (deln[vh.idx()] > singularity_surgery_tol) {
			toDeleteV.push_back(vh);
			for (auto fh : vh.faces()) toDeleteF.insert(fh);
		}
	}


	//eigen2ConnectedMatlab("deln", deln);

	//savePatches(getPath("fdel.obj"), m, toDeleteF);

	//std::cout << "toDeleteF size = " << toDeleteF.size() << std::endl;

#if 0
	std::set<OM::SmartEdgeHandle> banEdges = extract_boundary_edges(m, toDeleteF);

	//std::cout << "banEdges size = " << banEdges.size() << std::endl;

	//saveEdges(getPath("banEdges"), m, banEdges);

	int planeSeed = std::min_element(deln.data(), deln.data() + deln.size()) - deln.data();
	OM::SmartVertexHandle planeVh(planeSeed, &m);

	std::set<OM::SmartFaceHandle> isofaces = maximal_connected_face_region(m, planeVh, banEdges);
	//savePatches(getPath("isopatch.obj"), m, isofaces);

	std::vector<std::set<OM::SmartFaceHandle>> patches;
	while (!isofaces.empty()) { patches.emplace_back(extract_one_connected_component(isofaces)); };
#else
	std::vector<std::set<OM::SmartFaceHandle>> patches;
	while (!toDeleteF.empty()) { patches.emplace_back(extract_one_connected_component(toDeleteF)); };
#endif

	auto toDeletePatches = filter_singular_faces(patches);

	//savePatches(getPath("befexp.obj"), m, toDeleteF);
#if 1
	toDeleteF = expand_region(m, so, L, Av, toDeletePatches);
#endif
	//savePatches(getPath("aftexp.obj"), m, toDeleteF);


	return toDeleteF;
}

std::tuple<Eigen::MatrixX2d, Eigen::MatrixX3i> mesh_closed_curves(const std::vector<std::vector<Eigen::Vector2d>>& loops, double max_area);

Eigen::VectorXi point_indices(const std::vector<Eigen::Vector2d>& vset, const Eigen::MatrixX3d& query, double bb = 10) {
	PeriodicGridIndex ider(Eigen::Vector3d(-bb, -bb, -bb), Eigen::Vector3d(2 * bb, 2 * bb, 2 * bb), 1e-6);
	for (int i = 0; i < vset.size(); i++) {
		Eigen::Vector3d p; p << vset[i], 0;
		ider.insert(p);
	}
	Eigen::VectorXi idlist(query.rows());
	for (int i = 0; i < query.rows(); i++) {
		idlist[i] = ider.query(query.row(i));
	}
	return idlist;
}

Eigen::VectorXi point_indices(const Eigen::MatrixX3d& vset, const Eigen::MatrixX3d& query, double bb = 10) {
	PeriodicGridIndex ider(Eigen::Vector3d(-bb, -bb, -bb), Eigen::Vector3d(2 * bb, 2 * bb, 2 * bb), 1e-6);
	for (int i = 0; i < vset.rows(); i++) {
		Eigen::Vector3d p = vset.row(i).transpose();
		ider.insert(p);
	}
	Eigen::VectorXi idlist(query.rows());
	for (int i = 0; i < query.rows(); i++) {
		idlist[i] = ider.query(query.row(i));
	}
	return idlist;
}

Eigen::MatrixX3d fill_smooth_patch(const Eigen::MatrixX3d& V, const Eigen::MatrixX3i& F, const Eigen::VectorXi& Vfix, const std::vector<Eigen::Vector3d>& vloop) {

	Eigen::MatrixXi bedges;
	igl::boundary_facets(F, bedges);
	if (bedges.rows() > vloop.size() - 1) {
		std::cout << "mesh template has more boundary edges than target, potentially caused by edge splitting in delaunay mesh generation\n";
		throw std::runtime_error("unmatched bounndary");
	}

	Eigen::SparseMatrix<double> Lh; igl::cotmatrix(V, F, Lh);
	std::vector<bool> Vb(Vfix.rows(), false);
	for (int k = 0; k < Vfix.rows(); k++) Vb[k] = Vfix[k] != -1;
	for (int k = 0; k < Lh.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(Lh, k); it; ++it) {
			if (it.row() == it.col()) { if (Vb[it.row()]) it.valueRef() = 1; }
			else { if (Vb[it.row()]) it.valueRef() = 0; }
		}
	}
	//eigen2ConnectedMatlab("Lh", Lh);
	Eigen::SparseLU<Eigen::SparseMatrix<double>> so(Lh);
	Eigen::MatrixX3d bfix(V.rows(), 3); bfix.setZero();
	for (int k = 0; k < Vfix.rows(); k++) {
		if (Vfix[k] != -1) {
			bfix.row(k) = vloop[Vfix[k]].transpose();
		}
	}

	return so.solve(bfix);
}

std::tuple<
	PeriodSurfaceMesh,
	std::map<OM::SmartVertexHandle, OM::SmartVertexHandle>,
	std::vector<OM::SmartVertexHandle>
> extract_submesh(PeriodSurfaceMesh& m, const std::set<OM::SmartFaceHandle>& fhlist) {
	std::map<OM::SmartVertexHandle, OM::SmartVertexHandle> vh2sub;
	std::vector<OM::SmartVertexHandle> sub2vh;
	PeriodSurfaceMesh mnew;
	for (auto fh : fhlist) {
		OM::SmartVertexHandle fvnew[3];
		int counter = 0;
		for (auto vh : fh.vertices()) {
			if (!vh2sub.count(vh)) {
				auto p = m.point(vh);
				auto vhnew = mnew.add_vertex(p);
				vh2sub[vh] = vhnew;
				sub2vh.push_back(vh);
			}
			fvnew[counter] = vh2sub[vh];
			counter++;
		}
		mnew.add_face(fvnew[0], fvnew[1], fvnew[2]);
	}

	// sync vertices period
	auto p0 = mnew.point(vh2sub.begin()->second);
	std::vector<bool> syncd(mnew.n_vertices(), false); syncd[vh2sub.begin()->second.idx()] = true;
	std::queue<OM::SmartVertexHandle> fr;
	fr.push(vh2sub.begin()->second);
	while (!fr.empty()) {
		auto vh = fr.front(); fr.pop();
		for (auto vv : vh.vertices()) {
			if (syncd[vv.idx()]) continue;
			// sync period
			mnew.point(vv) = mnew.point(vh) + make_period(mnew.point(vv) - mnew.point(vh));
			syncd[vv.idx()] = true;
			fr.push(vv);
		}
	}

	return std::make_tuple(mnew, vh2sub, sub2vh);
}

auto localSmoothing(PeriodSurfaceMesh& m, const std::set<OM::SmartFaceHandle>& fhlist) {
	auto [mnew, vh2sub, sub2vh] = extract_submesh(m, fhlist);

	// get Vlist Flist
	auto [V, F] = mnew.getVFlist();

	//igl::write_triangle_mesh(getPath("mnewsync.obj"), V, F);

	// find constrained fixed dof 
	std::set<OM::SmartVertexHandle> vhfix;
	for (auto vh : mnew.vertices()) {
		if (vh.is_boundary()) {
			vhfix.insert(vh);
			for (auto vv : vh.vertices()) { vhfix.insert(vv); }
		}
	}
	Eigen::VectorXi fixdof(vhfix.size());
	Eigen::MatrixX3d Y(fixdof.rows(), 3);
	{ int counter = 0; for (auto vh : vhfix) { fixdof[counter] = vh.idx(); Y.row(counter) = toEigen(mnew.point(vh)).transpose(); counter++; } }

	// fairing 3 times
	for (int iter = 0; iter < 3; iter++) {
		// bilaplacian
		Eigen::SparseMatrix<double> L2;
		igl::cotmatrix(V, F, L2);
		Eigen::SparseMatrix<double>  M;
		igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
		L2 = (L2.transpose() * M.cwiseInverse() * L2).eval();

		// quadratic programming with fixed value
		igl::min_quad_with_fixed_data<double> data;
		Eigen::SparseMatrix<double> Aeq(0, mnew.n_vertices());
		if (!igl::min_quad_with_fixed_precompute(L2, fixdof, Aeq, true, data)) {
			igl::write_triangle_mesh(getPath("tris.obj"), V, F);
			std::cout << "min quad precompute failed\n"; 
			throw std::runtime_error("failed for min quad");
		}

		Eigen::MatrixXd B(mnew.n_vertices(), 3); B.setZero();
		Eigen::MatrixXd Beq(0, 3);
		if (!igl::min_quad_with_fixed_solve(data, B, Y, Beq, V)) { std::cout << "min quad solve failed\n"; }
	}

	//mnew.write(getPath("fillsmooth.obj"));

	// update solution to original mesh
	for (auto [vh, subvh] : vh2sub) {
		m.point(vh) = toOM(V.row(subvh.idx()));
	}
}

void post_smooth_patch(PeriodSurfaceMesh& m, const std::vector<OM::SmartVertexHandle>& vhfix, const std::vector<OM::SmartFaceHandle>& fhlist) {
	std::set<OM::SmartFaceHandle> fext(fhlist.begin(), fhlist.end());
	for (auto vh : vhfix) {
		for (auto fh : vh.faces()) { fext.insert(fh); }
	}
	localSmoothing(m, fext);
}

void relax_angle_distorsion(const std::vector<Eigen::Vector3d>& poly3, std::vector<Eigen::Vector2d>& poly2) {
	std::vector<double> angs(poly3.size() - 1);
	std::vector<double> lens3(poly3.size() - 1), lens2(poly3.size() - 1);
	int N = poly3.size() - 1;
	for (int k = 0; k < poly3.size() - 1; k++) {
		angs[k] = ((poly3[k + 1] - poly3[k]).normalized().dot((poly3[(k + N - 1) % N] - poly3[k]).normalized()));
		lens3[k] = (poly3[k + 1] - poly3[k]).norm();
		lens2[k] = (poly2[k + 1] - poly2[k]).norm();
	}
	double step = 0.01;
	double w_ang = 0.1, w_cfm = 0.5, w_len = 1;
	for (int iter = 0; iter < 20; iter++) {
		std::vector<Eigen::Vector2d> step_vec(poly2.size() - 1);
		for (int k = 0; k < poly2.size() - 1; k++) {
			using ADScalar = Eigen::AutoDiffScalar<Eigen::Vector2d>;
			Eigen::Vector<ADScalar, 2> vk;
			vk[0] = ADScalar(poly2[k][0], 2, 0); vk[1] = ADScalar(poly2[k][1], 2, 1);
			ADScalar obj = w_ang * Eigen::pow(((poly2[k + 1] - vk).normalized().dot((poly2[(k + N - 1) % N] - vk).normalized())) - angs[k], 2);
			obj += w_cfm * Eigen::pow((poly2[k + 1] - vk).norm() * lens3[(k + N - 1) % N] - (poly2[(k + N - 1) % N] - vk).norm() * lens3[k], 2);
			obj += w_len * Eigen::pow((poly2[k + 1] - vk).norm() - lens2[k], 2);
			step_vec[k] = -step * obj.derivatives();
		}
		for (int k = 0; k < poly2.size() - 1; k++) { poly2[k] += step_vec[k]; }
		*poly2.rbegin() = poly2[0];
	}
}

template<typename VLIST>
std::set<OM::SmartFaceHandle> extract_k_ring(int kmx, const VLIST& seedlist) {
	std::set<OM::SmartVertexHandle> vset(seedlist.begin(), seedlist.end());

	for (int iter = 0; iter < kmx; iter++) {
		std::set<OM::SmartVertexHandle> toAdd;
		for (auto vh : vset) { for (auto vv : vh.vertices()) { toAdd.insert(vv); } }
		for (auto vh : toAdd) vset.insert(vh);
	}

	std::set<OM::SmartFaceHandle> fset;
	for (auto vh : vset) { for (auto vf : vh.faces()) { fset.insert(vf); } }

	return fset;
}

auto select_small_patch(const Eigen::MatrixX3d& V, const Eigen::MatrixX3i& F, const Eigen::Vector3d& vseed) {
	Eigen::MatrixX2i edges;
	Eigen::SparseMatrix<int> A;
	igl::adjacency_matrix(F, A);
	igl::edges(A, edges);
	if (V.rows() + F.rows() - edges.rows() == 1) return std::make_tuple(V, F);
	Eigen::VectorXi Cid,Ks;
	igl::connected_components(A, Cid, Ks);
	int ktgt = -1;
	for (int i = 0; i < V.rows(); i++) {
		if ((V.row(i).transpose() - vseed).squaredNorm() < 1e-12) {
			ktgt = Cid[i];
			break;
		}
	}
	if (ktgt == -1) {
		std::cout << "cannot find the small mesh strips\n";
		throw std::runtime_error("cannot find the small mesh strips"); 
	}
	Eigen::MatrixX3i fsmall = F;
	int counter = 0;
	for (int i = 0; i < F.rows(); i++) {
		if (Cid[F(i, 0)] == ktgt) {
			fsmall.row(counter++) = F.row(i);
		}
	}
	fsmall.conservativeResize(counter, 3);
	Eigen::MatrixX3d NV;
	Eigen::MatrixX3i NF;
	Eigen::VectorXi I, J;
	igl::remove_unreferenced(V, fsmall, NV, NF, I, J);
	return std::make_tuple(NV, NF);
}

void fill_one_hole(PeriodSurfaceMesh& m, const std::vector<OM::SmartHalfedgeHandle>& heloop) {
#if 0
	// sync periodic vertices
	double total_len = 0;
	std::vector<OM::SmartVertexHandle> vhloop;
	std::vector<Eigen::Vector3d> vloop;
	vloop.emplace_back(toEigen(m.point(heloop.begin()->from())));
	for (auto he : heloop) {
		auto vec = make_period(m.calc_edge_vector(he));
		total_len += vec.norm();
		vloop.emplace_back(*vloop.rbegin() + toEigen(vec));
		vhloop.emplace_back(he.from());
	}

	//{ std::ofstream ofs(getPath("vloop"), std::ios::binary); ofs.write((const char*)vloop.data(), vloop.size() * sizeof(vloop[0])); }

	// flatten to R2
	std::vector<Eigen::Vector2d> vloop2;
	double ang = 0;
	vloop2.emplace_back(std::cos(ang), std::sin(ang));
	for (int i = 0; i < vloop.size() - 1; i++) {
		ang += (vloop[i + 1] - vloop[i]).norm() / total_len * 2 * M_PI;
		vloop2.emplace_back(std::cos(ang), std::sin(ang));
	}
	*vloop2.rbegin() = *vloop2.begin();

	//{ std::ofstream ofs(getPath("vloop2"), std::ios::binary); ofs.write((const char*)vloop2.data(), vloop2.size() * sizeof(vloop2[0])); }

	// relax angle distorsion
	relax_angle_distorsion(vloop, vloop2);

	//{ std::ofstream ofs(getPath("vloop2relx"), std::ios::binary); ofs.write((const char*)vloop2.data(), vloop2.size() * sizeof(vloop2[0])); }

	// 2D mesh 
	//double max_area = std::pow(total_len, 2) / 50;
	double max_area = 0.5;
	//double tgtlen = (M_PI * std::pow(total_len, 2) / 10);
	std::cout << "max_area = " << max_area << std::endl;
	auto [V2, F] = mesh_closed_curves({ vloop2 }, max_area);

	// copy to 3D
	Eigen::MatrixX3d V(V2.rows(), 3); V.setZero(); V.leftCols(2) = V2;

	//igl::write_triangle_mesh(getPath("hole.obj"), V, F);

	// reorient mesh
	{
		Eigen::Vector3d vp[3]; for (int k = 0; k < 3; k++) vp[k] = V.row(F(0, k)).transpose();
		if ((vp[1] - vp[0]).cross(vp[2] - vp[0])[2] < 0) { Eigen::VectorXi Fi = F.col(1); F.col(1) = F.col(2); F.col(2) = Fi; }
	}

	// find original indices in edge loop 
	auto Vid = point_indices(vloop2, V);

	// solve smoothed filling
	auto X = fill_smooth_patch(V, F, Vid, vloop);

	//igl::write_triangle_mesh(getPath("fill.obj"), X, F);

	// transfer edge loop indices to mesh indices
	for (int i = 0; i < Vid.size(); i++) {
		if (Vid[i] == -1) continue;
		if (Vid[i] > heloop.size()) {
			Vid[i] = heloop.begin()->from().idx();
		}
		else {
			Vid[i] = heloop[Vid[i]].from().idx();
		}
	}

	// append patch faces to mesh
	std::vector<OM::SmartFaceHandle> fnew;
	for (int i = 0; i < F.rows(); i++) {
		std::array<OM::VertexHandle, 3> fvh;
		for (int k = 0; k < 3; k++) {
			int vk = F(i, k);
			if (Vid[vk] == -1) {
				auto p = toOM(X.row(vk));
				Vid[vk] = m.add_vertex(p).idx();
			}
			fvh[k] = OM::VertexHandle(Vid[vk]);
		}
		fnew.emplace_back(m.add_face(fvh[0], fvh[1], fvh[2]));
	}

	// post smooth
	post_smooth_patch(m, vhloop, fnew);

	//m.saveUnitCell(getPath("mnew.obj"));

#else
	std::vector<OM::SmartVertexHandle> vseed;
	for (auto he : heloop) { vseed.emplace_back(he.from()); }

	// maybe filled by last oeprations
	//if (!vseed[0].is_boundary()) return;

	int hole_edge_num = heloop.size();
	
	auto patch = extract_k_ring(4, vseed);

	auto [mring, vh2sub, sub2vh] = extract_submesh(m, patch);
	
	auto [Vring, Fring] = mring.getVFlist();

	//igl::write_triangle_mesh(getPath("mring.obj"), Vring, Fring);

	auto [Vpatch, Fpatch] = patch_holes(Vring, Fring, { hole_edge_num });

	if (Vpatch.rows() == 0) return;

	try {
		std::tie(Vpatch, Fpatch) = select_small_patch(Vpatch, Fpatch, toEigen(mring.point(vh2sub[vseed[0]])));
	} catch (...) {
		igl::write_triangle_mesh(getPath("mring.obj"), Vring, Fring);
		igl::write_triangle_mesh(getPath("patch.obj"), Vpatch, Fpatch);
		std::cout << "seed = " << mring.point(vh2sub[vseed[0]]) << std::endl;
		std::cout << "cannot select small patch" << ", maybe filled by last oepration" << std::endl;
		return;
	}

	// get vertex index in submesh
	auto Vid = point_indices(Vring, Vpatch);

	// submesh to old mesh index
	for (int k = 0; k < Vid.size(); k++) { if (Vid[k] != -1) Vid[k] = sub2vh[Vid[k]].idx(); }

	// append new faces
	for (int fid = 0; fid < Fpatch.rows(); fid++) {
		OM::VertexHandle vhnew[3];
		for (int k = 0; k < 3; k++) {
			auto vk = Fpatch(fid, k);
			if (Vid[vk] == -1) { Vid[vk] = m.add_vertex(toOM(Vpatch.row(vk))).idx(); }
			vhnew[k] = OM::VertexHandle(Vid[vk]);
		}
		auto fhnew = m.add_face(vhnew[0], vhnew[1], vhnew[2]);
		if (!fhnew.is_valid()) {
			std::cout << "unoriented face\n";
			igl::write_triangle_mesh(getPath("mring_err.obj"), Vring, Fring);
			igl::write_triangle_mesh(getPath("patch_err.obj"), Vpatch, Fpatch);
			throw std::runtime_error("unoriented face");
		}
		patch.insert(fhnew);
	}

	localSmoothing(m, patch);

	//m.saveUnitCell(getPath("mfill.obj"));
#endif
}

void fill_holes(PeriodSurfaceMesh& m) {
	auto heloops = m.extract_boundary_loop();

	for (auto& heloop : heloops) {
		fill_one_hole(m, heloop);
	}

	for (auto vh : m.vertices()) {
		if (vh.is_boundary()) {
			std::cout << "Warning ! boundary still exist after filling holes, location = " << m.point(vh) << std::endl;
		}
	}
}

//std::array<bool, 3> connection_check(PeriodSurfaceMesh& m, const Eigen::VectorXi& Cid) {
//	auto [cotedge, Ev] = eval_cot_edges(m);
//	auto [Av, Nv] = eval_vertex_mass(m); Eigen::MatrixXd nmat = Eigen::MatrixXd::Map((const double*)Nv.data(), 3, m.n_vertices());
//	auto L = assemble_asym_cond_matrix(m, cotedge);
//	auto blist = assemble_asym_cond_vector(m, cotedge, Ev);
//	for (int k = 0; k < 3; k++) blist.col(k).array() -= blist.col(k).mean();
//
//	// solve system equation
//	Eigen::PardisoLDLT<Eigen::SparseMatrix<Real>> ldlt(L);
//	Eigen::MatrixX3d ulist = ldlt.solve(blist); ulist.row(0).setZero();
//	for (int k = 0; k < 3; k++) ulist.col(k).array() -= ulist.col(k).mean();
//
//	double As = Av.sum();
//	auto kAlist = eval_partial_asym_cond_matrix(m, blist, ulist, Cid);
//	
//	for (auto fh : m.faces()) {
//		auto fvh = m.getFaceVertexHandle(fh);
//		int cid = Cid[fvh[0].idx()];
//		bool rem = true;
//		for (int k = 0; k < 3; k++) {
//			rem = rem && kAlist[cid](0, 0) < 1e-3;
//		}
//	}
//}

void remove_unconnected_island(PeriodSurfaceMesh& m, const Eigen::VectorXi& Cid) {
#if 1
	m.period_shift();
	
#else
	std::vector<OM::SmartVertexHandle> vhseed(Cid.maxCoeff() + 1, OM::SmartVertexHandle(-1, &m));
	for (auto vh : m.vertices()) {
		int cid = Cid[vh.idx()];
		if (!vhseed[cid].is_valid()) { vhseed[cid] = vh; }
	}

	std::vector<bool> should_remove(vhseed.size(), true);
	for (int seedid = 0; seedid < vhseed.size(); seedid++) {
		auto v0h = vhseed[seedid];
		for (int axis = 0; axis < 3; axis++) {
			std::vector<double> vhdist(m.n_vertices(), 0);
			std::queue<std::pair<OM::SmartVertexHandle, double>> vfrt;
			vfrt.emplace(v0h, 0);
			double max_dist = -1;
			while (max_dist < 4 && !vfrt.empty()) {
				auto [fr, dist] = vfrt.front();
				vfrt.pop();
				if ((dist - vhdist[fr.idx()]) < -1e-6) continue;
				vhdist[fr.idx()] = dist;
				max_dist = std::max(max_dist, dist);
				for (auto vh : fr.vertices()) {
					auto vec = make_period(m.point(vh) - m.point(fr));
					double v1dist = dist + vec[axis];
					if (v1dist > vhdist[vh.idx()]) {
						vfrt.emplace(vh, v1dist);
						vhdist[vh.idx()] = v1dist;
						max_dist = std::max(max_dist, v1dist);
					}
				}
			}
			if (max_dist >= 4) {
				should_remove[seedid] = false;
				break;
			}
		}
	}

	for (auto vh : m.vertices()) {
		if (vh.deleted() || !vh.is_valid()) continue;
		auto cid = Cid[vh.idx()];
		if (should_remove[cid]) {
			m.delete_vertex(vh);
		}
	}

	m.garbage_collection();
#endif
}


void deleteisoisland(PeriodSurfaceMesh& m) {
	auto [V, Fper] = m.getVFlist();
	Eigen::SparseMatrix<double> A;
	//igl::adjacency_matrix(Fper, A);
	igl::facet_adjacency_matrix(Fper, A);
	Eigen::VectorXi Cid,Ks;
	igl::connected_components(A, Cid, Ks);
	// already connected
	if (Ks.rows() <= 1) return;
	std::cout << "Delete islands" << std::endl;
#if 0
	remove_unconnected_island(m, Cid);
#else
	int kMax = std::max_element(Ks.data(), Ks.data() + Ks.size()) - Ks.data();
	PeriodSurfaceMesh mnew;
	std::vector<OM::VertexHandle>  vh(V.rows(), OM::VertexHandle(-1));
	int thres = Ks[kMax] * std::min(island_cull, 1.);
	for (int fid = 0; fid < Fper.rows(); fid++) {
		Eigen::Vector3i fv = Fper.row(fid).transpose();
		bool remain = true;
		//for (int k = 0; k < 3; k++) remain = remain && Cid[fv[k]] == kMax;
		//remain = Ks[Cid[fv[0]]] > thres;
		remain = Ks[Cid[fid]] > thres;
		if (!remain) continue;
		OM::VertexHandle vhnew[3];
		for (int k = 0; k < 3; k++) {
			if (!vh[fv[k]].is_valid()) {
				vh[fv[k]] = mnew.add_vertex(toOM(V.row(fv[k])));
			}
			vhnew[k] = vh[fv[k]];
		}
		mnew.add_face(vhnew[0], vhnew[1], vhnew[2]);
	}
	mnew.saveUnitCell(getPath("mnew.obj"));
	//for (auto fh : mnew.faces()) {
	//	//std::cout << "fh valence = " << fh.valence() << std::endl;
	//	bool is_iso = true;
	//	for (auto he : fh.halfedges()) {
	//		if (!he.opp().is_boundary()) { is_iso = false; break; }
	//	}
	//	if (is_iso) mnew.delete_face(fh);
	//}
	//mnew.garbage_collection();
	//mnew.saveUnitCell(getPath("mnew.obj"));
	m = mnew;
#endif
}

bool msf::PeriodSurfaceMesh::surgery()
{
	auto toDeleteF = search_singular_strips(*this);
	if (toDeleteF.empty()) return false;
	
	//savePatches(getPath("fdel.obj"), *this, toDeleteF);

	for (auto fh : toDeleteF) delete_face(fh); 
	garbage_collection(); deleteisoisland(*this);
	//garbage_collection();
	
	//saveUnitCell(getPath("mdel.obj"));

	fill_holes(*this);
	
	deleteisoisland(*this);

	garbage_collection();

	return true;
}


void test_remove_island(std::string mfile) {
	PeriodSurfaceMesh m;
	m.read(mfile, false, false);
	m.mergePeriodVertices();
	m.delaunayRemesh(5, 0.1, 0.03, 1);
	m.surgery();
	m.saveUnitCell(getPath("rem.obj"));
}
