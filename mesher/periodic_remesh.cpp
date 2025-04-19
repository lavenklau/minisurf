#include "PeriodicMesher.h"

#include "geometry/bbox/BBox.h"
#include "CGAL/Simple_cartesian.h"
#include "CGAL/AABB_traits.h"
#include "CGAL/AABB_tree.h"
#include "CGAL/AABB_triangle_primitive.h"
#include "igl/writeOBJ.h"
#include "igl/extract_non_manifold_edge_curves.h"
#include "matlab/matlab_utils.h"
#include "cgal/mesh_intersection.h"
#include "igl/per_face_normals.h"
#include "igl/remove_unreferenced.h"
#include "cgal/cgal_utils.h"
#include "fundamental_forms.h"

std::string getPath(std::string);
extern msf::Real remesh_noise;
extern int rand_seed;
extern int debug_level;


typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef K::Triangle_3 Triangle;
typedef std::vector<Triangle>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;


//template<typename P> auto toOM(const P& p) { return OpenMesh::Vec3d(p[0], p[1], p[2]); }

template<typename P> auto toCGAL(const P& p) { return Point(p[0], p[1], p[2]); }

std::unique_ptr<Tree> generate_aabb_tree(const std::string& name, const msf::MeshSurface& mesh, msf::Real period = 2, msf::Real detach = 0.5) {
	static thread_local std::map<std::string, std::vector<Triangle>> triasmap;
	auto& trias = triasmap[name];
	trias.clear();
	Eigen::Vector3<msf::Real> border, Border;
	border.setConstant(period / 2 * 0.8);
	Border.setConstant(period / 2 * 1.2);
	msf::BBox bb(-border, border), BB(-Border, Border);
	msf::OM_Mesh m;
	// generate aabb tree
	for (auto& f : mesh.faces()) {
		std::array<OpenMesh::Vec3d, 3> Pts;
		int counter = 0;
		for (auto fv : mesh.fv_range(f)) {
			Pts[counter] = mesh.point(fv);
			counter++;
		}
		for (int i = 0; i < 5; i++) {
			auto vij = Pts[(i + 1) % 3] - Pts[i % 3];
			for (int k = 0; k < 3; k++) {
				if (vij[k] < -detach) {
					Pts[(i + 1) % 3][k] += period;
				} else if (vij[k] > detach) {
					Pts[i % 3][k] += period;
				}
			}
		}
		if (debug_level > 1) {
			auto v0 = m.add_vertex(Pts[0]);
			auto v1 = m.add_vertex(Pts[1]);
			auto v2 = m.add_vertex(Pts[2]);
			m.add_face(v0, v1, v2);
		}
		trias.emplace_back(Triangle(toCGAL(Pts[0]), toCGAL(Pts[1]), toCGAL(Pts[2])));
		Eigen::Vector3<msf::Real> tc = (msf::toEigen(Pts[0]) + msf::toEigen(Pts[1]) + msf::toEigen(Pts[2])) / 3;

		if (bb.isOut(tc)) {
			for (int x = -1; x <= 1; x++) {
				auto dupPts = Pts;
				for (int k = 0; k < 3; k++) dupPts[k][0] = Pts[k][0] + x * period;
				for (int y = -1; y <= 1; y++) {
					for (int k = 0; k < 3; k++) dupPts[k][1] = Pts[k][1] + y * period;
					for (int z = -1; z <= 1; z++) {
						for (int k = 0; k < 3; k++) dupPts[k][2] = Pts[k][2] + z * period;
						auto new_tc = tc;
						new_tc += Eigen::Vector3<msf::Real>(x, y, z) * period;
						if (!BB.isOut(new_tc)) {
							trias.emplace_back(Triangle(toCGAL(dupPts[0]), toCGAL(dupPts[1]), toCGAL(dupPts[2])));
						}
					}

				}
			}
		}
	}
	if (debug_level > 1) {
		OpenMesh::IO::write_mesh(m, getPath("aabb_mesh.obj"));
	}
	auto tree = std::make_unique<Tree>(trias.begin(), trias.end());
	return std::move(tree);
}

std::unique_ptr<Tree> generate_aabb_tree(const std::string& name, const Eigen::MatrixX3<msf::Real>& vlist, const Eigen::MatrixX3i& flist, msf::Real period = 2, msf::Real detach = 0.5) {
	static thread_local std::map<std::string, std::vector<Triangle>> triasmap;
	auto& trias = triasmap[name];
	trias.clear();
	// generate aabb tree
	for (int fi = 0; fi < flist.rows(); fi++) {
		std::array<Eigen::Vector3<msf::Real>, 3> Pts;
		for (int fv = 0; fv < 3; fv++) {
			Pts[fv] = vlist.row(flist(fi, fv)).transpose();
		}
		for (int i = 0; i < 5; i++) {
			auto vij = Pts[(i + 1) % 3] - Pts[i % 3];
			for (int k = 0; k < 3; k++) {
				if (vij[k] < - detach) {
					Pts[(i + 1) % 3][k] += period;
				} else if (vij[k] > detach) {
					Pts[i % 3][k] += period;
				}
			}
		}
		trias.emplace_back(Triangle(toCGAL(Pts[0]), toCGAL(Pts[1]), toCGAL(Pts[2])));
	}
	auto tree = std::make_unique<Tree>(trias.begin(), trias.end());
	return std::move(tree);
}

template<typename P>
P getClosePoint(Tree* tree, const P& p) {
	Point pq(p[0], p[1], p[2]);
	Point clst_point = tree->closest_point(pq);
	return  P(clst_point[0], clst_point[1], clst_point[2]);
}

bool msf::PeriodSurfaceMesh::is_cut_boundary(OM::SmartHalfedgeHandle he, const std::set<OM::SmartVertexHandle>& vcutset, bool strictOnSamePlane)
{
	bool boundary_edge = false;
	if ((vcutset.count(he.from()) && vcutset.count(he.to()))) {
		VertexFlag tag0, tag1;
		auto v0 = point(he.from()), v1 = point(he.to());
		tag0.set_period_boundary(1 - 1e-5, v0[0], v0[1], v0[2]);
		tag1.set_period_boundary(1 - 1e-5, v1[0], v1[1], v1[2]);
		if (!strictOnSamePlane) {
			boundary_edge = tag0.getPeriodFlag() & tag1.getPeriodFlag();
		} else {
			boundary_edge = tag0.getPeriodFlag() == tag1.getPeriodFlag();
		}
	}
	return boundary_edge;
}


template<typename P>
int period_degree(const P& p, double tol) {
	int deg = (std::abs(std::abs(p[0]) - 1) < tol) + (std::abs(std::abs(p[1]) - 1) < tol) + (std::abs(std::abs(p[2]) - 1) < tol);
	return deg;
}

std::vector<OpenMesh::SmartVertexHandle> msf::PeriodSurfaceMesh::periodic_remesh(int max_iter, const std::vector<OM::SmartVertexHandle>& vcut, Real tgt_len, int smooth_iter, int pert_iter, const MeshSurface* tgtmesh /*= nullptr*/)
{	
	if (has_boundary()) {
		std::cout << "\033[31m" << "Warning : periodically remeshing a mesh with boundary " << "\033[0m\n";
	}
	std::set<OM::SmartVertexHandle> vcutset(vcut.begin(), vcut.end());
	std::vector<HETag> hetags;

	if (debug_level > 1) savePeriodicMesh(getPath("rmshinput.obj"), vcutset);

	// deal with contradict cut set
	{
		for (auto fh : faces()) {
			auto fv = getFaceVertexHandle(fh);
			VertexFlag fvflg[3];
			int counter = 0;
			for (int k = 0; k < 3; k++) { counter += vcutset.count(fv[k]) && fv[k].is_valid() && !fv[k].deleted(); }
			double smin = 1e30, smax = -1e30;
			int kmin = -1, kmax = -1;
			if (counter != 3) continue;
			if (aspect_ratio(fh) < 100) continue;
#if 0
			std::cout << "Collapsing face ";
			for (int k = 0; k < 3; k++) {
				auto p = point(fv[k]);
				std::cout << p << " ";
				fvflg[k].set_period_boundary(0.9999, p[0], p[1], p[2]);
				double s = (toEigen(point(fv[k])) - Eigen::Vector3<Real>(-1, -1, -1)).cwiseAbs().minCoeff();
				if (s < smin) { smin = s; kmin = k; }
			}
			std::cout << std::endl;
			auto he1 = find_halfedge(fv[(kmin + 1) % 3], fv[kmin]);
			auto he2 = find_halfedge(fv[(kmin + 2) % 3], fv[kmin]);
			auto he2des = he1.next().opp();
			if (make_period(calc_edge_vector(he1)).norm() < make_period(calc_edge_vector(he2)).norm()) {
				if(is_collapse_ok(he1)) collapse(he1);
			} else {
				if (is_collapse_ok(he2)) collapse(he2);
			}
			auto vhdst = he2des.to();
			auto pmin = point(vhdst);
			VertexFlag vminflg;
			vminflg.set_period_boundary(fvflg[0].getPeriodFlag() | fvflg[1].getPeriodFlag() | fvflg[2].getPeriodFlag());
			//for (int k = 0; k < 3; k++) {
			//	if (vminflg.is_min_period(k)) { pmin[k] = -1; }
			//}
			//set_point(vhdst, pmin);
#else
			OM::SmartHalfedgeHandle he[3];
			he[0] = fh.halfedge(); he[1] = he[0].next(); he[2] = he[1].next();
			for (int i = 0; i < 3; i++) {
				Real s = make_period(calc_edge_vector(he[i])).norm();
				if (s < smin) { smin = s; kmin = i; }
			}
			auto p0 = toEigen(point(he[kmin].from()));
			auto p1 = toEigen(point(he[kmin].to()));
			//Eigen::Vector3<Real> v111(1, 1, 1);
			//if ((p0.cwiseAbs() - v111).cwiseAbs().minCoeff() < (p1.cwiseAbs() - v111).cwiseAbs().minCoeff()) {
			//	he[kmin] = he[kmin].opp();
			//}
			//std::cout << "deg from = " << period_degree(point(he[kmin].from()), 1e-5) << std::endl;
			//std::cout << "deg to   = " << period_degree(point(he[kmin].to()), 1e-5) << std::endl;
			if (period_degree(point(he[kmin].from()), 1e-5) > period_degree(point(he[kmin].to()), 1e-5)) {
				he[kmin] = he[kmin].opp();
			}
			std::cout << "Collapsing " << point(he[kmin].from()) << " to " << point(he[kmin].to()) << std::endl;
			if (is_collapse_ok(he[kmin])) {
				//std::cout << "face edge vec = \n" << face_edge_vector(fh) << std::endl;;
				//std::cout << "p = \n";
				//for (auto fv : fv_range(fh)) {
				//	std::cout << point(fv) << std::endl;
				//}
				//std::cout << "face aspect ratio = " << aspect_ratio(fh) << std::endl;
				collapse(he[kmin]); 
			}
#endif
		}
	}

	deduce_hetag(vcutset, hetags);

	Real edge_sum = 0;
	for (auto eh : edges()) {
		auto h0 = eh.halfedge(0);
		auto h0v = toEigen(calc_edge_vector(h0));
		if (hetags[h0.idx()].has_trans()) {
			h0v -= hetags[h0.idx()].torus_trans;
		}
		edge_sum += h0v.norm();
	}
	edge_sum /= n_edges();

	//Real lmin = std::min(edge_sum, tgt_len) * 0.8;
	//Real lmax = std::min(edge_sum, tgt_len) * 1.3333;
	Real lmin = tgt_len * 0.8;
	Real lmax = tgt_len * 1.3333;


	std::unique_ptr<Tree> tree;
	if (tgtmesh) {
		tree = generate_aabb_tree("remesh", *tgtmesh, 2, 0.5);
	} else {
		tree = generate_aabb_tree("remesh", *this, 2, 0.5);
	}

	auto [V, F] = getVFlist();

	bool output_step = debug_level > 1;
	//bool output_step = true;
	if (output_step) savePeriodicMesh(getPath("init.obj"), vcutset);

	for (int iter = 0; iter < max_iter; iter++) {
		//std::tie(V, F) = getVFlist();
		//igl::writeOBJ(getPath("split0.obj"), V, F);
		// split edge
		for (auto eh : edges()) {
			auto he = eh.halfedge(0);
			// do not split periodic boundary edges
			bool boundary_edge = is_cut_boundary(he, vcutset);
			// reverse halfedge direction, make sure it points to cutset
			if (vcutset.count(he.from())) he = eh.halfedge(1);
			auto v0 = point(he.from());
			auto v1 = point(he.to());
			auto hvec = v1 - v0;
			auto hvec_per = make_period(hvec, 2, 0.5);
			auto len = hvec_per.norm();
			if (len > lmax) {
				auto mid = v0 + hvec_per / 2;
				for (int k = 0; k < 3; k++) mid[k] = mid[k] < -1 - 1e-6 ? mid[k] + 2 : mid[k];
				auto vc = add_vertex(mid);
				split_edge(eh, vc);
				if (boundary_edge) vcutset.insert(vc);
			}
		}

		//hetags.resize(n_halfedges());

		//std::tie(V, F) = getVFlist();

		//igl::writeOBJ(getPath("split0.obj"), V, F);
		//OM::IO::write_mesh(*this, getPath("split0.obj"));

		deduce_hetag(vcutset, hetags);
		if (output_step) savePeriodicMesh(getPath("split.obj"), vcutset);
		if (output_step) { std::ofstream ofs(getPath("cutset")); for (auto h : vcutset) ofs << point(h) << std::endl; }

		// collapse edge
		for (auto he : halfedges()) {
			if (status(he).deleted() || !he.is_valid())continue;
			auto vh = he.from();
			auto [o, ring] = find1ring(vh, 0.5);
			auto [vhring, hering] = find1ring(vh);
			int id = 0;
			for (int i = 0; i < hering.size(); i++) { if (hering[i] == he) id = i; }
			auto v1 = ring[id];
			Real henorm = (v1 - o).norm();
			Real r_max = -1e30;
			for (int i = 0; i < ring.size(); i++) {
				if (i == id) continue;
				r_max = std::max((ring[i] - v1).norm(), r_max);
			}
			// do not create too long edge, but must collapse very short edge or sliver face on boundary
#if 0
			if (r_max > lmax && henorm > 5e-2 * lmin) {
				Real ar = aspect_ratio(he.face());
				if (ar < 50 && ar>0) continue;
			}
#else
			if (r_max > lmax && henorm > 5e-2 * lmin) { continue; }
#endif
			if (is_cut_boundary(he, vcutset, true) || !vcutset.count(he.from())) {
				if (is_collapse_ok(he)) {
					OM::Vec3d p = { 1,-0.103648,1 };
					if ((point(he.from()) - p).norm() < 1e-3 || (point(he.to()) - p).norm() < 1e-3) {
						std::cout << "he from = " << point(he.from()) << ", to = " << point(he.to()) << std::endl;
						std::cout << "is_cut_boundary ? " << is_cut_boundary(he, vcutset, true) << std::endl;
						std::cout << "cut from = " << vcutset.count(he.from()) << ", cut to = " << vcutset.count(he.to()) << std::endl;
					}
					vcutset.insert(he.from());
					collapse(he);

				}
			}
		}

		// remove redundant cutset points
		// vertices of one face all lies on boundary
		// leak some boundary nodes
		// We should consider the followinng cases:
		// 1. slider on corner
		// 2. two adjacent slier : 	https://s2.loli.net/2024/07/12/BCY1flNiuQRKZ7g.png
		if (1) {
			for (auto vh : vcutset) {
				if (vh.deleted() || !vh.is_valid()) continue;
				for (auto voh : vh.outgoing_halfedges()) {
					if (vcutset.count(voh.to())) {
						if (vcutset.count(voh.next().to())) {
							auto fh = voh.face();
							if (fh.deleted() || !fh.is_valid()) continue;
							if (aspect_ratio(fh) < 100) continue;
							int degmin = 1000;
							OM::SmartVertexHandle fvmin;
#if 0
							// the tip point of sliver face has degree 2, while two corners have degree 3
							for (auto fv : fh.vertices()) {
								int deg = 0;
								for (auto fvv : fv.vertices()) {
									if (vcutset.count(fvv)) { deg++; }
								}
								if (deg < degmin) { degmin = deg; fvmin = fv; }
							}
							// This condition is not always true 
#else
							Real hmax = -1e30;
							OM::SmartHalfedgeHandle fhmax;
							for (auto he : fh.halfedges()) {
								auto hv = make_period(calc_edge_vector(he));
								if (hv.norm() > hmax) {
									hmax = hv.norm(); fhmax = he;
								}
							}
							fvmin = fhmax.next().to();
#endif
							std::cout << "removinng " << point(fvmin) << " from cutset" << std::endl;
							vcutset.erase(fvmin);
							break;
						}
					}
				}
			}
		}
		//garbage_collection();
		deduce_hetag(vcutset, hetags);
		if (output_step) { std::ofstream ofs(getPath("cutset")); for (auto h : vcutset) ofs << point(h) << std::endl; }
		//OM::IO::write_mesh(*this, getPath("collapse0.obj"));
		//std::tie(V, F) = getVFlist();
		//igl::writeOBJ(getPath("collapse0.5.obj"), V, F);
		if (output_step) savePeriodicMesh(getPath("collapse.obj"), vcutset);

		// remove boundary sliver triangles

		// valence balance
		for (auto eh : edges()) {
			bool debug = false;
			if (status(eh).deleted() || !eh.is_valid()) continue;
			auto he = eh.halfedge(0);
			bool tri_boundary_edge = false;
			if (vcutset.count(he.from()) && vcutset.count(he.to())) {
				//continue;
				if (!(vcutset.count(he.next().to()) ^ vcutset.count(he.opp().next().to()))) continue;
				// For corrner edges, i.e., the opposite vertex is a bb-corner vertex or a bb-edge vertex
				if (is_cut_boundary(he, vcutset)) continue;
				////std::cout << "found tri-boundary edge" << point(he.from()) << ", " << point(he.to()) << std::endl;
				tri_boundary_edge = true;
			}
			else if (vcutset.count(he.next().to()) && vcutset.count(he.opp().next().to())) continue;
	
			int v_d0 = he.from().valence(), v_d1 = he.to().valence();
			int v_lef = he.next().to().valence(), v_rig = he.opp().next().to().valence();
			int old_valerr = std::abs(v_d0 - 6) + std::abs(v_d1 - 6) + std::abs(v_lef - 6) + std::abs(v_rig - 6);
			int new_valerr = std::abs(v_d0 - 1 - 6) + std::abs(v_d1 - 1 - 6) + std::abs(v_lef + 1 - 6) + std::abs(v_rig + 1 - 6);;
			/* valence based flip */
			if (old_valerr > new_valerr) {
				if (is_flip_ok(eh)) flip(eh);
			} /* degree based flip */
			else /*if (old_valerr == new_valerr)*/ {
				auto vh = he.from();
				auto [d0, d1, lef, rig] = find1diag(he, 0.8);
				Real deg[2] = {
					//angle(lef,d1,rig) + angle(rig,d0,lef),
					//angle(d0,lef,d1) + angle(d0,rig,d1)
					std::max(angle(lef, d1, rig) , angle(rig, d0, lef)),
					std::max(angle(d0, lef, d1)  , angle(d0, rig, d1))
				};
				if (std::abs(deg[1] - M_PI) < 0.1 && tri_boundary_edge) {
					if (is_flip_ok(eh)) { flip(eh); }
				}
				else if (deg[0] < deg[1] && old_valerr == new_valerr) {
					if (is_flip_ok(eh)) flip(eh);
				}
			}
		}

		//deduce_hetag(vcutset, hetags);
		if (output_step) savePeriodicMesh(getPath("valence.obj"), vcutset);
		//std::tie(V, F) = getVFlist();
		//savePeriodicMesh(getPath("valence.obj"), V, vcutset, hetags);
		//igl::writeOBJ(getPath("relax0.obj"), V, F);

		// tangential relaxation
		for (auto vh : vertices()) {
			if (!vh.is_valid() || status(vh).deleted()) { continue; }
			bool vclamp[3] = { false,false,false };
			if (vcutset.count(vh)) {
				auto p = point(vh);
				for (int k = 0; k < 3; k++) { vclamp[k] = (std::abs(p[k] + 1) < 1e-5); }
			}
			auto [o, ring] = find1ring(vh, 0.5);
			if (ring.size() == 0) continue;
			auto [vn, bary] = getVertexNormal(o, ring);
			CHECK_NAN(vn);
			bary = bary - (bary - o).dot(vn) * vn;
			CHECK_NAN(bary);
			for (int k = 0; k < 3; k++) { if (vclamp[k]) bary[k] = -1; }
			for (int k = 0; k < 3; k++) {
				if (bary[k] < -1) { bary[k] += 2; }
			}
			set_point(vh, OM::Vec3d(bary[0], bary[1], bary[2]));
		}
		//deduce_hetag(vcutset, hetags);
		//std::tie(V, F) = getVFlist();
		//igl::writeOBJ(getPath("relax0.obj"), V, F);
		if (output_step) savePeriodicMesh(getPath("relax.obj"), vcutset);

		// project to mesh
		for (auto vh : vertices()) {
			if (!vh.is_valid() || status(vh).deleted()) { continue; }
			if (is_isolated(vh)) continue;
			bool vclamp[3] = { false,false,false };
			if (vcutset.count(vh)) {
				auto p = point(vh);
				for (int k = 0; k < 3; k++) { vclamp[k] = std::abs(p[k] + 1) < 1e-6; }
			}
			auto p = point(vh);
			auto proj = getClosePoint(tree.get(), p);
			CHECK_NAN(toEigen(proj));
			int nclamp = vclamp[0] + vclamp[1] + vclamp[2];
			if (nclamp > 1) {
				auto nproj = p - proj;
				auto oldproj = proj; 
				for (int k = 0; k < 3; k++) {
					if (!vclamp[k]) {
						Real new_proj_k = (proj.dot(nproj) - (-1) * nproj[(k + 1) % 3] - (-1) * nproj[(k + 2) % 3]) / nproj[k];
						if (new_proj_k >= -1 && new_proj_k <= 1) { proj[k] = new_proj_k; }
						else {
							std::cout << "Warning : proj " << proj << " is pulled out of range\n";
							proj = oldproj; break;
						}
						break;
					}
				}
				for (int k = 0; k < 3; k++) { if (vclamp[k]) proj[k] = -1; }
			} else if (nclamp > 0) {
				//std::cout << "nclamp = " << nclamp << std::endl;
				//std::cout << "p    = " << p << std::endl;
				//std::cout << "proj = " << proj << std::endl;
				int kclamp = vclamp[0] ? 0 : (vclamp[1] ? 1 : 2);
				auto nproj = p - proj;
				Eigen::Matrix<double, 5, 5> A;
				A.setZero();
				A.block<3, 3>(0, 0).setIdentity();
				A.block<3, 1>(0, 3) = toEigen(nproj);
				A(kclamp, 4) = 1;
				A.block<1, 3>(3, 0) = toEigen(nproj).transpose();
				A(4, kclamp) = 1;
				Eigen::FullPivHouseholderQR<decltype(A)> qr(A);
				Eigen::Matrix<double, 5, 1> b;
				b << toEigen(p), nproj.dot(proj), -1;
				Eigen::Vector<double, 5> plam = A.householderQr().solve(b);
				auto newproj = toOM(plam.topRows(3));
				newproj[kclamp] = -1;
				auto oldproj = proj; oldproj[kclamp] = -1;
				for (int k = 0; k < 3; k++) {
					if (newproj[k] >= -1 && newproj[k] <= 1 && (std::abs)(newproj[k] - proj[k]) < 0.05) {
						proj[k] = newproj[k];
					} else {
						std::cout << "Warning : proj " << proj << " is pulled out of range\n";
						proj = oldproj; break;
					}
				}
				//std::cout << "dst  = " << proj << std::endl;
			}
			for (int k = 0; k < 3; k++) { if (proj[k] < -1) { proj[k] += 2; } }
			if ((p - proj).norm() < 0.5) { set_point(vh, proj); }
		}

		//std::tie(V, F) = getVFlist();

		//OM::IO::write_mesh(*this, getPath("project0.obj"));
		if (output_step) savePeriodicMesh(getPath("project.obj"), vcutset);

		//deduce_hetag(vcutset, hetags);
		//savePeriodicMesh(getPath("meshiter.obj"), vcutset);
	}

	if (has_boundary()) {
		std::cout << "\033[31m" << "Warning : remesh result contains boundary: " << "\033[0m\n";
		for (auto vh : vertices()) {
			if (vh.deleted() || !vh.is_valid()) continue;
			if (vh.is_boundary()) std::cout << point(vh) << std::endl;
		}
	}



	deduce_hetag(vcutset, hetags);

	// fix outrange vertex and clamp boundary vertices
	BBox bb(Eigen::Vector3<Real>(-1.0001, -1.0001, -1.0001), Eigen::Vector3<Real>(1.0001, 1.0001, 1.0001));
	for (auto vh : vertices()) {
		auto p = point(vh);
		for (int i = 0; i < 3; i++) {
			if (p[i] < -1.0001) p[i] += 2;
			if (p[i] > 1.0001) p[i] -= 2;
			if (std::abs(p[i] - 1) < 2e-5) p[i] = 1;
			if (std::abs(p[i] + 1) < 2e-5) p[i] = -1;
		}
		set_point(vh, p);
	}


	if (debug_level > 0) {
		savePeriodicMesh(getPath("beforesmooth.obj"), vcutset);
	}
	periodic_smooth(smooth_iter, 1);


	// fix sliver
	if (hasSliver(10)) {
		std::cout << "Warning : mesh contains sliver triangle" << std::endl;
		// fixing boundary sliver
		BBox bb(Eigen::Vector3<Real>(-0.9999, -0.9999, -0.9999), Eigen::Vector3<Real>(0.9999, 0.9999, 0.9999));
		BBox bbDomain(Eigen::Vector3<Real>(-1.0001, -1.0001, -1.0001), Eigen::Vector3<Real>(1.0001, 1.0001, 1.0001));
#if 1
		for (auto eh : edges()) {
			if (!eh.is_valid()) continue;
			if (eh.deleted()) continue;
			for (auto he : { eh.halfedge(0),eh.halfedge(1) }) {
				if (!bb.isOut(toEigen(point(he.from()))) || !bb.isOut(toEigen(point(he.to())))) continue;
				Real ar = aspect_ratio(he.face());
				if (ar > 100) {
					std::cout << "ar = " << ar << ", ";
					auto hemax = he;
					auto vhe = calc_edge_vector(he);
					Real max_len = make_period(vhe).norm();;
					auto vhe1 = calc_edge_vector(he.next());
					Real l1 = make_period(vhe1).norm();
					if (l1 > max_len) { max_len = l1; hemax = he.next(); }
					auto vhe2 = calc_edge_vector(he.next().next());
					Real l2 = make_period(vhe2).norm();
					if (l2 > max_len) { max_len = l2; hemax = he.next().next(); }
					auto pmid = point(hemax.from()) + make_period(calc_edge_vector(hemax)) / 2;
					if (bbDomain.isOut(toEigen(pmid))) pmid = point(hemax.to()) - make_period(calc_edge_vector(hemax)) / 2;
					auto vhmid = add_vertex(pmid);
					if (bb.isOut(toEigen(pmid))) { vcutset.insert(vhmid); }
					auto vhopp = hemax.next().to();
					split_edge(hemax.edge(), vhmid);
					auto hec = find_halfedge(vhopp, vhmid);
					if (is_collapse_ok(hec)) {
						std::cout << "Collapsing " << point(hec.from()) << " to " << point(hec.to()) << std::endl;
						collapse(hec);
					}
					break;
				}
			}
		}
		deduce_hetag(vcutset, hetags);
#endif
	}

	if (pert_iter) pertubMesh(pert_iter, remesh_noise, rand_seed);

	{
		//savePeriodicMesh(getPath("remesh.obj"), vcutset);
	}

	return std::vector<OM::SmartVertexHandle>(vcutset.begin(), vcutset.end());
}

void msf::PeriodSurfaceMesh::deduce_hetag(
	std::set<OM::SmartVertexHandle>& vcut, std::vector<HETag>& hetags, Real detach_distance/* = 0.5*/)
{
	std::unordered_set<Eigen::Vector3<Real>, Vector3Hash<Real>> vcutpos;
	for (auto vh : vcut) {
		if (!status(vh).deleted() && vh.is_valid()) { vcutpos.insert(toEigen(point(vh))); }
	}
	garbage_collection();

	std::set<OM::SmartVertexHandle> vcutnew;
	for (auto vh : vertices()) {
		if (status(vh).deleted() || !vh.is_valid()) { continue; }
		if (vcutpos.count(toEigen(point(vh)))) { vcutnew.insert(vh); }
	}
	vcut = vcutnew;

	hetags = std::vector<HETag>(n_halfedges());

	for (auto vh : vcutnew) {
		if (status(vh).deleted() || !vh.is_valid()) { continue; }
		for (auto voh : voh_range(vh)) {
			auto vec = calc_edge_vector(voh);
			for (int k = 0; k < 3; k++) {
				if (vec[k] < -detach_distance) {
					hetags[voh.idx()].torus_trans[k] = -2;
				} else if (vec[k] > detach_distance) {
					hetags[voh.idx()].torus_trans[k] = 2;
				}
			}
			hetags[voh.idx()].updateFlag();
		}
	}
}

Eigen::Vector3<msf::Real> msf::PeriodSurfaceMesh::bary1ring(
	const Eigen::Vector3<Real>& v, const std::vector<Eigen::Vector3<Real>>& ring)
{
	throw std::runtime_error("unimplemented");
}

std::vector<OpenMesh::SmartVertexHandle> msf::PeriodSurfaceMesh::periodic_remesh(int max_iter, std::string filename, Real tgt_len, int smoot_iter, int pter_iter /*= 0*/)
{
	read(filename, true, false);
	MeshSurface tgtmsh = *this;
	auto [V, F] = tgtmsh.getVFlist();
	std::tie(V, F) = cgal_remesh(V, F);
	read(V, F, true, false);
	mergePeriodEdges();
	auto vcut = mergePeriodVertices();
	auto vcutset = periodic_remesh(max_iter, vcut, tgt_len, smoot_iter, pter_iter, &tgtmsh);
	return vcutset;
}


void msf::PeriodSurfaceMesh::periodic_smooth(int max_iter, int type, bool clamp_period_boundary /*= true*/)
{
	// forward Euler MCF
	if (type == 0) {
		const Real step = 0.05;
		const Real eps = 1e-5;
		std::vector<OM::SmartVertexHandle> vhlist;
		for (auto vh : vertices()) vhlist.push_back(vh);
		std::vector<std::array<bool, 3>> vclamped(n_vertices());
		for (int i = 0; i < vhlist.size(); i++) {
			auto vh = vhlist[i]; if (status(vh).deleted() || !vh.is_valid()) continue;
			auto [o, ring] = find1ring(vh, 0.5);
			if (ring.empty()) continue;
			std::array<bool, 3> clamped = { false,false,false };
			for (int k = 0; k < 3; k++) { clamped[k] = o[k] < -1 + eps; }
			vclamped[vh.idx()] = clamped;
		}
		std::vector<Eigen::Vector3<Real>> Hvlist(n_vertices());
		for (int i = 0; i < vhlist.size(); i++) { Hvlist[i].setZero(); }

		for (int iter = 0; iter < max_iter; iter++) {
			// compute mean curvature normal direction
#pragma omp parallel for
			for (int i = 0; i < vhlist.size(); i++) {
				auto vh = vhlist[i]; if (status(vh).deleted() || !vh.is_valid()) continue;
				auto [o, ring] = find1ring(vh, 0.5);
				if (ring.empty()) continue;
				Eigen::Vector3<Real> Hv;
				meanH(o, ring, Hv);
				Hvlist[i] = Hv;
			}
#pragma omp parallel for
			for (int i = 0; i < vhlist.size(); i++) {
				auto vh = vhlist[i]; if (status(vh).deleted() || !vh.is_valid()) continue;
				Eigen::Vector3<Real> Hv = Hvlist[i];
				auto pnew = point(vh) + step * toOM(Hv);
				for (int k = 0; k < 3; k++) { if (vclamped[vh.idx()][k]) pnew[k] = -1; }
				set_point(vh, pnew);
			}
		}
	}
	else if (type == 1) {
		const Real step = 0.05;
		const Real eps = 1e-5;
		//std::vector<OM::Vec3d> olist, Hlist;
		for (int iter = 0; iter < max_iter; iter++) {
			for (auto vh : vertices()) {
				if (status(vh).deleted() || !vh.is_valid()) continue;
				auto [o, ring] = find1ring(vh, 0.5);
				if (ring.empty()) continue;
				Eigen::Vector3<Real> Hv;
				meanH(o, ring, Hv);
				//olist.push_back(toOM(o));
				//Hlist.push_back(toOM(Hv));
				bool clamped[3] = { false,false,false };
				for (int k = 0; k < 3; k++) { clamped[k] = o[k] < -1 + eps; }
				auto pnew = toOM(o) + step * toOM(Hv);
				if (clamp_period_boundary) {
					for (int k = 0; k < 3; k++) { if (clamped[k]) pnew[k] = -1; }
				}
				set_point(vh, pnew);
			}
		}
		//{
		//	std::ofstream ofs(getPath("Holist"));
		//	for (int i = 0; i < olist.size(); i++) {
		//		ofs << olist[i] << " " << Hlist[i] << std::endl;
		//	}
		//}
	}
	// backward Euler MCF
	else if (type == 2) {
		throw std::runtime_error("not implemeneted");
	}
}

bool msf::PeriodSurfaceMesh::hasSliver(Real angThres, Real period /*= 2*/, Real deduce /*= 0.5*/)
{
	Real cosThres = std::cos(angThres / 180 * M_PI);
	for (auto fh : faces()) {
		OM::Vec3d hv[3];
		auto hiter = fh.halfedge();
		hv[0] = calc_edge_vector(hiter);
		hiter = hiter.next();
		hv[1] = calc_edge_vector(hiter);
		hiter = hiter.next();
		hv[2] = calc_edge_vector(hiter);
		for (int i = 0; i < 3; i++) {
			hv[i] = make_period(hv[i], period, deduce).normalized();
		}
		for (int i = 0; i < 3; i++) {
			if (std::abs(hv[i].dot(hv[(i + 1) % 3])) > cosThres) {
				return true;
			}
		}
	}
	return false;
}

Eigen::Matrix3<msf::Real> msf::PeriodSurfaceMesh::face_edge_vector(OM::FaceHandle fh, msf::Real period /*= 2*/, msf::Real deduce /*= 0.5*/)
{
	Eigen::Matrix3<Real> hv;
	auto hiter = halfedge_handle(fh);
	hv.col(0) = make_period(toEigen(calc_edge_vector(hiter)), period, deduce);
	hiter = next_halfedge_handle(hiter);
	hv.col(1) = make_period(toEigen(calc_edge_vector(hiter)), period, deduce);
	hiter = next_halfedge_handle(hiter);
	hv.col(2) = make_period(toEigen(calc_edge_vector(hiter)), period, deduce);
	return hv;
}

msf::Real msf::PeriodSurfaceMesh::aspect_ratio(OM::FaceHandle fh, msf::Real period /*= 2*/, msf::Real deduce /*= 0.5*/)
{
	auto hv = face_edge_vector(fh,period,deduce);
	Real a = hv.col(0).norm(), b = hv.col(1).norm(), c = hv.col(2).norm();
	double s = (a + b + c) / 2.0;
	double AR = (a * b * c) / (8.0 * (s - a) * (s - b) * (s - c));
	return AR;
}



std::tuple<Eigen::MatrixX3<msf::Real>, Eigen::VectorX<msf::Real>> msf::PeriodSurfaceMesh::periodic_mean_curvature(void)
{
	Eigen::MatrixX3<Real> Hvlist(n_vertices(), 3);
	Eigen::VectorX<Real > Avlist(n_vertices());
	const Real eps = 1e-5;
	for (auto vh : vertices()) {
		if (status(vh).deleted() || !vh.is_valid()) continue;
		auto [o, ring] = find1ring(vh, 0.5);
		if (ring.empty()) continue;
		Eigen::Vector3<Real> Hv;
		meanH(o, ring, Hv);
		auto [Av, nv] = area1ring(o, ring, 2);
		Avlist[vh.idx()] = Av;
		Hvlist.row(vh.idx()) = Hv.transpose();
	}
	return { Hvlist, Avlist };

}

std::tuple<msf::Real,Eigen::Vector3<msf::Real>> msf::PeriodSurfaceMesh::area1ring(const Eigen::Vector3<Real>& v, const std::vector<Eigen::Vector3<Real>>& ring, int type)
{
	if (type == 1) {
		Real As = 0;
		Eigen::Vector3<Real> nv(0, 0, 0);
		for (int i = 0; i < ring.size(); i++) {
			Eigen::Vector3<Real> fv = (ring[i] - v).cross(ring[(i + 1) % ring.size()] - ring[i]);
			Real A = fv.norm() / 2;
			nv += A * fv.normalized();
			As += A;
		}
		return { As / 3, nv.normalized() };
	}
	else if (type == 2) {
		Real As = 0;
		Eigen::Vector3<Real> nv(0, 0, 0);
		for (int i = 0; i < ring.size(); i++) {
			Eigen::Vector3<Real> a = ring[i] - v;
			Eigen::Vector3<Real> b = ring[(i + 1) % ring.size()] - v;
			Real e_sq[] = { a.squaredNorm(),b.squaredNorm(),(a - b).squaredNorm() };
			if (e_sq[0] + e_sq[1] >= e_sq[2] && e_sq[1] + e_sq[2] >= e_sq[0] && e_sq[2] + e_sq[0] >= e_sq[1]) {
				// add Voronoi area
				Real cosA = (e_sq[1] + e_sq[2] - e_sq[0]) / 2 / std::sqrt(e_sq[1] * e_sq[2]);
				Real cosB = (e_sq[0] + e_sq[2] - e_sq[1]) / 2 / std::sqrt(e_sq[0] * e_sq[2]);
				Real cotA = cosA / std::sqrt(1 - cosA * cosA);
				Real cotB = cosB / std::sqrt(1 - cosB * cosB);
				// ToDo: check this
				Real A = cotA * e_sq[0] / 8 + cotB * e_sq[1] / 8;
				As += A;
				nv += A * (a.cross(b).normalized());
			} else if (e_sq[0] + e_sq[1] < e_sq[2]) {
				// half area if incident angle is obtuse
				Real A = a.cross(b).norm() / 4;
				As += A;
				nv += A * (a.cross(b).normalized());
			} else /* has obtuse angle but not the incident one */{
				// quater area if others
				Real A = a.cross(b).norm() / 8;
				As += A;
				nv += A * (a.cross(b).normalized());
			}
		}
		return { As, nv.normalized() };
	}
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


std::tuple<Eigen::MatrixX3<msf::Real>, Eigen::MatrixX3i, Eigen::VectorXi> msf::PeriodSurfaceMesh::weldMesh(
	const Eigen::MatrixX3<Real>& V1, const Eigen::MatrixX3i& F1, const Eigen::MatrixX3<Real>& V2, const Eigen::MatrixX3i& F2)
{
	auto tree1 = generate_aabb_tree("tree1", V1, F1, 2, 0.5);
	auto tree2 = generate_aabb_tree("tree2", V2, F2, 2, 0.5);
	auto vlist1 = V1, vlist2 = V2;
	auto flist1 = F1, flist2 = F2;
	mesh_split(vlist1, flist1, V2, F2);
	mesh_split(vlist2, flist2, V1, F1);
	Eigen::MatrixX3<Real> Vmerge(vlist1.rows() + vlist2.rows(), 3);
	Eigen::MatrixX3i Fmerge(flist1.rows() + flist2.rows(), 3);
	Vmerge << vlist1, vlist2;
	for (int i = 0; i < flist2.size(); i++) { flist2.data()[i] += vlist1.rows(); }
	Fmerge << flist1, flist2;

	// 
	const Real eps = 2e-3;
	std::tie(Vmerge, Fmerge) = remove_dup_vertices(Vmerge, Fmerge, 1e-6);
	// extract non manifold edges
	std::map<UDEdge, char> edge_refs = extractNonManifoldUDEdges(Fmerge);
	std::map<int, int> collapse_chain;
	std::vector<std::pair<Eigen::Vector3<Real>, Eigen::Vector3<Real>>> e_list;
	for (auto eref : edge_refs) {
		// find start and end
		int v0 = -1, v1 = -1;
		auto e = eref.first;
		v0 = e.first;
		while (collapse_chain.count(v0)) { v0 = collapse_chain[v0]; }
		v1 = e.second;
		while (collapse_chain.count(v1)) { v1 = collapse_chain[v1]; }
		Eigen::Vector3<Real> vec0 = Vmerge.row(v0).transpose();
		Eigen::Vector3<Real> vec1 = Vmerge.row(v1).transpose();
		//std::cout << vec0.transpose() << " " << vec1.transpose() << std::endl;
		if ((vec0 - vec1).squaredNorm() < eps * eps) {
			if (v0 > v1) std::swap(v0, v1);
			collapse_chain[v0] = v1;
			e_list.push_back({ vec0,vec1 });
		}
	}
	//{ std::ofstream ofs(getPath("e_list"), std::ios::binary);
	//ofs.write((char*)e_list.data(), sizeof(e_list[0])* e_list.size()); }
	std::cout << "find " << collapse_chain.size() << " length collapse chain" << std::endl;
	int counter = 0;
	for (int i = 0; i < Fmerge.rows(); i++) {
		for (int k = 0; k < 3; k++) {
			int vnew = Fmerge(i, k);
			while (collapse_chain.count(vnew)) { vnew = collapse_chain[vnew]; }
			Fmerge(i, k) = vnew;
		}
	}

	std::tie(Vmerge, Fmerge) = remove_dup_vertices(Vmerge, Fmerge);

#if 1
	auto udedgs = extractNonManifoldUDEdges(Fmerge);
	std::set<int> vcomplex;
	for (auto e : udedgs) {
		vcomplex.insert(e.first.first);
		vcomplex.insert(e.first.second);
	}
	std::vector<Eigen::Vector3<Real>> vappend;
	for (int i = 0; i < Fmerge.rows(); i++) {
		for (int k = 0; k < 3; k++) {
			int vid = Fmerge(i, k);
			if (vcomplex.count(vid)) {
				Fmerge(i, k) = Vmerge.rows() + vappend.size();
				vappend.push_back(Vmerge.row(vid).transpose());
			}
		}
	}
	
	int n_old = Vmerge.rows();
	Vmerge.conservativeResize(n_old + vappend.size(), 3);
	for (int i = 0; i < vappend.size(); i++) {
		Vmerge.row(n_old + i) = vappend[i].transpose();
	}

	Eigen::VectorXi I, J;
	Eigen::MatrixX3<Real> vnew;
	Eigen::MatrixX3i fnew;
	igl::remove_unreferenced(Vmerge, Fmerge, vnew, fnew, I);

	//igl::writeOBJ(getPath("mnew.obj"), vnew, fnew);
	auto components = connected_components(vnew, fnew);
	//{std::ofstream ofs(getPath("components"));
	//for (int i = 0; i < components.size(); i++) {
	//	for (int k = 0; k < components[i].size(); k++) { ofs << components[i][k] << " "; }
	//	ofs << std::endl;
	//}}
	std::vector<Eigen::MatrixX3<Real>> mvlist(components.size());
	std::vector<Eigen::MatrixX3i> mflist(components.size());
	int totalVRows = 0, totalFrows = 0;
	for (int i = 0; i < components.size(); i++) {
		mvlist[i] = vnew;
		mflist[i].resize(components[i].size(), 3);
		for (int k = 0; k < components[i].size(); k++) {
			mflist[i].row(k) = fnew.row(components[i][k]);
		}
		std::tie(mvlist[i], mflist[i]) = remove_dup_vertices(mvlist[i], mflist[i]);
		//igl::writeOBJ(getPath("mlist.obj"), mvlist[i], mflist[i]);
		totalVRows += mvlist[i].rows();
		totalFrows += mflist[i].rows();
	}
	vnew.resize(totalVRows, 3);
	fnew.resize(totalFrows, 3);
	int Vcur = 0, Fcur = 0;
	for (int i = 0; i < components.size(); i++) {
		vnew.block(Vcur, 0, mvlist[i].rows(), 3) = mvlist[i]; 
		mflist[i].array() += Vcur;
		fnew.block(Fcur, 0, mflist[i].rows(), 3) = mflist[i];
		Vcur += mvlist[i].rows();
		Fcur += mflist[i].rows();
	}
	
	
#else
	mesh_split_non_manifold_vertices(Vmerge, Fmerge);
#endif

	Eigen::VectorXi f2mid(fnew.rows(), 1);
	for (int i = 0; i < f2mid.size(); i++) {
		Eigen::Vector3i fv = fnew.row(i).transpose();
		Eigen::Vector3<Real> c(0, 0, 0);
		for (int k = 0; k < 3; k++) {
			c += vnew.row(fv[k]).transpose();
		}
		c /= 3;
		auto proj1 = getClosePoint(tree1.get(), c);
		auto proj2 = getClosePoint(tree2.get(), c);
		Real d1 = (proj1 - c).squaredNorm();
		Real d2 = (proj2 - c).squaredNorm();
		if (d1 < d2) {
			f2mid[i] = 1;
		} else {
			f2mid[i] = 2;
		}
	}
	
	return { vnew, fnew, f2mid };
}

void msf::PeriodSurfaceMesh::pertubMesh(int max_iter, Real strength, unsigned int seed)
{

	//Real edge_sum = 0;
	//for (auto eh : edges()) {
	//	auto h0 = eh.halfedge(0);
	//	auto h0v = toEigen(calc_edge_vector(h0));
	//	h0v = make_period(h0v);
	//	edge_sum += h0v.norm();
	//}
	//edge_sum /= n_edges();

	auto tree = generate_aabb_tree("remesh", *this, 2, 0.5);

	srand(seed);

	Eigen::VectorX<Real> wA(20);

	for (int iter = 0; iter < max_iter; iter++) {
		for (auto vh : vertices()) {
			auto p = point(vh);
			VertexFlag vflag;
			vflag.set_period_boundary(0.9999, p[0], p[1], p[2]);
			if (vflag.is_period_boundary()) continue;
			wA.setRandom();
			wA *= strength;
			auto [o, ring] = find1ring(vh, 0.8);
			std::vector<Eigen::Vector3<Real>> Ac(ring.size());
			Eigen::Vector3<Real> cring(0, 0, 0);
			Real w_s = 0;
			for (int i = 0; i < ring.size(); i++) {
				Ac[i] = (ring[(i + 1) % ring.size()] - o).cross(ring[i] - o);
				Real ai = Ac[i].norm() / 2;
				Eigen::Vector3<Real> c = (ring[(i + 1) % ring.size()] + ring[i] + o) / 3;
				cring += c * ai * (1 + wA[i]);
				w_s += (1 + wA[i]) * ai;
			}
			cring /= w_s;
			auto proj = getClosePoint(tree.get(), cring);

			bool invalid_pert = false;
			for (int i = 0; i < ring.size(); i++) {
				Eigen::Vector3<Real> newAc = (ring[(i + 1) % ring.size()] - proj).cross(ring[i] - proj);
				if (newAc.dot(Ac[i] / Ac[i].squaredNorm()) < 0.1) { invalid_pert = true; break; }
			}
			if (invalid_pert) continue;
			set_point(vh, toOM(proj));
		}
	}
}

// Adapted from source code of Repulsive Surface
//  https://github.com/icethrush/repulsive-surfaces

using namespace msf;

double msf::PeriodSurfaceMesh::periodEdgeLength(OM::EdgeHandle eh) const
{
	auto vec = calc_edge_vector(halfedge_handle(eh, 0));
	vec = make_period(vec, 2, 1);
	return vec.norm();
}

inline double edgeDihedralAngle(PeriodSurfaceMesh& m, OM::SmartEdgeHandle e) {
   // WARNING: Logic duplicated between cached and immediate version
	if (e.is_boundary() || !e.is_valid()) { return 0.; }

	auto [d0, d1, left, right] = m.find1diag(e.h0(), 1);

	Eigen::Vector3d N1 = (d1 - d0).cross(left - d0).normalized();
	Eigen::Vector3d N2 = (right - d0).cross(d1 - d0).normalized();
    auto pTail = d0;
	auto pTip = d1;
	auto edgeDir = (pTip - pTail).normalized();

	return atan2(edgeDir.dot(N1.cross(N2)), N1.dot(N2));
}


inline double vertexMeanCurvature(PeriodSurfaceMesh& m, OpenMesh::SmartVertexHandle v) {
   // WARNING: Logic duplicated between cached and immediate version
   double meanCurvature = 0.;
   for (auto he : v.outgoing_halfedges_ccw()) {
	   double len = m.periodEdgeLength(he.edge());
	   double alpha = edgeDihedralAngle(m, he.edge());
	   meanCurvature += alpha * len / 2.;
   }
   return meanCurvature/2.;
}

double getSmoothMeanCurvature(msf::PeriodSurfaceMesh& m, OM::SmartVertexHandle v)
{
	double A = m.period_vertex_dual_area(v);
	double S = vertexMeanCurvature(m, v);
	double K = S / A;
	return K;
}


double vertexAngleSum(PeriodSurfaceMesh& m, OM::SmartVertexHandle v) {
	double sum = 0;
	for (auto ch : v.incoming_halfedges()) {
		sum += m.period_sector_angle(ch);
	}
	return sum;
}

double vertexGaussianCurvature(PeriodSurfaceMesh& m, OM::SmartVertexHandle v) {
	return 2 * M_PI - vertexAngleSum(m, v);
}

double getSmoothGaussianCurvature(PeriodSurfaceMesh& m, OM::SmartVertexHandle v)
{
	double A = m.period_vertex_dual_area(v);
	double S = vertexGaussianCurvature(m, v);
	double K = S / A;
	return K;
}

// flatLength: specifies how long the target edge length should be in flat regions
// epsilon: controls how much variation in target length occurs due to curvature
double findGaussianTargetL(PeriodSurfaceMesh& m, OM::SmartEdgeHandle e, double flatLength, double epsilon)
{
	// Areas and curvatures are already required in main.cpp
	double averageK = 0;
	for (auto v : { e.vertex(0), e.vertex(1) }) {
		averageK += getSmoothGaussianCurvature(m, v);
	}
	averageK /= 2;
	double L = flatLength * epsilon / (fabs(sqrt(averageK)) + epsilon);
	return L;
	// return flatLength;
}

double findTotalCurvatureTargetL(PeriodSurfaceMesh& m, OM::SmartEdgeHandle e, double flatLength, double epsilon)
{
	// Areas and curvatures are already required in main.cpp
	double averageK = 0;
	for (auto v : { e.vertex(0), e.vertex(1) }) {
		auto [o, ring] = m.find1ring(v, 1);
		Compile1ring vring(o, ring);
		averageK += (4 * vring.H * vring.H - 2 * vring.K);
	}
	averageK /= 2;
	double L = flatLength * epsilon / (fabs(sqrt(averageK)) + epsilon);
	return L;
	// return flatLength;
}

double findMeanTargetL(PeriodSurfaceMesh& m, OpenMesh::SmartEdgeHandle e, double flatLength, double epsilon)
{
	// Areas and curvatures are already required in main.cpp
	double averageH = 0;
	auto he = e.h0();
#if 0
	for (auto vh : { he.from(), he.to() }) {
		auto [o, ring] = m.find1ring(vh, 1);
		auto [hv, A] = m.meanH(o, ring);
		averageH += hv.norm() / A;
	}
	averageH /= 2;
	double L = flatLength * epsilon / (fabs(averageH) + epsilon);
#else
	for (auto vh : { he.from(), he.to() }) {
		averageH += getSmoothMeanCurvature(m, vh);
	}
	averageH /= 2;
	double L = flatLength * epsilon / (fabs(averageH) + epsilon);
	return L;
#endif
	return L;
	// return flatLength;
}

template<typename Vector3>
inline double angle(const Vector3& u, const Vector3& v) {
	return std::acos(std::fmax(-1., std::fmin(1., u.normalized().dot(v.normalized()))));
}
template<typename Vector3>
inline double diamondAngle(Vector3 a, Vector3 b, Vector3 c, Vector3 d) // dihedral angle at edge a-b
{
	Vector3 n1 = (b - a).cross(c - a);
	Vector3 n2 = (b - d).cross(a - d);
	return M_PI - ::angle(n1, n2);
}

template<typename Vector3>
inline bool checkFoldover(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& x, double angle)
{
	return diamondAngle(a, b, c, x) < angle;
}

bool shouldCollapse(PeriodSurfaceMesh& m, OM::SmartEdgeHandle e)
{
	using Halfedge = OM::SmartHalfedgeHandle;
	using Vertex = OM::SmartVertexHandle;
	using Vector3 = OM::Vec3d;
	std::vector<Halfedge> toCheck;
	Vertex v1 = e.h0().from();
	Vertex v2 = e.h0().to();
	Vector3 midpoint = toOM(m.eval_period_edge(e.h0(), 0.5));
	// find (halfedge) link around the edge, starting with those surrounding v1
	Halfedge he = v1.halfedge();
	Halfedge st = he;
	do {
		he = he.next();
		if (he.from() != v2 && he.next().from() != v2) {
			toCheck.push_back(he);
		}
		he = he.next().opp();
	} while (he != st);
	// v2
	he = v2.halfedge();
	st = he;
	do {
		he = he.next();
		if (he.from() != v1 && he.next().from() != v1) {
			toCheck.push_back(he);
		}
		he = he.next().opp();
	} while (he != st);

	for (Halfedge he0 : toCheck) {
		Halfedge heT = he0.opp();
		Vertex v1 = heT.from();
		Vertex v2 = heT.next().from();
		Vertex v3 = heT.next().next().from();
		Vector3 a = m.point(v1);
		Vector3 b = m.point(v2);
		Vector3 c = m.point(v3);
		a = midpoint + make_period(a - midpoint);
		b = midpoint + make_period(b - midpoint);
		c = midpoint + make_period(c - midpoint);
		if (checkFoldover(a, b, c, midpoint, 2)) {
			// std::cout<<"prevented foldover"<<std::endl;
			return false;
		}
	}
	return true;
}

extern double delaunay_lengh_upp;
extern double delaunay_lengh_low;
bool adjustEdgeLengths(PeriodSurfaceMesh& m/*, const PeriodSurfaceMesh& mtgt*/, double flatLength, double epsilon, double minLength, bool curvatureAdaptive = true)
{
	bool didSplitOrCollapse = false;
	using Edge = typename OpenMesh::SmartEdgeHandle;
	using Halfedge = typename OpenMesh::SmartHalfedgeHandle;
	using Vector3 = Eigen::Vector3d;
	// queues of edges to CHECK to change
	std::vector<Edge> toSplit;
	std::vector<Edge> toCollapse;

	for (auto e : m.edges())
	{
		toSplit.push_back(e);
	}

	// actually do it
	//std::cerr << "Spliting..." << std::endl;
	while (!toSplit.empty())
	{
		Edge e = toSplit.back();
		toSplit.pop_back();
		double length_e = m.periodEdgeLength(e);
		//double threshold = (curvatureAdaptive) ? findMeanTargetL(m, e, flatLength, epsilon) : flatLength;
		//double threshold = (curvatureAdaptive) ? fabs(findGaussianTargetL(m, e, flatLength, epsilon)) : flatLength;
		double threshold = (curvatureAdaptive) ? fabs(findTotalCurvatureTargetL(m, e, flatLength, epsilon)) : flatLength;
		if (length_e > minLength && length_e > threshold * delaunay_lengh_upp)
		{
			didSplitOrCollapse = true;
			Vector3 newPos = m.eval_period_edge(e.h0(), 0.5);
			//Vector3 newPosOrig = edgeMidpoint(mesh, geometryOriginal, e);
			//Halfedge he = mesh->splitEdgeTriangular(e);
			//Vertex newV = he.vertex();
			//geometry->inputVertexPositions[newV] = newPos;
			//geometryOriginal->inputVertexPositions[newV] = newPosOrig;
			m.split(e, toOM(newPos));
		}
		else
		{
			toCollapse.push_back(e);
		}

	}
	//std::cerr << "Collapsing..." << std::endl;
	while (!toCollapse.empty())
	{
		Edge e = toCollapse.back();
		toCollapse.pop_back();
		if (!e.deleted() && e.is_valid()) // make sure it exists
		{
			//double threshold = (curvatureAdaptive) ? findMeanTargetL(m, e, flatLength, epsilon) : flatLength;
			//double threshold = (curvatureAdaptive) ? fabs(findGaussianTargetL(m, e, flatLength, epsilon)) : flatLength;
			double threshold = (curvatureAdaptive) ? fabs(findTotalCurvatureTargetL(m, e, flatLength, epsilon)) : flatLength;
			if (m.periodEdgeLength(e) < threshold * delaunay_lengh_low)
			{
				Vector3 newPos = m.eval_period_edge(e.h0(), 0.5);
				//Vector3 newPosOrig = edgeMidpoint(mesh, geometryOriginal, e);
				if (shouldCollapse(m, e)) {
					//auto vh = mesh->collapseEdgeTriangular(e);
					//if (v != Vertex()) {
					//	if (!v.isBoundary()) {
					//		geometry->inputVertexPositions[v] = newPos;
					//		//geometryOriginal->inputVertexPositions[v] = newPosOrig;
					//	}
					//}
					if (m.is_collapse_ok(e.h0())) {
						auto v1 = e.h0().to();
						m.collapse(e.h0());
						m.set_point(v1, toOM(newPos));
						didSplitOrCollapse = true;
					} else {
						std::cout << "invalid collapse " << m.point(e.h0().from()) << ", " << m.point(e.h0().to()) << std::endl;
					}
				}
			}
		}
	}

	m.garbage_collection();
	//mesh->validateConnectivity();
	////std::cerr<<"hereeee"<<std::endl;
	//mesh->compress();
	//geometry->refreshQuantities();
	return didSplitOrCollapse;
}

inline bool isDelaunay(PeriodSurfaceMesh& m, OM::SmartEdgeHandle e)
{
	auto he = e.h0();
	float angle1 = m.period_sector_angle(he.next());
	float angle2 = m.period_sector_angle(he.opp().next());
	return angle1 + angle2 <= M_PI;
}

inline bool isCrease(PeriodSurfaceMesh& m, OM::SmartEdgeHandle e) {
	auto [d0, d1, lef, rgt] = m.find1diag(e.h0(), 1);
	return angle((d1 - d0).cross(lef - d0).normalized(), (rgt - d0).cross(d1 - d0).normalized()) > M_PI / 3;
}

bool check_nan(const PeriodSurfaceMesh& m) {
	for (auto vh : m.vertices()) {
		auto p = m.point(vh);
		for (int k = 0; k < 3; k++) {
			if (std::isnan(p[k])) {
				std::cout << "vh = " << vh << ", p = " << p << std::endl;
				throw std::runtime_error("has nan");
				return false; 
			}
		}
	}
	return true;
}

void fixDelaunay(PeriodSurfaceMesh& m)
{
	//check_nan(m);

	using Edge = OM::SmartEdgeHandle;
	using Halfedge = OM::SmartHalfedgeHandle;
	// queue of edges to check if Delaunay
	std::queue<Edge> toCheck;
	// true if edge is currently in toCheck
	std::vector<bool> inQueue(m.n_edges() * 2, true); // todo: this maybe dangerous!
	// start with all edges
	for (Edge e : m.edges())
	{
		toCheck.push(e);
		inQueue.at(e.idx()) = true;
	}
	// counter and limit for number of flips
	int flipMax = 100 * m.n_vertices();
	int flipCnt = 0;
	bool has_invalid_flip = false;
	while (!toCheck.empty() && flipCnt < flipMax)
	{
		Edge e = toCheck.front();
		toCheck.pop();
		inQueue.at(e.idx()) = false;
		if (!e.is_valid() || e.deleted()) continue;
		// if not Delaunay, flip edge and enqueue the surrounding "diamond" edges (if not already)
		if (!e.is_boundary() && !isDelaunay(m, e) /*&& !isCrease(m, e)*/)
		{
			flipCnt++;
			Halfedge he = e.h0();
			Halfedge he1 = he.next();
			Halfedge he2 = he1.next();
			Halfedge he3 = he.opp().next();
			Halfedge he4 = he3.next();

			if (!inQueue.at(he1.edge().idx()))
			{
				toCheck.push(he1.edge());
				inQueue.at(he1.edge().idx()) = true;
			}
			if (!inQueue.at(he2.edge().idx()))
			{
				toCheck.push(he2.edge());
				inQueue.at(he2.edge().idx()) = true;
			}
			if (!inQueue.at(he3.edge().idx()))
			{
				toCheck.push(he3.edge());
				inQueue.at(he3.edge().idx()) = true;
			}
			if (!inQueue.at(he4.edge().idx()))
			{
				toCheck.push(he4.edge());
				inQueue.at(he4.edge().idx()) = true;
			}
			if (m.is_flip_ok(e)) {
				m.flip(e);
				//check_nan(m);
			} else {
				has_invalid_flip = true;
				//std::cout << "invalid flip " << m.point(e.h0().from()) << ", " << m.point(e.h0().to()) << std::endl;
				//m.period_shift();
				//m.savePeriodicMesh(getPath("fliperr.obj"), std::set<OM::SmartVertexHandle>{}, 1);
				//throw std::runtime_error("invalid flip");
			}
		}
	}
	if (has_invalid_flip) { std::cout << "invalid flip occurs" << std::endl; }
}

template<typename Vector3>
Vector3 findBarycenter(Vector3 p1, Vector3 p2, Vector3 p3)
{
	return (p1 + p2 + p3) / 3;
}

OM::Vec3d findBarycenter(PeriodSurfaceMesh& m, OM::SmartFaceHandle f)
{
	using Vector3 = OM::Vec3d;
	using Vertex = OM::SmartVertexHandle;
	// retrieve the face's vertices
	int index = 0;
	Vector3 p[3];
	for (Vertex v0 : f.vertices_ccw())
	{
		p[index] = m.point(v0);
		if (index) p[index] = p[0] + make_period(p[index] - p[0], 2, 1);
		index++;
	}
	return findBarycenter(p[0], p[1], p[2]);
}

template<typename Vector3>
Vector3 findCircumcenter(Vector3 p1, Vector3 p2, Vector3 p3)
{
	// barycentric coordinates of circumcenter
	double a = (p3 - p2).norm();
	double b = (p3 - p1).norm();
	double c = (p2 - p1).norm();
	double a2 = a * a;
	double b2 = b * b;
	double c2 = c * c;
	Vector3 O{ a2 * (b2 + c2 - a2), b2 * (c2 + a2 - b2), c2 * (a2 + b2 - c2) };
	// normalize to sum of 1
	O /= O[0] + O[1] + O[2];
	// change back to space
	return O[0] * p1 + O[1] * p2 + O[2] * p3;
}

auto findCircumcenter(PeriodSurfaceMesh& m, OM::SmartFaceHandle f, OM::SmartVertexHandle vbase)
{
	using Vector3 = OM::Vec3d;
	using Vertex = OM::SmartVertexHandle;
	// retrieve the face's vertices
	int index = 0;
	auto pbase = m.point(vbase);
	Vector3 p[3];
	for (Vertex v0 : f.vertices_ccw())
	{
		p[index] = m.point(v0);
		p[index] = pbase + make_period(p[index] - pbase, 2, 1);
		index++;
	}
	return findCircumcenter(p[0], p[1], p[2]);
}

template<typename Vector3>
inline Vector3 projectToPlane(Vector3 v, Vector3 norm)
{
	return v - norm * dot(norm, v);
}

auto vertexNormal(PeriodSurfaceMesh& m, OM::SmartVertexHandle v) {
	Eigen::Vector3d norm(0, 0, 0);
	auto [o, ring] = m.find1ring(v, 1);
	int N = ring.size();
	for (int k = 0; k < ring.size(); k++) {
		Eigen::Vector3d oh1 = ring[k] - o;
		Eigen::Vector3d oh2 = ring[(k + 1) % N] - o;
		Eigen::Vector3d nf = (oh1).cross(oh2);
		double ang = angle(oh1, oh2);
		norm += ang * nf.normalized();
	}
	//std::cerr<<norm<<std::endl;
	return toOM(norm.normalized());
}

void smoothByCircumcenter(PeriodSurfaceMesh& m)
{
	//geometry->requireFaceAreas();
	using Vector3 = OM::Vec3d;
	using Vertex = OM::SmartVertexHandle;
	using Face = OM::SmartFaceHandle;
	// smoothed vertex positions
	std::vector<Vector3> newVertexPosition(m.n_vertices());
	for (Vertex v : m.vertices())
	{
		newVertexPosition[v.idx()] = m.point(v); // default
		if (!v.is_boundary())
		{
			Vector3 updateDirection(0, 0, 0);
			//double totalD = 0;
			for (Face f : v.faces_ccw())
			{
				// add the center weighted by face area to the update direction
				Vector3 center;

				if (f.is_boundary())
				{
					center = findBarycenter(m, f);
				}
				else
				{
					center = findCircumcenter(m, f, v);
					//if (toEigen(center).hasNaN()) { std::cout << m.getFacePeriodVertex(f, 1) << std::endl; }
				}

				//double D = 1/findFaceTargetL(mesh, geometry, f, 1, 0.1);
				//D = D*D;
				updateDirection += m.period_face_area(f) * make_period(center - m.point(v));
				//totalD += geometry->faceArea(f) * D;
			}
			//std::cerr<<updateDirection<<std::endl;
			updateDirection /= (3 * m.period_vertex_dual_area(v));
			//updateDirection /= totalD;
			//std::cerr<<"  "<<updateDirection<<std::endl;
			// project update direction to tangent plane
			updateDirection = projectToPlane(updateDirection, vertexNormal(m, v));
			Vector3 newPos = m.point(v) + .5 * updateDirection;
			//if (toEigen(newPos).hasNaN()) { std::cout << m.point(v) << ", " << updateDirection << std::endl; throw std::runtime_error("has nan"); }
			newVertexPosition[v.idx()] = newPos;
		}
	}
	// update final vertices
	for (Vertex v : m.vertices())
	{
		m.point(v) = newVertexPosition[v.idx()];
	}
}

double estimate_average_length(PeriodSurfaceMesh& m) {
	double As = 0;
	for (auto fh : m.faces()) {
		double Af = m.period_face_area(fh);
		As += Af;
	}
	As /= m.n_faces();
	return std::sqrt(4 * As / std::sqrt(3));
}

void delete_degree3_faces(PeriodSurfaceMesh& m) {
	for (auto vh : m.vertices()) {
		if (vh.valence() == 3) {
			int counter = 0;
			OM::SmartVertexHandle vvh[3];
			for (auto vv : vh.vertices_ccw()) vvh[counter++] = vv;
			m.delete_vertex(vh);
			m.add_face(vvh[0], vvh[1], vvh[2]);
		}
	}
}

void msf::PeriodSurfaceMesh::delaunayRemesh(int inner_iter, double tgtlen /*= -1*/, double min_len /*= -1*/, double adaptive /*= 1e-6*/, int outer_iter /*= 1*/)
{
	period_shift();

	if (inner_iter <= 0) return;

	if (adaptive <= 0) adaptive = 1e-6;

	double flatLength;
	if (tgtlen < 0) {
		flatLength = estimate_average_length(*this);
	} else {
		flatLength = tgtlen;
	}

	if (min_len < 0 || min_len > flatLength) {
		min_len = flatLength / 4;
	}

	std::cout << "Delaunay target length = " << flatLength << std::endl;

	for (int i = 0; i < outer_iter; i++) {
		//std::cout << "fixing edges" << std::endl;
		if (!adjustEdgeLengths(*this, flatLength, 1. / adaptive, min_len)) break;
		//savePeriodicMesh("split.obj", std::set<OM::SmartVertexHandle>{}, 1);
		//std::cout << "flipping" << std::endl;
		fixDelaunay(*this);
		//delete_degree3_faces(*this);
		//savePeriodicMesh("flip.obj", std::set<OM::SmartVertexHandle>{}, 1);
		for (int j = 0; j < inner_iter; j++) {
			//std::cout << "smoothing" << std::endl;
			smoothByCircumcenter(*this);
			//std::cout << "flipping" << std::endl;
			fixDelaunay(*this);
			//delete_degree3_faces(*this);
		}
		//std::cout << "done" << std::endl;
	}
	delete_degree3_faces(*this);
	garbage_collection();
}

void test_delaunay_remesh(void) {
	std::string mfile = "D:/projects/minisurf/image/siggraph/adc-opt/P0.5-31.stl";
	PeriodSurfaceMesh m;
	m.readMergePeriodBoundary(mfile);
	m.pertubMesh(25, 1, 1);
	m.savePeriodicMesh("temp.obj", std::set<OM::SmartVertexHandle>{}, 1);
	m.delaunayRemesh(10, 0.04);
	m.split_unit_cell();
	m.savePeriodicMesh("dem.obj", std::set<OM::SmartVertexHandle>{}, 1);
}

void msf::PeriodSurfaceMesh::period_shift(int pm)
{
	for (auto vh : vertices()) {
		auto p = point(vh);
		for (int k = 0; k < 3; k++) {
			if (pm >= 0) {
				// shift right
				if (p[k] < -1 - 1e-5) p[k] += 2;
			} else {
				if (p[k] > 1 + 1e-5) p[k] -= 2;
			}
		}
		point(vh) = p;
	}
}

void msf::PeriodSurfaceMesh::period_shift(void)
{
	period_shift(1);
	period_shift(-1);
}

bool inbox(const OM::Vec3d& p) {
	for (int k = 0; k < 3; k++) {
		if (p[k] < -1 - 1e-5 || p[k]>1 + 1e-5) return false;
	}
	return true;
}

template<typename T>
void split_unit_cell(PeriodSurfaceMesh& m, std::vector<T>& vlist) {
	m.period_shift();
	//m.save("temp.obj");
	std::map<OM::SmartVertexHandle, T> vh2data;
	if (!vlist.empty()) { for (auto vh : m.vertices()) { vh2data[vh] = vlist[vh.idx()]; } }
	for (int iter = 0; iter < 1; iter++) {
		for (int axis = 0; axis < 3; axis++) {
			//std::vector<OM::EdgeHandle> toSplit;
			for (auto eh : m.edges()) {
				if (m.calc_edge_length(eh) > 1.) {
					auto he = eh.h0();
					auto p = m.point(he.from());
					auto p1 = m.point(he.to());
					p1 = p + make_period(p1 - p);
					if (fabs(fabs(p[axis]) - 1) < +1e-5) continue;
					if (fabs(fabs(p1[axis]) - 1) < +1e-5) continue;
					if (!((fabs(p1[axis]) <= 1) ^ (fabs(p[axis]) <= 1))) continue;
					double cut_p = p1[axis] > 1 ? 1 : -1;
					//toSplit.emplace_back(he.edge());
					double t = (cut_p - p[axis]) / (p1[axis] - p[axis]);
					OM::Vec3d c = t * p1 + (1 - t) * p;
					//if (!inbox(c)) {
					//	std::cout << "p0 = " << p << ", p1 = " << p1 << std::endl;
					//	std::cout << "c = " << c << std::endl;
					//}
					auto vhnew = m.split(he.edge(), c);
					if (vlist.size()) vh2data[vhnew] = vh2data[he.to()] * t + (1 - t) * vh2data[he.from()];
				}
			}
			m.period_shift();
		}
		//m.savePeriodicMesh("temp.obj", std::set<OM::SmartVertexHandle>{}, 1);
	}

	if (vlist.size()) {
		std::vector<T> vlistnew(m.n_vertices());
		for (auto [vh, p] : vh2data) { vlistnew[vh.idx()] = p; }
		vlist = vlistnew;
	}
	//m.period_shift();
	//m.save("temp.obj");
}

void msf::PeriodSurfaceMesh::split_unit_cell(std::vector<Eigen::Vector3d>& vlist)
{
	::split_unit_cell(*this, vlist);
}


void msf::PeriodSurfaceMesh::split_unit_cell(std::vector<double>& vlist)
{
	::split_unit_cell(*this, vlist);
}

void msf::PeriodSurfaceMesh::split_unit_cell(void)
{
	auto tmp = std::vector<int>{};
	::split_unit_cell(*this, tmp);
}

void msf::PeriodSurfaceMesh::clamp_period_boundary(double band /*= 2e-5*/)
{
	for (auto vh : vertices()) {
		auto o = point(vh);
		Eigen::Vector3d clm(0, 0, 0);
		for (int k = 0; k < 3; k++) {
			if (o[k] > 1 - band) clm[k] = 1;
			if (o[k] < -1 + band) clm[k] = -1;
		}
		if (clm.isZero()) continue;
		auto n = vertexNormal(*this, vh);
		// number of zero
		int nz = (clm[0] == 0) + (clm[1] == 0) + (clm[2] == 0);
		// no zero, clamp to vertex
		if (nz == 0) {
			set_point(vh, toOM(clm));
		}
		// only one zero, clamp to edge
		else if (nz == 1) {
			for (int k = 0; k < 3; k++) {
				if (clm[k] != 0) continue;
				double t = -n.dot(toOM(clm) - o) / n.dot(toOM(Eigen::Vector3d::Unit(k)));
				clm[k] = t; 
				set_point(vh, toOM(clm));
				break;
			}
		}
		// has two zero, clam to face
		else {
			for (int k = 0; k < 3; k++) {
				if (clm[k] == 0) continue;
				Eigen::Vector3d n1 = toEigen(n);
				Eigen::Vector3d n2 = Eigen::Vector3d::Unit(k).cross(n1);
				Eigen::Vector3d e1 = Eigen::Vector3d::Unit((k + 1) % 3);
				Eigen::Vector3d e2 = Eigen::Vector3d::Unit((k + 2) % 3);
				Eigen::Vector2d b; b << n1.dot(clm - toEigen(o)), n2.dot(clm - toEigen(o));
				Eigen::Matrix2d A;
				A << n1.dot(e1), n1.dot(e2),
					n2.dot(e1), n2.dot(e2);
				Eigen::Vector2d t = A.lu().solve(-b);
				clm += t[0] * e1 + t[1] * e2;
				set_point(vh, toOM(clm));
				break;
			}
		}
	}
}


