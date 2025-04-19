#include "skeletonMesher.h"
#include "Flag.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include "dir_utils.h"
#include "geometry/HashTypes.h"
#include "matlab/matlab_utils.h"
#include "CGAL/Simple_cartesian.h"
#include "CGAL/AABB_traits.h"
#include "CGAL/AABB_tree.h"
#include "CGAL/AABB_triangle_primitive.h"
#include "CGAL/AABB_segment_primitive.h"
#include "igl/writeOBJ.h"
#include "igl/writeSTL.h"
#include "igl/readOBJ.h"
#include "igl/read_triangle_mesh.h"
#include "igl/combine.h"
//#include "igl/orient_outward.h"
#include "igl/bfs_orient.h"
#include "igl/boundary_facets.h"
#include "vdb/openvdb_wrapper_t.h"
#include "igl/copyleft/marching_cubes.h"
#include "igl/extract_non_manifold_edge_curves.h"
#include "matlab/matlab_utils.h"
#include "cgal/mesh_intersection.h"
#include "igl/per_face_normals.h"
#include "igl/remove_unreferenced.h"
#include "triangle++/tpp_interface.hpp"
#include "geometry/HashTypes.h"
#include "PeriodicMesher.h"
#include "cgal/cgal_utils.h"
#include "reorient_mesh.h"

std::string getPath(std::string);
extern int debug_level;

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef K::Triangle_3 Triangle;
typedef std::vector<Triangle>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<K, Iterator> TriPrimitive;
typedef CGAL::AABB_traits<K, TriPrimitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> FaceTree;
typedef K::Segment_3 Segment;
typedef std::vector<Segment>::iterator SegIterator;
typedef CGAL::AABB_segment_primitive<K, SegIterator> SegPrimitive;
typedef CGAL::AABB_tree<CGAL::AABB_traits<K, SegPrimitive>> LineTree;
//typedef std::pair<Segment, OpenMesh::SmartHalfedgeHandle> Segment_with_attribute; // Segment with an integer attribute
//typedef std::vector<Segment_with_attribute>::iterator Segment_Iterator;

template<typename P> auto toCGAL(const P& p) { return ::Point(p[0], p[1], p[2]); }

SkeletonMesher::SkeletonMesher(const std::string& filename, int reso, Real period)
{
	read_skeleton(filename);
	buildSDF(reso, period);
}

bool SkeletonMesher::read_skeleton(const std::string& filename)
{
	try {
		auto ext = dir_utils::path2extension(filename);
		if (ext != ".sk") {
			Eigen::MatrixX3<Real> V;
			Eigen::MatrixX3i F;
			igl::read_triangle_mesh(filename, V, F);
			if (F.rows() == 0) throw std::runtime_error("invalid mesh");
			for (int i = 0; i < V.rows(); i++) {
				sk_vertex.emplace_back(V.row(i).transpose());
			}
			for (int i = 0; i < F.rows(); i++) {
				sk_tri.emplace_back(F.row(i).transpose());
			}
		} else {
			std::ifstream ifs(filename);
			char buf[1000];
			while (ifs.getline(buf, 1000)) {
				if (buf[0] == '#') continue;
				if (buf[0] == 'V') {
					Eigen::Vector3<Real> p;
					std::istringstream istr(buf + 1);
					istr >> p[0] >> p[1] >> p[2];
					sk_vertex.emplace_back(p);
				} else if (buf[0] == 'L') {
					Eigen::Vector2i lin;
					std::istringstream istr(buf + 1);
					istr >> lin[0] >> lin[1];
					sk_line.emplace_back(lin);
				} else if (buf[0] == 'T') {
					Eigen::Vector3i tri;
					std::istringstream istr(buf + 1);
					istr >> tri[0] >> tri[1] >> tri[2];
					sk_tri.emplace_back(tri);
				}
			}
		}
	}
	catch (...) {
		std::cout << "\033[31m" << "Skeleton file format error" << "\033[0m\n";
		return false;
	}
	return true;
}

void SkeletonMesher::buildSDF(int reso, Real period /*= 0*/)
{
	if (period != 0) {
		Eigen::Vector3<Real> dmin,dmax;
		dmin.setConstant(-period / 2 * 0.8);
		dmax.setConstant(period / 2 * 0.8);
		BBox bb_in(dmin, dmax), bb_out(dmin * 1.5, dmax * 1.5);
		// append periodic line elements
		int nl = sk_line.size();
		for (int i = 0; i < nl; i++) {
			auto p0 = sk_vertex[sk_line[i][0]];
			auto p1 = sk_vertex[sk_line[i][1]];
			if (bb_in.isOut(p0) || bb_in.isOut(p1)) {
				for (int j = 0; j < 27; j++) {
					if (j == 13) continue;
					Eigen::Vector3<Real> joff{ period * (j % 3 - 1), period * (j / 3 % 3 - 1), period * (j / 9 - 1) };
					sk_vertex.push_back(joff + p0);
					sk_vertex.push_back(joff + p1);
					sk_line.emplace_back(sk_vertex.size() - 2, sk_vertex.size() - 1);
				}
			}
		}
		int nt = sk_tri.size();
		// append periodic triangle elements
		for (int i = 0; i < nt; i++) {
			Eigen::Vector3<Real> p[3] = { sk_vertex[sk_tri[i][0]], sk_vertex[sk_tri[i][1]], sk_vertex[sk_tri[i][2]] };
			if (bb_in.isOut(p[0]) || bb_in.isOut(p[1]) || bb_in.isOut(p[2])) {
				for (int j = 0; j < 27; j++) {
					if (j == 13) continue;
					Eigen::Vector3<Real> joff{ period * (j % 3 - 1), period * (j / 3 % 3 - 1), period * (j / 9 - 1) };
					Eigen::Vector3<Real> pnew[3] = { joff + p[0], joff + p[1], joff + p[2] };
					if (bb_out.isOut(pnew[0]) && bb_out.isOut(pnew[1]) && bb_out.isOut(pnew[2])) {
						continue;
					}
					sk_vertex.push_back(pnew[0]); sk_vertex.push_back(pnew[1]); sk_vertex.push_back(pnew[2]);
					sk_tri.emplace_back(sk_vertex.size() - 3, sk_vertex.size() - 2, sk_vertex.size() - 1);
				}
			}
		}
	}

	// build AABB tree
	std::vector<Triangle> trilist;
	std::vector<Segment> linelist;
	for (int i = 0; i < sk_line.size(); i++) {
		auto vid = sk_line[i];
		Segment line(toCGAL(sk_vertex[vid[0]]), toCGAL(sk_vertex[vid[1]]));
		linelist.emplace_back(line);
	}
	for (int i = 0; i < sk_tri.size(); i++) {
		auto vid = sk_tri[i];
		Triangle tri(toCGAL(sk_vertex[vid[0]]), toCGAL(sk_vertex[vid[1]]), toCGAL(sk_vertex[vid[2]]));
		trilist.emplace_back(tri);
	}
	if (debug_level > 1) {
		Eigen::MatrixX3<Real> vlist(sk_vertex.size(), 3);
		Eigen::MatrixX3i flist(trilist.size(), 3);
		for (int i = 0; i < sk_vertex.size(); i++) {
			vlist.row(i) = sk_vertex[i].transpose();
		}
		for (int i = 0; i < sk_tri.size(); i++) {
			flist.row(i) = sk_tri[i].transpose();
		}
		igl::writeOBJ(getPath("trilist.obj"), vlist, flist);
	}
	FaceTree face_aabb(trilist.begin(), trilist.end());
	LineTree line_aabb(linelist.begin(), linelist.end());
	Eigen::Vector3<Real> dmin(-1, -1, -1), dmax(1, 1, 1);
	if (period != 0) {
		dmin.setConstant(-period / 2); dmax.setConstant(period / 2);
	}
	Eigen::Vector3<Real> d = dmax - dmin;
	int np = std::pow(reso + 1, 3);
	Eigen::VectorX<Real> gridDist(np);
	grid_pos.resize(3, np);
#pragma omp parallel for
	for (int zi = 0; zi < reso + 1; zi++) {
		Eigen::Vector3<Real> pi(0, 0, Real(zi) / reso);
		for (int yi = 0; yi < reso + 1; yi++) {
			pi[1] = Real(yi) / reso;
			for (int xi = 0; xi < reso + 1; xi++) {
				pi[0] = Real(xi) / reso;
				Eigen::Vector3<Real> p = pi.cwiseProduct(d) + dmin;
				double dist1 = 1e30, dist2 = 1e30;
				if (!linelist.empty()) {
					auto pclose = line_aabb.closest_point(toCGAL(p));
					dist1 = std::sqrt(CGAL::squared_distance(pclose, toCGAL(p)));
				}
				if (!trilist.empty()) {
					auto pclose = face_aabb.closest_point(toCGAL(p));
					dist2 = std::sqrt(CGAL::squared_distance(pclose, toCGAL(p)));
				}
				double dist = (std::min)(dist1, dist2);
				int pid = xi + yi * (reso + 1) + zi * std::pow(reso + 1, 2);
				gridDist[pid] = dist;
				grid_pos.col(pid) << p;
			}
		}
	}
	grid_pos = grid_pos.transpose().eval();
	grid_dist = gridDist;

	//eigen2ConnectedMatlab("gridpos", grid_pos);

	if (debug_level > 1) {
		int gs[3] = { reso + 1,reso + 1,reso + 1 };
		std::vector<float> gridvalues;
		std::transform(grid_dist.begin(), grid_dist.end(), std::back_inserter(gridvalues), [](const double& f) {return float(f); });
		openvdb_wrapper_t<float>::lexicalGrid2openVDBfile(getPath("vdbtest.vdb"), gs, gridvalues);
	}
}

extern std::string meshfile, meshfile_x;
extern std::vector<double> thick_list;
extern Real sample_reso;
extern Real remesh_tgtlen;
extern std::vector<double> remesh_tgtlen_list;
extern Real remesh_noise;
extern int rand_seed;
extern int remesh_iter;
extern int remesh_smooth_iter;
extern int remesh_smooth_type;
extern int remesh_pertub_iter;
extern int remesh_smooth_batch_size;
extern int log_detail;


std::tuple<Eigen::MatrixX3<Real>, Eigen::MatrixX3i > append_period_face(const Eigen::MatrixX3<Real>& vlist, const Eigen::MatrixX3i& flist);


std::tuple<Eigen::MatrixX3<msf::Real>, Eigen::MatrixX3i> removeDupVertices(const Eigen::MatrixX3<msf::Real>& v, const Eigen::MatrixX3i& f, msf::Real eps);

void clamp_vertices(Eigen::MatrixX3<Real>& V, Eigen::MatrixX3i& F);

std::tuple<Eigen::MatrixX3<Real>, Eigen::MatrixX3i> SkeletonMesher::extract_level_set(Real s, Real edge_len, bool closed)
{
	if (grid_dist.size() == 0) {
		throw std::runtime_error("the SDF values on Grid have not been evaluated");
	}
	int gs = std::pow(grid_dist.size() + 0.1, 0.33333333333);
	Eigen::MatrixX3<Real> V;
	Eigen::MatrixX3i F;

	// 1. extract level set surface
	std::cout << "Meshing isovalue level set at " << s << "...\n";
	igl::copyleft::marching_cubes(grid_dist, grid_pos,
		gs, gs, gs, s, V, F);
	std::cout << "finished\n";
	MeshSurface mtgt; mtgt.read(V, F);
	std::tie(V, F) = cgal_remesh(V, F);
	std::tie(V, F) = removeDupVertices(V, F, 2e-5);
	clamp_vertices(V, F);
	if (debug_level > 1) { igl::writeOBJ(getPath("mc.obj"), V, F); }
	
	// 2. remesh the mesh
	std::cout << "Remeshing isosurface...\n";
	PeriodSurfaceMesh tmesher;
	tmesher.read(V, F, true, false);
	auto vcut = tmesher.mergePeriodBoundary();
	vcut = tmesher.periodic_remesh(remesh_iter, vcut, edge_len, remesh_smooth_iter, 0, &mtgt);
	std::tie(V, F) = tmesher.getVFlist();
	Eigen::VectorXi new2old;
	std::tie(V, F, new2old) = tmesher.dupPeriodFaces(V, F);
	std::cout << "Finished\n";
	if (debug_level > 1) { igl::writeSTL(getPath("befapp.stl"), V, F); }
	
	// 3. attach period face
	if (closed) {
		std::cout << "Attaching period face...\n";
		std::tie(V, F) = append_period_face(V, F);
		// deal with inserted Steiner poinnts 
		// boundary_loop returns a list of vertices along the boundary loop(s) of the mesh.
		// boundary_facets returns a #E x 2 matrix of all the boundary edge indices.
		std::cout << "Finished\n";
	}


	return { V,F };
}

void clamp_vertices(Eigen::MatrixX3<Real>& V) {
	std::transform(V.data(), V.data()+V.size(), V.data(), [&](Real vi) {
		if (vi < -0.99999) return -1.;
		if (vi > 0.99999) return 1.;
		return vi;
		});
}

void clamp_vertices(Eigen::MatrixX3<Real>& V, Eigen::MatrixX3i& F) {
	Eigen::MatrixX2i e_nake;
	igl::boundary_facets(F, e_nake);
	std::sort(e_nake.data(), e_nake.data() + e_nake.size());
	int* last = std::unique(e_nake.data(), e_nake.data() + e_nake.size());

	std::vector<int> vbound(e_nake.data(), last);
	std::cout << "Clamping " << vbound.size() << " boundary vertices" << std::endl;
	std::transform(vbound.begin(), vbound.end(), vbound.begin(), [&](int vi) {
		for (int k = 0; k < 3; k++) {
			if (V(vi, k) < -1 + 2e-5) {
				V(vi, k) = -1; 
			} else if (V(vi, k) > 1 - 2e-5) {
				V(vi, k) = 1; 
			}
		}
		return vi;
		});
}

void SkeletonMesher::exportLelveSet(const std::vector<Real>& slist, const std::string& fileprefix)
{
	for (int i = 0; i < slist.size(); i++) {
		double s = slist[i];
		double elen = remesh_tgtlen;
		if (i < remesh_tgtlen_list.size()) { elen = remesh_tgtlen_list[i]; }
		auto [V, F] = extract_level_set(s, elen, true);
		clamp_vertices(V);
#if 0
		std::string filepath = fileprefix + "_" + std::to_string(s) + ".obj";
		igl::writeOBJ(filepath, V, F);
#else
		std::string filepath = fileprefix + "_" + std::to_string(s) + ".stl";
		igl::writeSTL(filepath, V, F, igl::FileEncoding::Binary);
#endif
	}
}

std::string getPath(std::string s);

void test_skeleton_mesher(void) {
	SkeletonMesher mesher(meshfile, sample_reso, 2);
	
	auto f = dir_utils::path2filename(meshfile);
	mesher.exportLelveSet(thick_list, getPath(f));
}

void test_close_mesh(void) {
	// 2. remesh the mesh
	std::cout << "Remeshing surface...\n";
	PeriodSurfaceMesh tmesher;
	tmesher.read(meshfile, true, false);
	auto [V, F] = tmesher.getVFlist();
	MeshSurface mtgt; mtgt.read(V, F);
	std::tie(V, F) = cgal_remesh(V, F);
	if (debug_level > 0) { igl::writeOBJ(getPath("cgalrm.obj"), V, F); }
	tmesher = PeriodSurfaceMesh();
	tmesher.read(V, F, true, false);
	auto vcut = tmesher.mergePeriodBoundary();
	vcut = tmesher.periodic_remesh(remesh_iter, vcut, remesh_tgtlen, remesh_smooth_iter, 0, &mtgt);
	std::tie(V, F) = tmesher.getVFlist();
	Eigen::VectorXi new2old;
	std::tie(V, F, new2old) = tmesher.dupPeriodFaces(V, F);
	std::cout << "Finished\n";

	// 3. attach period face
	std::cout << "Attaching period face...\n";
	std::tie(V, F) = append_period_face(V, F);
	std::cout << "Finished\n";
	clamp_vertices(V);
	auto f = dir_utils::path2filename(meshfile);
	std::string filepath = getPath(f) + "_close.stl";
	igl::writeSTL(filepath, V, F, igl::FileEncoding::Binary);
}


using namespace msf;

using namespace tpp;

void tricall_main(void) {
#if 0
	Eigen::MatrixX2<Real> V;
	Eigen::MatrixX2i E;
	V.resize(8, 2);
	V << 0, 0, 2, 0, 2, 2, 0, 2,
		0.5, 0.5, 1.5, 0.5, 1.5, 1.5, 0.5, 1.5;
	E.resize(8, 2);
	E << 0, 1, 1, 2, 2, 3, 3, 0,
		4, 7, 7, 6, 6, 5, 5, 4;
	Eigen::MatrixX2<Real> WV;
	Eigen::MatrixX3i WF;
	Eigen::MatrixX2i WE;
	Eigen::VectorXi J(V.rows());
/// @param[in] V  #V by 2 list of texture mesh vertices
/// @param[in] E  #E by 2 list of constraint edge indices into V
/// @param[in] flags  string of triangle flags should contain "-c" unless the
///     some subset of segments are known to enclose all other
///     points/segments.
/// @param[out] WV  #WV by 2 list of background mesh vertices 
/// @param[out] WF  #WF by 2 list of background mesh triangle indices into WV
/// @param[out] WE  #WE by 2 list of constraint edge indices into WV (might be smaller
///     than E because degenerate constraints have been removed)
/// @param[out] J  #V list of indices into WF/WE for each vertex in V
	igl::triangle::cdt(V, E, "-c", WV, WF, WE, J);
	Eigen::MatrixX3<Real> vlist(WV.rows(), 3); vlist.setZero();
	vlist.leftCols(2) = WV;
	igl::writeOBJ(getPath("cdt.obj"), WV, WF);
#else
    // prepare input
    std::vector<Delaunay::Point> delaunayInput;

    delaunayInput.push_back(Delaunay::Point(0, 0));
    delaunayInput.push_back(Delaunay::Point(2, 0));
    delaunayInput.push_back(Delaunay::Point(2, 2));
    delaunayInput.push_back(Delaunay::Point(0, 2));
    delaunayInput.push_back(Delaunay::Point(0.5, 0.5));
    delaunayInput.push_back(Delaunay::Point(1.5, 0.5));
    delaunayInput.push_back(Delaunay::Point(1.5, 1.5));
    delaunayInput.push_back(Delaunay::Point(0.5, 1.5));

    //  standard triangulation
    Delaunay trGenerator(delaunayInput);
    //  constrained Delaunay
    std::vector<Delaunay::Point> segments;
    std::vector<Delaunay::Point> holes;
	for (int i = 0; i < 4; i++) {
		segments.push_back(delaunayInput[i]);
		segments.push_back(delaunayInput[(i + 1) % 4]);
		//holes.push_back(delaunayInput[i + 4]);
		//holes.push_back(delaunayInput[(i + 1) % 4 + 4]);
	}
	for (int i = 0; i < 4; i++) {
		segments.push_back(delaunayInput[i + 4]);
		segments.push_back(delaunayInput[(i + 1) % 4 + 4]);
	}

	holes.emplace_back(Delaunay::Point(1.5, 1.5));
	trGenerator.setMaxArea(0.04f);
	trGenerator.useConvexHullWithSegments(true);
	trGenerator.setSegmentConstraint(segments);
	trGenerator.setHolesConstraint(holes);
	trGenerator.Triangulate(true);
	//trGenerator.TriangulateConf(true, Info);

	auto triCount = trGenerator.triangleCount();
	auto vCount = trGenerator.verticeCount();

	Eigen::MatrixX3d vlist(triCount * 3, 3);
	Eigen::MatrixX3i flist(triCount, 3);
    // iterate over triangles
	int counter = 0;
	for (FaceIterator fit = trGenerator.fbegin(); fit != trGenerator.fend(); ++fit) {
		// vertex id list
		Delaunay::Point tri_p[3];
        int vertexIdx1 = fit.Org(tri_p);
		int vertexIdx2 = fit.Dest(tri_p + 1);
		int vertexIdx3 = fit.Apex(tri_p + 2);
        // access data
		vlist.row(counter * 3) = Eigen::Vector3<Real>(tri_p[0][0], tri_p[0][1], 0).transpose();
		vlist.row(counter * 3 + 1) = Eigen::Vector3<Real>(tri_p[1][0], tri_p[1][1], 0).transpose();
		vlist.row(counter * 3 + 2) = Eigen::Vector3<Real>(tri_p[2][0], tri_p[2][1], 0).transpose();
		flist.row(counter) = Eigen::Vector3i(counter * 3, counter * 3 + 1, counter * 3 + 2).transpose();
		counter++;
    }
	
	igl::writeOBJ(getPath("cdt.obj"), vlist, flist);
#endif
}

//int box_edge_id(const Eigen::Vector3<Real>& p) {
//	VertexFlag vflag;
//	vflag.set_period_boundary(0.9999, p[0], p[1], p[2]);
//	if (!vflag.is_period_boundary()) return -1;
//	if (vflag.is_min_period(2)) {
//		if (vflag.is_min_period(0)) {
//
//		}
//	}
//}
Eigen::Vector3<Real> edge_vector(int period_id) {
	EdgeFlag eflag;
	eflag.set_period_id_type(period_id);
	Eigen::Vector3<Real> n(0, 0, 0);
	for (int i = 0; i < 3; i++) {
		if (!eflag.is_period_boundary(i)) {
			n[i] = 1;
		}
	}
	Eigen::Vector3<Real> ev = n;
	return ev;
}

inline double angle(const Eigen::Vector2<Real>& v1, const Eigen::Vector2<Real>& v2) {
	double ang = std::acos(v1.normalized().dot(v2.normalized()));
	Eigen::Matrix2<Real> v12;
	v12 << v1, v2;
	if (v12.determinant() < 0) {
		return -ang;
	} else {
		return ang;
	}
}

int winding_number(const Eigen::Vector2<Real>& p, const std::vector<Eigen::Vector2<Real>>& loop) {
	double ang_sum = 0;
	for (int i = 0; i < loop.size() - 1; i++) {
		double ang = angle(loop[i] - p, loop[i + 1] - p);
		//std::cout << ang << std::endl;
		ang_sum += ang;
	}
	double w = ang_sum / (2 * M_PI);
	return w > 0 ? (w + 0.5) : (w - 0.5);
}

std::vector<Eigen::Vector2<Real>> checkout_holes(const std::vector<std::vector<Eigen::Vector2<Real>>>& loops) {
	if (loops.size() <= 1) return {};

	std::vector<Eigen::Vector2<Real>> sample_points;
	{
		//std::ofstream ofs(getPath("hole_loops"));
		//for (int i = 0; i < loops.size(); i++) {
		//	for (int j = 0; j < loops[i].size(); j++) {
		//		ofs << loops[i][j][0] << " " << loops[i][j][1] << " ";
		//	}
		//	ofs << std::endl;
		//}
	}
	for (int i = 0; i < loops.size(); i++) {
		Eigen::Vector2<Real> v12 = loops[i][0] - loops[i][1];
		Eigen::Vector2<Real> n(-v12[1], v12[0]);
		Eigen::Vector2<Real> c = (loops[i][0] + loops[i][1]) / 2;
		Eigen::Vector2<Real> p[2] = { c - 0.02 * n, c + 0.02 * n };
		//std::cout << "p = " << p[0][0] << " " << p[0][1] << " ; " << p[1][0] << " " << p[1][1] << std::endl;
		int n_wind[2] = { 0 };
		for (int k = 0; k < loops.size(); k++) {
			n_wind[0] += winding_number(p[0], loops[k]);
			n_wind[1] += winding_number(p[1], loops[k]);
		}
		//std::cout << "nwind = " << n_wind[0] << ", " << n_wind[1] << std::endl;
		if (n_wind[0] == 0) { sample_points.push_back(p[0]); }
		if (n_wind[1] == 0) { sample_points.push_back(p[1]); }
	}
	return sample_points;
}

std::tuple<Eigen::MatrixX3<msf::Real>, Eigen::MatrixX3i> removeDupVertices(const Eigen::MatrixX3<msf::Real>& v, const Eigen::MatrixX3i& f, msf::Real eps /*= 1e-5*/);

std::tuple<Eigen::MatrixX2<Real>, Eigen::MatrixX3i> mesh_closed_curves(const std::vector<std::vector<Eigen::Vector2<Real>>>& loops, double max_area) {
	
	std::unordered_map<Eigen::Vector2<Real>, int> pid;
	for (int i = 0; i < loops.size(); i++) {
		for (int j = 0; j < loops[i].size(); j++) {
			if (pid.count(loops[i][j])) continue;
			pid[loops[i][j]] = pid.size();
		}
	}
	// Delaunay triangulation points
	std::vector<Delaunay::Point> delaunayInput(pid.size());
	// list points
	for (auto [p, id] : pid) {
		delaunayInput[id] = (Delaunay::Point(p[0], p[1]));
	}
    //  standard triangulation
    Delaunay trGenerator(delaunayInput);
    //  constrained Delaunay
    std::vector<Delaunay::Point> segments;
    std::vector<Delaunay::Point> holes;
	for (int i = 0; i < loops.size() ; i++) {
		for (int j = 0; j < loops[i].size() - 1; j++) {
			segments.push_back(Delaunay::Point(loops[i][j][0], loops[i][j][1]));
			segments.push_back(Delaunay::Point(loops[i][j + 1][0], loops[i][j + 1][1]));
		}
	}
	auto hole_points = checkout_holes(loops);
	for (int i = 0; i < hole_points.size(); i++) {
		holes.push_back(Delaunay::Point(hole_points[i][0], hole_points[i][1]));
	}
	//trGenerator.useConvexHullWithSegments(true);
	trGenerator.setSegmentConstraint(segments);
	trGenerator.setHolesConstraint(holes);
	trGenerator.setMaxArea(max_area);
	trGenerator.setMinAngle(30);
	trGenerator.Triangulate(true);

	int triCount = trGenerator.triangleCount();
	Eigen::MatrixX3d vlist(triCount * 3, 3);
	Eigen::MatrixX3i flist(triCount, 3);
	// iterate over triangles
	int counter = 0;
	for (FaceIterator fit = trGenerator.fbegin(); fit != trGenerator.fend(); ++fit) {
		// vertex id list
		Delaunay::Point tri_p[3];
		int vertexIdx1 = fit.Org(tri_p);
		int vertexIdx2 = fit.Dest(tri_p + 1);
		int vertexIdx3 = fit.Apex(tri_p + 2);
		// access data
		vlist.row(counter * 3) = Eigen::Vector3<Real>(tri_p[0][0], tri_p[0][1], 0).transpose();
		vlist.row(counter * 3 + 1) = Eigen::Vector3<Real>(tri_p[1][0], tri_p[1][1], 0).transpose();
		vlist.row(counter * 3 + 2) = Eigen::Vector3<Real>(tri_p[2][0], tri_p[2][1], 0).transpose();
		flist.row(counter) = Eigen::Vector3i(counter * 3, counter * 3 + 1, counter * 3 + 2).transpose();
		counter++;
	}
	
	std::tie(vlist, flist) = removeDupVertices(vlist, flist, 1e-5);

	//igl::writeOBJ(getPath("close.obj"), vlist, flist);

	Eigen::MatrixX2<Real> vlist_plane = vlist.leftCols(2);
	return { vlist_plane, flist };
}

bool is_eq(const Eigen::Vector3<Real>& v1, const Eigen::Vector3<Real>& v2, Real eps = 1e-5) {
	return (v1 - v2).squaredNorm() < eps * eps;
}

std::vector<std::vector<Eigen::Vector3<Real>>> extract_period_edges(MeshSurface& mesh) {
	double lavr = mesh.average_edge_length();
	//std::cout << "- mesh average length " << lavr << std::endl;
	// find cut point on boundary edge
	std::map<int, std::vector<std::pair<Eigen::Vector3<Real>, Eigen::Vector3<Real>>>> e2p;
	//std::vector<Eigen::Vector3<Real>> ev;
	//std::vector<Eigen::Vector3<Real>> nv;
	for (auto vh : mesh.vertices()) {
		if (!vh.is_boundary()) continue;
		auto p = toEigen(mesh.point(vh));
		VertexFlag vflag;
		vflag.set_period_boundary(0.9999, p[0], p[1], p[2]);
		int n_b = 0;
		for (int i = 0; i < 3; i++) { if (vflag.is_period_boundary(i)) n_b++; }
		if (n_b >= 2) {
			//ev.push_back(p);
			auto n = toEigen(mesh.normal(vh));
			//nv.push_back(n);
			int period_id = vflag.period_type_id();
			e2p[period_id].emplace_back(p, n);
		}
	}
	// find edge
	std::vector<std::pair<Eigen::Vector3<Real>, Eigen::Vector3<Real>>> elist;
	for (auto it : e2p) {
		auto ve = edge_vector(it.first);
		std::vector<std::pair<Real, int>> ticks;
		for (auto vn : it.second) {
			auto [v, n] = vn;
			ticks.emplace_back(v.dot(ve), n.dot(ve) > 0 ? 1 : -1);
		}
		auto pbase = it.second[0].first;
		auto tbase = ticks[0].first;
		std::sort(ticks.begin(), ticks.end(), [](auto& left, auto& right) {
			return left.first < right.first;
			});
		for (int i = 0; i < ticks.size() - 1; i++) {
			if (ticks[i].second < 0 && ticks[i + 1].second > 0) {
				Real dtick[2] = { ticks[i].first - tbase, ticks[i + 1].first - tbase };
				elist.emplace_back(dtick[0] * ve + pbase, dtick[1] * ve + pbase);
			}
		}
		if (ticks[0].second > 0) {
			Real dtick[2] = { -1 - tbase, ticks[0].first - tbase };
			elist.emplace_back(dtick[0] * ve + pbase, dtick[1] * ve + pbase);
		}
		if (ticks.back().second < 0) {
			Real dtick[2] = { ticks.back().first - tbase, 1 - tbase };
			elist.emplace_back(dtick[0] * ve + pbase, dtick[1] * ve + pbase);
		}
	}
	{
		//std::ofstream ofs(getPath("elist"));
		//for (int i = 0; i < elist.size(); i++) {
		//	ofs << elist[i].first.transpose() << " " << elist[i].second.transpose() << std::endl;
		//}
	}
	// refined edge
	int ne = elist.size();
	std::vector<std::vector<Eigen::Vector3<Real>>> refined_edges;
	for (int i = 0; i < ne; i++) {
		refined_edges.emplace_back();
		Eigen::Vector3<Real> dl = elist[i].second - elist[i].first;
		Real L = dl.norm();
		int ns = (std::max)((L / lavr) + 0.5, 1.);
		dl /= ns;
		for (int j = 0; j < ns + 1; j++) {
			Real t = Real(j) / ns;
			refined_edges.back().emplace_back(elist[i].first + j * dl);
		}
	}
	{
		//std::ofstream ofs(getPath("refined_edges"));
		//for (int i = 0; i < refined_edges.size(); i++) {
		//	for (int j = 0; j < refined_edges[i].size(); j++) {
		//		ofs << refined_edges[i][j].transpose() << " ";
		//	}
		//	ofs << std::endl;
		//}
	}
	// find edge loops on boundary faces
	std::set<OM::SmartHalfedgeHandle> loop_edges;
	std::vector<std::vector<OM::SmartHalfedgeHandle>> loops;
	for (auto he : mesh.halfedges()) {
		if (!he.is_boundary()) continue;
		if (loop_edges.count(he)) continue;
		auto p0 = mesh.point(he.from());
		auto p1 = mesh.point(he.to());
		VertexFlag v0; v0.set_period_boundary(0.9999, p0[0], p0[1], p0[2]);
		VertexFlag v1; v1.set_period_boundary(0.9999, p1[0], p1[1], p1[2]);
		if (!(v0.getPeriodFlag() & v1.getPeriodFlag())) continue;
		// find start edge
		auto he0 = he;
		do {
			p0 = mesh.point(he.from());
			v0 = VertexFlag();
			v0.set_period_boundary(0.9999, p0[0], p0[1], p0[2]);
			if (v0.is_period_edge()) break;
			he = he.prev();
		} while (he != he0);
		he0 = he;
		loops.emplace_back();
		do {
			loop_edges.insert(he);
			loops.back().push_back(he);
			p0 = mesh.point(he.to());
			//std::cout << mesh.point(he.from()) << " " << mesh.point(he.to()) << std::endl;
			v0 = VertexFlag();
			v0.set_period_boundary(0.9999, p0[0], p0[1], p0[2]);
			if (v0.is_period_edge()) break;
			he = he.next();
		} while (he != he0);
		//if (he == he0) { loops.back().erase(loops.back().end() - 1); }
	}
	auto check_boundary = [](VertexFlag vflag, const Eigen::Vector3<Real>& e_beg, const Eigen::Vector3<Real>& e_end) {
		VertexFlag flg[2];
		flg[0].set_period_boundary(0.9999, e_beg[0], e_beg[1], e_beg[2]);
		flg[1].set_period_boundary(0.9999, e_end[0], e_end[1], e_end[2]);
		return vflag.getPeriodFlag() & flg[0].getPeriodFlag() & flg[1].getPeriodFlag();
		};

	auto edge_flag = [](const Eigen::Vector3<Real>& e_beg, const Eigen::Vector3<Real>& e_end) {
		VertexFlag flg[2];
		flg[0].set_period_boundary(0.9999, e_beg[0], e_beg[1], e_beg[2]);
		flg[1].set_period_boundary(0.9999, e_end[0], e_end[1], e_end[2]);
		VertexFlag eflg;
		eflg.set_period_boundary(flg[0].getPeriodFlag() & flg[1].getPeriodFlag());
		return eflg;
		};
	auto find_next_edge = [&](VertexFlag vflag, const Eigen::Vector3<Real>& e_beg, const Eigen::Vector3<Real>& e_end) {
		std::pair<int, int> loop_id{ 0,0 };
		loop_id.first = 1;
		const double eps = 1e-6 * 1e-6;
		for (int i = 0; i < refined_edges.size(); i++) {
			auto e_p0 = *refined_edges[i].begin();
			auto e_p1 = *refined_edges[i].rbegin();
			//std::cout << "e_p0 = " << e_p0.transpose() << ", e_p1 = " << e_p1.transpose() << std::endl;
			if (!check_boundary(vflag, e_p0, e_p1)) continue;
			if (is_eq(*(refined_edges[i].begin() + 1), e_beg) && is_eq(*refined_edges[i].begin(), e_end)) continue;
			if (is_eq(*(refined_edges[i].rbegin() + 1), e_beg) && is_eq(*refined_edges[i].rbegin(), e_end)) continue;
			if (is_eq(e_p0, e_end)) {
				loop_id.second = i;
				return loop_id;
			}
			else if (is_eq(e_p1, e_end)) {
				loop_id.second = -(i + 1);
				return loop_id;
			}
		}
		loop_id.first = 2;
		for (int i = 0; i < loops.size(); i++) {
			auto e_p0 = toEigen(mesh.point(loops[i].begin()->from()));
			auto e_p1 = toEigen(mesh.point((loops[i].begin() + 1)->to()));
			// this may cause error
			if (check_boundary(vflag, e_p0, e_p1)) {
				if (is_eq(e_p0, e_end)) {
					loop_id.second = i;
					return loop_id;
				}
			}
		}
		std::cout << "\033[31mCannot found oriented edge loops\033[0m\n";
		throw std::runtime_error("failed found next  segments");
		};
	std::vector<std::vector<Eigen::Vector3<Real>>> loop_vlist;
	for (int i = 0; i < loops.size(); i++) {
		loop_vlist.emplace_back();
		auto p = toEigen(mesh.point(loops[i][0].from()));
		loop_vlist.back().push_back(p);
		for (int j = 0; j < loops[i].size(); j++) {
			loop_vlist[i].push_back(toEigen(mesh.point(loops[i][j].to())));
		}
	}
	{
		//std::ofstream ofs(getPath("loop_vlist"));
		//for (int i = 0; i < loop_vlist.size(); i++) {
		//	for (int j = 0; j < loop_vlist[i].size(); j++) {
		//		ofs << loop_vlist[i][j].transpose() << " ";
		//	}
		//	ofs << std::endl;
		//}
	}
	std::vector<bool> loop_merged(loops.size(), false);
	for (int i = 0; i < loops.size(); i++) {
		if (loops[i].begin()->from() == loops[i].rbegin()->to()) continue;
		if (loop_merged[i]) continue;
		// cut by boundary edge
		Eigen::Vector3<Real> vbeg = toEigen(mesh.point(loops[i].begin()->from()));
		Eigen::Vector3<Real> vend = toEigen(mesh.point(loops[i].begin()->to()));
		VertexFlag eflg = edge_flag(vbeg, vend);
		//std::cout << "flg = " << eflg._flag << std::endl;
		while (!is_eq(vbeg, vend)) {
			auto [strip_type, strip_id] = find_next_edge(eflg, *(loop_vlist[i].rbegin() + 1), *loop_vlist[i].rbegin());
			//std::cout << "strip_type = " << strip_type << ", " << "strip_id = " << strip_id << std::endl;
			if (strip_type == 1) {
				if (strip_id < 0) {
					loop_vlist[i].insert(loop_vlist[i].end(), refined_edges[-strip_id - 1].rbegin() + 1, refined_edges[-strip_id - 1].rend());
				}
				else {
					loop_vlist[i].insert(loop_vlist[i].end(), refined_edges[strip_id].begin() + 1, refined_edges[strip_id].end());
				}
			}
			else if (strip_type == 2) {
				if (strip_id < 0) {
					loop_vlist[i].insert(loop_vlist[i].end(), loop_vlist[-strip_id].rbegin() + 1, loop_vlist[-strip_id].rend());
					loop_merged[strip_id] = true;
				}
				else {
					loop_vlist[i].insert(loop_vlist[i].end(), loop_vlist[strip_id].begin() + 1, loop_vlist[strip_id].end());
					loop_merged[strip_id] = true;
				}
			}
			vbeg = *loop_vlist[i].begin();
			vend = *loop_vlist[i].rbegin();
		}
		loop_vlist[i][0] = loop_vlist[i].back();
		//std::cout << std::endl;
	}
	std::vector<std::vector<Eigen::Vector3<Real>>> edgeloops;
	for (int i = 0; i < loop_merged.size(); i++) {
		if (!loop_merged[i]) { edgeloops.push_back(loop_vlist[i]); }
	}
	return edgeloops;
}

std::tuple<std::vector<int>,
	std::vector<std::vector<Eigen::Vector2<Real>>> >
	flatten_loops(std::vector<std::vector<Eigen::Vector3<Real>>>& loops) {
	std::vector<int> atPlane(loops.size());
	std::vector<std::vector<Eigen::Vector2<Real>>> flattened;
	for (int i = 0; i < loops.size(); i++) {
		int counter = 0;
		int axis = 0;
		do {
			Eigen::Vector3<Real> e1 = loops[i][counter] - loops[i][0];
			Eigen::Vector3<Real> e2 = loops[i][(counter + 1) % loops[i].size()] - loops[i][0];
			Eigen::Vector3<Real> n = e1.cross(e2);
			n = n.cwiseAbs();
			if (n.norm() > 1e-6) {
				Eigen::Vector3<Real> nn = n.normalized();
				for (int k = 0; k < 3; k++) {
					if (nn[k] > 0.5) { axis = k + 1; break; }
				}
				if (axis == 0) {
					throw std::runtime_error("cannnot determine axis");
				}
				atPlane[i] = loops[i][0][axis - 1] < 0 ? -axis : axis;
				break;
			}
		} while (counter++ < loops[i].size());
		flattened.emplace_back();
		for (int j = 0; j < loops[i].size(); j++) {
			Eigen::Vector2<Real> pl;
			pl[0] = loops[i][j][(axis - 1 + 1) % 3];
			pl[1] = loops[i][j][(axis - 1 + 2) % 3];
			flattened.back().emplace_back(pl);
		}
	}
	return { atPlane, flattened };
}

std::tuple<std::vector<Eigen::MatrixX3<Real>>, std::vector<Eigen::MatrixX3i>> mesh_period_face(const Eigen::MatrixX3<Real>& vlist, const Eigen::MatrixX3i& flist) {
	MeshSurface mesh;
	mesh.read(vlist, flist);
	auto loops = extract_period_edges(mesh);
	{
		//std::ofstream ofs(getPath("loops"));
		//for (int i = 0; i < loops.size(); i++) {
		//	for (int j = 0; j < loops[i].size(); j++) {
		//		ofs << loops[i][j][0] << " " << loops[i][j][1] << " " << loops[i][j][2] << " ";
		//	}
		//	ofs << std::endl;
		//}
	}
	auto [axis, plsg] = flatten_loops(loops);
	{
		//std::ofstream ofs(getPath("plsg"));
		//for (int i = 0; i < plsg.size(); i++) {
		//	for (int j = 0; j < plsg[i].size(); j++) {
		//		ofs << plsg[i][j][0] << " " << plsg[i][j][1] << " ";
		//	}
		//	ofs << std::endl;
		//}
		//std::cout << "axis = " << axis[0] << axis[1] << axis[2] << axis[3] << axis[4] << axis[5] << std::endl;
	}
	std::vector<std::vector<std::vector<Eigen::Vector2<Real>>>> plane_loops(6);
	for (int i = 0; i < axis.size(); i++) {
		int axis_id = (axis[i] < 0 ? 0 : 3) + (std::abs(axis[i]) - 1);
		//std::cout << "plsg " << i << " = " << plsg[i].size() << std::endl;
		plane_loops[axis_id].emplace_back(plsg[i]);
		//std::cout << "pl " << axis_id << " = " << plane_loops[axis_id].size() << std::endl;
	}
	std::vector<std::tuple<Eigen::MatrixX2<Real>, Eigen::MatrixX3i>> period_meshes;
	for (int ax = 0; ax < 3; ax++) {
		auto [vlist, flist] = mesh_closed_curves(plane_loops[ax], std::pow(remesh_tgtlen, 2) / 2);
		period_meshes.emplace_back(vlist, flist);
	}
	std::tuple<std::vector<Eigen::MatrixX3<Real>>, std::vector<Eigen::MatrixX3i>> period_meshes_3;
	auto& [Vlist, Flist] = period_meshes_3;
	for (int i = 0; i < 6; i++) {
		int j = (i % 3 + 1) % 3;
		int k = (i % 3 + 2) % 3;
		auto& [v2, f2] = period_meshes[i % 3];
		Eigen::MatrixX3<Real> V(v2.rows(), 3);
		Eigen::MatrixX3i F = f2;
		V.col(i % 3).setConstant(i / 3 ? 1 : -1);
		V.col(j) = v2.col(0); V.col(k) = v2.col(1);
		Vlist.emplace_back(V);
		Flist.emplace_back(F);
	}
	return period_meshes_3;
}

template <typename T, typename Compare>
std::vector<std::size_t> sort_permutation(
	const std::vector<T>& vec, const Compare& compare)
{
	std::vector<std::size_t> p(vec.size());
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(),
		[&](std::size_t i, std::size_t j) { return compare(vec[i], vec[j]); });
	return p;
}
template <typename T>
std::vector<T> apply_permutation(
	const std::vector<T>& vec,
	const std::vector<std::size_t>& p)
{
	std::vector<T> sorted_vec(vec.size());
	std::transform(p.begin(), p.end(), sorted_vec.begin(),
		[&](std::size_t i) { return vec[i]; });
	return sorted_vec;
}


std::vector<Eigen::Vector3<Real>> generate_period_group(const Eigen::Vector3<Real>& v) {
	int on_boarder[3] = { 0 };
	for (int i = 0; i < 3; i++) {
		on_boarder[i] = v[i] <  -0.9999 || v[i]>0.9999;
	}
	std::vector<Eigen::Vector3<Real>> plist;
	for (Real xi : (on_boarder[0] ? std::vector<Real>{-1, 1} : std::vector<Real>{ v[0] })) {
		Eigen::Vector3<Real> p = v;
		p[0] = xi;
		for (Real yi : (on_boarder[1] ? std::vector<Real>{ -1, 1 } : std::vector<Real>{ v[1] })) {
			p[1] = yi;
			for (Real zi : (on_boarder[2] ? std::vector<Real>{ -1, 1 } : std::vector<Real>{ v[2] })) {
				p[2] = zi;
				//std::cout << "p = " << p.transpose() << std::endl;
				plist.push_back(p);
			}
		}
	}
	return plist;
}

std::tuple<Eigen::MatrixX3<Real>, Eigen::MatrixX3i > append_period_face(const Eigen::MatrixX3<Real>& vlist, const Eigen::MatrixX3i& flist) {
	auto [VV, FF] = mesh_period_face(vlist, flist);
	VV.push_back(vlist);
	FF.push_back(flist);
	Eigen::VectorXi vsizes, fsizes;
	Eigen::MatrixX3<Real> V;
	Eigen::MatrixX3i F;
	igl::combine(VV, FF, V, F);
	std::tie(V, F) = removeDupVertices(V, F, 1e-5);
	//igl::writeOBJ(getPath("beforemsplit.obj"), V, F);
	// split naked face
#if 1
	Eigen::MatrixX2i e_nake;
	igl::boundary_facets(F, e_nake);
	std::unordered_set<std::pair<Eigen::Vector3<Real>, int>> v_nake;
	for (int i = 0; i < e_nake.rows(); i++) {
		v_nake.insert(std::pair<Eigen::Vector3<Real>, int>(V.row(e_nake(i, 0)).transpose(), e_nake(i, 0)));
		v_nake.insert(std::pair<Eigen::Vector3<Real>, int>(V.row(e_nake(i, 1)).transpose(), e_nake(i, 1)));
	}
	// [[eid, vid],...]
	std::map<int, std::vector<int>> split_pairs;
	for (auto [p, vid] : v_nake) {
		for (int i = 0; i < e_nake.rows(); i++) {
			Eigen::Vector3<Real> a = V.row(e_nake(i, 0)).transpose();
			Eigen::Vector3<Real> b = V.row(e_nake(i, 1)).transpose();
			double t = (a - b).dot(p - b) / (a - b).squaredNorm();
			if (0.01 < t && t < 0.99 && is_eq(b + t * (a - b), p)) {
				split_pairs[i].emplace_back(vid);
			}
		}
	}
	std::map<int, std::vector<Eigen::Vector3i>> newFaces;
	std::vector<int> spliting_vertices;
	for (int fi = 0; fi < F.rows(); fi++) {
		Eigen::Vector3i fv = F.row(fi).transpose();
		int vk = -1;
		std::vector<int> vsplit;
		for (auto [eid, pmid] : split_pairs) {
			int ev[2] = { e_nake(eid,0), e_nake(eid,1) };
#if 1
			std::set<int> fvset{ fv[0],fv[1],fv[2] };
			fvset.erase(ev[0]); fvset.erase(ev[1]);
			if (fvset.size() > 1) continue;
			int vlast = *fvset.begin();
#else
			int vlast = fv[0] ^ fv[1] ^ fv[2] ^ ev[0] ^ ev[1];
			if (vlast != fv[0] + fv[1] + fv[2] - ev[0] - ev[1]) continue;
#endif
			for (int k = 0; k < 3; k++) {
				if (vlast == fv[k]) {
					vsplit = pmid;
					vk = k; break;
				}
			}
			if (vk > 0) break;
		}
		if (vk < 0) continue;
		Eigen::Vector3<Real> a = V.row(fv[(vk + 1) % 3]).transpose();
		Eigen::Vector3<Real> b = V.row(fv[(vk + 2) % 3]).transpose();
		std::vector<Real> tsplit;
		for (int k = 0; k < vsplit.size(); k++) {
			Real t = (V.row(vsplit[k]).transpose() - a).dot(b - a);
			tsplit.push_back(t);
		}
		auto p = sort_permutation(tsplit, std::less<Real>());
		apply_permutation(vsplit, p);
		spliting_vertices.insert(spliting_vertices.end(), vsplit.begin(), vsplit.end());
		vsplit.insert(vsplit.begin(), fv[(vk + 1) % 3]);
		vsplit.insert(vsplit.end(), fv[(vk + 2) % 3]);
		for (int i = 1; i < vsplit.size() - 1; i++) {
			Real s = i / (Real)(vsplit.size() - 1);
			V.row(vsplit[i]) = V.row(vsplit[0]) * (1 - s) + V.row(vsplit.back()) * s;
		}
		for (int i = 0; i < vsplit.size() - 1; i++) {
			newFaces[fi].emplace_back(fv[vk], vsplit[i], vsplit[i + 1]);
		}
	}
	int s_fold = F.rows();
	int n_fnew = s_fold;
	for (auto [fold, fnew] : newFaces) { n_fnew += fnew.size() - 1; }
	F.conservativeResize(n_fnew, 3);
	int counter = 0;
	for (auto [fold, fnew] : newFaces) {
		F.row(fold) = fnew[0].transpose();
		for (int i = 1; i < fnew.size(); i++) {
			F.row(counter + s_fold) = fnew[i].transpose();
			counter++;
		}
	}

	// split periodic vertices
	{
		//igl::writeOBJ(getPath("beforesplit.obj"), V, F);
		std::vector<Eigen::Vector3<Real>> new_split;
		PeriodicGridIndex indexer(Eigen::Vector3<Real>(-1.01, -1.01, -1.01), Eigen::Vector3<Real>(2.02, 2.02, 2.02));
		for (int i = 0; i < spliting_vertices.size(); i++) {
			Eigen::Vector3<Real> p = V.row(spliting_vertices[i]).transpose();
			//std::cout << "p = " << p.transpose() << std::endl;
			for (auto pp : generate_period_group(p)) {
				if (is_eq(pp, p)) continue;
				//new_split.push_back(pp);
				indexer.insert(pp);
			}
		}
		Eigen::Matrix3X<Real> pset = indexer.dumpPoints().transpose();
		for (int i = 0; i < pset.cols(); i++) {
			new_split.push_back(pset.col(i));
		}
		{
			//std::ofstream ofs(getPath("splitnew"), std::ios::binary);
			//for (int i = 0; i < new_split.size(); i++) { ofs.write((const char*)new_split[i].data(), sizeof(new_split[i])); }
			//ofs.close();
			//ofs.open(getPath("splitold"), std::ios::binary);
			//for (int i = 0; i < spliting_vertices.size(); i++) {
			//	Eigen::Vector3<Real> p = V.row(spliting_vertices[i]).transpose();
			//	ofs.write((const char*)p.data(), sizeof(p));
			//}
		}

		std::vector<int> period_edge_faces;
		std::map<int, std::vector<Eigen::Vector3i>> split_faces;
		for (int fid = 0; fid < F.rows(); fid++) {
			Eigen::Vector3i fv = F.row(fid).transpose();
			Eigen::Matrix3<Real> fvp;
			fvp << V.row(fv[0]).transpose(), V.row(fv[1]).transpose(), V.row(fv[2]).transpose();
			bool pbd[3] = { false,false,false };
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					pbd[i] = pbd[i] || (fvp(i, j) < -0.9999) || (fvp(i, j) > 0.9999);
				}
			}
			Eigen::Vector3<Real> n = (fvp.col(2) - fvp.col(0)).cross(fvp.col(1) - fvp.col(0));
			Eigen::Vector3<Real> t = -fvp.rowwise().sum() / 3 + n / std::sqrt(n.norm());
			int n_pbd = pbd[0] + pbd[1] + pbd[2];
			if (n_pbd < 2) { continue; }
			period_edge_faces.push_back(fid);
			Eigen::Matrix3<Real> A = fvp;
			A.colwise() += t;
			Eigen::FullPivHouseholderQR<Eigen::Matrix3<Real>> slv(A);
			std::vector<std::pair<Real, int>> e_split[3];
			for (int pid = 0; pid < new_split.size(); pid++) {
				auto pnew = new_split[pid];
				Eigen::Vector3<Real> c = slv.solve(pnew + t);
				if ((std::abs)(c.sum() - 1) > 1e-3) { continue; }
				bool is_out = false;
				bool is_in = false;
				for (int k = 0; k < 3; k++) {
					if (c[k] < -1e-3 || c[k]>1 + 1e-3) { is_out = true; }
					is_in = is_in || (c[k] > 1e-3 && c[k] < 1 - 1e-3);
				}
				if (is_out || !is_in) continue;
				//std::cout << "pnew = " << pnew.transpose() << ", " << "c = " << c.transpose() << std::endl;
				// these points must be in the interior of edges
				Eigen::Vector3<Real> ci = (c * 1e2).cast<int>().cast<Real>() / 1e2;
				for (int k = 0; k < 3; k++) {
					if (ci[k] == 0) { e_split[k].emplace_back(ci[(k + 2) % 3], pid); }
				}
			}
			// e_split[k][i] = (t, pid in new_split)
			if (e_split[0].empty() && e_split[1].empty() && e_split[2].empty()) { continue; }
			for (int k = 0; k < 3; k++) {
				if (e_split[k].empty()) continue;
				int new_fv[3] = { fv[0],fv[1],fv[2] };
				for (int kk = 0; kk < k; kk++) {
					if (e_split[kk].empty()) continue;
					new_fv[(kk + 1) % 3] = -e_split[kk].back().second - 1;
				}
				std::sort(e_split[k].begin(), e_split[k].end(), [](auto& p1, auto& p2) { return p1.first < p2.first; });
				split_faces[fid].emplace_back(new_fv[k], new_fv[(k + 1) % 3], -e_split[k][0].second - 1);
				for (int kk = 0; kk < e_split[k].size() - 1; kk++) {
					split_faces[fid].emplace_back(new_fv[k], -e_split[k][kk].second - 1, -e_split[k][kk + 1].second - 1);
				}
				split_faces[fid].emplace_back(new_fv[k], -e_split[k].back().second - 1, new_fv[(k + 2) % 3]);
			}
		}
		{
			//std::ofstream ofs(getPath("splitfaces"));
			//for (auto& [fid, fnew] : split_faces) {
			//	ofs << fid << " ";
			//	for (int j = 0; j < fnew.size(); j++) { ofs << fnew[j][0] << " " << fnew[j][1] << " " << fnew[j][2] << " "; }
			//	ofs << std::endl;
			//}
			//ofs.close();
			//ofs.open(getPath("edgefaces"));
			//for (int i = 0; i < period_edge_faces.size(); i++) {
			//	ofs << period_edge_faces[i] << " ";
			//}
			//ofs << std::endl;
		}
		// count inserted vertices
		std::map<int, int> inserted_vertices;
		n_fnew = 0;
		int nf_old = F.rows();
		int nrow_old = V.rows();
		for (auto& [fid, fnew] : split_faces) {
			n_fnew += fnew.size() - 1;
			for (auto& fv : fnew) {
				for (int k = 0; k < 3; k++) {
					if (fv[k] >= 0) continue;
					if (inserted_vertices.count(fv[k])) {
						fv[k] = inserted_vertices[fv[k]];
					} else {
						int vold = fv[k];
						fv[k] = inserted_vertices.size() + nrow_old;
						inserted_vertices[vold] = fv[k];
					}
				}
			}
		}

		V.conservativeResize(V.rows() + inserted_vertices.size(), 3);
		counter = 0;
		for (auto [split_id, vnew_id] : inserted_vertices) {
			V.row(vnew_id) = new_split[-split_id - 1].transpose();
		}
		F.conservativeResize(F.rows() + n_fnew, 3);
		counter = 0;
		for (const auto& [fid, fnew] : split_faces) {
			F.row(fid) = fnew[0].transpose();
			for (int k = 1; k < fnew.size(); k++) { F.row(nf_old + counter++) = fnew[k].transpose(); }
		}
	}
#else
#endif
	reorient_mesh(V, F);

	return { V, F };
}

void test_mesh_boundary(const std::string& meshfile) {
	Eigen::MatrixX3<Real> vlist;
	Eigen::MatrixX3i  flist;
	igl::readOBJ(meshfile, vlist, flist);
	auto [V, F] = append_period_face(vlist, flist);
	igl::writeOBJ(getPath("periodmesh.obj"), V, F);
}
