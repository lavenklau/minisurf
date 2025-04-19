#include "PeriodicMesher.h"
#include <unordered_set>

#include "CGAL/Simple_cartesian.h"
#include "CGAL/AABB_traits.h"
#include "CGAL/AABB_tree.h"
#include "CGAL/AABB_triangle_primitive.h"
#include "CGAL/AABB_segment_primitive.h"
#include "igl/writeOBJ.h"
#include "matlab/matlab_utils.h"
#include "fmt/core.h"
#include "igl/vertex_triangle_adjacency.h"
#include "igl/triangle_triangle_adjacency.h"

std::string getPath(std::string);

extern int debug_level;

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef K::Segment_3 Segment;
typedef std::pair<Segment, OpenMesh::SmartHalfedgeHandle> Segment_with_attribute; // Segment with an integer attribute
typedef std::vector<Segment_with_attribute>::iterator Segment_Iterator;

struct HandleSegment {
public:
	typedef Segment_Iterator Id;
	// the CGAL types returned
	typedef K::Point_3    Point;
	typedef Segment Datum;
private:
	Id m_it; // this is what the AABB tree will store internally
public:
	HandleSegment() {} // default constructor needed
	// the following constructor is the one that receives the iterators from the
	// iterator range given as input to the AABB_tree
	HandleSegment(Segment_Iterator a)
		: m_it(a) {}
	Id id() const { return m_it; }
	// on the fly conversion from the internal data
	// to the CGAL types
	Datum datum() const {
		return m_it->first; // assembles a triangle from three points
	}
	// returns one point which must be on the primitive
	Point reference_point() const {
		return m_it->first.point(0);
	}
};
typedef CGAL::AABB_traits<K,  HandleSegment> MyAABB_traits;
typedef CGAL::AABB_tree<MyAABB_traits> Tree;


std::string getPath(std::string);

template<typename P> auto toCGAL(const P& p) { return Point(p[0], p[1], p[2]); }

template<typename P1,typename P2>
msf::Real squaredDistance(const P1& v1, const P2& v2) {
	return 
		(v1[0] - v2[0]) * (v1[0] - v2[0]) +
		(v1[1] - v2[1]) * (v1[1] - v2[1]) +
		(v1[2] - v2[2]) * (v1[2] - v2[2]);
}

// default is 2^12 ~= 1e-6
//template<int PrecisionLevel = 22>
inline int isEnd(const Point p, const OpenMesh::Vec3d& v0, const OpenMesh::Vec3d& v1, msf::Real eps = 1e-6) {
	//constexpr msf::Real eps = msf::Real(1) / (1 << PrecisionLevel);
	int endid = -1;
	msf::Real d[2] = { squaredDistance(p,v0), squaredDistance(p,v1) };
	if (d[0] < eps*eps) return 0;
	if (d[1] < eps * eps) return 1;
	return -1;
}


std::unique_ptr<Tree> generate_aabb_tree(std::string name, msf::PeriodSurfaceMesh& mesh, std::vector<OpenMesh::SmartHalfedgeHandle>& helist) {
	static thread_local std::map<std::string, std::vector<Segment_with_attribute>> named_segments;
	named_segments[name].clear();
	auto& segments = named_segments[name];
	for (int i = 0; i < helist.size(); i++) {
		Point v[2] = {
			toCGAL(mesh.point(helist[i].from())),
			toCGAL(mesh.point(helist[i].to())),
		};
		Segment seg(v[0], v[1]);
		Segment_with_attribute heseg(seg, helist[i]);
		segments.push_back(heseg);
	}
	// Populate segments with your data
	// segments.push_back(std::make_pair(Segment(Point(x1, y1, z1), Point(x2, y2, z2)), attribute));
	auto tree = std::make_unique<Tree>(segments.begin(), segments.end());
	return std::move(tree);
}

void msf::PeriodSurfaceMesh::mergePeriodEdges(void)
{
#if 0
	std::vector<std::vector<OM::SmartHalfedgeHandle>> edgeChains;
	// (origin, normal)
	std::vector<std::pair<Eigen::Vector3<Real>, Eigen::Vector3<Real>>> chainPlane;
	std::vector<bool> passed(n_halfedges(), false);
	// checkout edge chain and normal
	for (auto he : halfedges()) {
		if (passed[he.idx()] == true) continue;
		passed[he.idx()] = true;
		// checkout edge chain
		auto hestart = he;
		edgeChains.emplace_back();
		chainPlane.emplace_back();
		bool normalTo[3] = { true, true, true };
		bool onRight[3] = { false,false,false };
		OM::Vec3d e[3] = { {1.,0.,0.},{0.,1.,0.},{0.,0.,1.} };
		Eigen::Vector3<Real> o(0, 0, 0);
		if (he.is_boundary()) {
			do {
				auto he_vec = calc_edge_vector(he).normalized();
				for (int i = 0; i < 3; i++) {
					if (point(he.from())[i] > 0.5) { onRight[i] = true; }
				}
				for (int i = 0; i < 3; i++) {
					if (!normalTo[i]) continue;
					if (std::abs(he_vec.dot(e[i])) > 1e-3) { normalTo[i] = false; }
				}
				edgeChains.rbegin()->push_back(he);
				he = hestart.next();
				passed[he.idx()] = true;
			} while (he != hestart);
			for (int i = 0; i < 3; i++) {
				if (normalTo[i]) {
					if (onRight[i]) { o = Eigen::Vector3<Real>::Unit(i); }
					chainPlane.emplace_back(std::pair<Eigen::Vector3<Real>, Eigen::Vector3<Real>>{o, toEigen(e[i]) });
				}
			}
		}
	}
	// classification 
#else
	
	// tag all boundary vertices
	std::map<OM::SmartVertexHandle, VertexFlag> vtag;
	Real border = 1 - 1e-6;
	for (auto vh : vertices()) {
		if (!vh.is_boundary()) continue;
		VertexFlag tag;
		auto p = point(vh);
		tag.set_period_boundary(border, p[0], p[1], p[2]);
		if (tag.is_period_boundary()) vtag[vh] = tag;
	}

	if (debug_level > 0) OM::IO::write_mesh(*this, getPath("beforecollapse.obj"));
	// collapse short boundary edges
#if 1
	for (auto eh : edges()) {
		if (eh.deleted() || !eh.is_valid()) continue;
		//if (!eh.is_boundary()) continue;
		auto he = eh.halfedge(0);;
		//if (he.from().is_boundary() || he.to().is_boundary());
		auto hv = calc_edge_vector(he);
		if (hv.norm() < 2e-5) {
			auto t0 = vtag[he.from()].getPeriodFlag();
			auto t1 = vtag[he.to()].getPeriodFlag();
			auto hflag = t0 & t1;
			if (t0 || t1) {
				if ((~hflag) & t0) {
					collapse(he.opp());
				} else {
					collapse(he);
				}
			}
		}
	}
	// remove 
	for (auto iter = vtag.begin(); iter != vtag.end();) {
		if (!iter->second.is_period_boundary()) {
			iter = vtag.erase(iter);
		} else if (iter->first.deleted() || !iter->first.is_valid()) {
			iter = vtag.erase(iter);
		} else {
			iter++;
		}
	}
#endif

	std::map<OM::SmartHalfedgeHandle, EdgeFlag> hflags;
	std::map<OM::SmartHalfedgeHandle, std::vector<OM::SmartHalfedgeHandle>> splitted_he;
	for (auto he : halfedges()) {
		if (he.is_boundary() && !he.deleted()) {
			EdgeFlag tag;
			if (vtag.count(he.from()) && vtag.count(he.to())) {
				tag.set_period_boundary(vtag[he.from()].getPeriodFlag() & vtag[he.to()].getPeriodFlag());
				hflags[he] = tag;
			}
		}
	}

	std::unique_ptr<Tree> aabb;
	{
		std::vector<OM::SmartHalfedgeHandle> hlist;
		for (auto h : hflags) {
			if (h.second.is_min_period()) {
				hlist.push_back(h.first);
				splitted_he[h.first].push_back(h.first);
			}
		}
		aabb = generate_aabb_tree("boundary edge", *this, hlist);
	}

	//auto heboundary_init = getBoundaryHalfedges();
	//for (auto he : hflags) { delete_face(he.first.opp().face()); }
	//auto heboundary_cut = getBoundaryHalfedges();
	//std::vector<OM::SmartHalfedgeHandle> heboundary_bold(heboundary_cut.size());
	//auto lastbold = std::set_difference(heboundary_init.begin(), heboundary_init.end(), heboundary_cut.begin(), heboundary_cut.end(), heboundary_bold.begin());
	//heboundary_bold.erase(lastbold, heboundary_bold.end());

	//std::map<OM::SmartHalfedgeHandle, int> group;
	//int cur_group = 1;
	//for (auto he2t : hflags) {
	//	auto he = he2t.first;
	//	if (group[he] != 0) continue;
	//	if (he2t.second.is_period_boundary()) {
	//		auto initflag = he2t.second.getPeriodFlag();
	//		do {
	//			group[he] = cur_group;
	//			he = he.next();
	//			if (!hflags.count(he)) { break; }
	//		} while (group[he] == 0);
	//		cur_group++;
	//	}
	//}

	if (debug_level > 0) OM::IO::write_mesh(*this, getPath("beforesplit.obj"));

	// project max boundary vertex to min boundary
	std::map<OM::SmartHalfedgeHandle, std::vector<OM::SmartVertexHandle>> splited;
	std::map<OM::SmartHalfedgeHandle, std::pair<OM::SmartHalfedgeHandle, OM::SmartHalfedgeHandle>> heprepos;
	bool debug = false;
	for (auto v : vtag) {
		if (v.first.deleted()) continue;
		auto tag = v.second;
		if (!tag.is_max_period()) { continue; }
		auto p = point(v.first);
		Eigen::Vector3<Real> trans(0, 0, 0);
		for (int k = 0; k < 3; k++) {
			if (tag.is_max_period(k)) { p[k] -= 2; trans[k] = 2; }
		}
		auto p_he = aabb->closest_point_and_primitive(toCGAL(p));
		auto p_close = p_he.first;
		auto he_close = p_he.second->second;
		// not that sometime there is no min point for periodic group, thus find wrong vertex, this is a workaround
		if (CGAL::squared_distance(p_close, toCGAL(p)) > 0.01 * 0.01) { std::cout << "Warning : distant projection occurred\n"; continue; }
		auto& hesplit = splitted_he[he_close];
		for (int k = 0; k < hesplit.size(); k++) {
			auto he = hesplit[k];
			//if (debug) { std::cout << point(he.from()) << ", " << point(he.to()) << std::endl; }
			auto isend = isEnd(p_close, point(he.from()), point(he.to()), 1e-5);
			if (isend == 0) {
				set_point(v.first, point(he.from()) + toOM(trans)); break;
			} else if (isend == 1) {
				set_point(v.first, point(he.to()) + toOM(trans)); break;
			} else {
				// split min boundary edge
				auto newvh = add_vertex(toOM(p_close));
				if (cosangle(toEigen(p_close), he) > 1 - 1e-4) { continue; }
				split(he.edge(), newvh);
				// set to new position
				auto newp = toOM(trans) + toOM(p_close);
				set_point(v.first, newp);
				if (point(he.to()) == toOM(p_close)) {
					hesplit.push_back(he.next()); break;
				} else {
					hesplit.push_back(he.prev()); break;
				}
			}
		}
	}
	//{
	//	std::ofstream ofs(getPath("heclose"));
	//	for (auto hes : splited) { ofs << point(hes.first.from()) << "  " << point(hes.first.to()) << std::endl; }
	//}

	// update vertex flag
	vtag.clear();
	for (auto vh : vertices()) {
		if (vh.deleted() || !vh.is_valid()) continue;
		if (!vh.is_boundary()) continue;
		VertexFlag tag;
		auto p = point(vh);
		tag.set_period_boundary(border, p[0], p[1], p[2]);
		if (tag.is_period_boundary()) vtag[vh] = tag;
	}

	if (debug_level > 0) OM::IO::write_mesh(*this, getPath("splitmin.obj"));

	// split max boundary edge to match min boundary edge
	splited.clear();
	splitted_he.clear();
	heprepos.clear();
	std::unique_ptr<Tree> aabb_max;
	{
		std::vector<OM::SmartHalfedgeHandle> hmax;
		for (auto h : hflags) {
			if (h.second.is_max_period()) {
				hmax.push_back(h.first);
				splitted_he[h.first].push_back(h.first);
			}
		}
		aabb_max = generate_aabb_tree("max boundary edge", *this, hmax);	
	}
	for (auto v : vtag) {
		if (v.first.deleted()) continue;
		auto tag = v.second;
		if (!tag.is_min_period()) { continue; }
		auto p = point(v.first);
		Eigen::Vector3<Real> trans(0, 0, 0);
		for (int k = 0; k < 3; k++) {
			if (tag.is_min_period(k)) { p[k] += 2; trans[k] = -2; }
		}
		auto p_he = aabb_max->closest_point_and_primitive(toCGAL(p));
		auto p_close = p_he.first;
		//if ((toOM(p_close) - OM::Vec3d{ 1, -0.249497, 0.40506 }).norm() < 1e-5) {
		//	std::cout << point(v.first) << " " << p << std::endl;
		//	auto he = p_he.second->second;
		//	std::cout << point(he.from()) << point(he.to()) << std::endl;
		//}
		// not that sometime there is no min point for periodic group, thus find wrong vertex, this is a workaround
		if (CGAL::squared_distance(p_close, toCGAL(p)) > 0.01 * 0.01) {
			std::cout << "Warning : distant projection occurred\n";
			continue; 
		}
		auto he_close = p_he.second->second;
		auto& hesplit = splitted_he[he_close];
		for (int k = 0; k < hesplit.size(); k++) {
			auto he = hesplit[k];
			auto isend = isEnd(p_close, point(he.from()), point(he.to()), 1e-5);
			if (isend == 0) {
				set_point(v.first, point(he.from()) + toOM(trans)); break;
			} else if (isend == 1) {
				set_point(v.first, point(he.to()) + toOM(trans)); break;
			} else {
				// split min boundary edge
				auto newvh = add_vertex(toOM(p_close));
				//auto he = he_close;
				//std::cout << "p = " << p[0] << " " << p[1] << " " << p[2] << std::endl;
				//std::cout << "pclose = " << p_close[0] << " " << p_close[1] << " " << p_close[2] << std::endl;
				//std::cout << "h0 = " << point(he_close.from()) << "; " << "h1 = " << point(he_close.to()) << std::endl;
				//fmt::print("d0 = {:e}, d1 = {:e}\n", (point(he_close.from()) - toOM(p_close)).norm(), (point(he_close.to()) - toOM(p_close)).norm());
				if (cosangle(toEigen(p_close), he) > 1 - 1e-4) {
					//std::cout << point(he.from()) << " " << point(he.to()) << ";";
					//he = he.next();
					continue;
				}
				//std::cout << std::endl;
				split(he.edge(), newvh);
				//splited[he_close].push_back(newvh);
				// set to new position
				auto newp = toOM(trans) + toOM(p_close);
				set_point(v.first, newp);
				//if (point(newvh) != point(he.to())) {
				//	std::cout << "newv   = " << point(newvh) << std::endl;
				//	std::cout << "he     = " << point(he.from()) << " ; " << point(he.to()) << std::endl;
				//	std::cout << "henext = " << point(he.next().from()) << " ; " << point(he.next().to()) << std::endl;
				//	std::cout << "heprev = " << point(he.prev().from()) << " ; " << point(he.prev().to()) << std::endl;
				//	throw std::runtime_error("topological error when splitting edge");
				//}
				if (point(he.to()) == toOM(p_close)) {
					hesplit.push_back(he.next()); break;
				} else {
					hesplit.push_back(he.prev()); break;
				}
			}
		}
	}


	garbage_collection();
	if (debug_level > 0)OM::IO::write_mesh(*this, getPath("mergeedge.obj"));
	//setPeriod();
#endif
}

std::vector<OpenMesh::SmartVertexHandle> msf::PeriodSurfaceMesh::mergePeriodVertices(void)
{
#if 0
	//{
	//	std::ofstream ofs(getPath("vper"));
	//	for (auto vh : vertices()) { if (_vflags[vh.idx()].is_period_boundary()) { ofs << point(vh) << std::endl; } }
	//}
	std::map<OM::SmartFaceHandle, std::array<OM::SmartVertexHandle, 3>> fold2new;
	//std::map<OM::SmartFaceHandle, Eigen::Vector3<Real>> vtrans4face;
	std::unordered_set<Eigen::Vector3<Real>, Vector3Hash<Real>> vPlist;
	std::map<OM::SmartVertexHandle, OM::SmartVertexHandle> toVmin;
	for (int i = 0; i < _vPeriodGroup.size(); i++) {
		auto grp = _vPeriodGroup[i];
		OM::SmartVertexHandle vminh;
		for (int k = 0; k < grp.size(); k++) {
			if (!_vflags[grp[k].idx()].is_max_period()) { vminh = OM::SmartVertexHandle(grp[k].idx(), this); break; }
		}
		// this case actually occurs, for mesh vertex happens to be the unit cell corner, but only 6 or 4 of them
		if (!vminh.is_valid()) {
			for (int k = 0; k < grp.size(); k++) {
				if (_vflags[grp[k].idx()].is_min_period()) { vminh = OM::SmartVertexHandle(grp[k].idx(), this); break; }
			}
		}
		for (int k = 0; k < grp.size(); k++) { toVmin[toSmart(grp[k])] = vminh; }
	}
	for (int i = 0; i < _vPeriodGroup.size(); i++) {
		auto grp = _vPeriodGroup[i];
		OM::SmartVertexHandle vminh;
		for (int k = 0; k < grp.size(); k++) {
			if (!_vflags[grp[k].idx()].is_max_period()) {
				vminh = OM::SmartVertexHandle(grp[k].idx(), this);
				break;
			}
		}
		// this case actually occurs, for mesh vertex happens to be the unit cell corner, but only 6 or 4 of them
		if (!vminh.is_valid()) {
			for (int k = 0; k < grp.size(); k++) {
				if (_vflags[grp[k].idx()].is_min_period()) { vminh = OM::SmartVertexHandle(grp[k].idx(), this); break; }
			}
		}
		if (!_vflags[vminh.idx()].is_min_period()) { std::cout << "Warning : no totally left period boundary vertices" << std::endl; }
		vPlist.insert(toEigen(point(vminh)));
		// delete old faces, add new
		for (int k = 0; k < grp.size(); k++) {
			if (grp[k] == vminh) continue;
			//Eigen::Vector3<Real> vtrans(0, 0, 0);
			//for (int j = 0; j < 3; j++) {
			//	if (_vflags[grp[k].idx()].is_max_period(j)) {
			//		vtrans[j] = 2;
			//	}
			//}
			for (auto vf : vf_range(grp[k])) {
				auto fv = getFaceVertexHandle(vf);
				for (int m = 0; m < 3; m++) {
					if (toVmin.count(fv[m])) {
						fv[m] = toVmin[fv[m]];
					}
				}
				//if (fv[2] != grp[k]) { printf("\033[31mInvalid face\033[0m\n"); }
				fold2new[vf] = { fv[0], fv[1], fv[2] };
				//vtrans4face[vf] = vtrans;
			}
		}
	}

	//{ std::ofstream ofs(getPath("fper")); for (auto f2n : fold2new) { ofs << f2n.first.idx() << std::endl; } }

	using vec_t = Eigen::Vector3<Real>;
	std::unordered_map<std::pair<vec_t, vec_t>, HETag, Vector3PairHash<Real>> he2tag;
	std::vector<OM::SmartHalfedgeHandle> hetrans;
	for (auto f2n : fold2new) {
		delete_face(f2n.first, true);
		auto fh = add_face(f2n.second[0], f2n.second[1], f2n.second[2]);
		if (!fh.is_valid()) {
			std::cout << "\033[31madd_face failed for " << f2n.second[0] << ", " << f2n.second[1] << ", " << f2n.second[2] << " \033[0m\n";
			throw std::runtime_error("failed to add face");
		}
		int overflow = 100;
		for (auto fhe : fh.halfedges()) {
			auto v0 = fhe.from(), v1 = fhe.to();
			HETag tag;
			for (int k = 0; k < 3; k++) {
				// if min boundary, it must be transfered
				tag.torus_trans[k] = (_vflags[v0.idx()].is_min_period(k) - _vflags[v1.idx()].is_min_period(k)) * 2;
			}
			tag.updateFlag();
			auto fhe_0 = toEigen(point(fhe.from())), fhe_1 = toEigen(point(fhe.to()));
			he2tag[{fhe_0, fhe_1}] = tag; 
			if (overflow-- < 0) {
				auto hit = fh.halfedge();
				std::cout << hit << " -> " << hit.next() << " -> " << hit.next().next() << std::endl;
				throw std::runtime_error("fh overflow, maybe bad face handle");
			}

		}
	}
	garbage_collection();

	//std::vector<HETag> hetags(n_halfedges());
	//// build hetag
	//for (auto he : halfedges()) {
	//	std::pair<vec_t, vec_t> he_01 = { toEigen(point(he.from())),  toEigen(point(he.to())) };
	//	if (he2tag.count(he_01)) {
	//		auto tag = he2tag[he_01];
	//		hetags[he.idx()] = tag;
	//		tag.torus_trans *= -1;
	//		hetags[he.opp().idx()] = tag;
	//	}
	//}
	// collect vcut
	std::vector<OM::SmartVertexHandle> vTlist;
	for (auto vh : vertices()) {
		if (vPlist.count(toEigen(point(vh)))) {
			vTlist.push_back(vh);
		}
	}
	//std::set<OM::SmartVertexHandle> vHset(vTlist.begin(), vTlist.end());
	//{ std::ofstream ofs(getPath("vright")); for (auto vH : vHset) { ofs << vH.idx() << std::endl; } }
	//for (int i = 0; i < hetags.size(); i++) {
	//	auto he = OM::SmartHalfedgeHandle(i, this);
	//	if (!hetags[i].has_trans()) continue;
	//	if (vHset.count(he.from())) vHset.erase(he.from());
	//	if (vHset.count(he.to())) vHset.erase(he.to());
	//}
	//std::cout << "vleft size = " << vHset.size() << std::endl;
	//{ std::ofstream ofs(getPath("vleft")); for (auto vH : vHset) { ofs << vH.idx() << std::endl; } }

	//OM::IO::write_mesh(*this, getPath("meshperiod.obj"));
	//savePeriodicMesh(getPath("mergevertx.obj"), vTlist);

	return vTlist;
#else
garbage_collection();
auto [vlist, flist] = getVFlist();
//std::tie(vlist, flist) = remove_dup_vertices(vlist, flist, 1e-5);
PeriodicGridIndex indexer(Eigen::Vector3<Real>(-1, -1, -1), Eigen::Vector3<Real>(2, 2, 2), 1e-6);
for (int i = 0; i < vlist.rows(); i++) {
	indexer.insert(vlist.row(i).transpose());
}
auto new_vlist = indexer.dumpPoints();
auto new_flist = flist;
int f_counter = 0;
for (int i = 0; i < flist.rows(); i++) {
	int fv[3];
	for (int k = 0; k < 3; k++) {
		fv[k] = indexer.query(vlist.row(flist(i, k)));
	}
	if (fv[0] == fv[1] || fv[1] == fv[2] || fv[0] == fv[2]) {
		continue;
	}
	new_flist.row(f_counter++) = Eigen::Vector3i(fv[0], fv[1], fv[2]).transpose();
}
new_flist.conservativeResize(f_counter, 3);
// orient face coherently
//Eigen::VectorXi VF, Ni;
Eigen::MatrixX3i FF;
//igl::vertex_triangle_adjacency(flist, new_vlist.rows(), VF, Ni);
igl::triangle_triangle_adjacency(new_flist, FF);
std::vector<bool> passed(new_flist.rows(), false);
std::queue<std::pair<int, int>> frontier;
frontier.emplace(0, 0);
passed[0] = true;
do {
	auto fr = frontier.front().second;
	auto fr_src = frontier.front().first;
	frontier.pop();
	if (fr != fr_src) {
		Eigen::Vector3i src = new_flist.row(fr_src).transpose();
		Eigen::Vector3i dst = new_flist.row(fr).transpose();
		int common_in_src[2], common_in_dst[2];
		int counter = 0;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (src[i] == dst[j]) {
					common_in_src[counter] = i; common_in_dst[counter] = j;
					counter++; break;
				}
			}
		}
		//common_in_dst[1] = (common_in_dst[1] + 3 + (common_in_dst[0] - common_in_src[0])) % 3;
		if (((common_in_dst[1] - common_in_dst[0]) - (common_in_src[1] - common_in_src[0])) % 3 == 0) {
			// this means wrong orientation
			std::swap(new_flist(fr, common_in_dst[0]), new_flist(fr, common_in_dst[1]));
		} 	
	}
	// add neighbor frontier 
	for (int i = 0; i < 3; i++) {
		try {
			if (passed.at(FF(fr, i))) continue;
		} catch (...) {
			printf("\033[31mOut of Range : FF(%d, %d) = %d\033[0m\n", fr, i, (int)FF(fr, i));
			std::cout << "F =\n" << new_vlist.row(new_flist(fr, 0)) << "\n" <<
				new_vlist.row(new_flist(fr, 1)) << "\n" <<
				new_vlist.row(new_flist(fr, 2)) << "\n";
		}
		passed[FF(fr, i)] = true;
		frontier.emplace(fr, FF(fr, i));
	}
} while (!frontier.empty());

OM_Mesh::clear();

for (int i = 0; i < new_vlist.rows(); i++) {
	add_vertex(toOM(new_vlist.row(i)));
}
for (int i = 0; i < new_flist.rows(); i++) {
	std::vector<OM::VertexHandle> vhlist;
	for (int k = 0; k < 3; k++) vhlist.emplace_back(new_flist(i, k));
	if (vhlist[0] == vhlist[1] || vhlist[1] == vhlist[2] || vhlist[0] == vhlist[2]) continue;
	auto fh = add_face(vhlist);
	if (!fh.is_valid()) {
		//OM::IO::write_mesh(*this, getPath("nonmanifold.obj"));
		std::cout << "face = " << point(vhlist[0]) << ";" << point(vhlist[1]) << ";" << point(vhlist[2]) << std::endl;
	}
}
request_normal_status();

//OM::IO::write_mesh(*this, getPath("mergevertex.obj"));

std::vector<OM::SmartVertexHandle> vTlist;
for (auto vh : vertices()) {
	VertexFlag vflag;
	auto p = point(vh);
	vflag.set_period_boundary(0.9999, p[0], p[1], p[2]);
	if (vflag.is_min_period()) {
		vTlist.push_back(vh);
	}
}
//savePeriodicMesh(getPath("mergevertex.obj"), vTlist);
return vTlist;
#endif
}

std::vector<OpenMesh::SmartVertexHandle> msf::PeriodSurfaceMesh::readMergePeriodBoundary(std::string fielanme)
{
	read(fielanme, true, false);
	mergePeriodEdges();
	return mergePeriodVertices();
}

std::vector<msf::OM::SmartVertexHandle> msf::PeriodSurfaceMesh::readMergePeriodBoundary(std::string fielanme, bool normalize)
{
	read(fielanme, normalize, false);
	mergePeriodEdges();
	return mergePeriodVertices();
}

std::set<msf::OM::SmartHalfedgeHandle> msf::PeriodSurfaceMesh::getBoundaryHalfedges(void)
{
	std::set<msf::OM::SmartHalfedgeHandle> helist;
	for (auto he : halfedges()) {
		if (status(he).deleted() || !he.is_valid()) continue;
		if (he.is_boundary()) { helist.insert(he); }
	}
	return helist;
}


