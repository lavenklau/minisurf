#define FMT_HEADER_ONLY
#include "PeriodicMesher.h"
#include "isosurface_generation.h"
#include "igl/writeOBJ.h"
#include <random>
#include <fmt/core.h>
#include <igl/connected_components.h>
#include <igl/adjacency_matrix.h>
#include <igl/AABB.h>
#include <igl/boundary_loop.h>
#include <igl/boundary_facets.h>
#include <igl/edges.h>
#include "matlab/matlab_utils.h"
#include "cgal/cgal_utils.h"


using namespace msf;

extern std::string getPath(std::string);
extern std::string out_mesh_type;

std::tuple<Eigen::MatrixX3<msf::Real>, Eigen::MatrixX3i> removeDupVertices(const Eigen::MatrixX3<msf::Real>& v, const Eigen::MatrixX3i& f, msf::Real eps);

void clamp_vertices(Eigen::MatrixX3<Real>& V, Eigen::MatrixX3i& F);
std::vector<msf::Real> genRandom(int number, Real low, Real upp) {
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::mt19937 gen(seed);
	// Define the range for the random number
	std::uniform_real_distribution<msf::Real> distrib(low, upp);
	std::vector<Real> seqRand;
	for (int i = 0; i < number; i++) { seqRand.push_back(distrib(gen)); }
	return seqRand;
}

msf::Real genRandom(Real low, Real upp) {
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::mt19937 gen(seed);
	// Define the range for the random number
	std::uniform_real_distribution<msf::Real> distrib(low, upp);
	return distrib(gen);
}

std::string GetCurrentDateTime() {
	// Get current time as time_point
	auto now = std::chrono::system_clock::now();
	// Convert time_point to time_t for converting to tm (broken-down time)
	std::time_t now_c = std::chrono::system_clock::to_time_t(now);
	// Convert to broken-down time
	std::tm now_tm = *std::localtime(&now_c);
	// Use stringstream to format the date and time
	std::stringstream ss;
	ss << std::put_time(&now_tm, "%Y%m%d-%H%M%S");
	return ss.str();
}

template<typename Scalar>
using Interval = std::pair<Scalar, Scalar>;

template<typename Scalar>
std::vector<Interval<Scalar>> mergeInterval(std::vector<Interval<Scalar>>& intervals) {
	std::vector<Interval<Scalar>> result;
	if (intervals.size() == 0) return result; // 区间集合为空直接返回
	// 排序的参数使用了lambda表达式
	std::sort(intervals.begin(), intervals.end(), [](const Interval<Scalar>& a, const Interval<Scalar>& b) {return a.first < b.first; });
	// 第一个区间就可以放进结果集里，后面如果重叠，在result上直接合并
	result.push_back(intervals[0]);
	for (int i = 1; i < intervals.size(); i++) {
		if (result.back().second >= intervals[i].first) { // 发现重叠区间
			// 合并区间，只更新右边界就好，因为result.back()的左边界一定是最小值，因为我们按照左边界排序的
			result.back().second = (std::max)(result.back().second, intervals[i].second);
		} else {
			result.push_back(intervals[i]); // 区间不重叠 
		}
	}
	return result;
}


std::pair<Eigen::Vector3<Real>, Eigen::Vector3i> sampleMeshTrigonometric(void) 
{
	
}

std::pair<Eigen::MatrixX3<Real>, Eigen::MatrixX3i>
getPeriodMesh(const Eigen::MatrixX3<Real>& V, const Eigen::MatrixX3i& F, Real period = 2, Real err = 1e-6)
{
	if (F.rows() == 0) return { V,F };
	// boundary vertices
	Eigen::MatrixX2i be;
	//igl::boundary_loop(F, bid);
	igl::boundary_facets(F, be);
	std::vector<int> bid(be.data(), be.data() + be.size());
	{ std::sort(bid.begin(), bid.end());
		auto uit = std::unique(bid.begin(), bid.end());
		bid.erase(uit, bid.end()); }
	//{ std::ofstream ofs(getPath("bid")); ofs << Eigen::VectorXi::Map(bid.data(), bid.size()); }
	std::map<int, VertexFlag> vflag;
	std::map<int, int> rootid;
	std::vector<std::pair<int, int>> period_pairs;
	for (int i = 0; i < bid.size(); i++) {
		for (int j = i + 1; j < bid.size(); j++) {
			Eigen::Vector3<Real> v12 = (V.row(bid[j]) - V.row(bid[i])).transpose();
			bool not_period_pair = false;
			for (int k = 0; k < 3; k++) {
				if (std::abs(v12[k]) > err && std::abs(v12[k]) < period - err) { not_period_pair = true; break; }
			}
			if (not_period_pair) continue;

			bool i_ge_j = true, j_ge_i = true;
			for (int k = 0; k < 3; k++) {
				if (v12[k] < -period + err) {
					//if (bid[i] == 23306) {
					//	std::cout << "j = " << bid[j] << ":" << V.row(bid[j]) << ", i = " << bid[i] << ":" << V.row(bid[i]) << std::endl;
					//}
					vflag[bid[i]].set_max_boundary(k);
					vflag[bid[j]].set_min_boundary(k);
					j_ge_i = false; 
				} else if (v12[k] > period - err) {
					//if (bid[j] == 23306) {
					//	std::cout << "j = " << bid[j] << ":" << V.row(bid[j]) << ", i = " << bid[i] << ":" << V.row(bid[i]) << std::endl;
					//}
					vflag[bid[i]].set_min_boundary(k);
					vflag[bid[j]].set_max_boundary(k);
					i_ge_j = false; 
				} else if (std::abs(v12[k]) < err) {}
				else {
					i_ge_j = j_ge_i = false;
				}
			}
			if (i_ge_j || j_ge_i) {
				period_pairs.emplace_back(bid[i], bid[j]);
			}
		}
	}
	for (int i = 0; i < period_pairs.size(); i++) {
		auto pr = period_pairs[i];
		auto iflg = vflag[pr.first].getPeriodFlag();
		auto jflg = vflag[pr.second].getPeriodFlag();
		//if (pr.first == 23306 || pr.second == 23306) {
		//	std::cout << "\033[32m i = " << pr.first << ":" << iflg << ", j = " << pr.second << ":" << jflg << "\033[0m" << std::endl;
		//}
		if (jflg > iflg) {
			rootid[pr.second] = pr.first;
		} else {
			rootid[pr.first] = pr.second;
		}
	}
	{
		//std::ofstream ofs(getPath("linkage"));
		//for (auto lk : rootid) { ofs << lk.first << " " << lk.second << std::endl; }
	}
	//std::cout << "============================================================" << std::endl;
	std::set<int> vmax, vmin;
	for (auto tag : vflag) {
		if (tag.second.is_max_period()) { vmax.insert(tag.first) ; }
		if (tag.second.is_min_period()) { vmin.insert(tag.first); }
	}
	std::vector<Eigen::Vector3<Real>> vper;
	std::map<int, int> perid;
	for (int i = 0; i < V.rows(); i++) {
		if (!vmax.count(i)) {
			perid[i] = vper.size();
			vper.push_back(V.row(i).transpose());
		}
	}
	for (int i = 0; i < V.rows(); i++) {
		if (vmax.count(i)) {
			int minid = i;
			while (rootid.count(minid)) {
				//std::cout << minid << "->";
				minid = rootid[minid]; 
			}
			if (vflag[minid].is_max_period()) { std::cout << "\033[31m" << "minid is max period " << minid << " = " << vflag[minid]._flag << "\033[0m\n"; }
			//std::cout << std::endl;
			perid[i] = perid[minid];
			//std::cout << "i = " << i << ", minid = " << minid << ", perid = " << perid[i] << std::endl;
		}
	}
	Eigen::MatrixX3<Real> Vnew(vper.size(), 3);
	for (int i = 0; i < vper.size(); i++) {
		Vnew.row(i) = vper[i].transpose();
	}
	Eigen::MatrixX3i Fnew(F.rows(), 3);
	for (int i = 0; i < F.size(); i++) {
		int oldid = F.data()[i];
		int newid = perid[oldid];
		Fnew.data()[i] = newid;
	}
	return { Vnew, Fnew };
}

bool is_t3_connected(const Eigen::MatrixX3<Real>& V, const Eigen::MatrixX3i& F, Real period = 1) {
	auto [Vper, Fper] = getPeriodMesh(V, F, 1);
	if (Fper.rows() == 0) return false;
	//PeriodSurfaceMesh::savePeriodicMesh(getPath("permesh-per.obj"), Vper, Fper, 1, 0.5);
	Eigen::SparseMatrix<Real> A;
	igl::adjacency_matrix(Fper, A);
	Eigen::VectorXi Cid,Ks;
	igl::connected_components(A, Cid, Ks);
	Eigen::MatrixX2i edgs;
	igl::edges(A, edgs);
	std::cout << Ks.rows() << " Components, " << edgs.rows() << " Edges, " << Vper.rows() << " Vertices" << std::endl;
	// *[axis][commponent][edge_interval_id]
	std::vector<std::vector<Interval<Real>>> edgeAxisIntervals[3];
	for (int i = 0; i < 3; i++) {
		// resize to number number of component
		edgeAxisIntervals[i].resize(Ks.rows());
	}
	for (int i = 0; i < edgs.rows(); i++) {
		int vi = edgs(i, 0), vj = edgs(i, 1);
		int c_id = Cid[vi];
		Eigen::Vector3<Real> pi = Vper.row(vi), pj = Vper.row(vj);
		Eigen::Vector3<Real> vec = (pj - pi).transpose();
		vec = make_period(vec, period, period * 0.5);
		for (int k = 0; k < 3; k++) {
			if (vec[k] > 0)
				edgeAxisIntervals[k][c_id].push_back(Interval<Real>(pi[k], pi[k] + vec[k]));
			else if (vec[k] < 0)
				edgeAxisIntervals[k][c_id].push_back(Interval<Real>(pj[k], pj[k] - vec[k]));
		}
	}
	std::vector<bool> connected(Ks.rows(), true);
	for (int k = 0; k < 3; k++) {
		for (int cid = 0; cid < edgeAxisIntervals[k].size(); cid++) {
			if (!connected[cid])  continue;
			edgeAxisIntervals[k][cid] = mergeInterval(edgeAxisIntervals[k][cid]);
			if (edgeAxisIntervals[k][cid][0].second - edgeAxisIntervals[k][cid][0].first < 1 - 1e-4) {
				connected[cid] = false;
			}
		}
	}
	for (int i = 0; i < Ks.rows(); i++) {
		if (!connected[i]) return false;
	}
	return true;
}

Eigen::MatrixX3<Real> generateSampleLocations(
	int reso, const std::pair<Eigen::Vector3<Real>, Eigen::Vector3<Real>>& bbox)
{
	int gs = reso + 1;
	int nsamples = gs * gs * gs;
	Eigen::MatrixX3<Real> sampleLocations(nsamples, 3);
	{
		size_t i = 0;
		for (size_t zi = 0; zi < gs; ++zi) {
			for (size_t yi = 0; yi < gs; ++yi) {
				for (size_t xi = 0; xi < gs; ++xi) {
					Eigen::Vector3<Real> p =
						bbox.first + (bbox.second - bbox.first).cwiseProduct(Eigen::Vector3<Real>(xi, yi, zi) / reso);
					sampleLocations.row(i) = p.transpose();
					++i;
				}
			}
		}
	}
	return sampleLocations;
}

// 2      :  rand choose 4 basis TPMS
// 10-14  :  manually choose (n-10)-th TPMS
//std::vector<std::function<Real(Real, Real, Real)>> sample_type_postprocessing(
//	int n, int sample_type, std::vector<double> c_cos[3], std::vector<double> c_sin[3],
//	std::vector<double> c_cos_cos[9], std::vector<double> c_sin_cos[9], std::vector<double> c_sin_sin[9], std::vector<int>& basislist,
//	Real randomness
//) 
//{
//	if (sample_type == 2) {
//		return std::vector<std::function<Real(Real, Real, Real)>>(n, [](Real, Real, Real) {return Real(0); });
//	} else if (sample_type >= 10 && sample_type < 14) {
//	}
//	else {
//		throw std::runtime_error("invalid sample type");
//	}
//	auto chooseBasis = genRandom(n, 0, 4);
//	basislist.resize(n);
//	//for (int k = 0; k < n; k++) basislist.push_back(chooseBasis[k]);
//	for (int i = 0; i < 3; i++) {
//		for (int k = 0; k < n; k++) {
//			c_cos[i][k] *= randomness;
//			c_sin[i][k] *= randomness;
//		}	
//		for (int j = 0; j < 3; j++) {
//			for (int k = 0; k < n; k++) {
//				c_cos_cos[i * 3 + j][k] *= randomness;
//				c_sin_cos[i * 3 + j][k] *= randomness;
//				c_sin_sin[i * 3 + j][k] *= randomness;
//			}
//		}
//	}
//
//	std::vector<std::function<Real(Real, Real, Real)>> postfunc;
//	for (int k = 0; k < n; k++) {
//		std::function<Real(Real, Real, Real)> postfix = [](Real x, Real y, Real z) ->Real { return 0; };
//		int basis = chooseBasis[k];
//		if (sample_type >= 10 && sample_type < 14) {
//			basis = sample_type - 10; 
//			chooseBasis[k] = basis;
//			basislist[k] = basis;
//		} 
//		if (basis == 0) {
//			// P surface
//			// cos 2pi x + cos 2pi y + cos 2pi z = 0
//			c_cos[0][k] += 1; c_cos[1][k] += 1; c_cos[2][k] += 1;
//		} else if (basis == 1) {
//			// G surface
//			// sin 2pix cos 2piy + sin 2piz cos 2pix + sin 2piy cos 2piz = 0
//			c_sin_cos[0 * 3 + 1][k] += 1; c_sin_cos[2 * 3 + 0][k] += 1; c_sin_cos[1 * 3 + 2][k] += 1;
//		} else if (basis == 2) {
//			// D surface
//			// cos 2pix cos 2piy cos 2piz - sin 2pix sin 2piy sin 2piz = 0
//			postfix = [](Real x, Real y, Real z) -> Real {
//				return std::cos(2 * M_PI * x) * std::cos(2 * M_PI * y) * std::cos(2 * M_PI * z) -
//					std::sin(2 * M_PI * x) * std::sin(2 * M_PI * y) * std::sin(2 * M_PI * z);
//				};
//		} else if (basis = 3) {
//			// IWP surface
//			// 2 ( cos 2pix cos 2piy + cos 2piy cos 2piz + cos 2piz cos 2pix ) - cos 2*2pix - cos 2*2piy - cos 2*2piz =0
//			c_cos_cos[0 * 3 + 1][k] += 2; c_cos_cos[1 * 3 + 2][k] += 2; c_cos_cos[2 * 3 + 0][k] += 2;
//			c_cos_cos[0 * 3 + 0][k] -= 2; c_cos_cos[1 * 3 + 1][k] -= 2; c_cos_cos[2 * 3 + 2][k] -= 2;
//			postfix = [](Real x, Real y, Real z) -> Real { return 3; };
//		}
//		postfunc.push_back(postfix);
//	}
//	return postfunc;
//}

std::function<Real(Real[3], Real[3])> sample_type_postprocessing(
	int n, int sample_type, double c_cos[3], double c_sin[3],
	double c_cos_cos[9], double c_sin_cos[9], double c_sin_sin[9], int& basis,
	Real randomness
) {
	if (sample_type == 2) {
		return [](Real[3], Real[3]) {return Real(0); };
	} else if (sample_type >= 10 && sample_type < 20) { }
	else {
		throw std::runtime_error("invalid sample type");
	}
	auto chooseBasis = genRandom(0, 4);
	//for (int k = 0; k < n; k++) basislist.push_back(chooseBasis[k]);
	for (int i = 0; i < 3; i++) {
		c_cos[i] *= randomness;
		c_sin[i] *= randomness;
		for (int j = 0; j < 3; j++) {
			c_cos_cos[i * 3 + j] *= randomness;
			c_sin_cos[i * 3 + j] *= randomness;
			c_sin_sin[i * 3 + j] *= randomness;
		}
	}

	std::function<Real(Real[3], Real[3])> postfunc;
	std::function<Real(Real[3], Real[3])> postfix = [](Real[3], Real[3]) ->Real { return 0; };
	basis = chooseBasis;
	if (sample_type >= 10 && sample_type < 20) {
		basis = sample_type - 10;
		chooseBasis = basis;
	}
	if (basis == 0) {
		// P surface
		// cos 2pi x + cos 2pi y + cos 2pi z = 0
		postfix = [](Real cosx[3], Real sinx[3]) -> Real {
			return cosx[0] + cosx[1] + cosx[2];
		};
	} else if (basis == 1) {
		// G surface
		// sin 2pix cos 2piy + sin 2piz cos 2pix + sin 2piy cos 2piz = 0
		postfix = [](Real cosx[3], Real sinx[3]) -> Real {
			return sinx[0] * cosx[1] + sinx[2] * cosx[0] + sinx[1] * cosx[2];
			};
	} else if (basis == 2) {
		// D surface
		// cos 2pix cos 2piy cos 2piz - sin 2pix sin 2piy sin 2piz = 0
		postfix = [](Real cosx[3], Real sinx[3]) -> Real {
			return cosx[0] * cosx[1] * cosx[2] - sinx[0] * sinx[1] * sinx[2];
			};
	} else if (basis == 3) {
		// IWP surface
		// 2 ( cos 2pix cos 2piy + cos 2piy cos 2piz + cos 2piz cos 2pix ) - cos 2*2pix - cos 2*2piy - cos 2*2piz =0
		postfix = [](Real cosx[3], Real sinx[3]) -> Real {
			Real val = 1.5;
			val += cosx[0] * cosx[1] + cosx[1] * cosx[2] + cosx[2] * cosx[0] -
				cosx[0] * cosx[0] - cosx[1] * cosx[1] - cosx[2] * cosx[2];
			return val; 
			};
	} else if (basis == 4) {
		// Schwartz D (diagmond)
		postfix = [](Real c[3], Real s[3]) -> Real {
			return s[0] * s[1] * s[2] + s[0] * c[1] * c[2] + c[0] * s[1] * c[2] + c[0] * c[1] * s[2];
			};
	} else if (basis == 5) {
		// Neovius
		postfix = [](Real cosx[3], Real sinx[3]) -> Real {
			return  cosx[0] + cosx[1] + cosx[2] + 4. / 3 * cosx[0] * cosx[1] * cosx[2];
			};
	}
	return postfix;
}

void rescale_vertices(Eigen::MatrixX3<Real>& plist) {
	Eigen::Vector3<Real> pmin, pmax;
	pmin.setConstant(1e30);
	pmax.setConstant(-1e30);
	for (int i = 0; i < plist.rows(); i++) {
		pmin = pmin.cwiseMin(plist.row(i).transpose());
		pmax = pmax.cwiseMax(plist.row(i).transpose());
	}
	Eigen::Vector3<Real> bmin(-1, -1, -1);
	Eigen::Vector3<Real> bmax(1, 1, 1);
	Eigen::Vector3<Real> s = (bmax - bmin).cwiseQuotient(pmax - pmin);
	for (int i = 0; i < plist.rows(); i++) {
		plist.row(i) = bmin.transpose() + (plist.row(i) - pmin.transpose()).cwiseProduct(s.transpose());
	}
}

void msf::PeriodSurfaceMesh::sampleMeshes(int number, int reso, int sample_type, Real randomness, int remshIter, int smooth_iter, Real edge_lengh)
{
	//int reso = 32;
	std::pair<Eigen::Vector3<Real>, Eigen::Vector3<Real>> bbox;
	bbox.first = Eigen::Vector3<Real>::Zero();
	bbox.second = Eigen::Vector3<Real>::Ones();
	auto p_sample = generateSampleLocations(reso, bbox);
	int n_sample = p_sample.rows();
	std::vector<Real> cosx[3], sinx[3];
	for (int k = 0; k < 3; k++) {
		cosx[k].resize(p_sample.rows()); sinx[k].resize(p_sample.rows());
		for (int i = 0; i < p_sample.rows(); i++) {
			Real pk = p_sample(i, k);
			//pk = (pk - 0.5) / (1 - 1e-6) + 0.5;
			pk = pk + 1e-6;
			cosx[k][i] = std::cos(2 * M_PI * pk);
			sinx[k][i] = std::sin(2 * M_PI * pk);
		}
	}
	double c_cos[3], c_sin[3], c_cos_cos[9], c_sin_cos[9], c_sin_sin[9];
	// get time string
	auto time_str = GetCurrentDateTime();
	
	int overflow = 1e5;
	// sample loop until generated enough
	for (int iter = 0; iter < number; ) {
		if (overflow-- < 0) break;
		double iso = 0;
		for (int i = 0; i < 3; i++) {
			c_cos[i] = genRandom(-1, 1);
			c_sin[i] = genRandom(-1, 1);
			for (int j = 0; j < 3; j++) {
				c_cos_cos[i * 3 + j] = genRandom( -1, 1);
				c_sin_cos[i * 3 + j] = genRandom( -1, 1);
				c_sin_sin[i * 3 + j] = genRandom( -1, 1);
			}
		}
		int basisType = -1;
		std::function<Real(Real[3], Real[3])> postfix;
		try {
			postfix = sample_type_postprocessing(number, sample_type, c_cos, c_sin, c_cos_cos, c_sin_cos, c_sin_sin, basisType, randomness);
		} catch (std::exception e) {
			std::cout << "invalid sample type : " << e.what() << std::endl;
			throw std::runtime_error("sample mesh failed");
		}
		// output random coefficient
		{
			std::ofstream ofs(getPath(fmt::format("{}-tpmscoeff", time_str)), std::ios::app);
			ofs << fmt::format("{}--m{}.obj :", time_str, iter);
			ofs << "\nsampleType = " << sample_type;
			ofs << "\nbasisType = " << basisType;
			ofs << "\ncos = " << c_cos[0] << " " << c_cos[1] << " " << c_cos[2];
			ofs << "\nsin = " << c_sin[0] << " " << c_sin[1] << " " << c_sin[2];
			ofs << "\ncoscos = "; for (int j = 0; j < 9; j++) { ofs << c_cos_cos[j] << " "; }
			ofs << "\nsincos = "; for (int j = 0; j < 9; j++) { ofs << c_sin_cos[j] << " "; }
			ofs << "\nsinsin = "; for (int j = 0; j < 9; j++) { ofs << c_sin_sin[j] << " "; }
			ofs << "\nisovalue = " << iso << std::endl;
			ofs.close();
		}
		std::cout << "Meshing sample " << iter << "..." << std::endl;
		Eigen::VectorX<Real> values(n_sample);
		values.setZero();
		Eigen::Matrix3X<Real> plist = p_sample.transpose();
		for (int n = 0; n < n_sample; n++) {
			Real val = 0;
			Real cx[3] = { cosx[0][n], cosx[1][n], cosx[2][n] };
			Real sx[3] = { sinx[0][n], sinx[1][n], sinx[2][n] };
			// weighted sum of basis function
			for (int i = 0; i < 3; i++) {
				val += c_cos[i] * cx[i];
				val += c_sin[i] * sx[i];
				for (int j = 0; j < 3; j++) {
					val += c_cos_cos[i * 3 + j] * cx[i] * cx[j];
					val += c_sin_cos[i * 3 + j] * sx[i] * cx[j];
					val += c_sin_sin[i * 3 + j] * sx[i] * sx[j];
				}
			}
			if (sample_type == 2 || (sample_type >= 10 && sample_type < 20)) {
				Eigen::Vector3<Real> p = plist.col(n);
				val += postfix(cx, sx);
			}
			values[n] = val - iso;
		}
		auto [V, F] = isosurf_mesh(reso, p_sample, values);
		//{ igl::writeOBJ(getPath("isomesh.obj"), V, F); }
		remove_dup_vertices(V, F);
		bool connected = is_t3_connected(V, F);
		std::cout << "Sample " << iter << " is " << (connected ? "Connected" : "Disconnected") << std::endl;
		if (connected) {
			rescale_vertices(V);
			MeshSurface mtgt; mtgt.read(V, F);
			std::tie(V, F) = cgal_remesh(V, F);
			std::tie(V, F) = removeDupVertices(V, F, 2e-5);
			clamp_vertices(V, F);
			std::cout << "Remeshing..." << std::endl;
			std::vector<OM::SmartVertexHandle> vcut;
			try {
				read(V, F, true, false);
				vcut = mergePeriodBoundary();
				periodic_remesh(remshIter, vcut, edge_lengh, smooth_iter, 0, &mtgt);
			} catch (std::exception e) {
				std::cout << "Error : " << e.what() << ", continuing" << std::endl;
				continue;
			}
			//igl::writeOBJ(getPath(fmt::format("{}--m{}.obj", time_str, iter)), V, F); 
			if (out_mesh_type.empty()) out_mesh_type = "obj";
			savePeriodicMesh(getPath(fmt::format("{}--m{}.{}", time_str, iter, out_mesh_type)), vcut);
			std::cout << "Remesh finished, " << n_vertices() << " vertices, " << n_faces() << " faces." << std::endl;
			iter++;
		}
	}
}

void msf::PeriodSurfaceMesh::savePeriodicMesh(std::string filename, const Eigen::MatrixX3<Real>& vlist, const Eigen::MatrixX3i& flist, Real period, Real detach /*= 0.5*/)
{
	auto newflist = flist;
	std::unordered_map<Eigen::Vector3<Real>, int, Vector3Hash<Real>> vnew2id;
	for (int i = 0; i < vlist.rows(); i++) {
		vnew2id[vlist.row(i).transpose()] = i;
	}
	std::vector<Eigen::Vector3<Real>> vappend;
	for (int i = 0; i < flist.rows(); i++) {
		Eigen::Vector3i fv = flist.row(i).transpose();

		Eigen::Vector3<Real> e[3];
		Eigen::Vector3<Real> vnew[3];
		bool has_period_face = false;
		for (int j = 0; j < 3; j++) { vnew[j] = vlist.row(fv[j]); }
		for (int j = 0; j < 3; j++) {
			e[j] = (vnew[(j + 1) % 3] - vnew[j]).transpose();
			for (int k = 0; k < 3; k++) {
				if (e[j][k] < -detach) {
					vnew[(j + 1) % 3][k] += period;
					has_period_face = true;
				} else if (e[j][k] > detach) {
					vnew[j][k] += period;
					has_period_face = true;
				}
			}
		}
		//. replace old face
		if (has_period_face) {
			for (int j = 0; j < 3; j++) {
				if (!vnew2id.count(vnew[j])) {
					fv[j] = vnew2id.size() + vappend.size();
					vappend.push_back(vnew[j]);
				} else {
					fv[j] = vnew2id[vnew[j]];
				}
			}
			//if (std::set<int>(fv.data(), fv.data() + 3).size() < 3) {
			//	std::cout << "\033[31mWarning duplicated vertices on face\033[0m " << fv.transpose() << " <- " << flist.row(i) << std::endl;
			//}
			newflist.row(i) = fv.transpose();
		}
	}
	int n_old = vlist.rows();
	auto newvlist = vlist;
	newvlist.conservativeResize(n_old + vappend.size(), 3);
	for (int i = 0; i < vappend.size(); i++) {
		newvlist.row(n_old + i) = vappend[i].transpose();
	}
	igl::writeOBJ(filename, newvlist, newflist);

}

bool msf::PeriodSurfaceMesh::read(const Eigen::MatrixX3<Real>& V, const Eigen::MatrixX3i& F, bool normalize /*= true*/, bool set_period /*= true*/, bool removeDup /*= true*/)
{
	Eigen::MatrixX3<Real> vlist = V;
	Eigen::MatrixX3i flist = F;
	if (removeDup) std::tie(vlist, flist) = remove_dup_vertices(V, F);
	OM_Mesh::clean(); OM_Mesh::clear(); OM_Mesh::reset_status();
	std::vector<OM::SmartVertexHandle> vhlist;
	for (int i = 0; i < vlist.rows(); i++) {
		vhlist.push_back(add_vertex(toOM(vlist.row(i))));
	}
	for (int i = 0; i < flist.rows(); i++) {
		//if (std::set<int>{flist(i, 0), flist(i, 1), flist(i, 2)}.size() == 2) {
		//	igl::writeOBJ(getPath("old.obj"), V, F);
		//	igl::writeOBJ(getPath("new.obj"), vlist, flist);
		//	std::cout << "invalid face " << flist.row(i) << std::endl;
		//	eigen2ConnectedMatlab("flist", flist);
		//}
		add_face(vhlist[flist(i, 0)], vhlist[flist(i, 1)], vhlist[flist(i, 2)]);
	}
	request_normal_status();
	updateBb();
	if (normalize)  normalizeBb();
	if (set_period) setPeriod();
	return true;
}



