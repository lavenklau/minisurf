#pragma once

#include <algorithm>
#include <set>
#define USE_NANOFLANN
#ifdef USE_NANOFLANN
#include "nanoflann.hpp"
#endif
#include <unordered_set>
#include <vector>
#include <Eigen/Eigen>
#include "HashTypes.h"
#pragma push_macro("_USE_MATH_DEFINES")
#define _USE_MATH_DEFINES
#include "PeriodicMesher.h"
// #include "utils/HashTypes.h"
#include "math.h"

#ifdef USE_NANOFLANN
namespace nanoflann {

	// 3d torus [-1,1]^3
	template <class T, class DataSource, typename _DistanceType = T>
	struct T3_Adaptor {
		typedef T ElementType;
		typedef _DistanceType DistanceType;

		const DataSource& data_source;

		T3_Adaptor(const DataSource& _data_source) : data_source(_data_source) {}

		inline DistanceType evalMetric(const T* a, const size_t b_idx, size_t size,
			DistanceType worst_dist = -1) const {
			DistanceType result = DistanceType();
			const T* last = a + size;
			const T* lastgroup = last - 3;
			size_t d = 0;

			/* Process 4 items with each loop for efficiency. */
			while (a < lastgroup) {
				const DistanceType diff0 = torusClamp(a[0] - data_source.kdtree_get_pt(b_idx, d++));
				const DistanceType diff1 = torusClamp(a[1] - data_source.kdtree_get_pt(b_idx, d++));
				const DistanceType diff2 = torusClamp(a[2] - data_source.kdtree_get_pt(b_idx, d++));
				const DistanceType diff3 = torusClamp(a[3] - data_source.kdtree_get_pt(b_idx, d++));
				result += diff0 * diff0 + diff1 * diff1 + diff2 * diff2 + diff3 * diff3;
				a += 4;
				if ((worst_dist > 0) && (result > worst_dist)) {
					return result;
				}
			}
			/* Process last 0-3 components.  Not needed for standard vector lengths. */
			while (a < last) {
				const DistanceType diff0 = torusClamp(*a++ - data_source.kdtree_get_pt(b_idx, d++));
				result += diff0 * diff0;
			}
			return result;
		}

		DistanceType torusClamp(DistanceType ab) const {
			ab += std::floor(abs(ab)) / 2 * 2 + 4;
			ab -= std::floor(ab) / 2 * 2;
			if (ab > 1) ab = 2 - ab;
			return ab;
		}

		/** Note: this assumes that input angles are already in the range [-pi,pi] */
		template <typename U, typename V>
		inline DistanceType accum_dist(const U a, const V b, const size_t) const {
			DistanceType result = DistanceType();
			result = b - a;
			// to [-1,1]^3
			result = torusClamp(result);
			return result * result;
		}
	};

	/** Metaprogramming helper traits class for the SO3_InnerProdQuat metric */
	struct metric_T3 : public Metric {
		template <class T, class DataSource> struct traits {
			typedef T3_Adaptor<T, DataSource> distance_t;
		};
	};

	template <typename T, int N = 3>
	struct PointCloud {
		struct Point { T x[N]; };

		using coord_t = T;  //!< The type of each coordinate

		std::vector<Point> pts;

		// Must return the number of data points
		inline size_t kdtree_get_point_count() const { return pts.size(); }

		// Returns the dim'th component of the idx'th point in the class:
		// Since this is inlined and the "dim" argument is typically an immediate
		// value, the
		//  "if/else's" are actually solved at compile time.
		inline T kdtree_get_pt(const size_t idx, const size_t dim) const
		{
			if (dim < N)
				return pts[idx].x[dim];
			else
				return pts[idx].x[N - 1];
		}

		// Optional bounding-box computation: return false to default to a standard
		// bbox computation loop.
		//   Return true if the BBOX was already computed by the class and returned
		//   in "bb" so it can be avoided to redo it again. Look at bb.size() to
		//   find out the expected dimensionality (e.g. 2 or 3 for point clouds)
		template <class BBOX>
		bool kdtree_get_bbox(BBOX& /* bb */) const {
			return false;
		}
	};

	template<typename T>
	Eigen::Matrix<T, 4, 1> embedTorusPoint(const T* p_data) {
		Eigen::Vector<T, 3> pcos, psin;
		//const T r1 = 1, r2 = 1, r3 = 1; // manifold is self intersected on R^4
		const T r1 = 8, r2 = 3, r3 = 1;
		for (int j = 0; j < 3; j++) {
			pcos[j] = cos(p_data[j] * M_PI); psin[j] = sin(p_data[j] * M_PI);
		}
		Eigen::Matrix<T, 4, 1> pemb;
		pemb[0] = pcos[0] * (r1 + pcos[1] * (r2 + r3 * pcos[2]));
		pemb[1] = psin[0] * (r1 + pcos[1] * (r2 + r3 * pcos[2]));
		pemb[2] =                 psin[1] * (r2 + r3 * pcos[2]);
		pemb[3] =                           (     r3 * psin[2]);
		return pemb;
	}

	template<typename T>
	PointCloud<T, 4> embedTorusPoints(const T* p_data, int n_p) {
		PointCloud<T, 4> pclouds;
		for (int i = 0; i < n_p * 3; i += 3) {
			Eigen::Matrix<T, 4, 1> p = embedTorusPoint(&p_data[i]);
			pclouds.pts.push_back(typename PointCloud<T, 4>::Point({ { p[0],p[1],p[2],p[3] } }));
		}
		return pclouds;
	}
};
#endif


template<typename Scalar>
bool operator<(const Eigen::Matrix<Scalar, 1, 1>& v_left, const Eigen::Matrix<Scalar, 1, 1>& v_right) {
	return v_left[0] < v_right[0];
}

template<typename Scalar, int N, std::enable_if_t< (N > 1), int> = 0 >
bool operator<(const Eigen::Matrix<Scalar, N, 1>& v_left, const Eigen::Matrix<Scalar, N, 1>& v_right) {
	bool les = false;
	if (v_left[N - 1] < v_right[N - 1]) return true;
	else if (v_left[N - 1] > v_right[N - 1])  return false;
	else {
		Eigen::Matrix<Scalar, N - 1, 1> v_left_head = v_left.template block<N - 1, 1>(0, 0);
		Eigen::Matrix<Scalar, N - 1, 1> v_right_head = v_right.template block<N - 1, 1>(0, 0);
		return v_left_head < v_right_head;
	}
}

namespace msf {

	template<typename T>
	T snapSmallToZero(T v, double epsi) {
		if (abs(v) < epsi) { return 0; }
		else return v;
	}

	//template<typename Vec3d>
	std::pair<std::vector<int>, std::vector<VertexFlag>> removeDupDof(std::vector<Eigen::Vector3d>& vlist, std::vector<std::array<int, 4>>& elist,
		std::unordered_set<std::array<int, 3>>& boundaryFaces, std::map<int, int>& vMuxMap, std::map<int, int>& vnew2old);

	struct torus_kd_tree_t {
		using kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor <
			nanoflann::L2_Adaptor<Real, nanoflann::PointCloud<Real, 4>>,
			nanoflann::PointCloud<Real, 4>, 4, size_t>;
		torus_kd_tree_t(Real x_t, Real y_t, Real z_t, Real eps = 1e-12) : _tol(eps), _period(x_t, y_t, z_t) {}
		void build(const Eigen::MatrixX3<Real>& pmat);
		int query(const Eigen::Vector3<Real>& p_torus);
	private:
		Real _tol;
		Eigen::Vector3<Real> _period;
		std::unique_ptr<kd_tree_t> _tree;
		//std::vector<Eigen::Vector3d> torusParamlist;
		nanoflann::PointCloud<msf::Real, 4> cloud;
	};
};

// counter-clock-wise vertex numbering
inline std::vector<std::array<int, 4>> splitBox2Tets(std::array<int, 8> bxCorner)
{
	std::vector<std::array<int, 4>> tets;
	tets.push_back({bxCorner[0], bxCorner[1], bxCorner[3], bxCorner[4]});
	tets.push_back({bxCorner[1], bxCorner[3], bxCorner[4], bxCorner[7]});
	tets.push_back({bxCorner[1], bxCorner[2], bxCorner[3], bxCorner[7]});
	tets.push_back({bxCorner[1], bxCorner[4], bxCorner[5], bxCorner[7]});
	tets.push_back({bxCorner[1], bxCorner[2], bxCorner[5], bxCorner[7]});
	tets.push_back({bxCorner[2], bxCorner[6], bxCorner[5], bxCorner[7]});
	return tets;
}

std::pair<std::vector<Eigen::Vector3d>, std::vector<std::array<int, 4>>> sampleTetmesh(int id);

#pragma pop_macro("_USE_MATH_DEFINES")
