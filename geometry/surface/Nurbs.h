#pragma once

#include "surface.h"
#include "Eigen/Eigen"
#include <limits>
#include <iostream>
#include <unsupported/Eigen/AutoDiff>

BEGIN_MINSURF_NAMESPACE

template<typename P, typename Scalar, typename T, int Order, int Level = 0>
struct DeBoor {
	static void eval(Scalar s, const T* _knots, const P* _cpts, int l, P p_k[]) {
		constexpr int n = Order;
		//constexpr int k = Order - Level;
		bool debug = s < 2.53 && s>2.51;
		P p_k_1[n + 1];
		for (int i = 0; i <= n; i++) {
			p_k_1[i] = _cpts[l - n + i]; if (debug) std::cout << p_k_1[i].transpose() << std::endl;
		}
		for (int k = 1; k <= n; k++) {
			for (int j = l - n + k; j < l + 1; j++) {
				Scalar alpha_k_j = (s - _knots[j]) / (_knots[j + n + 1 - k] - _knots[j]);
				if (debug) std::cout << "[" << k << ", " << j << "] = " << alpha_k_j << std::endl;
				int j_0 = j - (l - n + k);
				p_k_1[j_0] = (1 - alpha_k_j) * p_k_1[j_0] + alpha_k_j * p_k_1[j_0 + 1];
			}
		}
		p_k[0] = p_k_1[0];
		if(debug) std::cout << "p = " << p_k[0].transpose() << std::endl;
	}
};

//template<typename P, typename Scalar, typename T, int Order>
//struct DeBoor<P, Scalar, T, Order, Order> {
//	static void eval(Scalar s, const T* _knots, const P* _cpts, int l, P dst[]) {
//		constexpr int n = Order;
//		for (int i = l - n; i < l + 1; i++) {
//			dst[i] = _cpts[i];
//		}
//		return;
//	}
//};

template<typename P, typename Scalar, typename T>
void deBoor(Scalar s, int Order, int Level, const T* _knots, const P* _cpts, int l, P dst[]) {
	int k = Level;
	int n = Order;
	if (k == 0) {
		for (int i = l - n; i < l + 1; i++) {
			dst[i] = _cpts[i];
		}
		return;
	}
	std::vector<P> p_k_1(n - k + 1);
	deBoor<P, Scalar, T>(s, n, k - 1, _knots, _cpts, l, p_k_1.data());
	for (int i = l - n + k; i < l + 1; i++) {
		Scalar alpha_k_i = (s - _knots[i]) / (_knots[i + n + 1 - k] - _knots[i]);
		int ioff = i - (l - n + k);
		dst[i] = (1 - alpha_k_i) * p_k_1[ioff] + alpha_k_i * p_k_1[ioff + 1];
	}
	
}

template<typename Scalar, int Dim, int Order>
class BSpline {
	using P = Eigen::Vector<Scalar, Dim>;
	std::vector<Scalar> _knots;
	std::vector<P> _cpts;

	BSpline(const std::vector<Scalar>& knots, const std::vector<P>& ptlist) 
		: _knots(knots), _cpts(ptlist) { }

	P eval(Scalar t) const {
		int l = 0;
		for (int i = _knots.size() - 1; i >= 1; i--) {
			if (_knots[i - 1] <= t && t < _knots[i]) { l = i - 1; break; }
		}
		P p;
		DeBoor<P, Scalar, Scalar, Order>::eval(t, _knots.data(), _cpts.data(), l, &p);
		return p;
	}

	std::vector<P> diff(const std::vector<Scalar>& t) const {
		BSpline<Scalar, Dim, Order - 1> bdiff = diff();
		std::vector<P> tg(t.size());
		for (int i = 0; i < t.size(); i++) {
			tg[i] = bdiff.eval(t[i]);
		}
		return tg;
	}

	BSpline<Scalar, Dim, Order - 1> diff(void) const {
		std::vector<Scalar> knots(_knots.size() - 1);
		std::vector<Point> cpts(_cpts.size() - 1);
		for (int i = 0; i < knots.size(); i++) {
			knots[i] = _knots[i + 1];
			cpts[i] = Order / (_knots[i + Order + 1] - _knots[i + 1]) * (_cpts[i + 1] - _cpts[i]);
		}
		BSpline<Scalar, Dim, Order - 1> bdiff(knots, cpts);
		return bdiff;
	}
};

//inline void foo(void) {
//
//	Eigen::Vector4d wp;
//	wp.topRows<3>(0);
//}

template<typename Scalar, int Dim>
class BSpline<Scalar, Dim, -1> {
	using P = Eigen::Vector<Scalar, Dim>;
	std::vector<Scalar> _knots;
	std::vector<P> _cpts;
	int Order = -1;

public:
	BSpline(int order, const std::vector<Scalar>& knots, const std::vector<P>& ptlist)
		:Order(order), _knots(knots), _cpts(ptlist) { 
		if (_knots.size() - _cpts.size() != order + 1) {
			printf("Warning : Mismatch Bspline control %d points, %d knots, Order %d\n", _cpts.size(), _knots.size(), order);
		}
	}

	P eval(Scalar t) const {
		int l = 0;
		for (int i = _knots.size() - 1; i >= 1; i--) {
			if (_knots[i - 1] <= t && t < _knots[i]) { l = i - 1; break; }
		}
		P p;
		deBoor<P, Scalar, Scalar>(t, Order, Order, _knots.data(), _cpts.data(), l, &p);
		return p;
	}

	std::vector<P> diff(const std::vector<Scalar>& t) const {
		BSpline<Scalar, Dim, -1> bdiff = diff();
		std::vector<P> tg(t.size());
		for (int i = 0; i < t.size(); i++) {
			tg[i] = bdiff.eval(t[i]);
		}
		return tg;
	}

	BSpline<Scalar, Dim, -1> diff(void) const {
		std::vector<Scalar> knots(_knots.size() - 1);
		std::vector<Point> cpts(_cpts.size() - 1);
		for (int i = 0; i < knots.size(); i++) {
			knots[i] = _knots[i + 1];
			cpts[i] = Order / (_knots[i + Order + 1] - _knots[i + 1]) * (_cpts[i + 1] - _cpts[i]);
		}
		BSpline<Scalar, Dim, -1> bdiff(Order - 1, knots, cpts);
		return bdiff;
	}
};

template<typename Scalar, int Dim, int Order>
class NURBSSpline {
	using WP = Eigen::Vector<Scalar, Dim + 1>;
	using P = Eigen::Vector<Scalar, Dim >;
	std::vector<Scalar> _knots;
	std::vector<WP> _cpts;

	NURBSSpline<Scalar, Dim, Order - 1> diff(void) const {
		std::vector<Scalar> knots(_knots.size() - 1);
		std::vector<WP> cpts(_cpts.size() - 1);
		for (int i = 0; i < knots.size(); i++) {
			knots[i] = _knots[i + 1];
			cpts[i] = Order / (_knots[i + Order + 1] - _knots[i + 1]) * (_cpts[i + 1] - _cpts[i]);
		}
		NURBSSpline<Scalar, Dim, Order - 1> bdiff(knots, cpts);
		return bdiff;
	}

public:
	//using WP = Eigen::Vector<Scalar, Dim + 1>;
	//using P = Eigen::Vector<Scalar, Dim>;

	NURBSSpline(const std::vector<Scalar>& knots, const std::vector<P>& ptlist, const std::vector<Scalar>& weights)
		: _knots(knots) {
		if (ptlist.size() != weights.size()) {
			printf("Warning : Mismatch weights and control points\n");
		}
		if (_knots.size() - ptlist.size() != Order + 1) {
			printf("Warning : Mismatch Bspline control %d points, %d knots, Order %d\n", ptlist.size(), _knots.size(), Order);
		}
		for (int i = 0; i < ptlist.size(); i++) {
			WP wp;
			wp <<  weights[i] * ptlist[i], weights[i];
			_cpts.push_back(wp);
		}
	}

	NURBSSpline(const std::vector<Scalar>& knots, const std::vector<WP>& ptlist)
		: _knots(knots), _cpts(ptlist) { }

	std::pair<Scalar, Scalar> getKnotRange(void) {
		Scalar smin = std::numeric_limits<Scalar>::max();
		Scalar smax = std::numeric_limits<Scalar>::min();
		for (int i = 0; i < _knots.size(); i++) {
			if (smin > _knots[i]) smin = _knots[i];
			if (smax < _knots[i]) smax = _knots[i];
		}
		return { smin,smax };
	}

	P eval(Scalar t) const {
		int l = 0;
		for (int i = _knots.size() - 1; i >= 1; i--) {
			if (_knots[i - 1] <= t && t < _knots[i]) { l = i - 1; break; }
		}
		Eigen::Vector4d  wp;
		DeBoor<WP, Scalar, Scalar, Order>::eval(t, _knots.data(), _cpts.data(), l, &wp);
		return wp.block<3, 1>(0, 0) / wp[3];
	}
	WP eval_h(Scalar t) const {
		int l = 0;
		for (int i = 0; i < _knots.size() + 1; i++) {
			if (_knots[i] <= t && t < _knots[i + 1]) {
				l = i;
				break;
			}
		}
		Eigen::Vector4d  wp;
		DeBoor<WP, Scalar, Scalar, Order>::eval(t, _knots.data(), _cpts.data(), l, &wp);
		return wp;
	}

	std::vector<P> diff(const std::vector<Scalar>& t) const {
		NURBSSpline<Scalar, Dim, Order - 1> bdiff = diff();
		std::vector<P> tg(t.size());
		for (int i = 0; i < t.size(); i++) {
			WP p_h = eval_h(t[i]);
			WP tg_h = bdiff.eval_h(t[i]);
			tg[i] = tg_h.block<3, 1>(0, 0) / p_h[3] - p_h.block<3, 1>(0, 0) / pow(p_h[3], 2) * tg_h[3];
		}
		return tg;
	}
};


template<typename Scalar, int uOrder, int vOrder>
class NURBSSurface {
	Scalar _UDomain[2];
	Scalar _VDomain[2];

	std::vector<Scalar> _uknots, _vknots;

	using ADScalar = Eigen::AutoDiffScalar<Eigen::Vector2<Scalar>>;
	int _Ni, _Nj;

	using P = Eigen::Vector<Scalar, 3>;
	using WP = Eigen::Vector<Scalar, 4>;
	using V = P;


	std::vector<P> _Qpts;
	std::vector<P> _Qw;

private:
	std::tuple<int, int> locate(Scalar u, Scalar v) const {
		int ul = -1, vl = -1;
		for (int i = _uknots.size() - 1; i >= 1; i--) {
			if (_uknots[i - 1] <= u && u < _uknots[i]) { ul = i - 1; break; }
		}
		for (int i = _uknots.size() - 1; i >= 1; i--) {
			if (_uknots[i - 1] <= v && v < _uknots[i]) { vl = i - 1; break; }
		}
		return { ul,vl };
	}

public:
	P eval(Scalar u, Scalar v) const {
		// de Boor on u	
		auto [i, j] = locate(u, v);
		std::vector<WP> vpts(vOrder + 1);
		for (int vj = 0; vj < vOrder + 1; vj++) {
			DeBoor<WP, Scalar, Scalar, uOrder>::eval(u, _uknots.data(), _Qpts.data() + j * (_Ni - uOrder), i, &vpts[vj]);
		}
		WP result;
		DeBoor<WP, Scalar, Scalar, vOrder>::eval(v, _vknots.data() + i - vOrder, vpts.data(), 0, &result);
		return result.block<3, 1>(0, 0) / result[3];
	}

	// |qW| = |qpts| = |vknots| - order 
	NURBSSurface(const std::vector<Scalar>& uknots, const std::vector<Scalar>& vknots, const std::vector<Scalar>& qW, const std::vector<P>& qpts)
		: _uknots(uknots), _vknots(vknots), _Qw(qW), _Qpts(qpts) { 
		if (qW.size() != qpts.size()) {
			printf("Warning : Mismatch number of weights and control points\n");
		}
		if ((uknots.size() - (uOrder + 1)) * (vknots.size() - (vOrder + 1)) != qpts.size()) {
			printf("Warning : Mismatch surface knots and control points\n");
		}
	}

	auto tangent(Scalar u_, Scalar v_) const {
		auto [i, j] = locate(u_, v_);
		static thread_local std::vector<Eigen::Vector4<ADScalar>> vpts(vOrder + 1);
		vpts.resize(vOrder + 1);
		
		ADScalar u(u_, 2, 0), v(v_, 2, 1);
		for (int vj = 0; vj < vOrder + 1; vj++) {
			DeBoor<Eigen::Vector4<ADScalar>, ADScalar, Scalar, uOrder>::eval(u, _uknots.data(), _Qpts.data() + j * (_Ni - uOrder), i, &vpts[vj]);
		}
		Eigen::Vector4<ADScalar> result;
		DeBoor<Eigen::Vector4<ADScalar>, ADScalar, Scalar, vOrder>::eval(v, _vknots.data() + i - vOrder, vpts.data(), 0, &result);
		result.block<3, 1>(0, 0) /= result[3];
		Eigen::Matrix<Scalar, 3, 2> Tg;
		for (int i = 0; i < 3; i++) {
			Tg.rows(i) = result[i].derivatives().transpose();
		}
		return Tg;
	}

	auto metric(Scalar u_, Scalar v_) const {
		auto tg = tangent(u_, v_);
		Eigen::Matrix<Scalar, 2, 2> g;
		g(0, 0) = tg.col(0).dot(tg.col(0));
		g(0, 1) = tg.col(0).dot(tg.col(1));
		g(1, 1) = tg.col(1).dot(tg.col(1));
		g(1, 0) = g(0, 1);
		return g;
	}

	
};


END_MINSURF_NAMESPACE
