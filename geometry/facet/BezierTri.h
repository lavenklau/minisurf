#pragma once

#include "Facet.h"
#include <unsupported/Eigen/AutoDiff>
#include "facet/Quadrature.h"


BEGIN_MINSURF_NAMESPACE

inline int pid(int i, int j, int k) {
	int d = i + j + k;
	return (d + 1) * i - ((i * (i - 1)) / 2) + j;
}

inline std::array<int, 3> dip(int order, int k) {
	int at_i = 0, at_j = 0;
	for (int i = 0; i <= order; i++) {
		if (k > order - i) at_i = i + 1;
		else {
			at_j = k;
			break;
		}
		k -= order - i + 1;
	}
	return { at_i,at_j, order - at_i - at_j };
}

namespace detail {
	template<typename Scalar, int OrderEnd, int Order_, std::enable_if_t<OrderEnd == Order_, int> x = 0>
	inline Eigen::Vector3<Scalar> deCas_impl(const Point p[], int i, int j, int k, Scalar u, Scalar v) {
		return p[pid(i, j, k)];
	}
	template<typename Scalar, int OrderEnd, int Order_, std::enable_if_t<OrderEnd != Order_, int> x = 0>
	inline auto deCas_impl(const Point p[], int i, int j, int k, Scalar u, Scalar v) {
		return u * deCas_impl<Scalar, OrderEnd, Order_ + 1>(p, i + 1, j, k, u, v) + v * deCas_impl<Scalar, OrderEnd, Order_ + 1>(p, i, j + 1, k, u, v)
			+ (1 - u - v) * deCas_impl<Scalar, OrderEnd, Order_ + 1>(p, i, j, k + 1, u, v);
	}
	template<typename Scalar, int OrderEnd, int Order_, std::enable_if_t<OrderEnd == Order_, int> x = 0>
	inline auto d_deCas_impl(const Point p[], int i, int j, int k, Scalar u, Scalar v) {
		using V = Eigen::Vector3<Scalar>;
		return std::pair<V, V>{ V(0, 0, 0), V(0, 0, 0) };
	}
	template<typename Scalar, int OrderEnd, int Order_, std::enable_if_t<OrderEnd != Order_, int> x = 0>
	inline auto d_deCas_impl(const Point p[], int i, int j, int k, Scalar u, Scalar v) {
		std::pair<Eigen::Vector3<Scalar>, Eigen::Vector3<Scalar>> dp;
		dp.first.fill(0); dp.second.fill(0);
		std::pair<Eigen::Vector3<Scalar>, Eigen::Vector3<Scalar>> duv[3] = {
			d_deCas_impl<Scalar, OrderEnd, Order_ + 1>(p, i + 1, j, k, u, v),
			d_deCas_impl<Scalar, OrderEnd, Order_ + 1>(p, i , j + 1, k, u, v),
			d_deCas_impl<Scalar, OrderEnd, Order_ + 1>(p, i , j, k + 1, u, v) 
		};
		Eigen::Vector3<Scalar> b[3] = {
			deCas_impl<Scalar, OrderEnd, Order_ + 1>(p, i + 1, j, k, u, v),
			deCas_impl<Scalar, OrderEnd, Order_ + 1>(p, i , j + 1, k, u, v), 
			deCas_impl<Scalar, OrderEnd, Order_ + 1>(p, i , j, k + 1, u, v) };
		dp.first = b[0] + u * duv[0].first + v * duv[1].first - b[2] + (1 - u - v) * duv[2].first;
		dp.second = u * duv[0].second + b[1] + v * duv[1].second - b[2] + (1 - u - v) * duv[2].second;
		return dp;
	}

	template<typename T, typename Scalar>
	T deCas_impl(const T p[], int level, int i, int j, int k, Scalar u, Scalar v) {
		if (level == 0) return p[pid(i, j, k)];
		return u * deCas_impl(p, level - 1, i + 1, j, k) + v * deCas_impl(p, level - 1, i, j + 1, k) + (1 - u - v) * deCas_impl(p, level - 1, i, j, k + 1);
	}

	template<typename T, typename Scalar>
	T d_deCas_impl(const T p[], int level, int i, int j, int k, Scalar u, Scalar v, Scalar du, Scalar dv) {
		if (level == 0) return p[pid(i, j, k)];
		if (i == 0 && j == 0 && k == 0) {
			return level * (du * deCas_impl(p, level - 1, i + 1, j, k) + dv * deCas_impl(p, level - 1, i, j + 1, k) + (-du - dv) * deCas_impl(p, level - 1, i, j, k + 1));
		} else {
			return u * deCas_impl(p, level - 1, i + 1, j, k) + v * deCas_impl(p, level - 1, i, j + 1, k) + (1 - u - v) * deCas_impl(p, level - 1, i, j, k + 1);
		}
	}
};

template<typename Scalar, int Order>
class BezierTri
	: public Facet
{
	using Vector = Point;
public:
	static constexpr int n_cpts = (Order + 1) * (Order + 2) / 2;
private:
	Point _cpt[n_cpts];

public:
	auto& getPoint(int i, int j, int k) { return _cpt[pid(i, j, k)]; }
private:

	static constexpr int _fac(int k) {
		if (k == 0 || k == 1) {
			return 1; 
		} else {
			return k * _fac(k - 1);
		}
	}

	template<int Order_ = 0>
	static Point deCas(const Point p[], int i, int j, int k, Scalar u, Scalar v) {
		return detail::deCas_impl<Scalar, Order, Order_>(p, i, j, k, u, v);
	}
	template<int Order_ = 0>
	static std::pair<Point, Point> d_deCas(const Point p[], int i, int j, int k, Scalar u, Scalar v) {
		return detail::d_deCas_impl<Scalar, Order, Order_>(p, i, j, k, u, v);
	}

	// r_uu, r_uv, r_v
	template<int Order_ = 0>
	static std::tuple<Point, Point, Point> d2_deCas(const Point p[], int i, int j, int k, Scalar u, Scalar v) {
		using ADScalar = Eigen::AutoDiffScalar<Eigen::Vector2<Scalar>>;
		ADScalar u_ad = u, v_ad = v;
		u_ad.derivatives() = Eigen::Vector2<Scalar>::Unit(0);
		v_ad.derivatives() = Eigen::Vector2<Scalar>::Unit(1);
		auto dr = detail::d_deCas_impl<ADScalar, Order, Order_>(p, i, j, k, u_ad, v_ad);
		Eigen::Vector3<Scalar> r_uu, r_uv, r_vv;
		for (int i = 0; i < 3; i++) {
			r_uu[i] = dr.first[i].derivatives()[0];
			r_uv[i] = dr.first[i].derivatives()[1];
			r_vv[i] = dr.second[i].derivatives()[1];
		}
		return { r_uu, r_uv, r_vv };
	}


	static Scalar _B(int i, int j, int k, Scalar u, Scalar v) {
		if (Order != i + j + k)
			printf("issue with bezier poly degree: %d = %d + %d + %d\n", Order, i, j, k);
		if (i < 0 || j < 0 || k < 0)
			return 0.0;
		double b = _fac(Order);
		b *= pow(u, i);
		b *= pow(v, j);
		b *= pow(1.0 - u - v, k);
		b /= _fac(i);
		b /= _fac(j);
		b /= _fac(k);
		return b;
	}
	static constexpr Scalar small_eps = std::is_same_v<Scalar, float> ? 1e-5f : 1e-14;
	static std::tuple<Scalar, Scalar> _dB(int i, int j, int k, Scalar u, Scalar v) {
		if (Order != i + j + k)
			printf("issue with bezier poly degree: %d = %d + %d + %d\n", Order, i, j, k);
		if (i < 0 || j < 0 || k < 0)
			return 0.0;
		Scalar b = _fac(Order) / _fac(i) / _fac(j) / _fac(k);
		Scalar bu, bv;
		Scalar u_val = pow(u, i);
		Scalar du_val = i * pow(u, i - 1);
		b *= pow(u, i) * pow(v, j) * pow(1.0 - u - v, k);
		bu = (std::abs(u) < small_eps ? 0 : b * i / u) - (std::abs(1 - u - v) < small_eps ? 0 : b * k / (1 - u - v));
		bv = (std::abs(v) < small_eps ? 0 : b * j / v) - (std::abs(1 - u - v) < small_eps ? 0 : b * k / (1 - u - v));
		return { bu, bv };
	}

public:
	BezierTri(const Point p[]) {
		for (int i = 0; i < n_cpts; i++) { _cpt[i] = p[i]; }
	}
	BezierTri(void) = default;
	//BezierTri(const BezierTri<Scalar, Real>&) = default;
	Point eval(Scalar u, Scalar v) const {
		Point p(0, 0, 0);
#if 0
		for (int i = 0; i <= Order; i++) {
			for (int j = 0; j <= Order - i; j++) {
				int k = Order - i - j;
				int cpid = pid(i, j, k);
				p += _B(i, j, k, u, v) * _cpt[cpid];
			}
		}
#else
		p = deCas(_cpt, 0, 0, 0, u, v);
#endif
		return p;
	}
	std::tuple<Vector, Vector> TgEval(Scalar u, Scalar v) const {
#if 0
		Vector p_u(0, 0, 0), p_v(0, 0, 0);
		for (int i = 0; i <= Order; i++) {
			for (int j = 0; j <= Order - i; j++) {
				int k = Order - i - j;
				int cpid = pid(i, j, k);
				auto [du, dv] = _dB(i, j, k, u, v);
				p_u += du * _cpt[cpid];
				p_v += dv * _cpt[cpid];
			}
		}
#else
		auto [p_u, p_v] = d_deCas(_cpt, 0, 0, 0, u, v);
#endif
		return { p_u, p_v };
	}

	auto Hessian(Scalar u, Scalar v) const {
		return d2_deCas(_cpt, 0, 0, 0, u, v);
	}

	Vector normal(Scalar u, Scalar v) const {
		auto [p_u, p_v] = TgEval(u, v);
		Vector n = p_u.cross(p_v).normalized();
		return n;
	}

	Eigen::Matrix3<Scalar> frame(Scalar u, Scalar v) const {
		auto [p_u, p_v] = TgEval(u, v);
		Vector n = p_u.cross(p_v).normalized();
		Eigen::Matrix3<Scalar> fr;
		fr << p_u, p_v, n;
		return fr;
	}

	Real Quadrature(int order, std::function<Real(const Point& pt, const Vector& nrm)> func) {
		static Real UnitTriArea = 0.5;
		const Real (*q_uv)[2];
		const Real (*q_w);
		int q_size;
		std::tie(q_uv, q_w, q_size) = getTriQuadrature(order);
		Point p;
		Real s = 0;
		for (int i = 0; i < q_size; i++) {
			Real w = q_w[i];
			Real u = q_uv[i][0], v = q_uv[i][1];
			auto [tu, tv] = TgEval(u, v);
			Eigen::Matrix<Real, 3, 2> Tg; Tg << tu, tv;
			Eigen::Matrix2<Real> g = Tg.transpose() * Tg;
			Real J = (std::sqrt)(g.determinant());
			Vector nrm = normal(u, v);
			s += w * func(p, nrm) * J;
		}
		s *= UnitTriArea;
		return s;
	}
};

template<typename Scalar>
class BezierTri<Scalar, -1>
	: public Facet
{
	using Vector = Point;
private:
	static std::map<int, std::vector<Scalar>> _bcoeff;
	static std::map<std::pair<std::string, int>, std::vector<Eigen::MatrixX2<Real>>> _Tbuffer;
	static std::map<std::pair<std::string, int>, std::vector<Eigen::MatrixX3<Real>>> _Hbuffer;
	std::vector<Scalar>* _b;
	std::string qBuffer;
protected:
	std::vector<Point> _cpt;
	int _order;
public:
	auto& getPoint(int i, int j, int k) { return _cpt[pid(i, j, k)]; }
	auto& getPoint(int cid) { return _cpt[cid]; }
	auto& getPoint(int cid) const { return _cpt[cid]; }
	auto& getPoints(void) const { return _cpt; }
	int n_cpts(void) const { return(_order + 1) * (_order + 2) / 2; }
	static int n_cpts(int order) { return(order + 1) * (order + 2) / 2; }
	int getBezierOrder(void) const { return _order; }
	void periodSync(const BezierTri<Scalar, -1>& tri2, Scalar period = 2) {
		Eigen::Vector3<Real> d = getPoint(0) - tri2.getPoint(0);
		for (int i = 0; i < 3; i++) {
			if (d[i] < -period / 2.) { d[i] = period; }
			else if (d[i] > period / 2.) { d[i] = -period; }
			else { d[i] = 0; }
		}
		for (int i = 0; i < n_cpts(); i++) { _cpt[i] += d; }
	}
private:
	static constexpr int _fac(int k) {
		if (k == 0 || k == 1) {
			return 1;
		}
		else {
			return k * _fac(k - 1);
		}
	}
	static void initBernCoeff(int order) {
		if (!_bcoeff.count(order)) {
			_bcoeff[order].resize(n_cpts(order));
			std::vector<size_t> facs(order + 1);
			for (int i = 0; i <= order; i++) facs[i] = _fac(i);
			for (int i = 0; i <= order; i++) {
				for (int j = 0; j <= order - i; j++) {
					int k = order - i - j;
					size_t bi = facs[order] / facs[i] / facs[j] / facs[k];
					_bcoeff[order][pid(i, j, k)] = bi;
				}
			}
		}
	}
	void initBernCoeff(void) {
		initBernCoeff(_order);
		_b = &_bcoeff[_order];
	}

	static Scalar _duvw(int du, int dv, int i, int j, int k, Scalar u, Scalar v, Scalar w) {
		if (du == 0 && dv == 0) {
			Scalar ui = i <= 0 ? 1 : std::pow(u, i);
			Scalar vj = j <= 0 ? 1 : std::pow(v, j);
			Scalar wk = k <= 0 ? 1 : std::pow(w, k);
			return ui * vj * wk;
		}
		else if (du != 0) {
			return i * _duvw(du - 1, dv, i - 1, j, k, u, v, w) - k * _duvw(du - 1, dv, i, j, k - 1, u, v, w);
		}
		else if (dv != 0) {
			return j * _duvw(du, dv - 1, i, j - 1, k, u, v, w) - k * _duvw(du, dv - 1, i, j, k - 1, u, v, w);
		}
		return 0;
	}

	static Scalar _B(int i, int j, int k, Scalar u, Scalar v, const std::vector<Scalar>& bc) {
		Scalar b = bc[pid(i, j, k)];
		b *= pow(u, i); b *= pow(v, j); b *= pow(1.0 - u - v, k);
		return b;
	}
	static constexpr Scalar small_eps = std::is_same_v<Scalar, float> ? 1e-5f : 1e-14;
	static std::tuple<Scalar, Scalar> _dB(int i, int j, int k, Scalar u, Scalar v, const std::vector<Scalar>& bc) {
		Scalar b = bc[pid(i, j, k)];
		Scalar bu, bv;
		//Scalar u_val = pow(u, i);
		//Scalar du_val = i * pow(u, i - 1);
		//b *= pow(u, i) * pow(v, j) * pow(1.0 - u - v, k);
		//bu = (std::abs(u) < small_eps ? 0 : b * i / u) - (std::abs(1 - u - v) < small_eps ? 0 : b * k / (1 - u - v));
		//bv = (std::abs(v) < small_eps ? 0 : b * j / v) - (std::abs(1 - u - v) < small_eps ? 0 : b * k / (1 - u - v));
		bu = b * _duvw(1, 0, i, j, k, u, v, 1 - u - v);
		bv = b * _duvw(0, 1, i, j, k, u, v, 1 - u - v);
		return { bu, bv };
	}

	static Eigen::Matrix2<Scalar> _d2B(int i, int j, int k, Scalar u, Scalar v, const std::vector<Scalar>& bc) {
		Eigen::Matrix2<Scalar> H;
		Scalar b = bc[pid(i, j, k)];
		H(0, 0) = b * _duvw(2, 0, i, j, k, u, v, 1 - u - v);
		H(0, 1) = b * _duvw(1, 1, i, j, k, u, v, 1 - u - v);
		H(1, 1) = b * _duvw(0, 2, i, j, k, u, v, 1 - u - v);
		H(1, 0) = H(0, 1);
		return H;
	}

public:
	static void bufferSamples(const std::string& name, int order, const std::vector<Eigen::Vector2<Real>>& uvlist) {
		auto key = std::pair<std::string, int>(name, order);
		if (_Tbuffer.count(key) && _Hbuffer.count(key)) return;
		auto Hbasis = HessianBasis(order, uvlist);
		auto Tgbasis = TgBasis(order, uvlist);
		_Tbuffer[key] = Tgbasis;
		_Hbuffer[key] = Hbasis;
	}

	static auto& getTgBuffer(const std::string& name, int order) { return _Tbuffer.at(std::pair<std::string, int>(name, order)); }

	auto& getTgBuffer(void) const { return getTgBuffer(qBuffer, _order); }

	static auto& getHbuffer(const std::string& name, int order) { return _Hbuffer.at(std::pair<std::string, int>(name, order)); }

	auto& getHbuffer(void) const { return getHbuffer(qBuffer, _order); }

	void useBuffer(const std::string& qname) {
		if (_Tbuffer.count(std::pair<std::string, int>(qname, _order)) && _Hbuffer.count(std::pair<std::string, int>(qname, _order)))
			qBuffer = qname;
		else {
			throw std::runtime_error("unknown buffer");
		}
	}

	static auto& getBernCoeff(int order) {
		initBernCoeff(order);
		return _bcoeff[order];
	}

	BezierTri(const Point p[], int order)
		: _order(order)
	{
		_cpt.resize(n_cpts());
		for (int i = 0; i < n_cpts(); i++) { _cpt[i] = p[i]; }
		initBernCoeff();
	}
	BezierTri(void) = default;
	//BezierTri(const BezierTri<Scalar, Real>&) = default;
	Point eval(Scalar u, Scalar v) const {
		Point p(0, 0, 0);
		for (int i = 0; i <= _order; i++) {
			for (int j = 0; j <= _order - i; j++) {
				int k = _order - i - j;
				int cpid = pid(i, j, k);
				p += _B(i, j, k, u, v, *_b) * _cpt[cpid];
			}
		}
		return p;
	}
	std::tuple<Vector, Vector> TgEval(Scalar u, Scalar v) const {
		Vector p_u(0, 0, 0), p_v(0, 0, 0);
		for (int i = 0; i <= _order; i++) {
			for (int j = 0; j <= _order - i; j++) {
				int k = _order - i - j;
				int cpid = pid(i, j, k);
				auto [du, dv] = _dB(i, j, k, u, v, *_b);
				p_u += du * _cpt[cpid];
				p_v += dv * _cpt[cpid];
			}
		}
		return { p_u, p_v };
	}

	std::tuple<Vector, Vector> TgEval(int qid) const {
		auto& qb = getTgBuffer();
		Vector p_u(0, 0, 0), p_v(0, 0, 0);
		for (int i = 0; i < n_cpts(); i++) {
			p_u += _cpt[i] * qb[qid](i, 0);
			p_v += _cpt[i] * qb[qid](i, 1);
		}
		return { p_u, p_v };
	}

	static std::vector<Eigen::MatrixX2<Scalar>> TgBasis(int order, const std::vector<Eigen::Vector2<Scalar>>& ulist) {
		std::vector<Eigen::MatrixX2<Scalar>> basislist;
		for (int n = 0; n < ulist.size(); n++) {
			basislist.emplace_back(TgBasis(order, ulist[n]));
		}
		return basislist;
	}

	static Eigen::MatrixX2<Scalar> TgBasis(int order, const Eigen::Vector2<Scalar>& ulist) {
		int np = n_cpts(order);
		Eigen::MatrixX2<Scalar> gradBn(np, 2);
		for (int i = 0; i <= order; i++) {
			for (int j = 0; j <= order - i; j++) {
				int k = order - i - j;
				int cpid = pid(i, j, k);
				auto [du, dv] = _dB(i, j, k, ulist[0], ulist[1], getBernCoeff(order));
				gradBn(cpid, 0) = du; gradBn(cpid, 1) = dv;
			}
		}
		return gradBn;
	}


	std::vector<Eigen::MatrixX2<Scalar>> TgBasis(const std::vector<Eigen::Vector2<Scalar>>& ulist) const {
		return TgBasis(_order, ulist);
	}

	auto Hessian(Scalar u, Scalar v) const {
		using V = Eigen::Vector3<Scalar>;
		V r_uu(0, 0, 0), r_uv(0, 0, 0), r_vv(0, 0, 0);
		for (int i = 0; i <= _order; i++) {
			for (int j = 0; j <= _order - i; j++) {
				int k = _order - i - j;
				int cid = pid(i, j, k);
				auto Hijk = _d2B(i, j, k, u, v, *_b);
				r_uu += _cpt[cid] * Hijk(0, 0);
				r_uv += _cpt[cid] * Hijk(0, 1);
				r_vv += _cpt[cid] * Hijk(1, 1);
			}
		}
		return std::tuple<V, V, V>(r_uu, r_uv, r_vv);
	}

	static auto HessianBasis(int order, const std::vector<Eigen::Vector2<Scalar>>& ulist) {
		std::vector<Eigen::MatrixX3<Scalar>> basislist;
		for (int n = 0; n < ulist.size(); n++) {
			basislist.emplace_back(HessianBasis(order, ulist[n]));
		}
		return basislist;
	}

	static auto HessianBasis(int order, const Eigen::Vector2<Scalar>& ulist) {
		Eigen::MatrixX3<Scalar> hij(n_cpts(order), 3);
		for (int i = 0; i <= order; i++) {
			for (int j = 0; j <= order - i; j++) {
				int k = order - i - j; int cid = pid(i, j, k);
				auto H = _d2B(i, j, k, ulist[0], ulist[1], getBernCoeff(order));
				hij(cid, 0) = H(0, 0); hij(cid, 1) = H(0, 1); hij(cid, 2) = H(1, 1);
			}
		}
		return hij;
	}


	auto HessianBasis(const std::vector<Eigen::Vector2<Scalar>>& ulist) {
		return HessianBasis(_order, ulist); 
	}

	Eigen::Matrix2<Scalar> second_form(Scalar u, Scalar v) const {
		auto [r_uu, r_uv, r_vv] = Hessian(u, v);
		auto n = normal(u, v);
		Eigen::Matrix2<Scalar> bij;
		bij(0, 0) = r_uu.dot(n);
		bij(0, 1) = r_uv.dot(n);
		bij(1, 1) = r_vv.dot(n);
		bij(1, 0) = bij(0, 1);
		return bij;
	}

	Eigen::Matrix2<Scalar> first_form(Scalar u, Scalar v) const {
		auto [gi, gj] = TgEval(u, v);
		Eigen::Matrix2<Scalar> gij;
		gij << gi.dot(gi), gi.dot(gj), gj.dot(gi), gj.dot(gj);
		return gij;
	}

	Vector normal(Scalar u, Scalar v) const {
		auto [p_u, p_v] = TgEval(u, v);
		Vector n = p_u.cross(p_v).normalized();
		return n;
	}

	Eigen::Matrix3<Scalar> frame(Scalar u, Scalar v) const {
		auto [p_u, p_v] = TgEval(u, v);
		Vector n = p_u.cross(p_v).normalized();
		Eigen::Matrix3<Scalar> fr;
		fr << p_u, p_v, n;
		return fr;
	}

	Real QuadratureOnSurface(int order, std::function<Real(const Point& pt, const Vector& nrm)> func) {
		static Real UnitTriArea = 0.5;
		const Real (*q_uv)[2];
		const Real (*q_w);
		int q_size;
		std::tie(q_uv, q_w, q_size) = getTriQuadrature(order);
		Point p;
		Real s = 0;
		for (int i = 0; i < q_size; i++) {
			Real w = q_w[i];
			Real u = q_uv[i][0], v = q_uv[i][1];
			auto [tu, tv] = TgEval(u, v);
			Eigen::Matrix<Real, 3, 2> Tg; Tg << tu, tv;
			Eigen::Matrix2<Real> g = Tg.transpose() * Tg;
			Real J = (std::sqrt)(g.determinant());
			Vector nrm = normal(u, v);
			s += w * func(p, nrm) * J;
		}
		s *= UnitTriArea;
		return s;
	}
};

template<typename Scalar>
SELECTANY std::map<int, std::vector<Scalar>> msf::BezierTri<Scalar, -1>::_bcoeff;

template<typename Scalar>
SELECTANY std::map<std::pair<std::string, int>, std::vector<Eigen::MatrixX2<Real>>> msf::BezierTri<Scalar, -1>::_Tbuffer;
template<typename Scalar>
SELECTANY std::map<std::pair<std::string, int>, std::vector<Eigen::MatrixX3<Real>>> msf::BezierTri<Scalar, -1>::_Hbuffer;

END_MINSURF_NAMESPACE
