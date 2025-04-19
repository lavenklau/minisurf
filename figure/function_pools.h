#pragma once

#include <iostream>
#include <Eigen/Eigen>
#include <string>
#include <tuple>

namespace msf {

	namespace vfn {
		Eigen::Vector3d V1(const Eigen::Vector3d&);
		Eigen::Matrix3d dV1(const Eigen::Vector3d&);
	}

	namespace sfn {
		double f1(const Eigen::Vector3d&);
		Eigen::Vector3d df1(const Eigen::Vector3d&);

		double sphere(const Eigen::Vector3d&);
		Eigen::Vector3d dsphere(const Eigen::Vector3d&);
		Eigen::Matrix3d d2sphere(const Eigen::Vector3d&);

		double cylinder(const Eigen::Vector3d& );
		Eigen::Vector3d dcylinder(const Eigen::Vector3d&);
		Eigen::Matrix3d d2cylinder(const Eigen::Vector3d&);

		double catenoid(const Eigen::Vector3d&);
		Eigen::Vector3d dcatenoid(const Eigen::Vector3d&);
		Eigen::Matrix3d d2catenoid(const Eigen::Vector3d&);

		
		std::tuple<decltype(&sphere), decltype(&dsphere), decltype(&d2sphere)> get(std::string fname);

		namespace period {
			// Tri-P surface
			double tri_p(const Eigen::Vector3d&);
			Eigen::Vector3d dtri_p(const Eigen::Vector3d&);

			// distorted P surface
			double f2(const Eigen::Vector3d&);
			Eigen::Vector3d df2(const Eigen::Vector3d&);

			double f3(const Eigen::Vector3d&);
			Eigen::Vector3d df3(const Eigen::Vector3d&);

			double iwp(const Eigen::Vector3d&);
			Eigen::Vector3d diwp(const Eigen::Vector3d& p);

			double f4(const Eigen::Vector3d&);
			Eigen::Vector3d df4(const Eigen::Vector3d& p);

			double tri_p6(const Eigen::Vector3d& p);
			Eigen::Vector3d dtri_p6(const Eigen::Vector3d& p);

			std::tuple<decltype(&f2), decltype(&df2)> get(std::string fname);

			namespace D1 {
				double cos_quad(double t);
				double dcos_quad(double t);
				double d2cos_quad(double t);

				double f1(double t);
				double df1(double t);
				double d2f1(double t);
			}
		}

		Eigen::Matrix3d implicit_bform(const Eigen::Vector3d& df, const Eigen::Matrix3d& Hf);
	}
	struct E3Field {
		virtual Eigen::Vector3d at(const Eigen::Vector3d& p) const = 0;
		virtual Eigen::Matrix3d grad(const Eigen::Vector3d& p) const = 0;
	};


	struct AnalyticalSurface {
		virtual Eigen::Vector3d eval(double u, double v) const = 0;
		virtual Eigen::Matrix<double, 3, 2> tgnt(double u, double v) const = 0;
		virtual Eigen::Matrix<double, 3, 4> hess(double u, double v) const = 0;
		// return [urange, vrange]
		virtual std::tuple<Eigen::Vector2d, Eigen::Vector2d> domain() const = 0;
		virtual Eigen::Vector2d project(const Eigen::Vector3d& p) const {
			std::cout << "p = " << p.transpose() << std::endl;
			Eigen::Vector2d uv; uv.setRandom();
			auto [udom, vdom] = domain();
			uv[0] = (udom[1] - udom[0]) * uv[0] + udom[0];
			uv[1] = (vdom[1] - vdom[0]) * uv[1] + vdom[0];
			for (int iter = 0; iter < 20; iter++) {
				auto r = eval(uv[0], uv[1]);
				auto tg = tgnt(uv[0], uv[1]);
				Eigen::Vector2d d = tg.transpose() * (r - p);
				if (d.squaredNorm() < 1e-16) return uv;
				auto Hr = hess(uv[0], uv[1]);
				Eigen::Matrix2d H = (Hr.transpose() * (r - p)).reshaped(2, 2) + tg.transpose() * tg;
				uv -= H.lu().solve(d);
			}
			return uv;
		}
		Eigen::Matrix2d g(double u, double v) const { auto gtan = tgnt(u, v); return gtan.transpose() * gtan; }
		Eigen::Matrix2d Ig(double u, double v) const { return g(u, v).inverse(); }
		Eigen::Vector3d An(double u, double v) const { auto gtan = tgnt(u, v); return gtan.col(0).cross(gtan.col(1)) / 2; }
		Eigen::Vector3d n(double u, double v) const { return An(u, v).normalized(); }
		Eigen::Matrix2d b(double u, double v) const {
			auto d2r = hess(u, v);
			auto nv = n(u, v);
			return (nv.transpose() * d2r).reshaped(2, 2);
		}
		Eigen::Matrix2d Ib(double u, double v) const {
			auto b22 = b(u, v); auto ig = Ig(u, v);
			return ig * b22 * ig;
		}
		Eigen::Matrix<double, 3, 2> Itgnt(double u, double v) const {
			auto tg = tgnt(u, v);
			Eigen::Matrix<double, 3, 2> Itg = (tg.transpose() * tg).lu().solve(tg.transpose()).transpose();
			return Itg;
		}
		Eigen::Matrix<double, 4, 2> GM(double u, double v) const {
			auto H = hess(u, v);
			auto tg = tgnt(u, v);
			Eigen::Matrix<double, 3, 2> Itg = (tg.transpose() * tg).lu().solve(tg.transpose()).transpose();
			Eigen::Matrix<double, 4, 2> Gamma;
			Gamma.col(0) = (H.transpose() * Itg.col(0)).reshaped();
			Gamma.col(1) = (H.transpose() * Itg.col(1)).reshaped();
			return Gamma;
		}

		Eigen::Matrix2d gamma(const E3Field& V, double u, double v) const {
			auto p = eval(u, v);
			auto dp = tgnt(u, v);
			auto nv = n(u, v);
			auto dV = V.grad(p);
			auto Gm = GM(u, v);
			Eigen::Matrix2d dVtg = (dp.transpose() * dV * dp);
			dVtg = (dVtg + dVtg.transpose()).eval() / 2;
			auto bf = b(u, v);
			double V3 = V.at(p).dot(nv);
			return dVtg /*- bf * V3*/;
		}
	};


	struct Sphere
		: public AnalyticalSurface
	{

		virtual std::tuple<Eigen::Vector2d, Eigen::Vector2d> domain() const override;

		constexpr static double R = 0.5;
		virtual Eigen::Vector3d eval(double u, double v) const override {
			Eigen::Vector3d p;
			p[0] = R * std::cos(v) * std::cos(u);
			p[1] = R * std::cos(v) * std::sin(u);
			p[2] = R * std::sin(v);
			return p;
		}

		virtual Eigen::Matrix<double, 3, 2> tgnt(double u, double v) const override {
			Eigen::Matrix<double, 3, 2> Tg;
			Tg.col(0) = Eigen::Vector3d{ -R * std::cos(v) * std::sin(u), R * std::cos(v) * std::cos(u), 0 };
			Tg.col(1) = Eigen::Vector3d{ -R * std::cos(u) * std::sin(v), -R * std::sin(v) * std::sin(u), R * std::cos(v) };
			return Tg;
		}

		virtual Eigen::Vector2d project(const Eigen::Vector3d& p) const override;

		virtual Eigen::Matrix<double, 3, 4> hess(double u, double v) const override {
			double cu = std::cos(u), cv = std::cos(v), su = std::sin(u), sv = std::sin(v);
			Eigen::Matrix<double, 3, 4> Hr;
			Hr.row(0) = Eigen::Vector4d(-R * cv * cu, R * sv * su, R * sv * su, -R * cv * cu).transpose();
			Hr.row(1) = Eigen::Vector4d(-R * cv * su, -R * cu * sv, -R * cu * sv, -R * cv * su).transpose();
			Hr.row(2) = Eigen::Vector4d(0, 0, 0, -sv / 2);
			return Hr;
		}
	};


	struct V1Field : public E3Field
	{
		virtual Eigen::Vector3d at(const Eigen::Vector3d& p) const override {
			return vfn::V1(p);
		}

		virtual Eigen::Matrix3d grad(const Eigen::Vector3d& p) const override {
			return vfn::dV1(p);
		}
	};

	template<typename F>
	Eigen::Vector3d search_for_zero(F& f, const Eigen::Vector3d& x0, const Eigen::Vector3d& d, double eps, double tol = 1e-8) {
		double a = -eps, b = eps;
		double fx;
		Eigen::Vector3d n = d.normalized();
		Eigen::Vector3d xa = x0 - eps * n;
		Eigen::Vector3d xb = x0 + eps * n;
		double fa = f(xa), fb = f(xb);
		if (fa == 0) return xa;
		if (fb == 0) return xb;
		int iter = 0;
		Eigen::Vector3d x;
		do {
			double t = (a + b) / 2;
			x = x0 + n * t;
			fx = f(x);
			if (fx * fa < 0) {
				b = t; fb = fx;
			} else if (fx * fb < 0) {
				a = t; fa = fx;
			} else {
				return x;
			}
		} while ((b - a) > tol && iter++ < 100);
		return x;
	}
}