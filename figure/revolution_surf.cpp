#include "mesher/asymptotic_analysis.h"
#include "Eigen/PardisoSupport"
#include "material/materail.h"
#include "fmt/core.h"
#include "matlab/matlab_utils.h"
#include "function_pools.h"
#include "boost/algorithm/string.hpp"
#include "mesher/isosurface_generation.h"
#include "cgal/cgal_utils.h"
#include "igl/write_triangle_mesh.h"

using namespace msf;

inline double Cos(double t) { return std::cos(t); }
inline double Sin(double t) { return std::sin(t); }

inline double Cosh(double t) { return std::cosh(t); }
inline double Sinh(double t) { return std::sinh(t); }

inline double Power(double x, double e) { return std::pow(x, e); }
inline double Sqrt(double x) { return std::sqrt(x); }

#define EXPAND(p) double x=p[0],y=p[1],z=p[2];

#define Pi M_PI

#define Subscript(s,k) s[k-1]
#define Rule(src,dst) dst = src


void clamp_vertices(Eigen::MatrixX3d& V, Eigen::MatrixX3i& F);

template<typename F, typename DF>
inline void project_vertex_to_zero(PeriodSurfaceMesh& mesh, F flev, DF df) {
#pragma omp parallel for
	for (int vid = 0; vid < mesh.n_vertices(); vid++) {
		auto p = toEigen(mesh.points()[vid]);
		auto pz = search_for_zero(flev, p, df(p), 2e-2);
		auto vh = OM::VertexHandle(vid);
		mesh.set_point(vh, toOM(pz));
	}
	return;
}

struct RevoSurf 
{
	std::function<double(double)> frule, dfrule, d2frule;

	RevoSurf(std::function<double(double)> f_, std::function<double(double)> df_, std::function<double(double)> d2f_) 
		: frule(f_), dfrule(df_), d2frule(d2f_) {}

private:
	double _f(double x) const { return frule(x); }
	double _df(double x) const { return dfrule(x); }
	double _d2f(double x) const { return d2frule(x); }
public:
	Eigen::Vector3d eval(double x, double v) const {
		double y = _f(x) * std::cos(v);
		double z = _f(x) * std::sin(v);
		return { x,y,z };
	}

	Eigen::Matrix<double, 3, 2> tang(double x, double v) const {
		double f = _f(x), df = _df(x);
		Eigen::Vector3d g1 = { 1,df * std::cos(v), df * std::sin(v) };
		Eigen::Vector3d g2 = { 0, -f * std::sin(v), f * std::cos(v) };
		Eigen::Matrix<double, 3, 2> gtan;
		gtan << g1, g2;
		return gtan;
	}

	double area_form(double x) const {
		double f = _f(x), df = _df(x);
		return std::sqrt(f * f * (1 + df * df));
	}

	Eigen::Matrix2d bform(double x) const{
		double f = _f(x), df = _df(x), d2f = _d2f(x);
		Eigen::Matrix2d b; b.setZero();
		double A = std::sqrt(f * f * (1 + df * df));
		b(0, 0) = -f * d2f / A;
		b(1, 1) = f * f / A;
		return b;
	}

	Eigen::Matrix2d Chris(double x, const Eigen::Vector2d& ut) const {
		double f = _f(x), df = _df(x), d2f = _d2f(x);
		double y = ut[0], z = ut[1];
		Eigen::Matrix2d Chr;
		Chr << y * df * d2f / (1 + df * df), z* df / f, z* df / f, -y * f * df / (1 + df * df);
		return Chr;
	}

	Eigen::Matrix2d metric(double x) const {
		double df = _df(x), f = _f(x);
		return Eigen::Vector2d(1 + df * df, f * f).asDiagonal();
	}

	Eigen::Vector3d normal(double x, double v) const {
		double f = _f(x), df = _df(x), d2f = _d2f(x);
		Eigen::Vector3d n(f * df, -std::cos(v) * f, -f * std::sin(v));
		return n;
	}


	auto mesh(void) const {
		Eigen::Vector3d pl, pr;
		pl.setConstant(-1); pr.setConstant(1);
		auto flev = [=](double x, double y, double z) {return y * y + z * z - std::pow(frule(x), 2); };
		auto fgrad = [=](auto& p) {auto p1 = p; double f = _f(p[0]), df = _df(p[0]); p1.bottomRows(2).normalize(); return Eigen::Vector3d(df, -p1[1], -p1[2]); };
		auto [V, F] = isosurf_mesh(128, std::make_pair(pl, pr), flev, 0);
		igl::write_triangle_mesh(getPath("iglmesh.obj"), V, F);
		std::tie(V, F) = cgal_remesh(V, F);
		clamp_vertices(V, F);
		
		PeriodSurfaceMesh mesh;
		mesh.read(V, F, false, false, false);

		auto vcut = mesh.mergePeriodBoundary();

		mesh.periodic_remesh(20, vcut, 0.02, 0);

		project_vertex_to_zero(mesh, [=](auto& p) {return p[1] * p[1] + p[2] * p[2] - std::pow(frule(p[0]), 2); }, fgrad);

		mesh.period_shift();

		mesh.clamp_period_boundary(2e-4);

		mesh.garbage_collection();

		return mesh;
	}

	auto quad_matrix(double lam, double mu, double d2f, double df, double f) const {
		double 	Rule(f, v0), Rule(2 * mu, v1), Rule(lam + mu, v2), Rule(-lam, v3),
			Rule(df, v4), Rule(Power(v0, 2), v5), Rule(d2f, v6),
			Rule(Power(v4, 2), v7), Rule(Power(v4, 4), v8), Rule(1 / (v1 + lam), v9), Rule(1 + v7, v10),
			Rule(v7 * lam, v11), Rule(v8 * mu, v12), Rule(v8 * lam, v13), Rule(v0 * v1 * v6, v14), Rule(2 * v0 * v6 * lam, v15),
			Rule(v3 * v7, v16), Rule(v1 * v7, v17), Rule(2 * v11, v18), Rule(Power(v10, -2), v19),
			Rule(v5 * Power(v6, 2), v20), Rule(v20 * mu, v21), Rule(v20 * lam, v22), Rule(v10 * v5, v23),
			Rule(Sqrt(v23), v24), Rule(1 / Sqrt(v23), v25), Rule(Power(v23, -1.5), v26),
			Rule(v11 + v14 + v15 + lam, v27), Rule(v14 + v15 + v16 + v3, v28),
			Rule((4 * Pi * v28 * v9 * mu) / v10, v29), Rule((8 * Pi * v19 * v2 * v4 * (-1 + v20 - 2 * v7 - v8) * v9 * mu) / v0, v30),
			Rule(4 * Pi * v0 * v25 * v27 * v4 * v9 * mu, v31);
		double clist[] = { 8 * Pi * v19 * v2 * v24 * v9 * mu, 8 * Pi * v19 * v28 * v9 * mu, 8 * Pi * Power(v0, 3) * v26 * v27 * v4 * v9 * mu,
			16 * Pi * v2 * v25 * v5 * v9 * mu, 8 * Pi * v19 * v25 *
			(v12 + v13 + v17 + v18 + v2 + v21 + v22 + v0 * v16 * v6 + v0 * v3 * v6) * v9 * mu, v30, v29, v30,
			8 * Pi * v26 * v5 * v7 * v9 * (v12 + v13 + v17 + v18 + v2 + v21 + v22 + v0 * v11 * v6 + v0 * v6 * lam) * mu, v31,
			v29, v31, 8 * Pi * v2 * v24 * v9 * mu };
	
		double c0 = clist[0];
		Eigen::Vector3d c1 = Eigen::Vector3d::Map(clist + 1);
		Eigen::Matrix3d c2 = Eigen::Matrix3d::Map(clist + 4);
		return std::make_tuple(c0, c1, c2);
	}

	auto min_quadratic_energy(double lam, double mu, int n_steps) const {
		Eigen::SparseMatrix<double> C2(n_steps * 2, n_steps * 2);
		Eigen::VectorXd C1(n_steps * 2); C1.setZero();
		double C0 = 0;
		std::vector<Eigen::Triplet<double>> triplist;
		const double du = 2. / n_steps;
		for (int i = 0; i < n_steps; i++) {
			double ui = i * du - 1;
			double f = _f(ui), df = _df(ui), d2f = _d2f(ui);
			auto [c0, c1, c2] = quad_matrix(lam, mu, d2f, df, f);
			for (int k = 0; k < 2; k++) {
				for (int l = 0; l < 2; l++) {
					triplist.emplace_back(i * 2 + k, i * 2 + l, c2(k, l));
				}
			}
			int jnext = (i + 1) % n_steps;
			int jbef = (i + n_steps - 1) % n_steps;
			for (int k = 0; k < 2; k++) {
				triplist.emplace_back(i * 2 + k, jnext * 2 + 1, c2(k, 2) / 2 / du);
				triplist.emplace_back(jnext * 2 + 1, i * 2 + k, c2(k, 2) / 2 / du);
				triplist.emplace_back(i * 2 + k, jbef * 2 + 1, -c2(k, 2) / 2 / du);
				triplist.emplace_back(jbef * 2 + 1, i * 2 + k, -c2(k, 2) / 2 / du);
			}
			triplist.emplace_back(jnext * 2 + 1, jnext * 2 + 1, c2(2, 2) / 4 / du / du);
			triplist.emplace_back(jbef * 2 + 1, jbef * 2 + 1, c2(2, 2) / 4 / du / du);
			triplist.emplace_back(jbef * 2 + 1, jnext * 2 + 1, -c2(2, 2) / 4 / du / du);
			triplist.emplace_back(jnext * 2 + 1, jbef * 2 + 1, -c2(2, 2) / 4 / du / du);
			C1[i * 2] += c1[0]; C1[i * 2 + 1] += c1[1];
			C1[jnext * 2 + 1] += c1[2] / 2 / du; C1[jbef * 2 + 1] += -c1[2] / 2 / du;
			C0 += c0;
		}
		C2.setFromTriplets(triplist.begin(), triplist.end());
		C2 *= du ; C1 *= du; C0 *= du;
		Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> so(C2 * 2);
		Eigen::VectorXd x = -so.solve(C1);
		if (x.hasNaN()) { std::cout << "x is NaN!" << std::endl; throw std::runtime_error("NaN in x"); }
		return std::make_tuple(C0, C1, C2, x);
	}

	auto variation_min(double lam, double mu) const {
		constexpr int n_steps = 100000;
		double du = 2. / n_steps;
		//eigen2ConnectedMatlab("C2", C2);
		//eigen2ConnectedMatlab("C1", C1);
		auto [C0, C1, C2, x] = min_quadratic_energy(lam, mu, n_steps);
		double en = C0 + C1.dot(x) + x.dot(C2 * x);
		if (std::isnan(en)) { std::cout << "en is NaN!" << std::endl; throw std::runtime_error("NaN in en"); }
		std::vector<Eigen::Vector3d> plist, vlist;
		double v_sample = 0;
		for (int i = 0; i < n_steps; i++) {
			double ui = i * du - 1;
			auto p = eval(ui, v_sample);
			auto tan = tang(ui, v_sample);
			auto n = tan.col(0).cross(tan.col(1)).normalized();
			Eigen::Vector3d vp = x[i * 2] * n + x[i * 2 + 1] * tan.col(0);
			plist.emplace_back(p); vlist.emplace_back(vp);
		}
		return std::make_tuple(en, plist, vlist);
	}

	double itxraw(double lam, double mu, double d3f, double d2f, double df, double f, double vn, double dvn, double d2vn, double s, double tau, double dtau) const {
		double z = vn, dz = dvn, d2z = d2vn;
		return (4 * Pi * mu * (-4 * (lam + mu) * f * s * Power(1 + Power(df, 2), 3) * Sqrt(Power(f, 2) * (1 + Power(df, 2))) *
			(f * s - tau * df * Sqrt(Power(f, 2) * (1 + Power(df, 2)))) * (z * Sqrt(Power(f, 2) * (1 + Power(df, 2))) - f * df * (z * (1 + Power(df, 2)) +
					Sqrt(Power(f, 2) * (1 + Power(df, 2))) * dz) - Power(f, 2) * ((1 + Power(df, 2)) * dz + z * df * d2f)) + (2 * lam * s * Power(Power(f, 2) * (1 + Power(df, 2)), 2.5) *
				(z * Sqrt(Power(f, 2) * (1 + Power(df, 2))) - f * df * (z * (1 + Power(df, 2)) + Sqrt(Power(f, 2) * (1 + Power(df, 2))) * dz)\
					- Power(f, 2) * ((1 + Power(df, 2)) * dz + z * df * d2f)) * ((Power(Power(f, 2) * (1 + Power(df, 2)), 1.5) * dtau) / Power(f, 2) + f * s * d2f +
					Sqrt(Power(f, 2) * (1 + Power(df, 2))) * (1 + tau * df * d2f))) / Power(f, 2) + 2 * Power(f, 2) * z * (1 + Power(df, 2)) * ((lam + mu) * Power(1 + Power(df, 2), 2) *
				Power(f * s - tau * df * Sqrt(Power(f, 2) * (1 + Power(df, 2))), 2) * (3 + 3 * Power(df, 2) + f * d2f) - lam * f * (1 + Power(df, 2)) *
				(f * s - tau * df * Sqrt(Power(f, 2) * (1 + Power(df, 2)))) * (1 + Power(df, 2) - f * d2f) * ((Power(Power(f, 2) * (1 + Power(df, 2)), 1.5) *
					dtau) / Power(f, 2) + f * s * d2f + Sqrt(Power(f, 2) * (1 + Power(df, 2))) * (1 + tau * df * d2f)) - (lam + mu) * Power(f, 2) * (1 + Power(df, 2) +
					3 * f * d2f) * Power((Power(Power(f, 2) * (1 + Power(df, 2)), 1.5) * dtau) / Power(f, 2) + f * s * d2f + Sqrt(Power(f, 2) * (1 + Power(df, 2))) *
					(1 + tau * df * d2f), 2)) - 2 * lam * Power(f, 4) * Power(1 + Power(df, 2), 2) * (f * s - tau * df * Sqrt(Power(f, 2) * (1 + Power(df, 2)))) *
			(f * s * (df * dz * d2f + Power(df, 3) * dz * d2f + z * Power(d2f, 2) - d2z - 2 * Power(df, 2) * d2z - Power(df, 4) * d2z) + Sqrt(Power(f, 2) * (1 + Power(df, 2))) *
				(2 * Power(df, 3) * dz + df * (2 * dz + 3 * s * z * Power(d2f, 2)) - s * (dz * d2f + z * d3f) - s * Power(df, 2) * (dz * d2f + z * d3f))) +
			4 * (lam + mu) * Power(f, 5) * (1 + Power(df, 2)) * ((Power(Power(f, 2) * (1 + Power(df, 2)), 1.5) * dtau) / Power(f, 2) + f * s * d2f + Sqrt(Power(f, 2) * (1 + Power(df, 2))) *
				(1 + tau * df * d2f)) * (f * s * (df * dz * d2f + Power(df, 3) * dz * d2f + z * Power(d2f, 2) - d2z - 2 * Power(df, 2) * d2z - Power(df, 4) * d2z) +
				Sqrt(Power(f, 2) * (1 + Power(df, 2))) * (2 * Power(df, 3) * dz + df * (2 * dz + 3 * s * z * Power(d2f, 2)) - s * (dz * d2f + z * d3f) -
					s * Power(df, 2) * (dz * d2f + z * d3f))))) / ((lam + 2 * mu) * Power(f, 6) * Power(1 + Power(df, 2), 5));
	}

	double itx_new(double lam, double mu, double d2f, double df, double f, double vn, double dvn, double s, double ds, double tau, double dtau) const {
		double z = vn, dz = dvn;
		return (8 * Pi * mu * Sqrt(Power(f, 2) * (1 + Power(df, 2))) * (2 * (lam + mu) * s * z * Power(1 + Power(df, 2), 4) *
			(-(f * s) + tau * df * Sqrt(Power(f, 2) * (1 + Power(df, 2)))) * (Sqrt(Power(f, 2) * (1 + Power(df, 2))) - f * (df + Power(df, 3))) -
			lam * f * s * z * Power(1 + Power(df, 2), 3) * (-Sqrt(Power(f, 2) * (1 + Power(df, 2))) + f * (df + Power(df, 3))) *
			((Power(Power(f, 2) * (1 + Power(df, 2)), 1.5) * dtau) / Power(f, 2) + f * s * d2f + Sqrt(Power(f, 2) * (1 + Power(df, 2))) *
				(1 + tau * df * d2f)) + (f * z * Power(1 + Power(df, 2), 2) * ((lam + mu) * Power(1 + Power(df, 2), 2) *
					Power(f * s - tau * df * Sqrt(Power(f, 2) * (1 + Power(df, 2))), 2) * (3 + 3 * Power(df, 2) + f * d2f) - lam * f * (1 + Power(df, 2)) *
					(f * s - tau * df * Sqrt(Power(f, 2) * (1 + Power(df, 2)))) * (1 + Power(df, 2) - f * d2f) * ((Power(Power(f, 2) * (1 + Power(df, 2)), 1.5) *
						dtau) / Power(f, 2) + f * s * d2f + Sqrt(Power(f, 2) * (1 + Power(df, 2))) * (1 + tau * df * d2f)) - (lam + mu) * Power(f, 2) * (1 + Power(df, 2) +
							3 * f * d2f) * Power((Power(Power(f, 2) * (1 + Power(df, 2)), 1.5) * dtau) / Power(f, 2) + f * s * d2f + Sqrt(Power(f, 2) * (1 + Power(df, 2))) *
								(1 + tau * df * d2f), 2))) / Sqrt(Power(f, 2) * (1 + Power(df, 2))) - lam * (1 + Power(df, 2)) * (f * s - tau * df *
									Sqrt(Power(f, 2) * (1 + Power(df, 2)))) * (Power(1 + Power(df, 2), 2) * Power(Power(f, 2) * (1 + Power(df, 2)), 1.5) * ds * dz + s * z *
										Power(Power(f, 2) * (1 + Power(df, 2)), 1.5) * Power(d2f, 2) + Power(f, 3) * Power(1 + Power(df, 2), 2) * (2 * Power(df, 3) * dz +
											(z * ds - s * dz) * d2f + Power(df, 2) * (z * ds - s * dz) * d2f + df * (2 * dz +
												s * z * Power(d2f, 2)))) + 2 * (lam + mu) * f * ((Power(Power(f, 2) * (1 + Power(df, 2)), 1.5) * dtau) / Power(f, 2) + f * s * d2f + Sqrt(Power(f, 2) * (1 + Power(df, 2))) *
													(1 + tau * df * d2f)) * (Power(1 + Power(df, 2), 2) * Power(Power(f, 2) * (1 + Power(df, 2)), 1.5) * ds * dz + s * z *
														Power(Power(f, 2) * (1 + Power(df, 2)), 1.5) * Power(d2f, 2) + Power(f, 3) * Power(1 + Power(df, 2), 2) * (2 * Power(df, 3) * dz +
															(z * ds - s * dz) * d2f + Power(df, 2) * (z * ds - s * dz) * d2f + df * (2 * dz +
																s * z * Power(d2f, 2)))))) / ((lam + 2 * mu) * Power(f, 5) * Power(1 + Power(df, 2), 6));
	}


	double sens(int type, double lam, double mu, std::function<double(double)> vnfunc) {
		constexpr int n_steps = 100000;
		auto [C0, C1, C2, x] = min_quadratic_energy(lam, mu, n_steps);
		const double du = 2. / n_steps;
		double it = 0;
#pragma omp parallel for reduction(+:it)
		for (int i = 0; i < n_steps; i++) {
			double ui = i * du - 1;
			double s = x[i * 2], tau = x[i * 2 + 1];
			int inxt = (i + 1) % n_steps;
			int iprev = (i + n_steps - 1) % n_steps;
			double ds = (x[inxt * 2] - x[iprev * 2]) / 2 / du;
			double dtau = (x[inxt * 2 + 1] - x[iprev * 2 + 1]) / 2 / du;
			double f = frule(ui), df = dfrule(ui), d2f = d2frule(ui);
			double d3f = (d2frule(ui + du) - d2frule(ui - du)) / 2 / du;
			double vn = vnfunc(ui), dvn = (vnfunc(ui + du) - vnfunc(ui - du)) / 2 / du;
			double d2vn = (vnfunc(ui + du) + vnfunc(ui - du) - 2 * vnfunc(ui)) / du / du;
			double fi;
			if (type == 0) {
				fi = itxraw(lam, mu, d3f, d2f, df, f, vn, dvn, d2vn, s, tau, dtau) ;
			} else {
				fi = itx_new(lam, mu, d2f, df, f, vn, dvn, s, ds, tau, dtau);
			}
			it += fi;
		}
		return it * du;
	}

	double sens_test(int type, double lam, double mu, std::function<double(double)> vnfunc, std::function<double(double)> sfunc, std::function<double(double)> taufunc) {
		constexpr int n_steps = 100000;
		const double du = 2. / n_steps;
		double it = 0;
#pragma omp parallel for reduction(+:it)
		for (int i = 0; i < n_steps; i++) {
			double ui = i * du - 1;
			double s = sfunc(ui), tau = taufunc(ui);
			int inxt = (i + 1) % n_steps;
			int iprev = (i + n_steps - 1) % n_steps;
			double ds = (sfunc(ui + du) - sfunc(ui - du)) / 2 / du;
			double dtau = (taufunc(ui + du) - taufunc(ui - du)) / 2 / du;
			double f = frule(ui), df = dfrule(ui), d2f = d2frule(ui);
			double d3f = (d2frule(ui + du) - d2frule(ui - du)) / 2 / du;
			double vn = vnfunc(ui), dvn = (vnfunc(ui + du) - vnfunc(ui - du)) / 2 / du;
			double d2vn = (vnfunc(ui + du) + vnfunc(ui - du) - 2 * vnfunc(ui)) / du / du;
			double fi;
			if (type == 0) {
				fi = itxraw(lam, mu, d3f, d2f, df, f, vn, dvn, d2vn, s, tau, dtau) ;
			} else {
				fi = itx_new(lam, mu, d2f, df, f, vn, dvn, s, ds, tau, dtau);
				if (i == 30000) {
					std::cout << "ui = " << ui << "\n";
					std::cout << "fi = " << fi << std::endl;
					std::cout << "s = " << s << ", tau = " << tau << ", ds = " << ds << ", dtau = " << dtau << ", f = " << f << ", d3f = " << d3f << std::endl;
					std::cout << "vn = " << vn << ", d2vn = " << d2vn << std::endl;
				}
			}
			it += fi;
		}
		return it * du;
	}

	double area(void) const {
		constexpr int n_steps = 100000;
		double du = 2. / n_steps;
		double As = 0;
		for (int i = 0; i < n_steps; i++) {
			double A = area_form(du * i - 1);
			As += A * M_PI * 2 * du;
		}
		return As;
	}
};

Eigen::Matrix<double, 6, 6> eval_ads_matrix(PeriodSurfaceMesh& mesh, double lam0, double mu);

extern double eval_diameter(PeriodSurfaceMesh& m);

void test_revolution_surface(std::string obj) {

	auto [lam, mu] = mtl::lameCoeff(1, 0.3);
	double lam0 = 2 * lam * mu / (lam + 2 * mu);

	int N = 2;

	if (!obj.empty()) { N = *obj.rbegin() - '0'; }
	if (N > 10) { std::cout << "invalid objective" << std::endl; throw std::runtime_error(""); }

	auto f_ = [=](double x) { return sfn::period::D1::cos_quad(N * x); };
	auto df_ = [=](double x) { return N * sfn::period::D1::dcos_quad(N * x); };
	auto d2f_ = [=](double x) { return N * N * sfn::period::D1::d2cos_quad(N * x); };
	RevoSurf surf(f_, df_, d2f_);

	auto m = surf.mesh();
	m.saveUnitCell(getPath("revo.obj"));

	auto [en, plist, vlist] = surf.variation_min(lam, mu);

	double sol_semi_ana = en / surf.area();

	std::cout << "vari min = " << sol_semi_ana << std::endl;

	auto flev = [=](auto& p) {double x = p[0], y = p[1], z = p[2]; return y * y + z * z - std::pow(surf.frule(x), 2); };
	auto fgrad = [=](auto& p) {auto p1 = p; double f = f_(p[0]), df = df_(p[0]); p1.bottomRows(2).normalize(); return Eigen::Vector3d(df, -p1[1], -p1[2]); };


	std::ofstream logfile(getPath("errlog"), std::ios::app);
	logfile << "semi_ana = " << sol_semi_ana << std::endl;

	double h0 = 0.1;
	for (int iter = 0; iter < 15; iter++) {
		std::cout << "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = \n";
		double h = h0 * std::pow(0.8, iter);
		m.delaunayRemesh(5, h, h * 0.8, 1e-10, 5);
		
		project_vertex_to_zero(m, flev, fgrad);

		auto CA = eval_ads_matrix(m, lam0, mu);

		double hm = eval_diameter(m);

		std::cout << "hm = " << hm << std::endl;
		std::cout << "CA = \n" << CA << std::endl;

		logfile << "hm = " << hm << std::endl;
		logfile << "CA = \n" << CA << std::endl;
		logfile << "c11 = " << CA(0, 0) << std::endl;
	}
}

using namespace sfn::period;


extern std::vector<double> aux_number;
void sensitivity_test(std::string obj) {
	auto [lam, mu] = mtl::lameCoeff(1, 0.3);
	double lam0 = 2 * lam * mu / (lam + 2 * mu);

	int N = 1;
	auto f_ = [=](double x) { return sfn::period::D1::cos_quad(N * x); };
	auto df_ = [=](double x) { return N * sfn::period::D1::dcos_quad(N * x); };
	auto d2f_ = [=](double x) { return N * N * sfn::period::D1::d2cos_quad(N * x); };
	RevoSurf surf(f_, df_, d2f_);
	auto m = surf.mesh();
	m.saveUnitCell(getPath("revo.obj"));
	auto [en, plist, vlist] = surf.variation_min(lam, mu);
	double sol_semi_ana = en /*/ surf.area()*/;
	std::cout << "vari min = " << sol_semi_ana << std::endl;

	double h = 0.03;
	m.delaunayRemesh(5, h, h * 0.8, 1e-10, 5);
	m.saveUnitCell(getPath("revo.obj"));

	auto [frlist, vrings] = generate_face_frame_and_vring(m);
	auto [fm, um, Enmem, As] = asym_stif_fem(m, frlist, vrings, lam0, mu);
	auto CA = eval_ads_matrix(m, lam0, mu);
	//auto [v_sens, A_sens] = asym_stif_shape_sensitivity(m, frlist, vrings, um, fm, lam0, mu);
	std::cout << "num min = " << CA(0, 0) * surf.area() << std::endl;

	//{
	//	auto [Av, Nv] = eval_vertex_mass(m);
	//	std::ofstream ofs(getPath("Vn"), std::ios::binary);
	//	for (auto vh : m.vertices()) {
	//		auto p = toEigen(m.point(vh));
	//		Eigen::Vector3d vn = v_sens(vh.idx(), 0) * Nv[vh.idx()];
	//		ofs.write((const char*)p.data(), sizeof(p));
	//		ofs.write((const char*)vn.data(), sizeof(vn));
	//	}
	//}

	double scal = aux_number.at(0);

	auto vnfunc = [=](double x) { return  scal * (D1::f1(N * x) + 1);  };

	auto sfunc = [=](double x) {return std::cos(M_PI * x); };
	auto taufunc = [=](double x) {return std::sin(M_PI * x); };

	//double sens0_test = surf.sens_test(0, lam, mu, vnfunc, sfunc, taufunc);
	//double sens1_test = surf.sens_test(1, lam, mu, vnfunc, sfunc, taufunc);

	//std::cout << "sens0_test = " << sens0_test << ", sens1_test = " << sens1_test << std::endl;


	auto flev = [=](const auto& p) {double x = p[0], y = p[1], z = p[2]; return y * y + z * z - std::pow(surf.frule(x), 2); };
	auto fgrad = [=](const auto& p) {auto p1 = p; double f = f_(p[0]), df = df_(p[0]); p1.bottomRows(2).normalize(); return Eigen::Vector3d(df, -p1[1], -p1[2]); };
	{
		// num sens test
		auto vx = [=](double x) {return std::cos(M_PI * x); };
		auto vy = [=](double y) {return std::cos(M_PI * y); };
		Eigen::Matrix<double, -1, 6> um(m.n_vertices() * 3, 6); um.setZero();
		for (auto vh : m.vertices()) {
			auto p = m.point(vh); Eigen::Vector3d n = fgrad(toEigen(p)).normalized();
			double y = std::sqrt(p[1] * p[1] + p[2] * p[2]);
			Eigen::Vector3d uspace(vx(p[0]), vy(y), 0);
			Eigen::Vector3d up = rotateAlign(Eigen::Vector3d(0, y, 0).normalized(), Eigen::Vector3d(0, p[1], p[2]).normalized()) * uspace;
			um.block<3, 1>(vh.idx() * 3, 0) = up;
			if ((p - OM::Vec3d(0.300736, -0.498996, 0.41096)).norm() < 1e-2) {
				std::cout << "p = " << p << ", up = " << up.transpose() << std::endl;
				std::cout << "py = " << Eigen::Vector3d(p[0], y, 0).transpose() << std::endl;
				std::cout << "uspace = " << uspace.transpose() << std::endl;
				std::cout << "R = \n" << rotateAlign(Eigen::Vector3d(p[0], y, 0).normalized(), toEigen(p).normalized()) << std::endl;
			}
		}
		auto [v_sens, A_sens] = asym_stif_shape_sensitivity(m, frlist, vrings, um, fm, lam0, mu);

		double sens_num = 0;
		for (auto vh : m.vertices()) {
			auto p = m.point(vh);
			sens_num += vnfunc(p[0]) * v_sens(vh.idx(), 0);
		}
		eigen2ConnectedMatlab("vsens", v_sens);
		std::cout << "sens_num = " << sens_num << std::endl;
	}

	double sens0 = surf.sens(0, lam, mu, vnfunc);
	double sens1 = surf.sens(1, lam, mu, vnfunc);

	std::cout << "sens0 = " << sens0 << ", sens1 = " << sens1 << std::endl;


	//double sens_num = 0;
	//for (auto vh : m.vertices()) {
	//	auto p = m.point(vh);
	//	sens_num += vnfunc(p[0]) * v_sens(vh.idx(), 0);
	//}
	//std::cout << "sens_num = " << sens_num << std::endl;

	//for (auto vh : m.vertices()) {
	//	auto p = m.point(vh);
	//	Eigen::Vector3d n = fgrad(toEigen(p)).normalized();
	//	auto pnew = p + vnfunc(p[0]) * toOM(n);
	//	m.set_point(vh, pnew);
	//}
	//m.saveUnitCell(getPath("offset.obj"));
	//CA = eval_ads_matrix(m, lam0, mu);
	//std::cout << "num min offset = " << CA(0, 0) * surf.area() << std::endl;

	//auto f1_ = [=](double x) { return D1::cos_quad(x) + scal * (D1::f1(N * x) + 1); };
	//auto df1_ = [=](double x) { return D1::dcos_quad(x) + scal * N * D1::df1(N * x); };
	//auto d2f1_ = [=](double x) { return  D1::d2cos_quad(x) + scal * N * N * D1::d2f1(N * x); };
	//RevoSurf surf1(f1_, df1_, d2f1_);
	//auto m1 = surf1.mesh();
	//m1.saveUnitCell(getPath("revo1.obj"));
	//std::tie(en, plist, vlist) = surf1.variation_min(lam, mu);
	//double sol_semi_ana_1 = en /*/ surf1.area()*/;

	//std::cout << "Delta I_a = " << sol_semi_ana_1 - sol_semi_ana << std::endl;

	//auto flev1 = [=](auto& p) {double x = p[0], y = p[1], z = p[2]; return y * y + z * z - std::pow(surf1.frule(x), 2); };
	//auto fgrad1 = [=](auto& p) {auto p1 = p; double f = f1_(p[0]), df = df1_(p[0]); p1.bottomRows(2).normalize(); return Eigen::Vector3d(df, -p1[1], -p1[2]); };

	//project_vertex_to_zero(m, flev, fgrad);
	//m.saveUnitCell(getPath("dem.obj"));
	//CA = eval_asym_elastic_tensor(fm, um, Enmem, As);

	//auto [Av, Nv] = eval_vertex_mass(m);

	//Eigen::VectorXd Vn(m.n_vertices()); Vn.setZero();
	//for (auto vh : m.vertices()) {
	//	auto p = toEigen(m.point(vh));
	//	Eigen::Vector3d n = Nv[vh.idx()];
	//	auto p1 = search_for_zero(flev1, p, n, 0.1);
	//	Vn[vh.idx()] = (p1 - p).dot(n);
	//}

	//{
	//	std::ofstream ofs(getPath("Vn"), std::ios::binary);
	//	for (auto vh : m.vertices()) {
	//		auto p = toEigen(m.point(vh));
	//		Eigen::Vector3d vn = Vn[vh.idx()] * Nv[vh.idx()];
	//		ofs.write((const char*)p.data(), sizeof(p));
	//		ofs.write((const char*)vn.data(), sizeof(vn));
	//	}
	//}

	//std::tie(v_sens, A_sens) = asym_stif_shape_sensitivity(m, frlist, vrings, um, fm, lam0, mu);

	//std::cout << "delta I_a = " << Vn.dot(v_sens.col(0)) << std::endl;

	//double hm = eval_diameter(m);

	//std::cout << "hm = " << hm << std::endl;
	//std::cout << "CA = \n" << CA << std::endl;

}
