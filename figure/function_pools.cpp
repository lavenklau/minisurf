#define _USE_MATH_DEFINES
#include "function_pools.h"
#include <cmath>

inline double Cos(double t) { return std::cos(t); }
inline double Sin(double t) { return std::sin(t); }

inline double Cosh(double t) { return std::cosh(t); }
inline double Sinh(double t) { return std::sinh(t); }

inline double Power(double x, double e) { return std::pow(x, e); }

#define EXPAND(p) double x=p[0],y=p[1],z=p[2];

#define Pi M_PI

Eigen::Vector3d msf::vfn::V1(const Eigen::Vector3d& p)
{
	EXPAND(p);
	return{ Cos(x), Sin(x * y), Sin(x + z) };
}

Eigen::Matrix3d msf::vfn::dV1(const Eigen::Vector3d& p)
{
	EXPAND(p);
	Eigen::Matrix3d dv;
	dv << -Sin(x), 0, 0, y* Cos(x * y), x* Cos(x * y), 0, Cos(x + z), 0, Cos(x + z);
	return dv;
}

double msf::sfn::f1(const Eigen::Vector3d& p)
{
	EXPAND(p);
	return Sin(x * y + z);
}

Eigen::Vector3d msf::sfn::df1(const Eigen::Vector3d& p)
{
	EXPAND(p);
	return { y * Cos(x * y + z), x * Cos(x * y + z), Cos(x * y + z) };
}

double msf::sfn::sphere(const Eigen::Vector3d& p)
{
	EXPAND(p);
	return p.squaredNorm() - 0.5 * 0.5;
}

Eigen::Vector3d msf::sfn::dsphere(const Eigen::Vector3d& p)
{
	EXPAND(p);
	return 2 * p;
}

Eigen::Matrix3d msf::sfn::d2sphere(const Eigen::Vector3d&)
{
	return Eigen::Matrix3d::Identity() * 2;
}

double msf::sfn::cylinder(const Eigen::Vector3d& p)
{
	EXPAND(p);
	return x * x + y * y - 0.5 * 0.5;
}

Eigen::Vector3d msf::sfn::dcylinder(const Eigen::Vector3d& p)
{
	EXPAND(p);
	return Eigen::Vector3d(2 * x, 2 * y, 0);
}

Eigen::Matrix3d msf::sfn::d2cylinder(const Eigen::Vector3d&)
{
	Eigen::Matrix3d H; H.setZero();
	H(0, 0) = H(1, 1) = 2;
	return H;
}

double msf::sfn::catenoid(const Eigen::Vector3d& p)
{
	EXPAND(p);
	return Power(x, 2) + Power(y, 2) - 0.16000000000000003 * Power(Cosh(2.5 * z), 2);
}

Eigen::Vector3d msf::sfn::dcatenoid(const Eigen::Vector3d& p)
{
	EXPAND(p);
	return { 2 * x, 2 * y, -0.8 * Cosh(2.5 * z) * Sinh(2.5 * z) };
}

Eigen::Matrix3d msf::sfn::d2catenoid(const Eigen::Vector3d& p)
{
	EXPAND(p);
	Eigen::Matrix3d H; H.setZero();
	H(0, 0) = H(1, 1) = 2;
	H(2, 2) = -2. * Power(Cosh(2.5 * z), 2) - 2. * Power(Sinh(2.5 * z), 2);
	return H;
}

std::tuple<decltype(&msf::sfn::sphere), decltype(&msf::sfn::dsphere), decltype(&msf::sfn::d2sphere)> msf::sfn::get(std::string fname)
{
	if (fname == "sphere") {
		return std::make_tuple(sphere, dsphere, d2sphere);
	}
	else if (fname == "cylinder") {
		return std::make_tuple(cylinder, dcylinder, d2cylinder);
	}
	else if (fname == "catenoid") {
		return std::make_tuple(catenoid, dcatenoid, d2catenoid);
	}
	throw std::runtime_error(std::string("unknown function name ") + fname);
}

double msf::sfn::period::tri_p(const Eigen::Vector3d& p)
{
	EXPAND(p);
	return Cos(Pi * x) + Cos(Pi * y) + Cos(Pi * z);
}

Eigen::Vector3d msf::sfn::period::dtri_p(const Eigen::Vector3d& p)
{
	EXPAND(p);
	return { -Pi * Sin(Pi * x), -Pi * Sin(Pi * y), -Pi * Sin(Pi * z) };
}

double msf::sfn::period::f2(const Eigen::Vector3d& p)
{
	EXPAND(p);
	return 1.2 * Cos(Pi * x) + 0.9 * Cos(Pi * y) + 1.5 * Cos(Pi * z);
}

Eigen::Vector3d msf::sfn::period::df2(const Eigen::Vector3d& p)
{
	EXPAND(p);
	return { -3.7699111843077517 * Sin(Pi * x), -2.827433388230814 * Sin(Pi * y), -4.71238898038469 * Sin(Pi * z) };
}

double msf::sfn::period::f3(const Eigen::Vector3d& p)
{
	EXPAND(p);
	return 1.2 * Cos(Pi * x) + 0.9 * Cos(Pi * y) + 0.3 * Cos(2 * Pi * (-x + y)) + 0.3 * Cos(2 * Pi * (x + y)) +
		1.5 * Cos(Pi * z);
}

Eigen::Vector3d msf::sfn::period::df3(const Eigen::Vector3d& p)
{
	EXPAND(p);
	return { -3.7699111843077517 * Sin(Pi * x) + 1.8849555921538759 * Sin(2 * Pi * (-x + y)) -
	1.8849555921538759 * Sin(2 * Pi * (x + y)),
   -2.827433388230814 * Sin(Pi * y) - 1.8849555921538759 * Sin(2 * Pi * (-x + y)) -
	1.8849555921538759 * Sin(2 * Pi * (x + y)),-4.71238898038469 * Sin(Pi * z) };
}

double msf::sfn::period::iwp(const Eigen::Vector3d& p)
{
	EXPAND(p);
	return -Cos(2 * Pi * x) + Cos(Pi * x) * Cos(Pi * y) - Cos(2 * Pi * y) + Cos(Pi * x) * Cos(Pi * z) +
		Cos(Pi * y) * Cos(Pi * z) - Cos(2 * Pi * z);
}

Eigen::Vector3d msf::sfn::period::diwp(const Eigen::Vector3d& p)
{
	EXPAND(p);
	return { -(Pi * Cos(Pi * y) * Sin(Pi * x)) - Pi * Cos(Pi * z) * Sin(Pi * x) + 2 * Pi * Sin(2 * Pi * x),
		-(Pi * Cos(Pi * x) * Sin(Pi * y)) - Pi * Cos(Pi * z) * Sin(Pi * y) + 2 * Pi * Sin(2 * Pi * y),
		-(Pi * Cos(Pi * x) * Sin(Pi * z)) - Pi * Cos(Pi * y) * Sin(Pi * z) + 2 * Pi * Sin(2 * Pi * z) };
}

double msf::sfn::period::f4(const Eigen::Vector3d& p)
{
	EXPAND(p);
	return 0.2 + Cos(Pi * x) + Cos(Pi * y + Pi * z) + Cos(Pi * (x + y)) * Sin(Pi * z);
}

Eigen::Vector3d msf::sfn::period::df4(const Eigen::Vector3d& p)
{
	EXPAND(p);
	return { -(Pi * Sin(Pi * x)) - Pi * Sin(Pi * (x + y)) * Sin(Pi * z),
   -(Pi * Sin(Pi * (x + y)) * Sin(Pi * z)) - Pi * Sin(Pi * y + Pi * z),
   Pi * Cos(Pi * (x + y)) * Cos(Pi * z) - Pi * Sin(Pi * y + Pi * z) };
}

double msf::sfn::period::tri_p6(const Eigen::Vector3d& p)
{
	EXPAND(p);
	return Cos(6 * Pi * x) + Cos(6 * Pi * y) + Cos(6 * Pi * z);
}

Eigen::Vector3d msf::sfn::period::dtri_p6(const Eigen::Vector3d& p)
{
	EXPAND(p);
	return { -6 * Pi * Sin(6 * Pi * x), -6 * Pi * Sin(6 * Pi * y), -6 * Pi * Sin(6 * Pi * z) };
}

std::tuple<decltype(&msf::sfn::period::f2), decltype(&msf::sfn::period::df2)>  msf::sfn::period::get(std::string fname)
{
	if (fname == "f2") {
		return std::make_tuple(msf::sfn::period::f2, msf::sfn::period::df2);
	}
	else if (fname == "f3") {
		return std::make_tuple(msf::sfn::period::f3, msf::sfn::period::df3);
	}
	else if (fname == "f4") {
		return std::make_tuple(msf::sfn::period::f4, msf::sfn::period::df4);
	}
	else if (fname == "tri_p") {
		return std::make_tuple(msf::sfn::period::tri_p, msf::sfn::period::dtri_p);
	}
	else if (fname == "tri_p6") {
		return std::make_tuple(msf::sfn::period::tri_p6, msf::sfn::period::dtri_p6);
	}

	throw std::runtime_error(std::string("no function with name ") + fname);
}

std::tuple<Eigen::Vector2d, Eigen::Vector2d> msf::Sphere::domain() const
{
	return std::make_tuple(Eigen::Vector2d(0, 2 * M_PI), Eigen::Vector2d(-M_PI / 2, M_PI / 2));
}

Eigen::Vector2d msf::Sphere::project(const Eigen::Vector3d& p) const
{
	//std::cout << "p = " << p.transpose() << '\t';
	double u = std::atan(p[1] / p[0]);
	if (p[0] < 0) {
		u += M_PI;
	}
	else if (p[1] < 0) {
		u += M_PI * 2;
	}
	double v = std::asin(p[2] / R);
	Eigen::Vector2d uv{ u,v };
	//std::cout << "uv = " << uv.transpose() << std::endl;
	return { u,v };
}


#define Subscript(s,k) s[k-1]
#define Rule(src,dst) dst = src


Eigen::Matrix3d msf::sfn::implicit_bform(const Eigen::Vector3d& df, const Eigen::Matrix3d& Hf) {
	double s[9] = {
		df[0],df[1],df[2],
		Hf(0,0),Hf(0,1),Hf(0,2),Hf(1,1),Hf(1,2),Hf(2,2)
	};

	double Rule(Subscript(s, 9), v0), Rule(Subscript(s, 2), v1), Rule(Subscript(s, 1), v2),
		Rule(Subscript(s, 8), v3), Rule(Subscript(s, 3), v4), Rule(Subscript(s, 7), v5),
		Rule(Subscript(s, 6), v6), Rule(Subscript(s, 5), v7), Rule(Subscript(s, 4), v8),
		Rule(Power(v1, 4), v9), Rule(Power(v1, 2), v10), Rule(Power(v2, 2), v11), Rule(Power(v2, 4), v12),
		Rule(Power(v1, 3), v13), Rule(Power(v4, 2), v14), Rule(Power(v2, 3), v15), Rule(Power(v4, 3), v16),
		Rule(Power(v4, 4), v17), Rule(2 * v1 * v11 * v3 * v4, v18), Rule(2 * v10 * v2 * v4 * v6, v19),
		Rule(2 * v1 * v14 * v2 * v7, v20), Rule(Power(v10 + v11 + v14, -2.5), v21),
		Rule(v21 * (-(v10 * v11 * v3) - v12 * v3 - 2 * v10 * v14 * v3 - v11 * v14 * v3 + v0 * v1 * v11 * v4 +
			v0 * v13 * v4 + v1 * v16 * v5 + v1 * v11 * v4 * v5 + v1 * v15 * v6 + v13 * v2 * v6 - v1 * v14 * v2 * v6 +
			v16 * v2 * v7 + v15 * v4 * v7 - v10 * v2 * v4 * v7 - v1 * v11 * v4 * v8), v22),
		Rule(v21 * (v1 * v15 * v3 + v13 * v2 * v3 - v1 * v14 * v2 * v3 + v0 * v15 * v4 + v0 * v10 * v2 * v4 -
			v10 * v2 * v4 * v5 - v10 * v11 * v6 - v10 * v14 * v6 - 2 * v11 * v14 * v6 + v1 * v16 * v7 - v1 * v11 * v4 * v7 +
			v13 * v4 * v7 + v16 * v2 * v8 + v10 * v2 * v4 * v8 - v6 * v9), v23),
		Rule(v21 * (-(v0 * v1 * v14 * v2) + v16 * v2 * v3 + v15 * v3 * v4 - v10 * v2 * v3 * v4 + v1 * v15 * v5 +
			v1 * v14 * v2 * v5 + v1 * v16 * v6 - v1 * v11 * v4 * v6 + v13 * v4 * v6 - 2 * v10 * v11 * v7 - v10 * v14 * v7 -
			v11 * v14 * v7 - v17 * v7 + v13 * v2 * v8 + v1 * v14 * v2 * v8), v24);

	double bform[3][3] = { {v21 * (-(v0 * v11 * v14) + v19 + v20 - 2 * v1 * v11 * v3 * v4 - v10 * v11 * v5 + 2 * v16 * v2 * v6 +
		2 * v13 * v2 * v7 - 2 * v10 * v14 * v8 - v17 * v8 - v8 * v9), v24, v23},
		{v24, v21 * (-(v0 * v10 * v14) + v18 + v20 + 2 * v1 * v16 * v3 - v12 * v5 - 2 * v11 * v14 * v5 -
			v17 * v5 - 2 * v10 * v2 * v4 * v6 + 2 * v1 * v15 * v7 - v10 * v11 * v8), v22},
		{v23, v22, v21 * (-2 * v0 * v10 * v11 - v0 * v12 + v18 + v19 + 2 * v13 * v3 * v4 - v10 * v14 * v5 +
			2 * v15 * v4 * v6 - 2 * v1 * v14 * v2 * v7 - v11 * v14 * v8 - v0 * v9)} };

	return Eigen::Matrix3d::Map(bform[0]);
}


double msf::sfn::period::D1::cos_quad(double t)
{
	const double quad = 1. / 4;
	return (2 + std::cos(M_PI * t)) * quad;
}

double msf::sfn::period::D1::dcos_quad(double t)
{
	const double quad = 1. / 4;
	return -quad * M_PI * std::sin(M_PI * t);
}

double msf::sfn::period::D1::d2cos_quad(double t)
{
	const double quad = 1. / 4;
	return -quad * M_PI * M_PI * std::cos(M_PI * t);
}

double msf::sfn::period::D1::f1(double t)
{
	return std::sin(M_PI * t);
}

double msf::sfn::period::D1::df1(double t)
{
	return M_PI * std::cos(M_PI * t);
}

double msf::sfn::period::D1::d2f1(double t)
{
	return -M_PI * M_PI * std::sin(M_PI * t);
}
