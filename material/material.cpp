#include <iostream>
#include "materail.h"
#include "cmath"

#define Power(x,y) std::pow(x,y)
#define Rule(src,dst) dst = src


template<typename Mat>
inline auto Subscript(const Mat& m, int i, int j) { return m(i - 1, j - 1); }

auto compliance_matrix(double lam, double mu) {
	Eigen::Matrix<double, 6, 6> S; S.setZero();
	S.block<3, 3>(0, 0).setConstant(-lam / (6 * lam * mu + 4 * mu * mu));
	S.block<3, 3>(0, 0).diagonal().setConstant((lam + mu) / (3 * lam * mu + 2 * mu * mu));
	S.block<3, 3>(3, 3).diagonal().setConstant(1. / mu);
	return S;
}

std::tuple<Eigen::Matrix<double, 6, 6>, Eigen::Matrix<double, 6, 6>> dS(double lam, double mu) {
	Eigen::Matrix<double, 6, 6> dSdlam, dSdmu;
	dSdlam.setZero(); dSdmu.setZero();
	dSdlam.block<3, 3>(0, 0).setConstant(-1. / std::pow(3 * lam + 2 * mu, 2));

	double ds11dmu = -(3 * std::pow(lam, 2) + 4 * lam * mu + 2 * std::pow(mu, 2)) / std::pow(mu, 2) / std::pow(3 * lam + 2 * mu, 2);
	double ds12dmu = lam * (6 * lam + 8 * mu) / std::pow(6 * lam * mu + 4 * mu * mu, 2);
	double ds44dmu = -1. / std::pow(mu, 2);
	dSdmu.block<3, 3>(0, 0).setConstant(ds12dmu);
	dSdmu.block<3, 3>(0, 0).diagonal().setConstant(ds11dmu);
	dSdmu.block<3, 3>(3, 3).diagonal().setConstant(ds44dmu);

	return std::make_tuple(dSdlam, dSdmu);
}


std::tuple<Eigen::Matrix<double, 6, 6>, Eigen::Matrix<double, 6, 6>, Eigen::Matrix<double, 6, 6>> d2S(double lam, double mu) {
	Eigen::Matrix<double, 6, 6> d2Sdlam2, d2Sdmu2, d2Sdlamdmu;
	d2Sdlam2.setZero(); d2Sdmu2.setZero(); d2Sdlamdmu.setZero();

	d2Sdlam2.block<3, 3>(0, 0).setConstant(6. / std::pow(3 * lam + 2 * mu, 3));

	double a = (2 * (9 * Power(lam, 3) + 18 * Power(lam, 2) * mu + 12 * lam * Power(mu, 2) + 4 * Power(mu, 3))) /
		(Power(mu, 3) * Power(3 * lam + 2 * mu, 3));
	double b = (-3 * lam * (3 * Power(lam, 2) + 6 * lam * mu + 4 * Power(mu, 2))) / (Power(mu, 3) * Power(3 * lam + 2 * mu, 3));
	d2Sdmu2.block<3, 3>(0, 0).setConstant(b);
	d2Sdmu2.block<3, 3>(0, 0).diagonal().setConstant(a);
	d2Sdmu2.block<3, 3>(3, 3).diagonal().setConstant(2. / Power(mu, 3));
	
	a = 4 / Power(3 * lam + 2 * mu, 3);
	d2Sdlamdmu.block<3, 3>(0, 0).setConstant(a);

	return std::make_tuple(d2Sdlam2, d2Sdmu2, d2Sdlamdmu);
}


auto hessian_err(const Eigen::Matrix<double, 6, 6>& s, double& lam, double& mu) {
	double Rule(Power(mu, 5), v0), Rule(Power(mu, 4), v1), Rule(Power(mu, 3), v2), Rule(Power(lam, 2), v3),
		Rule(Power(mu, 2), v4), Rule(Power(lam, 3), v5), Rule(Power(lam, 4), v6), Rule(Subscript(s, 6, 6), v7),
		Rule(Subscript(s, 5, 5), v8), Rule(Subscript(s, 4, 4), v9), Rule(Subscript(s, 3, 3), v10),
		Rule(Subscript(s, 3, 2), v11), Rule(Subscript(s, 3, 1), v12), Rule(Subscript(s, 2, 2), v13),
		Rule(Subscript(s, 2, 1), v14), Rule(Subscript(s, 1, 1), v15), Rule(Power(3 * lam + 2 * mu, -4), v16),
		Rule(-9 + 6 * v10 * lam + 12 * v11 * lam + 12 * v12 * lam + 6 * v13 * lam + 12 * v14 * lam + 6 * v15 * lam + 4 * v10 * mu + 
			8 * v11 * mu + 8 * v12 * mu + 4 * v13 * mu + 8 * v14 * mu + 4 * v15 * mu, v17), Rule(-4 * v16 * v17, v18);
	double H[2][2] = { { -6 * v16 * v17,v18},{v18,(v16 *
	   (1224 * v1 - 32 * v0 * v10 - 32 * v0 * v13 - 32 * v0 * v15 - 288 * v10 * v2 * v3 + 288 * v11 * v2 * v3 +
		 288 * v12 * v2 * v3 - 288 * v13 * v2 * v3 + 288 * v14 * v2 * v3 - 288 * v15 * v2 * v3 + 16200 * v3 * v4 -
		 288 * v10 * v4 * v5 + 288 * v11 * v4 * v5 + 288 * v12 * v4 * v5 - 288 * v13 * v4 * v5 + 288 * v14 * v4 * v5 -
		 288 * v15 * v4 * v5 + 6075 * v6 - 256 * v0 * v7 - 3456 * v2 * v3 * v7 - 3456 * v4 * v5 * v7 - 256 * v0 * v8 -
		 3456 * v2 * v3 * v8 - 3456 * v4 * v5 * v8 - 256 * v0 * v9 - 3456 * v2 * v3 * v9 - 3456 * v4 * v5 * v9 -
		 144 * v1 * v10 * lam + 96 * v1 * v11 * lam + 96 * v1 * v12 * lam - 144 * v1 * v13 * lam + 96 * v1 * v14 * lam -
		 144 * v1 * v15 * lam + 7200 * v2 * lam - 1536 * v1 * v7 * lam - 1536 * v1 * v8 * lam - 1536 * v1 * v9 * lam +
		 16200 * v5 * mu - 108 * v10 * v6 * mu + 108 * v11 * v6 * mu + 108 * v12 * v6 * mu - 108 * v13 * v6 * mu +
		 108 * v14 * v6 * mu - 108 * v15 * v6 * mu - 1296 * v6 * v7 * mu - 1296 * v6 * v8 * mu - 1296 * v6 * v9 * mu)) /
	 Power(mu,4)} };
	return Eigen::Matrix2d::Map(H[0]).eval();
}

auto newton_setp_s(const Eigen::Matrix<double, 6, 6>& c, double& lam, double& mu) {
	double Rule(Power(mu, 2), v0), Rule(Power(lam, 2), v1), Rule(Subscript(c, 6, 6), v2),
		Rule(Subscript(c, 5, 5), v3), Rule(Subscript(c, 4, 4), v4), Rule(Subscript(c, 3, 3), v5),
		Rule(Subscript(c, 3, 2), v6), Rule(Subscript(c, 3, 1), v7), Rule(Subscript(c, 2, 2), v8),
		Rule(Subscript(c, 2, 1), v9), Rule(Subscript(c, 1, 1), v10), Rule(Power(v5, 2), v11),
		Rule(Power(v6, 2), v12), Rule(Power(v7, 2), v13), Rule(Power(v8, 2), v14), Rule(Power(v9, 2), v15),
		Rule(Power(v10, 2), v16), Rule(4 * v5 * mu, v17), Rule(4 * v8 * mu, v18), Rule(4 * v10 * mu, v19),
		Rule(1 / (-225 + v17 + v18 + v19 + 48 * v2 * mu + 48 * v3 * mu + 48 * v4 * mu - 4 * v6 * mu - 4 * v7 * mu - 4 * v9 * mu),
			v20);

	double dlamdmu[2] = { (v20 * (-96 * v0 * v10 - 675 * v1 * v10 + 48 * v0 * v2 + 48 * v0 * v3 + 48 * v0 * v4 - 96 * v0 * v5 -
		675 * v1 * v5 - 204 * v0 * v6 - 1350 * v1 * v6 - 204 * v0 * v7 - 1350 * v1 * v7 - 96 * v0 * v8 -
		675 * v1 * v8 - 204 * v0 * v9 - 1350 * v1 * v9 + 675 * lam + 8 * v0 * v11 * lam - 16 * v0 * v12 * lam -
		16 * v0 * v13 * lam + 8 * v0 * v14 * lam - 16 * v0 * v15 * lam + 8 * v0 * v16 * lam + 96 * v0 * v10 * v2 * lam +
		96 * v0 * v10 * v3 * lam + 96 * v0 * v10 * v4 * lam + 16 * v0 * v10 * v5 * lam + 96 * v0 * v2 * v5 * lam + 96 * v0 * v3 * v5 * lam +
		96 * v0 * v4 * v5 * lam + 8 * v0 * v10 * v6 * lam + 192 * v0 * v2 * v6 * lam + 192 * v0 * v3 * v6 * lam + 192 * v0 * v4 * v6 * lam +
		8 * v0 * v5 * v6 * lam + 8 * v0 * v10 * v7 * lam + 192 * v0 * v2 * v7 * lam + 192 * v0 * v3 * v7 * lam + 192 * v0 * v4 * v7 * lam +
		8 * v0 * v5 * v7 * lam - 32 * v0 * v6 * v7 * lam + 16 * v0 * v10 * v8 * lam + 96 * v0 * v2 * v8 * lam + 96 * v0 * v3 * v8 * lam +
		96 * v0 * v4 * v8 * lam + 16 * v0 * v5 * v8 * lam + 8 * v0 * v6 * v8 * lam + 8 * v0 * v7 * v8 * lam + 8 * v0 * v10 * v9 * lam +
		192 * v0 * v2 * v9 * lam + 192 * v0 * v3 * v9 * lam + 192 * v0 * v4 * v9 * lam + 8 * v0 * v5 * v9 * lam - 32 * v0 * v6 * v9 * lam -
		32 * v0 * v7 * v9 * lam + 8 * v0 * v8 * v9 * lam + 12 * v1 * v11 * mu - 24 * v1 * v12 * mu - 24 * v1 * v13 * mu +
		12 * v1 * v14 * mu - 24 * v1 * v15 * mu + 12 * v1 * v16 * mu + 144 * v1 * v10 * v2 * mu + 144 * v1 * v10 * v3 * mu +
		144 * v1 * v10 * v4 * mu + 24 * v1 * v10 * v5 * mu + 144 * v1 * v2 * v5 * mu + 144 * v1 * v3 * v5 * mu +
		144 * v1 * v4 * v5 * mu + 12 * v1 * v10 * v6 * mu + 288 * v1 * v2 * v6 * mu + 288 * v1 * v3 * v6 * mu +
		288 * v1 * v4 * v6 * mu + 12 * v1 * v5 * v6 * mu + 12 * v1 * v10 * v7 * mu + 288 * v1 * v2 * v7 * mu + 288 * v1 * v3 * v7 * mu +
		288 * v1 * v4 * v7 * mu + 12 * v1 * v5 * v7 * mu - 48 * v1 * v6 * v7 * mu + 24 * v1 * v10 * v8 * mu + 144 * v1 * v2 * v8 * mu +
		144 * v1 * v3 * v8 * mu + 144 * v1 * v4 * v8 * mu + 24 * v1 * v5 * v8 * mu + 12 * v1 * v6 * v8 * mu + 12 * v1 * v7 * v8 * mu +
		12 * v1 * v10 * v9 * mu + 288 * v1 * v2 * v9 * mu + 288 * v1 * v3 * v9 * mu + 288 * v1 * v4 * v9 * mu + 12 * v1 * v5 * v9 * mu -
		48 * v1 * v6 * v9 * mu - 48 * v1 * v7 * v9 * mu + 12 * v1 * v8 * v9 * mu - 612 * v10 * lam * mu - 144 * v2 * lam * mu -
		144 * v3 * lam * mu - 144 * v4 * lam * mu - 612 * v5 * lam * mu - 1188 * v6 * lam * mu - 1188 * v7 * lam * mu - 612 * v8 * lam * mu -
		1188 * v9 * lam * mu)) /
	(-9 + v17 + v18 + v19 + 6 * v10 * lam + 6 * v5 * lam + 12 * v6 * lam + 12 * v7 * lam + 6 * v8 * lam + 12 * v9 * lam +
	  8 * v6 * mu + 8 * v7 * mu + 8 * v9 * mu),v20 * mu *
	(-75 + 2 * v10 * mu + 24 * v2 * mu + 24 * v3 * mu + 24 * v4 * mu + 2 * v5 * mu - 2 * v6 * mu - 2 * v7 * mu + 2 * v8 * mu -
	  2 * v9 * mu) };
	return std::make_tuple(dlamdmu[0], dlamdmu[1]);
}

double s_err(const Eigen::Matrix<double, 6, 6>& c, double& lam, double& mu) {
	double Rule(Power(mu, 4), v0), Rule(Power(mu, 3), v1), Rule(Power(mu, 2), v2), Rule(Power(lam, 2), v3),
		Rule(Subscript(c, 6, 6), v4), Rule(Subscript(c, 5, 5), v5), Rule(Subscript(c, 4, 4), v6),
		Rule(Subscript(c, 3, 3), v7), Rule(Subscript(c, 3, 2), v8), Rule(Subscript(c, 3, 1), v9),
		Rule(Subscript(c, 2, 2), v10), Rule(Subscript(c, 2, 1), v11), Rule(Subscript(c, 1, 1), v12),
		Rule(Power(v4, 2), v13), Rule(Power(Subscript(c, 6, 5), 2), v14),
		Rule(Power(Subscript(c, 6, 4), 2), v15), Rule(Power(Subscript(c, 6, 3), 2), v16),
		Rule(Power(Subscript(c, 6, 2), 2), v17), Rule(Power(Subscript(c, 6, 1), 2), v18),
		Rule(Power(v5, 2), v19), Rule(Power(Subscript(c, 5, 4), 2), v20),
		Rule(Power(Subscript(c, 5, 3), 2), v21), Rule(Power(Subscript(c, 5, 2), 2), v22),
		Rule(Power(Subscript(c, 5, 1), 2), v23), Rule(Power(v6, 2), v24),
		Rule(Power(Subscript(c, 4, 3), 2), v25), Rule(Power(Subscript(c, 4, 2), 2), v26),
		Rule(Power(Subscript(c, 4, 1), 2), v27), Rule(Power(v7, 2), v28), Rule(Power(v8, 2), v29),
		Rule(Power(v9, 2), v30), Rule(Power(v10, 2), v31), Rule(Power(v11, 2), v32),
		Rule(Power(v12, 2), v33);
	
	double serr = (-8 * v1 * v10 - 8 * v1 * v12 + 32 * v0 * v13 + 64 * v0 * v14 + 64 * v0 * v15 + 32 * v0 * v16 + 32 * v0 * v17 +
		32 * v0 * v18 + 32 * v0 * v19 + 102 * v2 + 64 * v0 * v20 + 32 * v0 * v21 + 32 * v0 * v22 + 32 * v0 * v23 +
		32 * v0 * v24 + 32 * v0 * v25 + 32 * v0 * v26 + 32 * v0 * v27 + 8 * v0 * v28 + 16 * v0 * v29 + 225 * v3 +
		72 * v13 * v2 * v3 + 144 * v14 * v2 * v3 + 144 * v15 * v2 * v3 + 72 * v16 * v2 * v3 + 72 * v17 * v2 * v3 +
		72 * v18 * v2 * v3 + 72 * v19 * v2 * v3 + 144 * v2 * v20 * v3 + 72 * v2 * v21 * v3 + 72 * v2 * v22 * v3 +
		72 * v2 * v23 * v3 + 72 * v2 * v24 * v3 + 72 * v2 * v25 * v3 + 72 * v2 * v26 * v3 + 72 * v2 * v27 * v3 +
		18 * v2 * v28 * v3 + 36 * v2 * v29 * v3 + 16 * v0 * v30 + 36 * v2 * v3 * v30 + 8 * v0 * v31 + 18 * v2 * v3 * v31 +
		16 * v0 * v32 + 36 * v2 * v3 * v32 + 8 * v0 * v33 + 18 * v2 * v3 * v33 - 64 * v1 * v4 - 64 * v1 * v5 - 64 * v1 * v6 -
		8 * v1 * v7 + 96 * v1 * v13 * lam + 192 * v1 * v14 * lam + 192 * v1 * v15 * lam + 96 * v1 * v16 * lam + 96 * v1 * v17 * lam +
		96 * v1 * v18 * lam + 96 * v1 * v19 * lam - 20 * v10 * v2 * lam + 8 * v11 * v2 * lam - 20 * v12 * v2 * lam + 192 * v1 * v20 * lam +
		96 * v1 * v21 * lam + 96 * v1 * v22 * lam + 96 * v1 * v23 * lam + 96 * v1 * v24 * lam + 96 * v1 * v25 * lam + 96 * v1 * v26 * lam +
		96 * v1 * v27 * lam + 24 * v1 * v28 * lam + 48 * v1 * v29 * lam + 48 * v1 * v30 * lam + 24 * v1 * v31 * lam + 48 * v1 * v32 * lam +
		24 * v1 * v33 * lam - 192 * v2 * v4 * lam - 192 * v2 * v5 * lam - 192 * v2 * v6 * lam - 20 * v2 * v7 * lam + 8 * v2 * v8 * lam +
		8 * v2 * v9 * lam - 12 * v10 * v3 * mu + 12 * v11 * v3 * mu - 12 * v12 * v3 * mu - 144 * v3 * v4 * mu - 144 * v3 * v5 * mu -
		144 * v3 * v6 * mu - 12 * v3 * v7 * mu + 12 * v3 * v8 * mu + 12 * v3 * v9 * mu + 300 * lam * mu) /
		(2. * Power(mu, 2) * Power(3 * lam + 2 * mu, 2));
	return serr;
}

Eigen::Matrix<double, 6, 6> msf::mtl::isotropic_projection_s(const Eigen::Matrix<double, 6, 6>& S, double& lam, double& mu) {
	for (int iter = 0; iter < 10; iter++) {
		{double err = s_err(S, lam, mu); std::cout << "err = " << err << std::endl; }
		std::cout << "H =\n" << hessian_err(S, lam, mu);
		auto [dlam, dmu] = newton_setp_s(S, lam, mu);
		lam += dlam;
		mu += dmu;
	}
	return compliance_matrix(lam, mu);
}