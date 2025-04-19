#include <iostream>
#include <fstream>
#include "BezierTri.h"

using namespace msf;

void test_bezier_tri(void) {
	Point cpt[BezierTri<Real, 3>::n_cpts];

	auto z = [](Real x, Real y) {
		//return std::sqrt(std::sin(2 * 3.14 * (x - 1. / 3)) + std::pow(y - 1. / 3, 2));
		return (std::sin(2 * 3.14 * (x - 1. / 3)) + std::pow(y - 1. / 3, 2));
	};
	for (int i = 0; i <= 3; i++) {
		for (int j = 0; j <= 3 - i; j++) {
			int k = 3 - i - j;
			Real u = i / 3., v = j / 3.;
			cpt[pid(i, j, k)] = Eigen::Vector3<Real>(u, v, z(u, v));
		}
	}

	msf::BezierTri<Real, 3> tri(cpt);
	
	std::vector<Point> samples;
	std::vector<std::pair<Point, Point>> sampleTgs;
	for (Real u = 0; u < 1; u += 0.03) {
		for (Real v = 0; v < 1 - u; v += 0.03) {
			Real s = 1 - u - v;
			samples.emplace_back(tri.eval(u, v));
			auto [tu, tv] = tri.TgEval(u, v);
			sampleTgs.emplace_back(std::pair<Point, Point> { tu, tv });
		}
	}
	
	std::ofstream ofs("sampleBez");
	for (int i = 0; i < samples.size(); i++) {
		ofs << samples[i].transpose() << std::endl;
	}
	ofs.close();

	ofs.open("sampleTg");
	for (int i = 0; i < sampleTgs.size(); i++) {
		ofs << sampleTgs[i].first.transpose() << " " << sampleTgs[i].second.transpose() << std::endl;
	}
	ofs.close();
}

void test_bezier_hessian(void) {
	Point cpt[BezierTri<Real, 3>::n_cpts];

	auto z = [](Real x, Real y) {
		return (std::sin(2 * 3.14 * (x - 1. / 3)) + std::pow(y - 1. / 3, 2));
	};
	for (int i = 0; i <= 3; i++) {
		for (int j = 0; j <= 3 - i; j++) {
			int k = 3 - i - j;
			Real u = i / 3., v = j / 3.;
			cpt[pid(i, j, k)] = Eigen::Vector3<Real>(u, v, z(u, v));
		}
	}
	msf::BezierTri<Real, 3> tri(cpt);

	auto dr2 = tri.Hessian(0.5, 0.5);
}
