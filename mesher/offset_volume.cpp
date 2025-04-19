#include "Config.h"
#include <vector>
#include "igl/read_triangle_mesh.h"
#include "igl/per_vertex_normals.h"	

using namespace msf;

inline Real Power(Real x, int p) { return (std::pow)(x, p); }

inline Real getElement(const Eigen::Matrix3<Real>& m, int row, int col) {
	return m(row - 1, col - 1);
}

// one row for one verte
Real tri_shell_volume(const Eigen::Matrix3<Real>& p, const Eigen::Matrix3<Real>& n, Real eps) {
	(-(Power(eps, 3) * getElement(n, 1, 3) * getElement(n, 2, 2) * getElement(n, 3, 1)) +
		Power(eps, 3) * getElement(n, 1, 2) * getElement(n, 2, 3) * getElement(n, 3, 1) +
		Power(eps, 3) * getElement(n, 1, 3) * getElement(n, 2, 1) * getElement(n, 3, 2) -
		Power(eps, 3) * getElement(n, 1, 1) * getElement(n, 2, 3) * getElement(n, 3, 2) -
		Power(eps, 3) * getElement(n, 1, 2) * getElement(n, 2, 1) * getElement(n, 3, 3) +
		Power(eps, 3) * getElement(n, 1, 1) * getElement(n, 2, 2) * getElement(n, 3, 3) -
		eps * getElement(n, 1, 3) * getElement(p, 1, 2) * getElement(p, 2, 1) -
		eps * getElement(n, 2, 3) * getElement(p, 1, 2) * getElement(p, 2, 1) -
		eps * getElement(n, 3, 3) * getElement(p, 1, 2) * getElement(p, 2, 1) +
		eps * getElement(n, 1, 2) * getElement(p, 1, 3) * getElement(p, 2, 1) +
		eps * getElement(n, 2, 2) * getElement(p, 1, 3) * getElement(p, 2, 1) +
		eps * getElement(n, 3, 2) * getElement(p, 1, 3) * getElement(p, 2, 1) +
		eps * getElement(n, 1, 3) * getElement(p, 1, 1) * getElement(p, 2, 2) +
		eps * getElement(n, 2, 3) * getElement(p, 1, 1) * getElement(p, 2, 2) +
		eps * getElement(n, 3, 3) * getElement(p, 1, 1) * getElement(p, 2, 2) -
		eps * getElement(n, 1, 1) * getElement(p, 1, 3) * getElement(p, 2, 2) -
		eps * getElement(n, 2, 1) * getElement(p, 1, 3) * getElement(p, 2, 2) -
		eps * getElement(n, 3, 1) * getElement(p, 1, 3) * getElement(p, 2, 2) -
		eps * getElement(n, 1, 2) * getElement(p, 1, 1) * getElement(p, 2, 3) -
		eps * getElement(n, 2, 2) * getElement(p, 1, 1) * getElement(p, 2, 3) -
		eps * getElement(n, 3, 2) * getElement(p, 1, 1) * getElement(p, 2, 3) +
		eps * getElement(n, 1, 1) * getElement(p, 1, 2) * getElement(p, 2, 3) +
		eps * getElement(n, 2, 1) * getElement(p, 1, 2) * getElement(p, 2, 3) +
		eps * getElement(n, 3, 1) * getElement(p, 1, 2) * getElement(p, 2, 3) +
		eps * getElement(n, 1, 3) * getElement(p, 1, 2) * getElement(p, 3, 1) +
		eps * getElement(n, 2, 3) * getElement(p, 1, 2) * getElement(p, 3, 1) +
		eps * getElement(n, 3, 3) * getElement(p, 1, 2) * getElement(p, 3, 1) -
		eps * getElement(n, 1, 2) * getElement(p, 1, 3) * getElement(p, 3, 1) -
		eps * getElement(n, 2, 2) * getElement(p, 1, 3) * getElement(p, 3, 1) -
		eps * getElement(n, 3, 2) * getElement(p, 1, 3) * getElement(p, 3, 1) -
		eps * getElement(n, 1, 3) * getElement(p, 2, 2) * getElement(p, 3, 1) -
		eps * getElement(n, 2, 3) * getElement(p, 2, 2) * getElement(p, 3, 1) -
		eps * getElement(n, 3, 3) * getElement(p, 2, 2) * getElement(p, 3, 1) +
		eps * getElement(n, 1, 2) * getElement(p, 2, 3) * getElement(p, 3, 1) +
		eps * getElement(n, 2, 2) * getElement(p, 2, 3) * getElement(p, 3, 1) +
		eps * getElement(n, 3, 2) * getElement(p, 2, 3) * getElement(p, 3, 1) -
		eps * getElement(n, 1, 3) * getElement(p, 1, 1) * getElement(p, 3, 2) -
		eps * getElement(n, 2, 3) * getElement(p, 1, 1) * getElement(p, 3, 2) -
		eps * getElement(n, 3, 3) * getElement(p, 1, 1) * getElement(p, 3, 2) +
		eps * getElement(n, 1, 1) * getElement(p, 1, 3) * getElement(p, 3, 2) +
		eps * getElement(n, 2, 1) * getElement(p, 1, 3) * getElement(p, 3, 2) +
		eps * getElement(n, 3, 1) * getElement(p, 1, 3) * getElement(p, 3, 2) +
		eps * getElement(n, 1, 3) * getElement(p, 2, 1) * getElement(p, 3, 2) +
		eps * getElement(n, 2, 3) * getElement(p, 2, 1) * getElement(p, 3, 2) +
		eps * getElement(n, 3, 3) * getElement(p, 2, 1) * getElement(p, 3, 2) -
		eps * getElement(n, 1, 1) * getElement(p, 2, 3) * getElement(p, 3, 2) -
		eps * getElement(n, 2, 1) * getElement(p, 2, 3) * getElement(p, 3, 2) -
		eps * getElement(n, 3, 1) * getElement(p, 2, 3) * getElement(p, 3, 2) +
		eps * getElement(n, 1, 2) * getElement(p, 1, 1) * getElement(p, 3, 3) +
		eps * getElement(n, 2, 2) * getElement(p, 1, 1) * getElement(p, 3, 3) +
		eps * getElement(n, 3, 2) * getElement(p, 1, 1) * getElement(p, 3, 3) -
		eps * getElement(n, 1, 1) * getElement(p, 1, 2) * getElement(p, 3, 3) -
		eps * getElement(n, 2, 1) * getElement(p, 1, 2) * getElement(p, 3, 3) -
		eps * getElement(n, 3, 1) * getElement(p, 1, 2) * getElement(p, 3, 3) -
		eps * getElement(n, 1, 2) * getElement(p, 2, 1) * getElement(p, 3, 3) -
		eps * getElement(n, 2, 2) * getElement(p, 2, 1) * getElement(p, 3, 3) -
		eps * getElement(n, 3, 2) * getElement(p, 2, 1) * getElement(p, 3, 3) +
		eps * getElement(n, 1, 1) * getElement(p, 2, 2) * getElement(p, 3, 3) +
		eps * getElement(n, 2, 1) * getElement(p, 2, 2) * getElement(p, 3, 3) +
		eps * getElement(n, 3, 1) * getElement(p, 2, 2) * getElement(p, 3, 3)) / 3;
}

void tri_shell_volume(const Eigen::Matrix3<Real>& p, const Eigen::Matrix3<Real>& n, Real& c_eps3, Real& c_eps) {
	c_eps3 = ((-(getElement(n, 1, 3) * getElement(n, 2, 2) * getElement(n, 3, 1)) +
		getElement(n, 1, 2) * getElement(n, 2, 3) * getElement(n, 3, 1) +
		getElement(n, 1, 3) * getElement(n, 2, 1) * getElement(n, 3, 2) -
		getElement(n, 1, 1) * getElement(n, 2, 3) * getElement(n, 3, 2) -
		getElement(n, 1, 2) * getElement(n, 2, 1) * getElement(n, 3, 3) +
		getElement(n, 1, 1) * getElement(n, 2, 2) * getElement(n, 3, 3))) / 3.;

	c_eps = (-(getElement(n, 1, 3) * getElement(p, 1, 2) * getElement(p, 2, 1)) -
		getElement(n, 2, 3) * getElement(p, 1, 2) * getElement(p, 2, 1) -
		getElement(n, 3, 3) * getElement(p, 1, 2) * getElement(p, 2, 1) +
		getElement(n, 1, 2) * getElement(p, 1, 3) * getElement(p, 2, 1) +
		getElement(n, 2, 2) * getElement(p, 1, 3) * getElement(p, 2, 1) +
		getElement(n, 3, 2) * getElement(p, 1, 3) * getElement(p, 2, 1) +
		getElement(n, 1, 3) * getElement(p, 1, 1) * getElement(p, 2, 2) +
		getElement(n, 2, 3) * getElement(p, 1, 1) * getElement(p, 2, 2) +
		getElement(n, 3, 3) * getElement(p, 1, 1) * getElement(p, 2, 2) -
		getElement(n, 1, 1) * getElement(p, 1, 3) * getElement(p, 2, 2) -
		getElement(n, 2, 1) * getElement(p, 1, 3) * getElement(p, 2, 2) -
		getElement(n, 3, 1) * getElement(p, 1, 3) * getElement(p, 2, 2) -
		getElement(n, 1, 2) * getElement(p, 1, 1) * getElement(p, 2, 3) -
		getElement(n, 2, 2) * getElement(p, 1, 1) * getElement(p, 2, 3) -
		getElement(n, 3, 2) * getElement(p, 1, 1) * getElement(p, 2, 3) +
		getElement(n, 1, 1) * getElement(p, 1, 2) * getElement(p, 2, 3) +
		getElement(n, 2, 1) * getElement(p, 1, 2) * getElement(p, 2, 3) +
		getElement(n, 3, 1) * getElement(p, 1, 2) * getElement(p, 2, 3) +
		getElement(n, 1, 3) * getElement(p, 1, 2) * getElement(p, 3, 1) +
		getElement(n, 2, 3) * getElement(p, 1, 2) * getElement(p, 3, 1) +
		getElement(n, 3, 3) * getElement(p, 1, 2) * getElement(p, 3, 1) -
		getElement(n, 1, 2) * getElement(p, 1, 3) * getElement(p, 3, 1) -
		getElement(n, 2, 2) * getElement(p, 1, 3) * getElement(p, 3, 1) -
		getElement(n, 3, 2) * getElement(p, 1, 3) * getElement(p, 3, 1) -
		getElement(n, 1, 3) * getElement(p, 2, 2) * getElement(p, 3, 1) -
		getElement(n, 2, 3) * getElement(p, 2, 2) * getElement(p, 3, 1) -
		getElement(n, 3, 3) * getElement(p, 2, 2) * getElement(p, 3, 1) +
		getElement(n, 1, 2) * getElement(p, 2, 3) * getElement(p, 3, 1) +
		getElement(n, 2, 2) * getElement(p, 2, 3) * getElement(p, 3, 1) +
		getElement(n, 3, 2) * getElement(p, 2, 3) * getElement(p, 3, 1) -
		getElement(n, 1, 3) * getElement(p, 1, 1) * getElement(p, 3, 2) -
		getElement(n, 2, 3) * getElement(p, 1, 1) * getElement(p, 3, 2) -
		getElement(n, 3, 3) * getElement(p, 1, 1) * getElement(p, 3, 2) +
		getElement(n, 1, 1) * getElement(p, 1, 3) * getElement(p, 3, 2) +
		getElement(n, 2, 1) * getElement(p, 1, 3) * getElement(p, 3, 2) +
		getElement(n, 3, 1) * getElement(p, 1, 3) * getElement(p, 3, 2) +
		getElement(n, 1, 3) * getElement(p, 2, 1) * getElement(p, 3, 2) +
		getElement(n, 2, 3) * getElement(p, 2, 1) * getElement(p, 3, 2) +
		getElement(n, 3, 3) * getElement(p, 2, 1) * getElement(p, 3, 2) -
		getElement(n, 1, 1) * getElement(p, 2, 3) * getElement(p, 3, 2) -
		getElement(n, 2, 1) * getElement(p, 2, 3) * getElement(p, 3, 2) -
		getElement(n, 3, 1) * getElement(p, 2, 3) * getElement(p, 3, 2) +
		getElement(n, 1, 2) * getElement(p, 1, 1) * getElement(p, 3, 3) +
		getElement(n, 2, 2) * getElement(p, 1, 1) * getElement(p, 3, 3) +
		getElement(n, 3, 2) * getElement(p, 1, 1) * getElement(p, 3, 3) -
		getElement(n, 1, 1) * getElement(p, 1, 2) * getElement(p, 3, 3) -
		getElement(n, 2, 1) * getElement(p, 1, 2) * getElement(p, 3, 3) -
		getElement(n, 3, 1) * getElement(p, 1, 2) * getElement(p, 3, 3) -
		getElement(n, 1, 2) * getElement(p, 2, 1) * getElement(p, 3, 3) -
		getElement(n, 2, 2) * getElement(p, 2, 1) * getElement(p, 3, 3) -
		getElement(n, 3, 2) * getElement(p, 2, 1) * getElement(p, 3, 3) +
		getElement(n, 1, 1) * getElement(p, 2, 2) * getElement(p, 3, 3) +
		getElement(n, 2, 1) * getElement(p, 2, 2) * getElement(p, 3, 3) +
		getElement(n, 3, 1) * getElement(p, 2, 2) * getElement(p, 3, 3)) / 3;
}

std::vector<Real> offset_volume(std::string mfile, const std::vector<Real>& half_thick) {
	Eigen::MatrixX3<Real> V;
	Eigen::MatrixX3i F;
	igl::read_triangle_mesh(mfile, V, F);
	Eigen::MatrixX3<Real> N;
	igl::per_vertex_normals(V, F, N);
	std::vector<Real> c_eps3, c_eps;
	for (int i = 0; i < F.rows(); i++) {
		Eigen::Matrix3<Real> Fv;
		Eigen::Matrix3<Real> Nv;
		for (int k = 0; k < 3; k++) {
			Fv.row(k) = V.row(F(i, k));
			Nv.row(k) = N.row(F(i, k));
		}
		Real eps3, eps;
		tri_shell_volume(Fv, Nv, eps3, eps);
		c_eps3.push_back(eps3);
		c_eps.push_back(eps);
	}
	std::vector<Real> volist;
	for (int i = 0; i < half_thick.size(); i++) {
		Real t = half_thick[i];
		Real t3 = (std::pow)(t, 3);
		Real vol = 0;
		for (int k = 0; k < c_eps3.size(); k++) {
			vol += c_eps3[k] * t3 + c_eps[k] * t;
		}
		volist.push_back(vol);
	}
	for (int i = 0; i < volist.size(); i++) {
		std::cout << volist[i] << " ";
	}
	std::cout << std::endl;
	return volist;
}

