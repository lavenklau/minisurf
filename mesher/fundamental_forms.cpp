#include <iostream>
#include "fundamental_forms.h"
#include <Eigen/Eigen>
#include <vector>
#include "flow/triangle_quadrature.h"
#include "PeriodicMesher.h"

#define _USE_MATH_DEFINES
#include "math.h"

//#define _USE_VERTEX_MACRO_STRAIN_PROJECTION
#define _USE_VERTEX_DISPLACEMENT_PROJECTION
#define _USE_SECOND_FUND_FORM


void Compile1ring::compile1ring(const Eigen::Vector3d& o, const std::vector<Eigen::Vector3d>& ring)
{
	// mean curvature vector
	Eigen::Vector3d Hv(0, 0, 0);

	double theta_sum = 0; // angle-based summation of normal

	std::vector<double> esq_inc(ring.size());
	std::vector<double> esq_opp(ring.size());
	std::vector<double> cot_alpha(ring.size());
	std::vector<double> cot_beta(ring.size());
	std::vector<double> theta(ring.size());
	for (int i = 0; i < ring.size(); i++) {
		esq_inc[i] = (ring[i] - o).squaredNorm();
		esq_opp[i] = (ring[(i + 1) % ring.size()] - ring[i]).squaredNorm();
	}
	for (int i = 0; i < ring.size(); i++) {
		double e_sq[] = { esq_inc[i], esq_inc[(i + 1) % ring.size()], esq_opp[i] };
		double cosA = (e_sq[0] + e_sq[2] - e_sq[1]) / 2 / std::sqrt(e_sq[0] * e_sq[2]);
		double cosB = (e_sq[1] + e_sq[2] - e_sq[0]) / 2 / std::sqrt(e_sq[1] * e_sq[2]);
		double cotA = cosA / std::sqrt(1 - cosA * cosA);
		double cotB = cosB / std::sqrt(1 - cosB * cosB);
		cot_alpha[(i + 1) % ring.size()] = cotA;
		cot_beta[i] = cotB;

		double theta_i = std::acos((e_sq[0] + e_sq[1] - e_sq[2]) / 2 / std::sqrt(e_sq[0] * e_sq[1]));
		theta[i] = theta_i;
	}


	for (int i = 0; i < ring.size(); i++) {
		Eigen::Vector3<double> a = ring[i] - o;
		Eigen::Vector3<double> b = ring[(i + 1) % ring.size()] - o;
		double e_sq[] = { esq_inc[i], esq_inc[(i + 1) % ring.size()], esq_opp[i] };
		double theta_i = theta[i];
		theta_sum += theta_i;
		Eigen::Vector3d axb = a.cross(b);
		mass += axb.norm() / 2 / 3;
		nv += axb.normalized() * theta_i;
		Hv += (cot_alpha[i] + cot_beta[i]) / 2 * a;
		double cotA = cot_alpha[(i + 1) % ring.size()];
		double cotB = cot_beta[i];

		if (e_sq[0] + e_sq[1] >= e_sq[2] && e_sq[1] + e_sq[2] >= e_sq[0] && e_sq[2] + e_sq[0] >= e_sq[1]) {
			// add Voronoi area
			double A = cotA * e_sq[1] / 8 + cotB * e_sq[0] / 8;
			As += A;
		}
		else if (e_sq[0] + e_sq[1] < e_sq[2]) {
			// half area if incident angle is obtuse
			double A = axb.norm() / 4;
			As += A;
		}
		else /* has obtuse angle but not the incident one */ {
			// quater area if others
			double A = axb.norm() / 8;
			As += A;
		}
	}

	nv.normalize();

	K = (2 * M_PI - theta_sum) / As;

	Lx = Hv;

	H = Hv.norm() / As / 2;

	if (Hv.dot(nv) < 0) H *= -1;
}

//void Compile1ring::compileHalfring(const Eigen::Vector3d& o, const std::vector<Eigen::Vector3d>& ring)
//{
//	
//}

inline double dot2(const Eigen::Matrix2d& m1, const Eigen::Matrix2d& m2) { return m1.cwiseProduct(m2).sum(); }

std::tuple<double, double, double> fundamental_forms(const Eigen::Vector3d& v12, const Compile1ring& v1, const Compile1ring& v2)
{
	double I = v12.squaredNorm();
	double III = std::pow(std::acos(v1.nv.dot(v2.nv)), 2);
	double Kavg = (v1.K + v2.K) / 2;
	double Havg = (v1.H + v2.H);
	double II = (III + Kavg * I) /*/ Havg*/;
	return { I, II, III };
}
 
// from stretch of edges to the strain of face
// {eps11, eps22, 2 eps12} = B .{l1, l2, l3}
Eigen::Matrix<double, 3, 3> strain_matrix_edge_stretch(const Eigen::Matrix3d& tri, const Eigen::Vector3d& e1, const Eigen::Vector3d& e2) {
	Eigen::Vector2d v12;
	Eigen::Matrix<double, 3, 2> fram;
	fram << e1, e2;
	Eigen::Matrix3d tri1;
	tri1 << tri.col(1), tri.col(2), tri.col(0);
	Eigen::Matrix<double, 2, 3> d = fram.transpose() * (tri1 - tri);
	Eigen::Matrix3d A;
	A <<
		d(0, 0)* d(0, 0), d(1, 0)* d(1, 0), d(0, 0)* d(1, 0),
		d(0, 1)* d(0, 1), d(1, 1)* d(1, 1), d(0, 1)* d(1, 1),
		d(0, 2)* d(0, 2), d(1, 2)* d(1, 2), d(0, 2)* d(1, 2);
	Eigen::Matrix3d B = A.inverse();
	return B;
}

// from offset of edges to strain of face
// {eps11, eps22, 2 eps12} = B.{u12, u23, u30}
Eigen::Matrix<double, 3, 6> strain_matrix_edge_offset(const Eigen::Matrix3d& tri, const Eigen::Vector3d& e1, const Eigen::Vector3d& e2) {
	Eigen::Matrix2d eij;
	eij <<
		(tri.col(0) - tri.col(2)).dot(e1), (tri.col(1) - tri.col(2)).dot(e1),
		(tri.col(0) - tri.col(2)).dot(e2), (tri.col(1) - tri.col(2)).dot(e2);
	Eigen::Matrix3d ugrad;
	Eigen::Matrix<double, 3, 2> C;
	C << 1, 1,
		-1, 0,
		0, -1;
	Eigen::Matrix<double, 3, 2> S = C * eij.inverse();
	Eigen::Matrix<double, 3, 6> B;
	B <<
		S(0, 0), 0, S(1, 0), 0, S(2, 0), 0,
		0, S(0, 1), 0, S(1, 1), 0, S(2, 1),
		S(0, 1), S(0, 0), S(1, 1), S(1, 0), S(2, 1), S(2, 0);
	return B;
}

Eigen::Vector3d second_fundamental_form(const Eigen::Matrix3d& tri, const Eigen::Matrix3d& nv) {
	double be12, be23, be31;
	be12 = -(nv.col(1) - nv.col(0)).dot(tri.col(1) - tri.col(0));
	be23 = -(nv.col(2) - nv.col(1)).dot(tri.col(2) - tri.col(1));
	be31 = -(nv.col(0) - nv.col(2)).dot(tri.col(0) - tri.col(2));
	return { be12, be23, be31 };
}

Eigen::Vector3d second_fundamental_form(const Eigen::Matrix3d& tri, const Compile1ring& v1, const Compile1ring& v2, const Compile1ring& v3) {
#if 0
	double be12, be23, be31;
	double x1, x2;
	std::tie(x1, be12, x2) = fundamental_forms(tri.col(1) - tri.col(0), v1, v2);
	std::tie(x1, be23, x2) = fundamental_forms(tri.col(2) - tri.col(1), v2, v3);
	std::tie(x1, be31, x2) = fundamental_forms(tri.col(0) - tri.col(2), v3, v1);
	return { be12, be23, be31 };
#else
	double be12, be23, be31;
	be12 = -(v2.nv - v1.nv).dot(tri.col(1) - tri.col(0));
	be23 = -(v3.nv - v2.nv).dot(tri.col(2) - tri.col(1));
	be31 = -(v1.nv - v3.nv).dot(tri.col(0) - tri.col(2));
	return { be12, be23, be31 };
#endif
}


Eigen::Matrix3d second_fundamental_form(const Eigen::Matrix3d& tri, const Eigen::Matrix3d& fr, const Compile1ring v[3])
{
	auto be = second_fundamental_form(tri, v[0], v[1], v[2]);
	auto Bn = strain_matrix_edge_stretch(tri, fr.col(0), fr.col(1));
	return fr.leftCols(2) * fromvoigt((Bn * be).eval()) * fr.leftCols(2).transpose();
}

Eigen::Matrix3d second_fundamental_form(const Eigen::Matrix3d& tri, const Eigen::Matrix3d& fr, const Eigen::Matrix3d& nv)
{
	auto be = second_fundamental_form(tri, nv);
	auto Bn = strain_matrix_edge_stretch(tri, fr.col(0), fr.col(1));
	return fr.leftCols(2) * fromvoigt((Bn * be).eval()) * fr.leftCols(2).transpose();
}

Eigen::Matrix2d second_fundamental_form_2d(const Eigen::Matrix3d& tri, const Eigen::Matrix3d& fr, const Eigen::Matrix3d& nv)
{
	auto be = second_fundamental_form(tri, nv);
	auto Bn = strain_matrix_edge_stretch(tri, fr.col(0), fr.col(1));
	return  fromvoigt((Bn * be).eval());
}

std::tuple<
	Eigen::Matrix<double, 9, 9>, Eigen::Vector<double, 9>,
	Eigen::Vector3d, Eigen::Vector3d
> membrane_element_stif_matrix_vector(
	double lam0, double mu,
	const Eigen::Matrix3d& tri, const Compile1ring v[3],
	const Eigen::Matrix3d& epsm
) {
	//Eigen::Matrix3d Kt, Kn;
	Eigen::Vector3d e1 = (tri.col(1) - tri.col(0)).normalized();
	Eigen::Vector3d e2 = (tri.col(1) - tri.col(0)).cross(tri.col(2) - tri.col(0)).cross(e1).normalized();
	auto Bn = strain_matrix_edge_stretch(tri, e1, e2);

	auto be = second_fundamental_form(tri, v[0], v[1], v[2]);

	Eigen::Matrix3d D;
	D << lam0 + 2 * mu, lam0, 0,
		lam0, lam0 + 2 * mu, 0,
		0, 0, mu;

	Eigen::Matrix<double, 3, 9> B;

	auto Bt = strain_matrix_edge_offset(tri, e1, e2);

	B.block<3, 6>(0, 3) = Bt;

	static auto qlist = qd::tri::symq::get(2);

	Eigen::Matrix<double, 3, 2> P; P << e1, e2;
	Eigen::Matrix2d pEps = P.transpose() * epsm * P;
	Eigen::Vector3d pm(pEps(0, 0), pEps(1, 1), 2 * pEps(0, 1));

	Eigen::Vector<double, 9> f; f.setZero();
	Eigen::Matrix<double, 9, 9> K; K.setZero();
	double A = (tri.col(1) - tri.col(0)).cross(tri.col(2) - tri.col(0)).norm() / 2;

	for (int i = 0; i < qlist.size(); i++) {
		Eigen::Vector3d q(qlist[i].c[0], qlist[i].c[1], qlist[i].c[2]);
		double w = qlist[i].w;

		B.block<3, 3>(0, 0) = -Bn * be * q.transpose();

		K += B.transpose() * (D * w) * B;

		f += -B.transpose() * (D * w) * pm;
	}
	K *= A;
	f *= A;

	return std::make_tuple(K, f, e1, e2);
}

Eigen::Matrix3d rotateAlign(const Eigen::Vector3d& vfrom, const Eigen::Vector3d& vto) {
	//Eigen::Vector3d v = vfrom.cross(vto);
	Eigen::Vector3d v = vfrom.cross(vto - vfrom);
	double s = v.squaredNorm();
	double c = vfrom.dot(vto);
	Eigen::Matrix3d vx;
	vx <<
		0, -v[2], v[1],
		v[2], 0, -v[0],
		-v[1], v[0], 0;
	Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
	if (s < 1e-28) { return R; }
#if 0
	R += vx + vx * vx * (1 - c) / s;
#else
	R += vx + vx * vx * 1 / (1 + c);
#endif
	return R;
}

Eigen::Matrix<double, 4, 9> tgnt_displacement_grad_matrix(
	const Eigen::Matrix3d& tri, const Compile1ring v[3],
	const Eigen::Vector3d& n, const Eigen::Matrix<double, 3, 2>& P, bool vertex_proj = true
) {
	Eigen::Matrix2d V;
	V << P.transpose() * (tri.col(1) - tri.col(0)), P.transpose()* (tri.col(2) - tri.col(0));
	V = V.inverse().eval();
	Eigen::Matrix<double, 3, 2> S;
	S << -1, -1, 1, 0, 0, 1;
	S.applyOnTheRight(V);

	Eigen::Matrix3d R[3];
	for (int i = 0; i < 3; i++) {
		R[i] = rotateAlign(v[i].nv, n);
	}
	
	Eigen::Matrix<double, 4, 9> dudx; dudx.setZero();
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {
					dudx(i + j * 2, k + l * 3) += P(k, i) * S(l, j);
				}
			}
		}
	}

#ifdef _USE_VERTEX_DISPLACEMENT_PROJECTION
	if (vertex_proj) {
		for (int k = 0; k < 3; k++) {
			dudx.middleCols(k * 3, 3).applyOnTheRight(R[k]);
		}
	}
#endif

	return dudx;
}

Eigen::Matrix<double, 3, 9>  tgnt_membrane_strain_displacement_matrix(const Eigen::Matrix3d& tri, const Compile1ring v[3], const Eigen::Vector3d& n, const Eigen::Matrix<double, 3, 2>& P) {
	auto dudx = tgnt_displacement_grad_matrix(tri, v, n, P);
	Eigen::Matrix<double, 3, 9> B0;
	// inplane strain part
	B0 << dudx.row(0), dudx.row(3), dudx.row(1) + dudx.row(2);
	return B0;
}

Eigen::Matrix<double, 3, 9> tgnt_membrane_strain_displacement_matrix_novp(const Eigen::Matrix3d& tri, const Compile1ring v[3], const Eigen::Vector3d& n, const Eigen::Matrix<double, 3, 2>& P)
{
	auto dudx = tgnt_displacement_grad_matrix(tri, v, n, P, false);
	Eigen::Matrix<double, 3, 9> B0;
	// inplane strain part
	B0 << dudx.row(0), dudx.row(3), dudx.row(1) + dudx.row(2);
	return B0;
}

Eigen::Matrix<double, 3, 2> vertex_projection(const Eigen::Matrix<double, 3, 2>& Pf, const Eigen::Vector3d& n, const Compile1ring& v) {
	auto R = rotateAlign(v.nv, n);
	return R.transpose() * Pf;
}

Eigen::Matrix3d face_frame(const Eigen::Matrix3d& tri) {
	Eigen::Vector3d n = (tri.col(1) - tri.col(0)).cross(tri.col(2) - tri.col(0)) / 2;

	Eigen::Vector3d e1 = (tri.col(1) - tri.col(0)).normalized();
	Eigen::Vector3d e2 = n.cross(e1).normalized();
	Eigen::Matrix3d fr;
	fr << e1, e2, n;
	return fr;
}

void project_macro_strain(const Eigen::Matrix3d& fr, const Compile1ring v[3], Eigen::Matrix<double, 3, 6>& pmf, Eigen::Matrix<double, 3, 6> pmvlist[3]) {
	Eigen::Vector3d n = fr.col(2).normalized();

	Eigen::Matrix<double, 3, 2> Pv[3]; 

	for (int k = 0; k < 3; k++) {
		auto R = rotateAlign(v[k].nv, n);
		Pv[k] = R.transpose() * fr.leftCols(2);
	}

	for (int i = 0; i < 3; i++) {
		// face projection
		Eigen::Matrix3d epsm; epsm.setZero();
		epsm(i, i) = 1;
		Eigen::Matrix2d peps = fr.leftCols(2).transpose() * epsm * fr.leftCols(2);
		Eigen::Vector3d pm(peps(0, 0), peps(1, 1), 2 * peps(0, 1));
		pmf.col(i) = pm;
		// vertex projection
		for (int k = 0; k < 3; k++) {
			pmvlist[k].col(i) = voigt((Pv[k].transpose() * epsm * Pv[k]).eval());
		}

		// face projection
		epsm.setZero();
		int j = (i + 1) % 3, k = (i + 2) % 3;
		epsm(j, k) = epsm(k, j) = 0.5;
		peps = fr.leftCols(2).transpose() * epsm * fr.leftCols(2);
		pm = Eigen::Vector3d(peps(0, 0), peps(1, 1), 2 * peps(0, 1));
		pmf.col(i + 3) = pm;
		// vertex projection
		for (int k = 0; k < 3; k++) {
			pmvlist[k].col(i + 3) = voigt((Pv[k].transpose() * epsm * Pv[k]).eval());
		}
	}

}

std::tuple<
	Eigen::Matrix<double, 9, 9>, Eigen::Matrix<double, 9, 6>
> membrane_element_vertex_stif_matrix_vector(
	double lam0, double mu,
	const Eigen::Matrix3d& tri, const Eigen::Matrix3d& fr, const Compile1ring v[3]
) {
	Eigen::Vector3d n = fr.col(2);
	double A = n.norm();
	n /= A;

	Eigen::Vector3d e1 = fr.col(0);
	Eigen::Vector3d e2 = fr.col(1);

	Eigen::Matrix<double, 3, 2> P; P << e1, e2;
	Eigen::Matrix<double, 3, 6> pmf;
	Eigen::Matrix<double, 3, 6> pmvlist[3];

	project_macro_strain(fr, v, pmf, pmvlist);
	
	
	Eigen::Matrix<double, 9, 9> K; K.setZero();
	Eigen::Matrix<double, 9, 6> f; f.setZero();

	Eigen::Matrix3d D;
	D <<
		lam0 + 2 * mu, lam0, 0,
		lam0, lam0 + 2 * mu, 0,
		0, 0, mu;


	Eigen::Matrix<double, 3, 9> B0 = tgnt_membrane_strain_displacement_matrix(tri, v,n, P);
	// inplane strain part
	
	static auto qlist = qd::tri::symq::get(2);

	auto Bn = strain_matrix_edge_stretch(tri, e1, e2);
	auto be = second_fundamental_form(tri, v[0], v[1], v[2]);

	for (int i = 0; i < qlist.size(); i++) {
		Eigen::Vector3d q(qlist[i].c[0], qlist[i].c[1], qlist[i].c[2]);
		double w = qlist[i].w;
		auto B = B0;

		Eigen::Vector<double, 9> qn;
#if 0
		for (int j = 0; j < 3; j++) { qn.block<3, 1>(j * 3, 0) = q[j] * n; }
#else
		for (int j = 0; j < 3; j++) { qn.block<3, 1>(j * 3, 0) = q[j] * v[j].nv; }
#endif

#ifdef _USE_SECOND_FUND_FORM
		B -= Bn * be * qn.transpose();
#endif

		K += B.transpose() * (D * w) * B;

#ifndef _USE_VERTEX_MACRO_STRAIN_PROJECTION
		f += -B.transpose() * (D * w) * pmf;
#else
		f += -B.transpose() * (D * w) * (pmvlist[0] * q[0] + pmvlist[1] * q[1] + pmvlist[2] * q[2]);
#endif
	}
	K *= A; f *= A;
	return std::make_tuple(K, f);
}

std::tuple<
	Eigen::Matrix<double, 9, 9>, Eigen::Matrix<double, 9, 6>
> membrane_element_vertex_stif_matrix_vector(
	double lam0, double mu,
	const Eigen::Matrix3d& tri, const Eigen::Matrix3d& fr, const Compile1ring v[3], const Eigen::Matrix3d bv[3]
) {
	Eigen::Vector3d n = fr.col(2);
	double A = n.norm();
	n /= A;

	Eigen::Vector3d e1 = fr.col(0), e2 = fr.col(1);

	Eigen::Matrix<double, 3, 2> P; P << e1, e2;
	Eigen::Matrix<double, 3, 6> pmf;
	Eigen::Matrix<double, 3, 6> pmvlist[3];
	Eigen::Matrix2d bv_e[3];

	for (int k = 0; k < 3; k++) bv_e[k] = P.transpose() * bv[k] * P;

	project_macro_strain(fr, v, pmf, pmvlist);


	Eigen::Matrix<double, 9, 9> K; K.setZero();
	Eigen::Matrix<double, 9, 6> f; f.setZero();

	Eigen::Matrix3d D;
	D <<
		lam0 + 2 * mu, lam0, 0,
		lam0, lam0 + 2 * mu, 0,
		0, 0, mu;


	Eigen::Matrix<double, 3, 9> B0 = tgnt_membrane_strain_displacement_matrix(tri, v, n, fr.leftCols(2));
	// inplane strain part

	static auto qlist = qd::tri::symq::get(3);

	for (int i = 0; i < qlist.size(); i++) {
		Eigen::Vector3d q(qlist[i].c[0], qlist[i].c[1], qlist[i].c[2]);
		double w = qlist[i].w;
		auto B = B0;

		Eigen::Vector<double, 9> qn;
#if 0
		for (int j = 0; j < 3; j++) { qn.block<3, 1>(j * 3, 0) = q[j] * n; }
#else
		for (int j = 0; j < 3; j++) { qn.block<3, 1>(j * 3, 0) = q[j] * v[j].nv; }
#endif

#ifdef _USE_SECOND_FUND_FORM
		B -= voigt((bv_e[0] * q[0] + bv_e[1] * q[1] + bv_e[2] * q[2]).eval()) * qn.transpose();
#endif

		K += B.transpose() * (D * w) * B;

#ifndef _USE_VERTEX_MACRO_STRAIN_PROJECTION
		f += -B.transpose() * (D * w) * pmf;
#else
		f += -B.transpose() * (D * w) * (pmvlist[0] * q[0] + pmvlist[1] * q[1] + pmvlist[2] * q[2]);
#endif
	}
	K *= A; f *= A;
	return std::make_tuple(K, f);

}

std::tuple<Eigen::Matrix<double, 6, 6>, double> membrane_strain_energy(double lam0, double mu, const Eigen::Matrix3d& tri, const Eigen::Matrix3d& fr, const Compile1ring v[3])
{
	Eigen::Vector3d Av = (tri.col(1) - tri.col(0)).cross(tri.col(2) - tri.col(0))/2;
	double A = Av.norm();
	Eigen::Vector3d n = Av / A;

	Eigen::Matrix3d P = Eigen::Matrix3d::Identity() - n * n.transpose();
	
	Eigen::Matrix<double, 6, 6> EM; EM.setZero();

#ifndef _USE_VERTEX_MACRO_STRAIN_PROJECTION
	for (int i = 0; i < 6; i++) {
		Eigen::Matrix3d Ei = P * fromvoigt(Eigen::Vector<double, 6>::Unit(i).eval()) * P;
		for (int j = i; j < 6; j++) {
			Eigen::Matrix3d Ej = P * fromvoigt(Eigen::Vector<double, 6>::Unit(j).eval()) * P;
			double tr_i = Ei.trace();
			double tr_j = Ej.trace();
			double f2 = Ei.reshaped().dot(Ej.reshaped());
			double E = lam0 * tr_i * tr_j + 2 * mu * f2;
			EM(i, j) = EM(j, i) = E;
		}
	}
#else
	Eigen::Matrix2d Ev[3][6];
	Eigen::Matrix<double, 3, 2> Pv[3];
	for (int k = 0; k < 3; k++) {
		Pv[k] = vertex_projection(fr.leftCols(2), fr.col(2).normalized(), v[k]);
		for (int i = 0; i < 6; i++) {
			Ev[k][i] = Pv[k].transpose() * fromvoigt(Eigen::Vector<double, 6>::Unit(i).eval()) * Pv[k];
		}
	}
	static auto qlist = qd::tri::symq::get(2);
	for (auto qp : qlist) {
		Eigen::Vector3d q(qp.c[0], qp.c[1], qp.c[2]);
		double w = qp.w;
		for (int i = 0; i < 6; i++) {
			Eigen::Matrix2d Ei = Ev[0][i] * q[0] + Ev[1][i] * q[1] + Ev[2][i] * q[2];
			double tr_i = Ei.trace();
			for (int j = i; j < 6; j++) {
				Eigen::Matrix2d Ej = Ev[0][j] * q[0] + Ev[1][j] * q[1] + Ev[2][j] * q[2];
				double tr_j = Ej.trace();
				double f2 = dot2(Ei, Ej);
				double E = lam0 * tr_i * tr_j + 2 * mu * f2;
				EM(i, j) += E * w;
				EM(j, i) = EM(i, j);
			}
		}
	}
#endif

	//{
	//	Eigen::Vector3d e1 = n.cross(tri.col(1) - tri.col(0)).normalized();
	//	Eigen::Vector3d e2 = n.cross(e1).normalized();
	//	std::cout << "e1.e2 = " << e1.dot(e2) << "; e1.n = " << e1.dot(n) << "; e2.n = " << e2.dot(n) << std::endl;
	//	Eigen::Matrix3d D;
	//	D << lam0 + 2 * mu, lam0, 0,
	//		lam0, lam0 + 2 * mu, 0,
	//		0, 0, mu;
	//	Eigen::Matrix<double, 3, 2> P1;
	//	P1 << e1, e2;
	//	Eigen::Matrix2d p_eps = P1.transpose() * epsm * P1;
	//	Eigen::Vector3d p(p_eps(0, 0), p_eps(1, 1), 2 * p_eps(0, 1));
	//	double E1 = p.dot(D * p);
	//	std::cout << "E = " << E << ";  E1 = " << E1 << std::endl;
	//}

	return { EM * A , A };
}

double membrane_strain_energy(double lam0, double mu, const Eigen::Matrix3d& tri, const Eigen::Matrix3d& fr, const Compile1ring v[3], const Eigen::Vector<double, 9>& um, const Eigen::Matrix3d& eps_M)
{
	Eigen::Vector3d Av = fr.col(2);
	double A = Av.norm();
	Eigen::Vector3d n = Av / A;

	//Eigen::Matrix3d P = Eigen::Matrix3d::Identity() - n * n.transpose();
	Eigen::Matrix2d eps_e = fr.leftCols(2).transpose() * eps_M * fr.leftCols(2);
	
	double EM = 0;

	auto B0 = tgnt_membrane_strain_displacement_matrix(tri, v, n, fr.leftCols(2));
	auto be = second_fundamental_form(tri, v[0], v[1], v[2]);

	auto Bn = strain_matrix_edge_stretch(tri, fr.col(0), fr.col(1));

	static auto qlist = qd::tri::symq::get(2);

	for (auto qp : qlist) {
		Eigen::Vector3d q(qp.c[0], qp.c[1], qp.c[2]);
		double w = qp.w;
		auto B = B0;
		Eigen::Vector<double, 9> qn;
		for (int j = 0; j < 3; j++) { qn.block<3, 1>(j * 3, 0) = q[j] * v[j].nv; }
		Eigen::Matrix<double, 3, 9> Bb = Bn * be * qn.transpose();
		B -= Bb;
		Eigen::Matrix2d gam = fromvoigt((B * um).eval());
		gam += eps_e;
		EM += (lam0 * gam.trace() * gam.trace() + 2 * mu * dot2(gam, gam)) * w;
	}
	EM *= A;
	return EM;
}

Eigen::Vector<double, 6> voigt(const Eigen::Matrix3d& eps)
{
	Eigen::Vector<double, 6> veps;
	for (int i = 0; i < 3; i++) {
		veps[i] = eps(i, i);
		veps[i + 3] = 2 * eps((i + 1) % 3, (i + 2) % 3);
	}
	return veps;
}

Eigen::Vector<double, 3> voigt(const Eigen::Matrix2d& eps)
{
	Eigen::Vector<double, 3> veps;
	for (int i = 0; i < 2; i++) {
		veps[i] = eps(i, i);
	}
	veps[2] = eps(0, 1) * 2;
	return veps;
}

Eigen::Vector<double, 6> voigt_stress(const Eigen::Matrix3d& eps)
{
	Eigen::Vector<double, 6> veps;
	for (int i = 0; i < 3; i++) {
		veps[i] = eps(i, i);
		veps[i + 3] = eps((i + 1) % 3, (i + 2) % 3);
	}
	return veps;
}

Eigen::Matrix3d fromvoigt(const Eigen::Vector<double, 6>& eps)
{
	Eigen::Matrix3d E;
	for (int i = 0; i < 3; i++) {
		E(i, i) = eps[i];
		int j = (i + 1) % 3, k = (i + 2) % 3;
		E(j, k) = E(k, j) = eps[i + 3] / 2;
	}
	return E;
}


Eigen::Matrix2d fromvoigt(const Eigen::Vector<double, 3>& eps)
{
	Eigen::Matrix2d E;
	E(0, 0) = eps[0];
	E(1, 1) = eps[1];
	E(0, 1) = eps[2] / 2;
	E(1, 0) = E(0, 1);
	return E;
}


Eigen::Vector3d area_shape_derivative(
	const Eigen::Matrix3d& tri, const Eigen::Matrix3d& fr,
	const Compile1ring v[3]
) {
	auto Bn = strain_matrix_edge_stretch(tri, fr.col(0), fr.col(1));
	auto be = second_fundamental_form(tri, v[0], v[1], v[2]);
	Eigen::Vector3d bform = (Bn * be).eval();
	double H2 = bform[0] + bform[1];
	Eigen::Vector3d vn(H2, H2, H2);

	double A = fr.col(2).norm();

	return -vn * A / 3;
}

Eigen::Vector3d circumcenter(const Eigen::Matrix3d& tri) {
	Eigen::Vector3d a = tri.col(1) - tri.col(0);
	Eigen::Vector3d b = tri.col(2) - tri.col(0);
	Eigen::Vector3d An = a.cross(b);
	Eigen::Vector3d p = (a.squaredNorm() * b - b.squaredNorm() * a).cross(An) / (2 * An.squaredNorm()) + tri.col(0);
	return p;
}

// sdiv = circulation on {e12, e23, e31}
Eigen::Matrix<double, 3, 21> asym_stif_element_shape_derivative(
	const Eigen::Matrix3d& tri, const Eigen::Matrix3d& fr,
	const Compile1ring v[3],
	double lam0, double mu,
	const Eigen::Matrix<double, 9, 6>& um,
	Eigen::VectorXd* usr_ptr
) {
	auto& usr = *usr_ptr;
	Eigen::Vector3d pc = circumcenter(tri);

	Eigen::Matrix<double, 3, 21> Vn;
	Vn.setZero();
	static auto qlist = qd::tri::symq::get(4);

	Eigen::Vector3d n = fr.col(2);
	double A = n.norm(); n /= A;

	Eigen::Vector3d e1 = fr.col(0);
	Eigen::Vector3d e2 = fr.col(1);

	Eigen::Matrix<double, 3, 2> face_fr; face_fr << e1, e2;

	Eigen::Matrix<double, 2, 3> tri2 = face_fr.transpose() * (tri.colwise() - pc);

	Eigen::Matrix<double, 3, 9> B0 = tgnt_membrane_strain_displacement_matrix(tri, v, n, face_fr);
	auto Bn = strain_matrix_edge_stretch(tri, e1, e2);
	auto be = second_fundamental_form(tri, v[0], v[1], v[2]);

	Eigen::Matrix2d Em[6] = {};
	Eigen::Vector2d PEmn[6] = {};
	Eigen::Vector2d dn[3] = {};
	Eigen::Vector2d dV[3] = {};
	Eigen::Matrix2d dEm[3][6];
	Eigen::Vector<double, 6> utm[6];
	Eigen::Vector3d unm[6];
	Eigen::Vector2d dun[6];
	Eigen::Matrix2d dut[6];

	for (int i = 0; i < 6; i++) {
		Eigen::Matrix3d Eij = fromvoigt(Eigen::Vector<double, 6>::Unit(i).eval());
		Em[i] = face_fr.transpose() * Eij * face_fr;
		PEmn[i] = face_fr.transpose() * Eij * n;
		utm[i] = (um.col(i).reshaped(3, 3).transpose() * face_fr).transpose().reshaped();
		for (int j = 0; j < 3; j++) {
			unm[i][j] = um.block<3, 1>(j * 3, i).dot(v[j].nv);
		}
	}

	for (int i = 0; i < 3; i++) {
		Eigen::Matrix<double, 3, 2> dVp;
		dVp << tri.col(i) - tri.col((i + 1) % 3),
			tri.col(i) - tri.col((i + 2) % 3);
		dV[i] = (dVp.transpose() * face_fr).lu().solve(Eigen::Vector2d::Ones());
		dn[i] = dV[i];
		for (int j = 0; j < 6; j++) {
			dEm[i][j] = dn[i] * PEmn[j].transpose();
			dEm[i][j] = (dEm[i][j] + dEm[i][j].transpose()).eval();
		}
	}
	for (int i = 0; i < 6; i++) {
		dun[i] = dV[0] * unm[i][0] + dV[1] * unm[i][1] + dV[2] * unm[i][2];
		dut[i] = (tgnt_displacement_grad_matrix(tri, v, n, face_fr) * um.col(i)).reshaped(2, 2);
	}

#ifdef _USE_VERTEX_MACRO_STRAIN_PROJECTION
	Eigen::Vector2d PEmnvlist[3][6] = {};
	for (int k = 0; k < 3; k++) {
		auto Pv = vertex_projection(face_fr, n, v[k]);
		for (int i = 0; i < 6; i++) {
			PEmnvlist[k][i] = Pv.transpose() * fromvoigt(Eigen::Vector<double, 6>::Unit(i).eval()) * v[k].nv;
		}
	}
#endif

	//if (dun[0].cwiseAbs().maxCoeff() > 1e30) {
	//	std::cout << "dun[0] = " << dun[0].transpose() << std::endl;
	//}

	//bool debug = (tri.rowwise().mean() - Eigen::Vector3d(0.300736, -0.498996, 0.41096)).norm() < 1e-2;

	// derivative of metric
	for (auto qp : qlist) {
		Eigen::Vector3d q(qp.c[0], qp.c[1], qp.c[2]);
		double w = qp.w;

		auto B = B0;

		Eigen::Vector<double, 9> qn;

		for (int j = 0; j < 3; j++) { qn.block<3, 1>(j * 3, 0) = q[j] * v[j].nv; }

		Eigen::Matrix<double, 3, 9> Bb = Bn * be * qn.transpose();
		// negative?
		Eigen::Matrix2d bform = fromvoigt((Bn * be).eval());
		Eigen::Matrix2d cform = bform * bform;


		B -= Bb;

		//Eigen::Matrix2d budV[3][6];
		//Eigen::Matrix2d dVub[3][6];
		//Eigen::Matrix2d dbdt[3];
		//for (int k = 0; k < 3; k++) {
		//	dbdt[k] = -cform * q[k];
		//	for (int j = 0; j < 6; j++) {
		//		budV[k][j] = bform * utm[j].reshaped(2, 3) * q * dV[k].transpose();
		//		budV[k][j] = (budV[k][j] + budV[k][j].transpose()).eval() / 2;
		//		dVub[k][j] = dV[k].dot(utm[j].block<2, 1>(k * 2, 0)) * bform;
		//	}
		//}

#ifdef _USE_VERTEX_MACRO_STRAIN_PROJECTION
		Eigen::Matrix2d dEm[3][6]; // overwrite outer scope
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 6; j++) {
				Eigen::Vector2d pEMn = { 0,0 };
				for (int k = 0; k < 3; k++) pEMn += PEmnvlist[k][j] * q[k];
				dEm[i][j] = dn[i] * pEMn.transpose();
				dEm[i][j] = (dEm[i][j] + dEm[i][j].transpose()).eval();
			}
		}
#endif

		Eigen::Vector3d vn = q;

		int counter = 0;
		for (int i = 0; i < 6; i++) {
			double ui_3 = unm[i].dot(q);
			Eigen::Vector2d ui_t = utm[i].reshaped(2, 3) * q;
			for (int j = i; j < 6; j++) {
				double uj_3 = unm[j].dot(q);
				Eigen::Vector2d uj_t = utm[j].reshaped(2, 3) * q;
				Eigen::Matrix2d gam_i = fromvoigt((B * um.col(i)).eval()) + Em[i], gam_j = fromvoigt((B * um.col(j)).eval()) + Em[j];
				double Bw = (lam0 * (gam_i.trace() * dot2(gam_j, bform) + gam_j.trace() * dot2(gam_i, bform))
					+ 4 * mu * dot2((gam_i * gam_j).eval(), bform));
				double Aw = bform.trace() / 2 * (lam0 * gam_i.trace() * gam_j.trace() + 2 * mu * dot2(gam_i, gam_j));
				{
					if (i == 0 && j == 0 && usr_ptr) {
						usr[0] += dot2(gam_i, gam_j) * A * w;
						usr[1] += dot2(bform, bform) * A * w;
						usr[2] += dot2(ui_3 * bform, uj_3 * bform) * A * w;
						usr[3] += dot2(fromvoigt((B0 * um.col(i)).eval()), fromvoigt((B0 * um.col(j)).eval())) * A * w;
						usr[4] += dot2(Em[i], Em[j]) * A * w;
						usr[5] += dot2(fromvoigt((B * um.col(i)).eval()), fromvoigt((B * um.col(j)).eval())) * A * w;
					}
					//if (i == 0 && j == 0 && debug) {
//#pragma omp critical
						//{
						//	auto pp = (tri * q).eval();
						//	std::cout << "=============================================================================================\n";
						//	std::cout << "p = " << (tri * q).transpose() << std::endl;
						//	std::cout << "u3 = " << ui_3 << ", ut = " << (fr.leftCols(2) * ui_t).transpose() << std::endl;
						//	std::cout << "unorm = " << std::sqrt(ui_t.squaredNorm() + ui_3 * ui_3) << std::endl;
						//	std::cout << "Bw = " << Bw << ", Aw = " << Aw << std::endl;
						//	std::cout << "um = \n" << um.col(0).reshaped(3, 3) << std::endl;
						//	std::cout << "v = " << std::atan2(pp[2], pp[1]) << std::endl;
						//	std::cout << "gam_i = \n" << fr.leftCols(2) * fromvoigt((B * um.col(i)).eval()) * fr.leftCols(2).transpose() << std::endl;
						//	std::cout << "gam_i norm = " << fromvoigt((B * um.col(i)).eval()).squaredNorm() << std::endl;
						//	std::cout << "bform = \n" << fr.leftCols(2) * bform * fr.leftCols(2).transpose() << std::endl;
						//	std::cout << "Em[0] = \n" << fr.leftCols(2) * Em[0] * fr.leftCols(2).transpose() << std::endl;
						//	std::cout << "=============================================================================================\n";
						//}
					//}
				}
				//Eigen::Matrix3d deps; deps.setZero();
				for (int k = 0; k < 3; k++) {
					Eigen::Matrix2d dgam_i, dgam_j;
					dgam_i = bform * ui_t * dV[k].transpose() - vn[k] * bform * dut[i] + dun[i] * dV[k].transpose() + ui_3 * vn[k] * cform;
					dgam_j = bform * uj_t * dV[k].transpose() - vn[k] * bform * dut[j] + dun[j] * dV[k].transpose() + uj_3 * vn[k] * cform;
					dgam_i = (dgam_i + dgam_i.transpose()).eval() / 2; dgam_j = (dgam_j + dgam_j.transpose()).eval() / 2;
					double C =
						lam0 * gam_i.trace() * dgam_j.trace() + 2 * mu * dot2(gam_i, dgam_j) +
						lam0 * gam_j.trace() * dgam_i.trace() + 2 * mu * dot2(gam_j, dgam_i);
					double D =
						lam0 * gam_i.trace() * dEm[k][j].trace() + 2 * mu * dot2(gam_i, dEm[k][j]) +
						lam0 * gam_j.trace() * dEm[k][i].trace() + 2 * mu * dot2(gam_j, dEm[k][i]);
					Vn(k, counter) += 2 * vn[k] * Bw * w;
					Vn(k, counter) += C * w;
					Vn(k, counter) += D * w;
					Vn(k, counter) -= 2 * vn[k] * Aw * w;
					//Vn(k, counter) += vn[k] * dot2(gam_i, gam_j) * w;
					//if (i == 0 && j == 0) { usr[6] += dot2(dEm[k][0], dEm[k][0]) * w; }
				}
				//if (debug && i == 0 && j == 0) {
				//	Eigen::Vector3d qpos = tri * q;
				//	Eigen::Vector3d vn_vertex;
				//	Eigen::Matrix2d dEmq; dEmq.setZero();
				//	Eigen::Matrix2d dgam_0; dgam_0.setZero();
				//	for (int vi = 0; vi < 3; vi++) {
				//		double f = tri.col(vi).bottomRows(2).norm();
				//		double theta = std::atan2(tri(2, vi), tri(1, vi));
				//		double x = tri(0, vi);
				//		Eigen::Vector3d n(-M_PI / 4 * std::sin(M_PI * x), -tri(1, vi), -tri(2, vi)); n.bottomRows(2).normalize();
				//		n.normalize();
				//		Eigen::Vector3d ui(std::cos(M_PI * x), std::cos(M_PI * f) * std::cos(theta), std::cos(M_PI * f) * std::sin(theta));
				//		double vn_i = 0.03 * (std::sin(M_PI * x) + 1);
				//		vn_vertex[vi] = vn_i;
				//		dEmq += dEm[vi][0] * vn_i;
				//		dgam_0 += vn_i * (bform * ui_t * dV[vi].transpose() - q[vi] * bform * dut[i] + dun[i] * dV[vi].transpose() + ui_3 * q[vi] * cform);
				//	}
				//	double theta = std::atan2(qpos[2], qpos[1]);
				//	std::cout << "=============================================================================================\n";
				//	std::cout << "qpos = " << qpos.transpose() << std::endl;
				//	std::cout << "vn_vertex = " << vn_vertex.transpose() << std::endl;
				//	std::cout << "theta = " << theta << std::endl;
				//	std::cout << "dEmq = \n" << face_fr * dEmq * face_fr.transpose() << std::endl;
				//	std::cout << "dgam_0 = \n" << face_fr * dgam_0 * face_fr.transpose() << std::endl;
				//	std::cout << "budvn = \n" << face_fr * (bform * ui_t * (dV[0] * vn_vertex[0] + dV[1] * vn_vertex[1] + dV[2] * vn_vertex[2]).transpose()) * face_fr.transpose() << std::endl;
				//	std::cout << "vnbDut = \n" << face_fr * (bform * dut[0] * (q[0] * vn_vertex[0] + q[1] * vn_vertex[1] + q[2] * vn_vertex[2])) * face_fr.transpose() << std::endl;
				//	std::cout << "dundvn = \n" << face_fr * (dun[i] * (dV[0] * vn_vertex[0] + dV[1] * vn_vertex[1] + dV[2] * vn_vertex[2]).transpose()) * face_fr.transpose() << std::endl;
				//	std::cout << "u3vnc = \n" << face_fr * ui_3 * cform * q.dot(vn_vertex) * face_fr.transpose() << std::endl;
				//	std::cout << "=============================================================================================\n";
				//}
				counter++;
			}
		}
	}

	Vn *=  A;

	{
		//if (std::abs(Vn.cwiseAbs().maxCoeff()) > 1e30) {
		//	std::cout << "* * * * * * * * * * * * * * * * ERROR * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n";
		//	std::cout << "err tri = \n" << tri << std::endl;
		//	std::cout << "fr = \n" << fr << std::endl;
		//	std::cout << "um = \n" << um.col(0).reshaped(3, 3) << std::endl;
		//	std::cout << "dun = " << dun[0].transpose() << std::endl;
		//	std::cout << " compute dun = " << (dV[0] * unm[0][0] + dV[1] * unm[0][1] + dV[2] * unm[0][2]).eval().transpose() << std::endl;
		//	std::cout << "unm = " << unm[0].transpose() << std::endl;
		//	std::cout << "dV = " << dV[0].transpose() << std::endl;
		//	std::cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n";
		//}
	}

	{
		//Eigen::Vector3d tgt = { 0.145963, 0.227324, 0.420735 };
		//if ((pc - tgt).norm() < 5e-3 && usr_ptr) {
		//	auto eps1 = fromvoigt((B0 * um.col(0)).eval());
		//	Eigen::Matrix3d dem; dem.setZero();
		//	for (int k = 0; k < 3; k++) {
		//		dem += face_fr * (dEm[k][0] * usr[7 + k]) * face_fr.transpose();
		//	}
		//	std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
		//	std::cout << "tgt/A = " << usr[3] / A << ",\n tri = \n" << tri << std::endl;
		//	std::cout << "u1 =\n" << um.col(0).reshaped(3, 3) << std::endl;
		//	std::cout << "eps1 =\n" << face_fr * eps1 * face_fr.transpose() << std::endl;
		//	//std::cout << "dem =\n" << dem << std::endl;
		//	std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
		//}
	}
	// derivative of membrane strain
	

	return Vn;
}

using namespace msf;

void solve_global_smooth_bform(PeriodSurfaceMesh& m, const std::vector<Eigen::Matrix3d>& frlist, const std::vector<Compile1ring>& vrings) {
	Eigen::VectorXd blist(m.n_faces() * 3);
	std::vector<std::array<Eigen::Matrix2d, 3>> fts(m.n_faces());
	std::vector<bool> edge_passed(m.n_edges(), false);
	for (auto fh : m.faces()) {
		for (auto he : m.fh_ccw_range(fh)) {
			
		}
	}
}
