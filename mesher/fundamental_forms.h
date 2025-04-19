#pragma once
#include <Eigen/Eigen>

struct Compile1ring {
public:
	// local average region
	double As = 0;
	// approximation of vertex normal
	Eigen::Vector3<double> nv = { 0,0,0 };
	Eigen::Vector3<double> Lx = { 0,0,0 };
	// mean & Gauss curvature
	double H = 0, K = 0;
	double mass = 0;

	Compile1ring(void) = default;
private:
	void compile1ring(const Eigen::Vector3d& o, const std::vector<Eigen::Vector3d>& ring);
	//void compileHalfring(const Eigen::Vector3d& o, const std::vector<Eigen::Vector3d>& ring);

public:
	Compile1ring(const Eigen::Vector3d& o, const std::vector<Eigen::Vector3d>& ring) { compile1ring(o, ring); }
};

std::tuple<double, double, double> fundamental_forms(const Eigen::Vector3d& v12, const Compile1ring& v1, const Compile1ring& v2);
Eigen::Matrix<double, 3, 3> strain_matrix_edge_stretch(const Eigen::Matrix3d& tri, const Eigen::Vector3d& e1, const Eigen::Vector3d& e2);
Eigen::Vector3d second_fundamental_form(const Eigen::Matrix3d& tri, const Compile1ring& v1, const Compile1ring& v2, const Compile1ring& v3);
Eigen::Matrix3d second_fundamental_form(const Eigen::Matrix3d& tri, const Eigen::Matrix3d& fr, const Compile1ring v[3]);
Eigen::Matrix3d second_fundamental_form(const Eigen::Matrix3d& tri, const Eigen::Matrix3d& fr, const Eigen::Matrix3d& nv);
Eigen::Vector3d second_fundamental_form(const Eigen::Matrix3d& tri, const Eigen::Matrix3d& nv);
Eigen::Matrix2d second_fundamental_form_2d(const Eigen::Matrix3d& tri, const Eigen::Matrix3d& fr, const Eigen::Matrix3d& nv);
std::tuple< Eigen::Matrix<double, 9, 9>, Eigen::Vector<double, 9>, Eigen::Vector3d, Eigen::Vector3d >
membrane_element_stif_matrix_vector(double lam0, double mu, const Eigen::Matrix3d& tri, const Compile1ring v[3], const Eigen::Matrix3d& epsm);

Eigen::Matrix<double, 3, 9>  tgnt_membrane_strain_displacement_matrix(const Eigen::Matrix3d& tri, const Compile1ring v[3], const Eigen::Vector3d& n, const Eigen::Matrix<double, 3, 2>& P);
Eigen::Matrix<double, 3, 9>  tgnt_membrane_strain_displacement_matrix_novp(const Eigen::Matrix3d& tri, const Compile1ring v[3], const Eigen::Vector3d& n, const Eigen::Matrix<double, 3, 2>& P);
Eigen::Matrix3d face_frame(const Eigen::Matrix3d& tri);
std::tuple< Eigen::Matrix<double, 9, 9>, Eigen::Matrix<double, 9, 6> > membrane_element_vertex_stif_matrix_vector(double lam0, double mu, const Eigen::Matrix3d& tri, const Eigen::Matrix3d& fr, const Compile1ring v[3]);
std::tuple< Eigen::Matrix<double, 9, 9>, Eigen::Matrix<double, 9, 6> > membrane_element_vertex_stif_matrix_vector(double lam0, double mu, const Eigen::Matrix3d& tri, const Eigen::Matrix3d& fr, const Compile1ring v[3], const Eigen::Matrix3d bv[3]);
std::tuple<Eigen::Matrix<double, 6, 6>, double> membrane_strain_energy(double lam0, double mu, const Eigen::Matrix3d& tri, const Eigen::Matrix3d& fr, const Compile1ring v[3]);

double membrane_strain_energy(double lam0, double mu, const Eigen::Matrix3d& tri, const Eigen::Matrix3d& fr, const Compile1ring v[3], const Eigen::Vector<double, 9>& um, const Eigen::Matrix3d& eps_M);
Eigen::Vector<double, 6> voigt(const Eigen::Matrix3d& eps);
Eigen::Vector<double, 6> voigt_stress(const Eigen::Matrix3d& eps);
Eigen::Matrix3d fromvoigt(const Eigen::Vector<double, 6>& eps);
Eigen::Vector3d area_shape_derivative(const Eigen::Matrix3d& tri, const Eigen::Matrix3d& fr, const Compile1ring v[3]);
Eigen::Matrix<double, 3, 21> asym_stif_element_shape_derivative(
	const Eigen::Matrix3d& tri, const Eigen::Matrix3d& fr, const Compile1ring v[3], double lam0, double mu, const Eigen::Matrix<double, 9, 6>& um,
	Eigen::VectorXd* usr
);
Eigen::Vector<double, 3> voigt(const Eigen::Matrix2d& eps);
Eigen::Matrix2d fromvoigt(const Eigen::Vector<double, 3>& eps);

Eigen::Matrix3d rotateAlign(const Eigen::Vector3d& vfrom, const Eigen::Vector3d& vto);
