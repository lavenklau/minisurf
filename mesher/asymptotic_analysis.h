#pragma once

#include <fstream>
#include "PeriodicMesher.h"
#include <Eigen/Eigen>
#include <unsupported/Eigen/AutoDiff>
#include <tuple>
#include "fundamental_forms.h"
#include <array>

std::tuple<Eigen::Matrix<double, -1, 21>, Eigen::VectorXd> asym_stif_shape_sensitivity(msf::PeriodSurfaceMesh& m, const std::vector<Eigen::Matrix3d>& fr, const std::vector<Compile1ring>& vrings, const Eigen::Matrix<double, -1, 6>& um, const Eigen::Matrix<double, -1, 6>& fm, double lam0, double mu);

std::tuple<
	Eigen::Matrix<double, -1, 6>,
	Eigen::Matrix<double, -1, 6>,
	Eigen::Matrix<double, 6, 6>,
	double
> asym_stif_fem(msf::PeriodSurfaceMesh& m, const std::vector<Eigen::Matrix3d>& frlist, const std::vector<Compile1ring>& vrings, double lam0, double mu);

std::tuple<
	Eigen::SparseMatrix<double>,
	Eigen::Matrix<double, -1, 6>,
	Eigen::Matrix<double, 6, 6>,
	double
> assembly_stif_fem(msf::PeriodSurfaceMesh& m, const std::vector<Eigen::Matrix3d>& frlist, const std::vector<Compile1ring>& vrings, double lam0, double mu);

std::tuple<
	Eigen::Matrix<double, -1, 6>,
	Eigen::Matrix<double, -1, 6>,
	Eigen::Matrix<double, 6, 6>,
	double 
> asym_stif_fem(msf::PeriodSurfaceMesh& m, const std::vector<Eigen::Matrix3d>& frlist, const std::vector<Compile1ring>& vrings, double lam0, double mu, const std::vector<Eigen::Matrix3d>& bvlist);

Eigen::Matrix<double, 6, 6> eval_asym_elastic_tensor(const Eigen::Matrix<double, -1, 6>& fm, const Eigen::Matrix<double, -1, 6>& um, const Eigen::Matrix<double, 6, 6>& Enmem, double As);

Eigen::Matrix3d eval_asym_cond_matrix(msf::PeriodSurfaceMesh& m, const Eigen::Matrix<double, -1, 3>& blist, const Eigen::Matrix<double, -1, 3>& ulist, double As);

std::vector<Eigen::Matrix3d> eval_partial_asym_cond_matrix(msf::PeriodSurfaceMesh& m, const Eigen::Matrix<double, -1, 3>& blist, const Eigen::Matrix<double, -1, 3>& ulist, const Eigen::VectorXi& Cid);

Eigen::SparseMatrix<double> asym_stif_fem_matrix(msf::PeriodSurfaceMesh& m,const  std::vector<Eigen::Matrix3d>& frlist, const std::vector<Compile1ring>& vrings, double lam0, double mu, const std::vector<Eigen::Matrix3d>& bvlist);

Eigen::SparseMatrix<double> asym_stif_fem_matrix(msf::PeriodSurfaceMesh& m, const  std::vector<Eigen::Matrix3d>& frlist, const std::vector<Compile1ring>& vrings, double lam0, double mu);


std::tuple<Eigen::VectorXd, std::vector<Eigen::Vector3d>> eval_vertex_mass(msf::PeriodSurfaceMesh& m);

template<typename T>
auto retrieve_element_vertex_data(msf::PeriodSurfaceMesh& m, OpenMesh::FaceHandle fh, const std::vector<T>& datalist) {
	std::array<T, 3> e_data;
	auto fvh = m.getFaceVertexHandle(fh);
	for (int k = 0; k < 3; k++) e_data[k] = datalist[fvh[k].idx()];
	return e_data;
}

Eigen::MatrixX3d willmore_H2_gradient(msf::PeriodSurfaceMesh& m, const Eigen::SparseMatrix<double>& L, Eigen::VectorXd Av, std::vector<Eigen::Vector3d>& Nv);

Eigen::MatrixX3d willmore_H2_gradient(msf::PeriodSurfaceMesh& m);

std::tuple<std::vector<double>, std::vector<Eigen::Vector3<double>>> eval_cot_edges(msf::PeriodSurfaceMesh& m);

Eigen::VectorXd mean_curvature_vector(const msf::PeriodSurfaceMesh& m, const std::vector<double>& cote, const std::vector<Eigen::Vector3<double>>& evec);

Eigen::Matrix<double, 2, 3> scalar_gradient_matrix(const Eigen::Matrix3d& tri, const Eigen::Matrix3d& fr);

Eigen::MatrixX3d assemble_asym_cond_vector(msf::PeriodSurfaceMesh& m, const std::vector<double>& cotedge, const std::vector<Eigen::Vector3<double>>& Ev);

//Eigen::Matrix<double, -1, 6> asym_cond_sensitivity(msf::PeriodSurfaceMesh& m, const Eigen::MatrixX3d& ulist);

std::tuple<Eigen::Matrix<double, -1, 6>, Eigen::VectorXd> asym_cond_sensitivity(msf::PeriodSurfaceMesh& m, const std::vector<Compile1ring>& vrings, const Eigen::MatrixX3d& ulist);

Eigen::SparseMatrix<double> assemble_asym_cond_matrix(msf::PeriodSurfaceMesh& m, const std::vector<double>& cotedge);

Eigen::AutoDiffScalar<Eigen::Vector<double, 6>> adc_objective(std::string type, const Eigen::Matrix3d& kA);

Eigen::AutoDiffScalar<Eigen::Vector<double, 21>> ads_objective(std::string type, const Eigen::Matrix<double, 6, 6>& CA);

double search_ads_step(bool enable, msf::PeriodSurfaceMesh& m, const Eigen::VectorXd& step_vector, double obj_last, double gTp, double max_step,
	std::function<Eigen::AutoDiffScalar<Eigen::Vector<double, 21>>(const Eigen::Matrix<double, 6, 6>&)> objfunc);

double search_adc_step(bool enable, msf::PeriodSurfaceMesh& m, const Eigen::VectorXd& step_vector, double obj_last, double gTp, double max_step,
	std::function<Eigen::AutoDiffScalar<Eigen::Vector<double, 6>>(const Eigen::Matrix<double, 3, 3>&)> objfunc);

double expect_descent(const std::vector<Eigen::Vector3d>& nlist, const Eigen::VectorXd& step_vector, const Eigen::VectorXd& Gn, const Eigen::VectorXd& Av);

Eigen::Vector<double, 6> extract_cond_matrix_coeff(const Eigen::Matrix3d& kA);

Eigen::Vector<double, 21> extract_elastic_matrix_coeff(const Eigen::Matrix<double, 6, 6>& CA);

std::tuple<std::array<std::pair<int, int>, 21>, Eigen::Matrix<int, 6, 6>> elastic_matrix_flatthen_map(void);

std::string getPath(std::string s);

std::tuple<std::vector<Eigen::Matrix3d>, std::vector<Compile1ring>> generate_face_frame_and_vring(const msf::PeriodSurfaceMesh& m);
std::vector<Compile1ring> generate_vring(const msf::PeriodSurfaceMesh& m);
std::vector<Eigen::Matrix3d> generate_face_frame(const msf::PeriodSurfaceMesh& m);

double maximal_unflip_step(msf::PeriodSurfaceMesh& m, const Eigen::VectorXd& step_vector, double max_step, double angle_thres = M_PI_2);

void log_sensitivity(std::string filename, msf::PeriodSurfaceMesh& m, const Eigen::VectorXd& dfdvn, const Eigen::MatrixX3d& g_willmore);

Eigen::VectorXd assemble_mean_curvature_flow(const std::vector<Compile1ring>& vrings);

std::ofstream logger(void);

struct CheckConverge {
	constexpr static int his_size = 50;
private:
	int iter = 0;
	Eigen::Vector<double, his_size> last_objs;
	Eigen::Vector<double, his_size> last_steps;
	double tol = 1e-6;
	double stepTol = 2e-3;
	double c0 = 1;
	//double k_c = 20;
	double s_c = 0.5;
	auto cid(int it) const { return (it + 2 * his_size) % his_size; }
	// [k,b,err]
	std::tuple<double, double, double> lin_regr(void) const;
	bool is_converge() const;
public:
	CheckConverge(double obj_tol_, double step_tol, double c0_, double sc_)
		: tol(obj_tol_), stepTol(step_tol), c0(c0_)/*, k_c(kc_)*/, s_c(sc_), iter(0) {
		for (int i = 0; i < last_objs.size(); i++) { last_objs[i] = 1e30 * i; last_steps[i] = 1e30; }
	}
	double estimate_next_step(double tmax);
	double estimate_precondition(double cmax);
	bool operator()(double obj, double step);
};

extern int log_detail;

extern double max_time_step;

extern  double weight_willmore;

extern bool disable_line_search;
