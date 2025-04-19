#define FMT_HEADER_ONLY
//#include "mesher/MinsurfMesher.h"
#include "mesher/PeriodicMesher.h"
//#include "solver/HomoEnergyOpter.h"
//#include "mesher/PeriodicBezierMesh.h"
#include "CLI/CLI.hpp"
#include "fstream"
//#include "solver/ShellFEM.h"
#include "matlab/matlab_utils.h"
#include "mesher/dir_utils.h"
#include "fmt/core.h"
#include "cgal/homotopy_utils.h"
#include "cgal/mesh_intersection.h"
#include "igl/read_triangle_mesh.h"
#include "igl/writeOBJ.h"
#include "igl/principal_curvature.h"
//#include "reebhantu_utils.h"
#include <random>

using namespace msf;

std::string linefile;
std::string meshfile, meshfile_x;
std::string jobName;
std::string objective;
std::string boundaryFile;
std::string outdir;
double ext_ratio = 0.2;
bool rescale = false;
double converge_tol = 1e-3;
double converge_kc = 0;
double converge_sc = 0;
double step_tol = 1e-3; // step_tol * max_time_step is real step toleration
int max_iter = 1001;
double precondition_strength = 1;
int delaunay_remesh_outer_iter = 1;
int delaunay_remesh_inner_iter = 5;
double island_cull = 0.333;
double delaunay_lengh_upp = 1.5;
double delaunay_lengh_low = 0.5;
double singularity_surgery_tol = 25;
int surgery_type = 2;
bool asym_no_remesh = false;
bool asym_no_surgery = false;
double max_time_step = 0.1;
double weight_willmore = 0;
double weight_area = 0;
bool disable_line_search = false;
double thick = 0.001;
std::vector<double> thick_list;
double fourier_base_decay = 2;
double fourier_regul = 5;
std::string out_mesh_type = "";
int n_sample = 10;
int sample_type = 1;
Real sample_randness = 0.1;
Real sample_reso = 32;
Real remesh_tgtlen = 0.05;
std::vector<double> remesh_tgtlen_list;
Real remesh_noise = 0.05;
int rand_seed = -1;
int remesh_iter = 5;
int remesh_smooth_iter = 5;
int remesh_smooth_type = 0;
int remesh_pertub_iter = 0;
int remesh_smooth_batch_size = 100;
int log_detail = 0;
int mesh_order = 0;
int mesh_quadr_order = 10;
int mesh_area_pre_min_iter = 0;
int mesh_gk = 0;
bool force_manifold = false;
std::vector<Real> aux_number;
int debug_level = -1;

bool is_mesh_file(const std::string& filename) {
	auto ext = dir_utils::path2extension(filename);
	return ext == ".obj" || ext == ".stl" || ext == ".off";
}

std::vector<std::vector<double>> readPolyLine(std::string filename) {
	std::ifstream ifs(filename);
	if (!ifs) {
		printf("\033[31mCannot open file %s\033[0m\n", filename.c_str());
		throw std::runtime_error("cannot open file");
	}

	char buf[1000];
	std::vector<std::vector<double>> coords;
	while (ifs.getline(buf, 1000)) {
		if (buf[0] == 'L') {
			coords.emplace_back();
		}
		else {
			coords.rbegin()->push_back(std::stod(buf));
		}
	}
	ifs.close();

	// check if polyline is legal
	for (int i = 0; i < coords.size(); i++) {
		if (coords[i].size() % 3 != 0) {
			printf("\033[31m coordinates should be multiple of 3\033[0m\n");
			throw std::runtime_error("illegal format");
		}
		if (coords[i].size() <= 3) {
			printf("\033[31m polyline must contains at least three points\033[0m\n");
			throw std::runtime_error("illegal format");
		}
		else {
			int sbeg = 0, ebeg = coords[i].size() - 3;
			for (int j = 0; j < 3; j++) {
				if (coords[i][sbeg + j] != coords[i][ebeg + j]) {
					printf("\033[31mlast coords should be same as start\033[0m\n");
					throw std::runtime_error("last coords should be same as start");
				}
			}
		}
	}

	return coords;
}

std::string getPath(std::string s) {
	std::filesystem::path pat(outdir);
	if (!std::filesystem::exists(pat)) {
		std::cout << "creating directory " << pat << std::endl;
		bool suc = std::filesystem::create_directories(pat);
		if (!suc) {
			std::cout << "Failed" << std::endl;
			throw std::runtime_error("failed to create directories");
		}
	}
	return (pat / s).string();
}

static std::string getTimeString(void) {
	// Get current time as time_point
	auto now = std::chrono::system_clock::now();
	// Convert time_point to time_t for converting to tm (broken-down time)
	std::time_t now_c = std::chrono::system_clock::to_time_t(now);
	// Convert to broken-down time
	std::tm now_tm = *std::localtime(&now_c);
	// Use stringstream to format the date and time
	std::stringstream ss;
	ss << std::put_time(&now_tm, "[%Y-%m-%d %H:%M:%S]");
	return ss.str();
}
std::ofstream getLog() {
	std::ofstream ofs(getPath("minsurf.log"), std::ios::app);
	ofs << getTimeString() << ">>";
	return ofs;
}
std::ofstream getLog(const std::string& logfile) {
	std::ofstream ofs(getPath(logfile), std::ios::app);
	ofs << getTimeString() << ">>";
	return ofs;
}

void parse(int argc, char** argv) {
	CLI::App app;
	app.add_option("--line", linefile, "file that contains polyline coordinates information");
	app.add_option("--max-iter", max_iter, "allowed max iteration");
	app.add_option("--conv-tol", converge_tol, "allowed max iteration");
	app.add_option("--conv-kc", converge_kc, "precondition estimate parameter kc");
	app.add_option("--conv-sc", converge_sc, "precondition estimate parameter sc");
	app.add_option("--step-tol", step_tol, "minimal step tolerence for converge check");
	app.add_option("--surgery-tol", singularity_surgery_tol, "threshold of singularity measure for surgery");
	app.add_option("--surgery-type", surgery_type, "type of singularity measure");
	app.add_option("--island-cull", island_cull, "threshold for culling island (faces ratio to maximal connected region)");
	app.add_option("--prec", precondition_strength, "allowed max iteration");
	app.add_option("--no-remesh", asym_no_remesh, "disable remesh process");
	app.add_option("--derem-out", delaunay_remesh_outer_iter, "delaunay remesh outer iteration");
	app.add_option("--derem-inn", delaunay_remesh_inner_iter, "delaunay remesh outer iteration");
	app.add_option("--derem-elon", delaunay_lengh_upp, "delaunay remesh edge split threshold");
	app.add_option("--derem-esht", delaunay_lengh_low, "delaunay remesh edge collapse threshold");
	app.add_option("--no-surgery", asym_no_surgery, "allowed max iteration");
	app.add_option("--dt", max_time_step, "maximal time step");
	app.add_option("--ww", weight_willmore, "weight of willmore energy");
	app.add_option("--wa", weight_area, "weight of mean curvature flow");
	app.add_option("--ls", disable_line_search, "maximal time step");
	app.add_option("-m", meshfile, "periodic mesh file to be optimized");
	app.add_option("--mx", meshfile_x, "additional periodic mesh file to be optimized");
	app.add_option("-j", jobName, "job name");
	app.add_option("--obj", objective, "objective for optimization");
	app.add_option("--thick", thick, "shell thickness");
	app.add_option("-s", rescale, "rescale model to unit cube");
	app.add_option("--bd", boundaryFile, "boundary condition file");
	app.add_option("-o", outdir, "output directory");
	app.add_option("-g", debug_level, "debug level for logging");
	app.add_option("--ext-ratio", ext_ratio, "periodic extension ratio");
	app.add_option("--fbd", fourier_base_decay, "fourier base decay");
	app.add_option("--freg", fourier_regul, "function regularity");
	app.add_option("--nsample", n_sample, "number of samples");
	app.add_option("--sample-type", sample_type, "sample type of tpms");
	app.add_option("--sample-rand", sample_randness, "randomness[0-1] of coefficient");
	app.add_option("--sample-reso", sample_reso, "resolution of sampled SDF");
	app.add_option("--remesh-len", remesh_tgtlen, "target length for remeshed edge");
	app.add_option("--remesh-len-list", remesh_tgtlen_list, "list of target length for remeshed edge");
	app.add_option("--remesh-iter", remesh_iter, "remesh iteration");
	app.add_option("--remesh-smthiter", remesh_smooth_iter, "remesh smooth iteration");
	app.add_option("--remesh-smthtype", remesh_smooth_type, "remesh smooth type");
	app.add_option("--remesh-prtiter", remesh_pertub_iter, "remesh perturb iteration");
	app.add_option("--remesh-smth-szbatch", remesh_smooth_batch_size, "smooth round of one iteration");
	app.add_option("--remesh-noise", remesh_noise, "randomly perturbation strength");
	app.add_option("--morder", mesh_order, "order of mesh to represent the geometry");
	app.add_option("--mqorder", mesh_quadr_order, "order of mesh to represent the geometry");
	app.add_option("--mgk", mesh_gk, "order of geometry smoothness");
	app.add_option("--maminiter", mesh_area_pre_min_iter, "pre minimize area iteration");
	app.add_option("--rand-seed", rand_seed, "rand seed");
	app.add_option("--log-detail", log_detail, "log detail level, 0 = no additional log");
	app.add_option("--thick-list", thick_list, "list of thickness for test");
	app.add_option("--fmanifold", force_manifold, "force input transfered to manifold");
	app.add_option("--aux", aux_number, "auxiliary number for debugging");
	app.add_option("--omt", out_mesh_type, "type of output mesh file");

	app.parse(argc, argv);

	if (debug_level >= 0) {
		std::ofstream ofs(getPath("config"));
		ofs << app.config_to_str();
	}
}


void test_period_remesh(void) {
	PeriodSurfaceMesh mesher;
	auto vcut = mesher.readMergePeriodBoundary(meshfile, false);
	auto cutset = mesher.periodic_remesh(remesh_iter, vcut, remesh_tgtlen, remesh_smooth_iter, remesh_pertub_iter);
	auto fn = dir_utils::path2filename(meshfile);
	std::string src_ext = std::string(".");
	if (!out_mesh_type.empty()) src_ext += out_mesh_type;
	else src_ext = dir_utils::path2extension(meshfile);
	mesher.savePeriodicMesh(getPath(fn + ".rmsh" + src_ext), cutset);
}

void test_sample_mesh(void) {
	PeriodSurfaceMesh mesher;
	mesher.sampleMeshes(n_sample, sample_reso, sample_type, sample_randness, remesh_iter, remesh_smooth_iter, remesh_tgtlen);
}

inline Real hs_bound(Real E, Real nu) {
	double lam = E * nu / (1 + nu) / (1 - 2 * nu);
	double mu = E / 2 / (1 + nu);
	double lam0 = 2 * lam * mu / (lam + 2 * mu);
	return 4 * (lam0 + mu) / 9;
}



extern void test_gmsh(void);
extern void test_bezier_tri(void);
extern void test_nurbs(void);
extern void test_abel_jacobi(std::string mfile);
void test_skeleton_mesher(void);
void test_close_mesh(void);
void test_flow(std::string meshfile);
//extern "C" {
	void tricall_main(void); 
//}
void test_mesh_boundary(const std::string& meshfile);
int cgal_remesh_main(void);

void test_membraine_energy(std::string meshfile);

void test_orient_mesh(void);

void test_period_mesh_genus(std::string mfile);

std::vector<Real> offset_volume(std::string mfile, const std::vector<Real>& half_thick);

//void create_template_map(void);

void test_conjugate_surf(std::string meshfile);

void asymp_directional_cond(std::string meshfile);

void show_asym_stif(std::string mfile);

void test_asymptotic_stiffness(std::string meshfile);

void test_asymptotic_conductivity(std::string meshfile);

void test_plane_stress_element(void);

void test_grad_precondition(void);

void test_second_fundamental_form(std::string meshfile);

void test_energy_min(std::string mfile);

void test_energy_accuracy(void);

void test_refinement(void);

void test_bform_perturbation(void);

void test_delaunay_remesh(void);

void test_willmore_energy(void);

void test_inextensible_displacement(std::string meshfile);

void abalation_ads_willmore_remesh(void);

void tailor_ads(std::string meshfile, std::string obj);

void test_surgery(void);

void test_gamma_accuracy(void);

void test_plane_stress_element(std::string mfile);


void test_asym_convergence(std::string obj);

void test_complex_topology(std::string obj_type);


void test_stiffness_matrix_precondition();

void plot_precondition_gradient(std::string mfile, std::string obj_type);

void optimal_cond(std::string meshfile, std::string obj_type);

void test_revolution_surface(std::string obj);

void sensitivity_test(std::string obj);

void test_remove_island(std::string mfile);

void switch_user_routine(std::string name) {
	if (name == "gradprec") {
		//test_plane_stress_element();
		test_grad_precondition();
	}
	else if (name == "bform") {
		test_second_fundamental_form(meshfile);
	}
	else if (name == "emmin") {
		test_energy_min(meshfile);
	}
	else if (name == "gammacc") {
		test_gamma_accuracy();
	}
	else if (name == "eplane") {
		test_plane_stress_element(meshfile);
	}
	else if (name == "enacc") {
		test_energy_accuracy();
	}
	else if (name == "testrefine") {
		test_refinement();
	}
	else if (name == "bpert") {
		test_bform_perturbation();
	}
	else if (name == "testdem") {
		test_delaunay_remesh();
	}
	else if (name == "tailads") {
		tailor_ads(meshfile, objective);
	}
	else if (name == "surj") {
		test_surgery();
	}
	else if (name == "hconv") {
		test_asym_convergence(objective);
	}
	else if (name == "cplxtop") {
		test_complex_topology(objective);
	}
	else if (name == "ksprec") {
		test_stiffness_matrix_precondition();
	}
	else if (name == "figrad") {
		plot_precondition_gradient(meshfile, objective);
	}
	else if (name == "optcond") {
		optimal_cond(meshfile, objective);
	}
	else if (name == "revo") {
		test_revolution_surface(objective);
	}
	else if (name == "testsens") {
		sensitivity_test(objective);
	}
	else if (name == "island") {
		test_remove_island(meshfile);
	}
	else {
		std::cout << "Invalid task name" << std::endl;
		return;
	}
}

void check_periodic_mesh(void) {
	PeriodSurfaceMesh m;
	m.read(meshfile, true, false);
	int n_oldedges = m.n_edges();
	m.mergePeriodEdges();
	int n_newedges = m.n_edges();
	if (n_oldedges != n_newedges) {
		std::cout << "Non periodic mesh!" << std::endl;
		throw std::runtime_error("non-periodic");
	}
}

extern void test_symmetrization();


int main(int argc, char** argv) {
	//create_template_map();
	// parse command line option
	parse(argc, argv);
	//test_symmetrization();

	if (jobName == "helloworld") {
		std::cout << "Hello World!" << std::endl;
	}
	else if (jobName == "periodremesh") {
		test_period_remesh();
	}
	else if (jobName == "checkperiod") {
		check_periodic_mesh();
	}
	else if (jobName == "samplemesh") {
		test_sample_mesh();
	}
	else if (jobName == "thicksk") {
		test_skeleton_mesher();
	}
	else if (jobName == "closemesh") {
		test_close_mesh();
		//cgal_remesh_main();
	}
	else if (jobName == "reorient") {
		test_orient_mesh();
	}
	else if (jobName == "offsetvol") {
		offset_volume(meshfile, thick_list);
	}
	else if (jobName == "genus") {
		test_period_mesh_genus(meshfile);
	}
	else if (jobName == "conj") {
		test_conjugate_surf(meshfile);
	}
	else if (jobName == "ads") {
		show_asym_stif(meshfile);
	}
	else if (jobName == "ass") {
		test_asymptotic_stiffness(meshfile);
	}
	else {
		switch_user_routine(jobName);
	}
	return 0;
}
