#include "asymptotic_analysis.h"
#include "matlab/matlab_utils.h"

using namespace msf;

std::tuple<double, double, double> CheckConverge::lin_regr(void) const
{
	Eigen::Vector<double, his_size - 1> xtick;
	Eigen::Vector<double, his_size - 1> ytick;
	double cur_time = 0;
	for (int i = 0; i < his_size - 1; i++) {
		xtick[his_size - 2 - i] = cur_time;
		ytick[his_size - 2 - i] = last_objs[cid(iter - 1 - i)];
		cur_time -= last_steps[cid(iter - i - 2)];
	}
	//eigen2ConnectedMatlab("xk", xtick);
	//eigen2ConnectedMatlab("yk", ytick);
	int n = xtick.size();
	double k = (n * xtick.dot(ytick) - xtick.sum() * ytick.sum()) / (n * xtick.squaredNorm() - std::pow(xtick.sum(), 2));
	double x_tick_mean = xtick.mean();
	double y_tick_mean = ytick.mean();
	double b = y_tick_mean - k * x_tick_mean;
	double err = (k * xtick + b * Eigen::Vector<double, his_size - 1>::Ones() - ytick).cwiseAbs().maxCoeff();
	return std::make_tuple(k, b, err);
}

bool CheckConverge::is_converge() const
{
	bool small_step = true;
	bool obj_conv = true;
	for (int i = 0; i < 5; i++) {
		small_step = small_step && last_steps[cid(iter - 1 - i)] <= stepTol;
	}
#if 0
	for (int i = 0; i < 5; i++) {
		obj_conv = obj_conv &&
			std::abs((last_objs[cid(iter - i)] - last_objs[cid(iter - i - 1)])
				/ (0.02 + std::abs(last_objs[cid(iter - i)]))) <= tol;
	}
#else
	auto [k, b, err] = lin_regr();
	obj_conv = (iter > his_size) && k < tol && err < tol * 2;
#endif
	return obj_conv || small_step;
}

double CheckConverge::estimate_next_step(double tmax)
{
	double dt = 0;
	for (int k = 0; k < 10; k++) dt += last_steps[cid(iter - 1 - k)];
	return std::min(tmax, 2 * dt / 10);
}

double CheckConverge::estimate_precondition(double cmax)
{
	double c = c0;
#if 0
	if (iter < his_size) {
		return c0;
	}
	else if (iter > his_size) {
		double slope = std::abs(lin_regr());
		c = c0 * std::pow(slope, s_c);
	}
	else {
		//  correction
		double slope = std::abs(lin_regr());
		c = c0;
		c0 = c / std::pow(slope, s_c);
		//eigen2ConnectedMatlab("obhis", last_objs);
		//eigen2ConnectedMatlab("tihis", last_steps);
		std::cout << "Update c0 = " << c0 << ", slope = " << slope << std::endl;
	}
#else
	if (iter < his_size) {
		return c0;
	}
	else {
		auto [k, b, err] = lin_regr();
		double slope = std::abs(k);
		c = c0 * std::pow(slope, s_c);
	}
#endif
	return std::min(c, cmax);
}

bool CheckConverge::operator()(double obj, double step)
{
	last_objs[cid(iter)] = obj;
	last_steps[cid(iter)] = step;
	iter++;
	bool conv = is_converge();
	return conv;
}


std::tuple<Eigen::VectorXd, std::vector<Eigen::Vector3d>> eval_vertex_mass(msf::PeriodSurfaceMesh& m) {
#if 0
	Real As = 0;
	Eigen::VectorXd Av(m.n_vertices());
	std::vector<Eigen::Vector3d> Nv(m.n_vertices());
	for (auto vh : m.vertices()) {
		auto [o, ring] = m.find1ring(vh, 0.5);
		auto [a, nv] = m.area1ring(o, ring, 1);
		Av[vh.idx()] = a;
		As += a;
		Nv[vh.idx()] = nv;
	}
#else
	Eigen::VectorXd  Av(m.n_vertices()); Av.setZero();
	std::vector<Eigen::Vector3d> Nv(m.n_vertices());
	std::vector<Eigen::Vector3d> Anf(m.n_faces());
	for (auto fh : m.faces()) {
		auto tri = m.getFacePeriodVertex(fh, 1);
		Anf[fh.idx()] = (tri.col(1) - tri.col(0)).cross(tri.col(2) - tri.col(0)) / 2;
	}
	std::vector<OM::SmartVertexHandle> vhlist;
	for (auto vh : m.vertices()) { vhlist.emplace_back(vh); }
#pragma omp parallel for
	for (int vid = 0; vid < vhlist.size(); vid++) {
		auto vh = vhlist[vid];
		Eigen::Vector3d n(0, 0, 0);
		double As = 0;
		for (auto ih : vh.incoming_halfedges()) {
			auto nf = Anf[ih.face().idx()];
			double Af = nf.norm();
			As += Af;
			n += m.period_sector_angle(ih) * nf / Af;
		}
		Av[vh.idx()] = As / 3;
		Nv[vh.idx()] = n.normalized();
	}
#endif
	return std::make_tuple(Av, Nv);
}

void log_sensitivity(std::string filename, msf::PeriodSurfaceMesh& m, const Eigen::VectorXd& dfdvn, const Eigen::MatrixX3d& g_willmore)
{
	auto [Av, Nv] = eval_vertex_mass(m);

	std::ofstream ofs(filename, std::ios::binary);
	for (auto vh : m.vertices()) {
		auto p = m.point(vh); ofs.write((const char*)p.data(), sizeof(p));
		Eigen::Vector3d df = dfdvn[vh.idx()] * Nv[vh.idx()];
		ofs.write((const char*)df.data(), sizeof(df));
		Eigen::Vector3d dw = g_willmore.row(vh.idx()).transpose();
		ofs.write((const char*)dw.data(), sizeof(dw));
	}
}

double expect_descent(const std::vector<Eigen::Vector3d>& nlist, const Eigen::VectorXd& step_vector, const Eigen::VectorXd& Gn, const Eigen::VectorXd& Av)
{
	double gTp = 0;
	for (int i = 0; i < nlist.size(); i++) {
		auto nv = nlist[i];
		double vn = nv.dot(step_vector.block<3, 1>(i * 3, 0));
		gTp += Gn[i] * vn * Av[i];
	}
	return gTp;
}

double maximal_unflip_step(msf::PeriodSurfaceMesh& m, const Eigen::VectorXd& step_vector, double max_step, double angle_thres /*= M_PI_2*/) {
	double cos_thres = std::cos(angle_thres);
	for (auto fh : m.faces()) {
		double step = max_step;
		auto V = m.getFacePeriodVertex(fh, 1);
		Eigen::Vector3d n0 = (V.col(1) - V.col(0)).cross(V.col(2) - V.col(0)).normalized();
		auto fvh = m.getFaceVertexHandle(fh);
		Eigen::Matrix3d dV;
		dV << step_vector.block<3, 1>(fvh[0].idx() * 3, 0),
			step_vector.block<3, 1>(fvh[1].idx() * 3, 0),
			step_vector.block<3, 1>(fvh[2].idx() * 3, 0);
		do {
			Eigen::Matrix3d V1 = V + step * dV;
			Eigen::Vector3d n1 = (V1.col(1) - V1.col(0)).cross(V1.col(2) - V1.col(0)).normalized();
			if (n1.dot(n0) > cos_thres) {
				max_step = step; break;
			}
			step *= 0.7;;
		} while (step > 1e-10);
	}
	return max_step;
}


