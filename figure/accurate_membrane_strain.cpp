//#include "solver/ShellFEM.h"
#include "mesher/fundamental_forms.h"
#include "mesher/PeriodicMesher.h"
#include <Eigen/PardisoSupport>
#include <unsupported/Eigen/AutoDiff>
#include "matlab/matlab_utils.h"
#include "material/materail.h"
#include "function_pools.h"
#include "mesher/dir_utils.h"
#include "mesher/asymptotic_analysis.h"
#include "mesher/isosurface_generation.h"
#include "cgal/cgal_utils.h"
#include <fstream>
#include "igl/per_vertex_normals.h"
#include "igl/write_triangle_mesh.h"


using namespace msf;

extern std::vector<Real> aux_number;

std::tuple<Eigen::MatrixX3<msf::Real>, Eigen::MatrixX3i> removeDupVertices(const Eigen::MatrixX3<msf::Real>& v, const Eigen::MatrixX3i& f, msf::Real eps /*= 1e-5*/);

Eigen::Matrix<double, 3, 2> principle_curvature(const Eigen::Matrix3d& fr, const Eigen::Matrix3d& bform) {
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> so(fr.leftCols(2).transpose() * bform * fr.leftCols(2));
	Eigen::Vector2d d = so.eigenvalues().cwiseAbs();
	Eigen::Matrix<double, 3, 2> curvature_vector = fr.leftCols(2) * so.eigenvectors() * d.asDiagonal();
	return curvature_vector;
}

void test_gamma_accuracy(void) {
	double Y = 1, nu = 0.3;
	auto [lam, mu] = mtl::lameCoeff(Y, nu);
	double lam0 = 2 * lam * mu / (lam + 2 * mu);
	std::cout << "lam = " << lam << ", lam0 = " << lam0 << ", mu = " << mu << std::endl;


#if 0
	PeriodSurfaceMesh m;
	std::string mfile = "D:/projects/minisurf/temp/testformat/ball.stl";
	m.read(mfile, false, false, false);
	int ndof = m.n_vertices() * 3;
	std::vector<Compile1ring> vrings;
	for (auto vh : m.vertices()) { auto [o, ring] = m.find1ring(vh, 1); vrings.emplace_back(o, ring); }

	Sphere Ball;
	V1Field V1;

	std::vector<Eigen::Vector3d> ulist(m.n_vertices());
	for (auto vh : m.vertices()) {
		auto p = toEigen(m.point(vh));
		auto v = V1.at(p);
		ulist[vh.idx()] = v;
	}

	std::vector<Eigen::Vector3d> fp(m.n_faces());
	std::vector<Eigen::Vector2d> fuv(m.n_faces());
	std::vector<OM::SmartFaceHandle> fhlist(m.n_faces());
	std::vector<Eigen::Matrix3d> gam_ana(m.n_faces());
	std::vector<Eigen::Matrix3d> gam_raw(m.n_faces());
	std::vector<Eigen::Matrix3d> gam_bf(m.n_faces());
	std::vector<Eigen::Matrix3d> gam_vpbf(m.n_faces());

	for (auto fh : m.faces()) {
		fhlist[fh.idx()] = fh;
		auto tri = m.getFacePeriodVertex(fh, 1);
		fp[fh.idx()] = tri.rowwise().mean();
		auto uv = Ball.project(fp[fh.idx()]);
		fuv[fh.idx()] = uv;
		fp[fh.idx()] = Ball.eval(uv[0], uv[1]);
		auto Itg = Ball.Itgnt(uv[0], uv[1]);
		gam_ana[fh.idx()] = Itg * Ball.gamma(V1, uv[0], uv[1]) * Itg.transpose();
	}

	{
		double u = 1, v = 1;
		auto Itg = Ball.Itgnt(u, v);
		auto gm = Ball.gamma(V1, u, v);
		auto p = Ball.eval(u, v);
		std::cout << "Itg =\n" << Itg << std::endl;
		std::cout << "gamma =\n" << gm << std::endl;
		std::cout << "p =\n" << p.transpose() << std::endl;
		std::cout << "p =\n" << Ball.n(u, v).transpose() << std::endl;
		std::cout << "b =\n" << Ball.b(u, v) << std::endl;
		std::cout << "V3 = \n" << V1.at(p).transpose() << std::endl;
	}

	{
		std::ofstream ofs("fp", std::ios::binary);
		ofs.write((const char*)fp.data(), fp.size() * sizeof(Eigen::Vector3d));
	}


#pragma omp parallel for
	for (int fid = 0; fid < fhlist.size(); fid++) {
		auto fh = fhlist[fid];
		auto fvh = m.getFaceVertexHandle(fh);
		auto tri = m.getFacePeriodVertex(fh,1);
		auto fr = face_frame(tri);
		auto n = fr.col(2).normalized();
		Eigen::Vector<double, 9> uf;
		Compile1ring v[3];
		for (int k = 0; k < 3; k++) uf.block<3, 1>(k * 3, 0) = ulist[fvh[k].idx()];
		for (int k = 0; k < 3; k++) v[k] = vrings[fvh[k].idx()];
		auto B0 = tgnt_membrane_strain_displacement_matrix_novp(tri, v, n, fr.leftCols(2));
		auto B1 = tgnt_membrane_strain_displacement_matrix(tri, v, n, fr.leftCols(2));
		auto bf = second_fundamental_form(tri, v[0], v[1], v[2]);
		auto Bn = strain_matrix_edge_stretch(tri, fr.col(0), fr.col(1));
		double u3 = (uf.middleRows(0, 3).dot(v[0].nv) + uf.middleRows(3, 3).dot(v[1].nv) + uf.middleRows(6, 3).dot(v[2].nv)) / 3;

		auto raw = fromvoigt((B0 * uf).eval());
		gam_raw[fh.idx()] = fr.leftCols(2) * raw * fr.leftCols(2).transpose();

		auto raw_bf = fromvoigt((B0 * uf - u3 * Bn * bf).eval());
		gam_bf[fh.idx()] = fr.leftCols(2) * raw_bf * fr.leftCols(2).transpose();

		auto vpbf = fromvoigt((B1 * uf - u3 * Bn * bf).eval());
		gam_vpbf[fh.idx()] = fr.leftCols(2) * vpbf * fr.leftCols(2).transpose();
	}

	{
		std::cout << "uv = " << fuv[123].transpose() << std::endl;
		std::cout << "gam = \n" << gam_ana[123] << std::endl;
	}

	{
		std::ofstream ofs("err", std::ios::binary);
		for (int fid = 0; fid < fhlist.size(); fid++) {
			double er0 = std::sqrt((gam_raw[fid] - gam_ana[fid]).squaredNorm());
			double er1 = std::sqrt((gam_bf[fid] - gam_ana[fid]).squaredNorm());
			double er2 = std::sqrt((gam_vpbf[fid] - gam_ana[fid]).squaredNorm());
			ofs.write((const char*)&er0, sizeof(double));
			ofs.write((const char*)&er1, sizeof(double));
			ofs.write((const char*)&er2, sizeof(double));
		}
	}
#else

	V1Field V1;

	for (auto sn : { "sphere", "cylinder", "catenoid" }) {
		auto [f, df, d2f] = sfn::get(sn);
		Eigen::Vector3d pl, pr;
		pl.setConstant(-1); pr.setConstant(1);
		Eigen::MatrixX3d V;
		Eigen::MatrixX3i F;
		std::tie(V, F) = isosurf_mesh(128, std::make_pair(pl, pr), [=](double x, double y, double z) {return f({ x,y,z }); }, 0);
		std::tie(V, F) = cgal_remesh(V, F);

		std::tie(V, F) = removeDupVertices(V, F, 1e-5);
		
		// project mesh vertices to sphere
		for (int vid = 0; vid < V.rows(); vid++) {
			Eigen::Vector3d x0 = V.row(vid).transpose().eval();
			Eigen::Vector3d d = df(x0);
			V.row(vid) = search_for_zero(f, x0, d, 2e-2).transpose();
		}

		Eigen::MatrixX3d Nv;
		igl::per_vertex_normals(V, F, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE, Nv);
		for (int i = 0; i < Nv.rows(); i++) { if (Nv.row(i).hasNaN()) { std::cout << "Vertex " << V.row(i) << " no normal" << std::endl; } }

		std::vector<double> err_plane, err_bf, err_vtbf;
		std::vector<Eigen::Matrix<double, 3, 2>> vcur_exact, vcur_approx;

		for (int fid = 0; fid < F.rows(); fid++) {
			Eigen::Vector3i fv(F(fid, 0), F(fid, 1), F(fid, 2));
			
			Eigen::Matrix3d tri; tri << V.row(fv[0]), V.row(fv[1]), V.row(fv[2]);
			tri.transposeInPlace();

			Eigen::Matrix3d U; U << V1.at(tri.col(0)), V1.at(tri.col(1)), V1.at(tri.col(2));

			Eigen::Vector3d n = (tri.col(1) - tri.col(0)).cross(tri.col(2) - tri.col(0)).normalized();

			auto fr = face_frame(tri);

			Eigen::Vector3d c = tri.rowwise().mean();

			Eigen::Vector3d g = df(c);

			c = search_for_zero(f, c, g, 2e-2); g = df(c);

			Eigen::Matrix3d Vgrad = V1.grad(c);

			Eigen::Matrix3d Hf = d2f(c);

			Eigen::Matrix3d bform_exact = sfn::implicit_bform(g, Hf);

			vcur_exact.push_back(principle_curvature(fr, bform_exact));

			// exact membrane strain
			Eigen::Matrix3d P = Eigen::Matrix3d::Identity() - g.normalized() * g.normalized().transpose();
			Eigen::Matrix3d gam_exact = P * Vgrad * P;

			Eigen::Matrix3d Nvf;
			for (int i = 0; i < 3; i++) { Nvf.col(i) = Nv.row(fv[i]).transpose(); }
			auto bform_approx = second_fundamental_form(tri, fr, Nvf);

			vcur_approx.push_back(principle_curvature(fr, bform_approx));

			// plane element
			Eigen::Matrix<double, 2, 3> u2 = fr.leftCols(2).transpose() * U;
			Eigen::Matrix2d du; du << u2.col(1) - u2.col(0), u2.col(2) - u2.col(0);
			Eigen::Matrix2d dv; dv << fr.leftCols(2).transpose() * (tri.col(1) - tri.col(0)),
				fr.leftCols(2).transpose()* (tri.col(2) - tri.col(0));
			Eigen::Matrix2d dudv = du * dv.inverse();
			Eigen::Matrix3d gam_plane = fr.leftCols(2) * (dudv + dudv.transpose()) / 2 * fr.leftCols(2).transpose();
			err_plane.push_back((gam_plane - gam_exact).norm());

			// plane element + approx bform
			double u3 = 0;
			for (int i = 0; i < 3; i++) { u3 += U.col(i).dot(Nv.row(fv[i]).transpose()) / 3; }
			//if (n.dot(g) < 0) u3 *= -1;
			Eigen::Matrix3d gam_bform = gam_plane - u3 * bform_approx;
			err_bf.push_back((gam_bform - gam_exact).norm());

			// vertex projection + approx bform
			for (int i = 0; i < 3; i++) u2.col(i) = fr.leftCols(2).transpose() * rotateAlign(Nv.row(fv[i]).transpose(), n) * U.col(i);
			du << u2.col(1) - u2.col(0), u2.col(2) - u2.col(0);
			dudv = du * dv.inverse();
			Eigen::Matrix3d gam_vtproj_bform = fr.leftCols(2) * (dudv + dudv.transpose()) / 2 * fr.leftCols(2).transpose() - u3 * bform_approx;
			err_vtbf.push_back((gam_vtproj_bform - gam_exact).norm());
		}

		igl::write_triangle_mesh(getPath(std::string(sn) + ".obj"), V, F);
		{ std::ofstream ofs(getPath(std::string(sn) + ".fcur"), std::ios::binary); ofs.write((const char*)vcur_exact.data(), vcur_exact.size() * sizeof(vcur_exact[0])); }
		{ std::ofstream ofs(getPath(std::string(sn) + ".fcurapprox"), std::ios::binary); ofs.write((const char*)vcur_approx.data(), vcur_approx.size() * sizeof(vcur_approx[0])); }
		{ std::ofstream ofs(getPath(std::string(sn) + ".err_pl"), std::ios::binary); ofs.write((const char*)err_plane.data(), err_plane.size() * sizeof(double)); }
		{ std::ofstream ofs(getPath(std::string(sn) + ".err_bf"), std::ios::binary); ofs.write((const char*)err_bf.data(), err_bf.size() * sizeof(double)); }
		{ std::ofstream ofs(getPath(std::string(sn) + ".err_vtbf"), std::ios::binary); ofs.write((const char*)err_vtbf.data(), err_vtbf.size() * sizeof(double)); }
	}
#endif
}

void test_energy_accuracy(void) {
	double Y = 1, nu = 0.3;
	auto [lam, mu] = mtl::lameCoeff(Y, nu);
	double lam0 = 2 * lam * mu / (lam + 2 * mu);
	std::cout << "lam = " << lam << ", lam0 = " << lam0 << ", mu = " << mu << std::endl;

	std::vector<std::string> mlist{
		"D:/projects/minisurf/image/siggraph/refine-mesh/ball0.stl",
		"D:/projects/minisurf/image/siggraph/refine-mesh/ball1.stl"
	};

	Eigen::Matrix3d eps_M; eps_M.setZero(); eps_M(0, 0) = 1;

	for (auto mfile : mlist) {
		PeriodSurfaceMesh m;
		m.read(mfile, false, false, false);
		int ndof = m.n_vertices() * 3;
		std::vector<Compile1ring> vrings;
		std::vector<Eigen::Matrix3d> frlist(m.n_faces());
		std::vector<OM::SmartFaceHandle> fhlist(m.n_faces());
		for (auto vh : m.vertices()) { auto [o, ring] = m.find1ring(vh, 1); vrings.emplace_back(o, ring); }
		for (auto fh : m.faces()) {
			auto tri = m.getFacePeriodVertex(fh, 1);
			frlist[fh.idx()] = face_frame(tri);
			auto fvh = m.getFaceVertexHandle(fh);
			fhlist[fh.idx()] = fh;
		}

		Sphere Ball; V1Field V1;
		std::vector<Eigen::Vector3d> ulist(m.n_vertices());
		for (auto vh : m.vertices()) {
			auto p = toEigen(m.point(vh));
			auto v = V1.at(p);
			ulist[vh.idx()] = v;
		}

		double Es = 0;
		for (int fid = 0; fid < fhlist.size(); fid++) {
			auto fh = fhlist[fid];
			auto fvh = m.getFaceVertexHandle(fh);
			auto tri = m.getFacePeriodVertex(fh, 1);
			auto fr = frlist[fh.idx()];
			auto n = fr.col(2).normalized();
			Eigen::Vector<double, 9> uf;
			Compile1ring v[3];
			for (int k = 0; k < 3; k++) uf.block<3, 1>(k * 3, 0) = ulist[fvh[k].idx()];
			for (int k = 0; k < 3; k++) v[k] = vrings[fvh[k].idx()];

			double ee = membrane_strain_energy(lam0, mu, tri, fr, v, uf, eps_M);
			Es += ee;
		}
		std::cout << dir_utils::path2filename(mfile) << ",  E = " << Es << std::endl;
	}

}


auto plane_element_matrix_vector(const Eigen::Matrix3d& tri, const Eigen::Matrix3d& fr, double lam0, double mu) {
	Eigen::Matrix2d V;
	V << tri.col(1) - tri.col(0), tri.col(2) - tri.col(0);
	Eigen::Matrix<double, 3, 2> S;
	S << -1, -1, 1, 0, 0, 1;
	Eigen::Matrix<double, 3, 2> Rgt = S * V.inverse();
	Eigen::Matrix<double, 2, 3> Lef = fr.leftCols(2).transpose();

	Eigen::Matrix<double, 3, 9> Bt; Bt.setZero();

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			Bt(0, i * 3 + j) += Lef(0, j) * Rgt(i, 0);
			Bt(1, i * 3 + j) += Lef(1, j) * Rgt(i, 1);
			Bt(2, i * 3 + j) += Lef(0, j) * Rgt(i, 1) + Lef(1, j) * Rgt(i, 0);
		}
	}

	Eigen::Matrix<double, 3, 3> D; D.setZero();
	D.block<2, 2>(0, 0).setConstant(lam0); 
	D.block<2, 2>(0, 0).diagonal() += Eigen::Vector2d(2 * mu, 2 * mu);
	D(2, 2) = mu;

	double Ae = fr.col(2).norm();

	Eigen::Matrix<double, 9, 9> Ke = Bt.transpose() * D * Bt * Ae;

	Eigen::Matrix<double, 6, 6> Em; Em.setZero();

	Eigen::Matrix<double, 9, 6> fe;
	for (int i = 0; i < 6; i++) {
		auto eps_i = fromvoigt(Eigen::Vector<double, 6>::Unit(i).eval());
		Eigen::Matrix2d E_i = fr.leftCols(2).transpose() * eps_i * fr.leftCols(2);
		fe.col(i) = Bt.transpose() * D * voigt(E_i) * Ae;
		for (int j = i; j < 6; j++) {
			auto eps_j = fromvoigt(Eigen::Vector<double, 6>::Unit(j).eval());
			Eigen::Matrix2d E_j = fr.leftCols(2).transpose() * eps_j * fr.leftCols(2);
			Em(i, j) = voigt(E_i).dot(D * voigt(E_j)) * Ae;
			Em(j, i) = Em(i, j);
		}
	}

	return std::make_tuple(Ke, fe, Em, Ae);
}



auto plane_element_asym_stif_matrix_vector(PeriodSurfaceMesh& m, double lam0, double mu) {
	Eigen::SparseMatrix<double> K(m.n_vertices() * 3, m.n_vertices() * 3);
	Eigen::Matrix<double, -1, 6> f(m.n_vertices() * 3, 6); f.setZero();

	Eigen::Matrix<double, 6, 6> Em; Em.setZero();

	double As = 0;

	std::vector<Eigen::Triplet<double>> triplist;

	for (auto fh : m.faces()) {
		auto tri = m.getFacePeriodVertex(fh, 1);
		auto fvh = m.getFaceVertexHandle(fh);
		auto fr = face_frame(tri);
		auto [Ke, fe, Ee, Ae] = plane_element_matrix_vector(tri, fr, lam0, mu);
		Em += Ee; As += Ae;
		for (int i = 0; i < 9; i++) {
			int vi = fvh[i / 3].idx();
			for (int j = 0; j < 9; j++) {
				int vj = fvh[j / 3].idx();
				triplist.emplace_back(vi * 3 + i % 3, vj * 3 + j % 3, Ke(i, j));
			}
			f.row(vi * 3 + i % 3) += fe.row(i);
		}
	}
	K.setFromTriplets(triplist.begin(), triplist.end());
	return std::make_tuple(K, f, Em, As);
}

Eigen::VectorXd remove_translation(const Eigen::VectorXd& v, const Eigen::VectorXd& m) {
	Eigen::MatrixX3d V = v.reshaped(3, v.size() / 3).transpose();
	V.rowwise() -= (m.transpose() * V / m.sum()).eval();
	return V.transpose().reshaped();
}

auto plane_element_asym_fem(PeriodSurfaceMesh& m, double lam0, double mu) {
	auto [K, f, Em, As] = plane_element_asym_stif_matrix_vector(m, lam0, mu);

	auto [Av, Nv] = eval_vertex_mass(m);

	for (int i = 0; i < 6; i++) { 
		Eigen::MatrixX3d Fi = f.col(i).reshaped(3, f.rows() / 3).transpose();
		for (int k = 0; k < 3; k++) { Fi.col(k).array() -= Fi.col(k).mean(); }
		f.col(i) = Fi.transpose().reshaped();
	}

	//eigen2ConnectedMatlab("fm", f);

	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> so(K);

	Eigen::Matrix<double, -1, 6> u = so.solve(f);

	for (int i = 0; i < 6; i++) u.col(i) = remove_translation(u.col(i), Av);

	//eigen2ConnectedMatlab("um", u);

	auto CA = eval_asym_elastic_tensor(f, u, Em, As);

	std::cout << "plane Em = \n" << Em << std::endl;

	return std::make_tuple(CA, u);
}

auto my_element_asym_fem(PeriodSurfaceMesh& m, double lam0, double mu) {
	auto [frlist, vrings] = generate_face_frame_and_vring(m);

	auto [fm, um, Enmem, As] = asym_stif_fem(m, frlist, vrings, lam0, mu);

	Eigen::Matrix<double, 6, 6> CA;
	for (int i = 0; i < 6; i++) {
		for (int j = i; j < 6; j++) {
			CA(i, j) = (Enmem(i, j) - um.col(i).dot(fm.col(j))) / As;
			CA(j, i) = CA(i, j);
		}
	}

	std::cout << "plane Em = \n" << Enmem << std::endl;

	return std::make_tuple(CA, um);
}

void save_displacement_mesh(std::string filepth, const PeriodSurfaceMesh& m, const Eigen::VectorXd& v) {
	Eigen::VectorXd nv = v.reshaped(3, v.size() / 3).colwise().norm().reshaped();
	double vmax = nv.maxCoeff();
	auto m_new = m;
	for (auto vh : m_new.vertices()) {
		auto p = m_new.point(vh) + aux_number.at(0) * toOM(v.block<3, 1>(vh.idx() * 3, 0)) / vmax;
		m_new.set_point(vh, p);
	}
	m_new.saveUnitCell(filepth);
}

void test_plane_stress_element(std::string mfile) {
	double Y = 1, nu = 0.3;
	auto [lam, mu] = mtl::lameCoeff(Y, nu);
	double lam0 = 2 * lam * mu / (lam + 2 * mu);
	std::cout << "lam = " << lam << ", lam0 = " << lam0 << ", mu = " << mu << std::endl;
	PeriodSurfaceMesh m;
	//std::string mfile = "D:/projects/minisurf/image/siggraph/accu_gamma/P0.2-2.stl";
	//std::string mfile = "D:/projects/minisurf/data/Pcell.off";
	m.readMergePeriodBoundary(mfile);
	int ndof = m.n_vertices() * 3;

	auto [CA_my, u_my] = my_element_asym_fem(m, lam0, mu);
	auto [CA_pl, u_pl] = plane_element_asym_fem(m, lam0, mu);

	std::cout << "CA_my = \n" << CA_my << std::endl;
	std::cout << "CA_pl = \n" << CA_pl << std::endl;


	// hydrostatic loadding element
	/// save vertex position
	{
		std::ofstream ofs(getPath("vpos"), std::ios::binary);
		for (auto vh : m.vertices()) { auto p = m.point(vh); ofs.write((const char*)p.data(), sizeof(p)); }
	}
	/// Our element
	Eigen::VectorXd u_b(u_my.rows()); u_b.setZero();
	for (int i = 0; i < 3; i++) { u_b += u_my.col(i) / 3; }
	{
		std::ofstream ofs(getPath("ub_ours"), std::ios::binary);
		ofs.write((const char*)u_b.data(), u_b.size() * sizeof(double));
	}
	save_displacement_mesh(getPath("m_ours.obj"), m, u_b);
	/// plane element
	u_b.setZero(); for (int i = 0; i < 3; i++) { u_b += u_pl.col(i) / 3; }
	{
		std::ofstream ofs(getPath("ub_pl"), std::ios::binary);
		ofs.write((const char*)u_b.data(), u_b.size() * sizeof(double));
	}
	save_displacement_mesh(getPath("m_pl.obj"), m, u_b);
}

void test_asym_stif_shape_sensitive(void) {
	double Y = 1, nu = 0.3;

	auto [lam, mu] = mtl::lameCoeff(Y, nu);
	double lam0 = 2 * lam * mu / (lam + 2 * mu);

	std::cout << "lam = " << lam << ", lam0 = " << lam0 << ", mu = " << mu << std::endl;

	PeriodSurfaceMesh m;
	std::string mfile = "D:/projects/minisurf/temp/testformat/ball.stl";
	m.read(mfile, false, false, false);
	int ndof = m.n_vertices() * 3;

	std::vector<Compile1ring> vrings;
	for (auto vh : m.vertices()) { auto [o, ring] = m.find1ring(vh, 1); vrings.emplace_back(o, ring); }

	//auto L =  m.getLaplacian();
	//Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> so(L);

	Eigen::Matrix<double, -1, 6> um(ndof, 6);
	Eigen::VectorXd vn(m.n_vertices());
	um.setRandom(); um /= um.norm();


	// predefined vector fields
	for (auto vh : m.vertices()) {
		auto p = m.point(vh);
		Eigen::Vector3d up(std::cos(p[0]), std::sin(p[0] * p[1]), std::sin(p[0] + p[2]));
		um.block<3, 1>(vh.idx() * 3, 0) = up;
		vn[vh.idx()] = std::sin(p[0] * p[1] + p[2]);
	}

	{
		m.write("ballread.obj");
		std::ofstream ofs("vnball", std::ios::binary);
		ofs.write((const char*)vn.data(), vn.size() * sizeof(double)); ofs.close();
		ofs.open("umball", std::ios::binary);
		ofs.write((const char*)um.col(0).data(), sizeof(double) * um.rows()); ofs.close();
		ofs.open("nball", std::ios::binary);
		for (int k = 0; k < vrings.size(); k++) { ofs.write((const char*)vrings[k].nv.data(), sizeof(vrings[k].nv)); }
	}

	std::vector<Eigen::Matrix3d> frlist(m.n_faces());

	Eigen::Matrix<double, -1, 21> Vn_en(m.n_vertices(), 21); Vn_en.setZero();
	Eigen::VectorXd Vn_A(m.n_vertices()); Vn_A.setZero();

	eigen2ConnectedMatlab("um", um);

	Eigen::Matrix3d Em; Em.setZero(); Em(0, 0) = 1;

	double Int_em = 0;

	std::vector<OM::SmartFaceHandle> fhlist(m.faces_begin(), m.faces_end());

	Eigen::VectorXd usr(10); usr.setZero();

#pragma omp parallel for
	for (int fid = 0; fid < fhlist.size(); fid++) {
		auto fh = fhlist[fid];
		auto tri = m.getFacePeriodVertex(fh, 1);
		auto fvh = m.getFaceVertexHandle(fh);
		Compile1ring v[3];
		for (int i = 0; i < 3; i++) v[i] = vrings[fvh[i].idx()];

		Eigen::Vector3d An = (tri.col(1) - tri.col(0)).cross(tri.col(2) - tri.col(0)) / 2;
		Eigen::Matrix3d fr;
		fr.col(0) = (tri.col(1) - tri.col(0)).normalized();
		fr.col(1) = An.cross(fr.col(0)).normalized();
		fr.col(2) = An;

		Eigen::Matrix<double, 9, 6> ue;
		for (int k = 0; k < 3; k++) {
			ue.middleRows(k * 3, 3) = um.middleRows(fvh[k].idx() * 3, 3);
		}


		Eigen::VectorXd usr_e(10); usr_e.setZero(); 
		usr_e[7] = vn[fvh[0].idx()]; usr_e[8] = vn[fvh[1].idx()]; usr_e[9] = vn[fvh[2].idx()];
		// integral part (nominator)
		auto Vne = asym_stif_element_shape_derivative(tri, fr, v, lam0, mu, ue, &usr_e);

		// area part (denominator)
		//As += fr[fh.idx()].col(2).norm();
		auto dAdv = area_shape_derivative(tri, fr, v);

		// accumulation
#pragma omp critical 
		{
			for (int k = 0; k < 3; k++) { Vn_en.row(fvh[k].idx()) += Vne.row(k); }
			// test part1
			Int_em += (fr.leftCols(2).transpose() * Em * fr.leftCols(2)).squaredNorm() * An.norm();
			for (int i = 0; i < 3; i++) { Vn_A[fvh[i].idx()] += dAdv[i]; }
			usr += usr_e;
		}
	}

	//eigen2ConnectedMatlab("Vn_en", Vn_en);
	//eigen2ConnectedMatlab("Vn_A", Vn_A);

	double dE = vn.dot(Vn_en.col(0));
	double dA = vn.dot(Vn_A);
	std::cout << "Int_em = " << Int_em << std::endl;
	std::cout << "dE = " << dE << std::endl;
	std::cout << "dA = " << dA << std::endl;
	std::cout << "usr = \n" << usr << std::endl;
}

Eigen::Vector<double, 9> retrieve_element_vector(PeriodSurfaceMesh& m, OM::SmartFaceHandle fh, const Eigen::VectorXd& u) {
	auto fvh = m.getFaceVertexHandle(fh);
	Eigen::Vector<double, 9> ue;
	for (int k = 0; k < 3; k++) { ue.block<3, 1>(k * 3, 0) = u.block<3, 1>(fvh[k].idx() * 3, 0); }
	return ue;
}

Eigen::Matrix3d retrieve_element_bform(PeriodSurfaceMesh& m, OM::SmartFaceHandle fh, const std::vector<Eigen::Matrix3d>& frlist, const std::vector<Compile1ring>& vrings) {
	Compile1ring v[3]; auto fvh = m.getFaceVertexHandle(fh);
	for (int k = 0; k < 3; k++) v[k] = vrings[fvh[k].idx()];
	return second_fundamental_form(m.getFacePeriodVertex(fh, 1), frlist[fh.idx()], v);
}

double calc_element_membrane_energy(
	double lam0, double mu,
	PeriodSurfaceMesh& m, OM::SmartFaceHandle fh, const Eigen::Matrix3d& fr,
	const std::vector<Compile1ring>& vrings,
	const Eigen::VectorXd& u, const Eigen::Matrix3d& eps_M
) {
	auto fvh = m.getFaceVertexHandle(fh);
	auto tri = m.getFacePeriodVertex(fh, 1);
	Compile1ring v[3];
	for (int k = 0; k < 3; k++) v[k] = vrings[fvh[k].idx()];
	Eigen::Vector<double, 9> ue;
	for (int k = 0; k < 3; k++) { ue.block<3, 1>(k * 3, 0) = u.block<3, 1>(fvh[k].idx() * 3, 0); }
	double ee = membrane_strain_energy(lam0, mu, tri, fr, v, ue, eps_M);
	return ee;
}

std::vector<double> read_double_number(std::string dfile) {
	std::ifstream ifs(dfile, std::ifstream::ate | std::ifstream::binary);
	size_t flen = ifs.tellg();
	ifs.seekg(std::ifstream::beg);
	std::vector<double> flist(flen / sizeof(double));
	ifs.read((char*)flist.data(), flen);
	return flist;
}


void test_refinement(void) {
	double Y = 1, nu = 0.3;
	auto [lam, mu] = mtl::lameCoeff(Y, nu);
	double lam0 = 2 * lam * mu / (lam + 2 * mu);
	std::cout << "lam = " << lam << ", lam0 = " << lam0 << ", mu = " << mu << std::endl;
	std::vector<std::string> mfilelist{
		//"D:/projects/minisurf/image/siggraph/refine-mesh/curtube03.stl",
		//"D:/projects/minisurf/image/siggraph/refine-mesh/curtube02.stl",
		//"D:/projects/minisurf/image/siggraph/refine-mesh/curtube01.stl"
		"D:/projects/minisurf/image/siggraph/refine-mesh/psub030.stl",
		"D:/projects/minisurf/image/siggraph/refine-mesh/psub015.stl"
	};
	std::vector<std::string> bvfilelist{
		"D:/projects/minisurf/image/siggraph/refine-mesh/psub030.bv",
		"D:/projects/minisurf/image/siggraph/refine-mesh/psub015.bv"
	};
	Eigen::Matrix3d eps_M; eps_M.setZero(); eps_M(0, 0) = 1;

	std::vector<PeriodSurfaceMesh> mlist;
	std::vector<std::vector<Eigen::Matrix3d>> bvlist;
	std::vector<PeriodicGridIndex> bvIndx;
	std::map<OM::SmartVertexHandle, OM::SmartVertexHandle> v2fine;
	std::map<OM::SmartEdgeHandle, OM::SmartVertexHandle> e2fine;
	std::map<OM::SmartFaceHandle, std::array<OM::SmartFaceHandle, 4>> f2fine;
	std::vector<Eigen::VectorXd> ulist;
	for (auto bfile : bvfilelist) {
		PeriodicGridIndex indx(Eigen::Vector3d(-1, -1, -1), Eigen::Vector3d(2, 2, 2), 1e-3);
		bvlist.emplace_back();
		bvIndx.emplace_back(indx);
		auto flist = read_double_number(bfile);
		for (int k = 0; k < flist.size(); k += 12) {
			Eigen::Vector3d p(flist[k], flist[k + 1], flist[k + 2]);
			Eigen::Matrix3d m;
			std::copy(&flist[k + 3], &flist[k + 12], m.data());
			bvlist.rbegin()->emplace_back(m);
			bvIndx.rbegin()->insert(p);
			int qid = bvIndx.rbegin()->query(p);
			bvlist.rbegin()->at(qid) = *bvlist.rbegin()->rbegin();
		}
	}
	for (auto mfile : mfilelist) {
		PeriodSurfaceMesh m;
		m.readMergePeriodBoundary(mfile, false);
		mlist.emplace_back(m);
	}
	// update bvlist
	for (int i = 0; i < bvlist.size(); i++) {
		std::vector<Eigen::Matrix3d> bvnew(mlist[i].n_vertices());
		for (auto vh : mlist[i].vertices()) {
			auto p = mlist[i].point(vh);
			int qid = bvIndx[i].query(p);
			if (qid == -1) { std::cout << "cannot find " << p << std::endl; }
			bvnew[vh.idx()] = bvlist[i][qid];
		}
		bvlist[i] = bvnew;
	}
	{
		std::ofstream ofs("bvread", std::ios::binary);

		std::vector<Compile1ring> vrings;
		for (auto vh : mlist[1].vertices()) { auto [o, ring] = mlist[1].find1ring(vh, 1); vrings.emplace_back(o, ring); }
		Eigen::Vector3d vr; vr.setRandom();
		Eigen::VectorXd nrm_err(mlist[1].n_vertices());
		for (auto vh : mlist[1].vertices()) {
			Eigen::Matrix<double, 3, 2> P;
			P.col(0) = vrings[vh.idx()].nv.cross(vr).normalized();
			P.col(1) = vrings[vh.idx()].nv.cross(P.col(0)).normalized();

			auto p = mlist[1].point(vh);
			ofs.write((const char*)p.data(), sizeof(p));
			auto vfr = bvlist[1][vh.idx()];
			nrm_err[vh.idx()] = (vfr * vrings[vh.idx()].nv).norm();
			if ((vfr - vfr.transpose()).norm() > 1e-12) { std::cout << "vfr not symmetric " << std::endl; }
			Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> so(P.transpose() * vfr * P);
			Eigen::Vector2d d = so.eigenvalues().cwiseAbs();
			Eigen::Matrix<double, 3, 2> V = P * so.eigenvectors() * d.asDiagonal();
			ofs.write((const char*)V.col(0).data(), sizeof(vr));
			ofs.write((const char*)V.col(1).data(), sizeof(vr));
			//ofs.write((const char*)vfr.col(0).data(), sizeof(vr));
			//ofs.write((const char*)vfr.col(1).data(), sizeof(vr));
		}
		std::sort(nrm_err.begin(), nrm_err.end());
		std::cout << "nrm_err min = " << nrm_err[0] << ", max = " << nrm_err[nrm_err.rows() - 1] << std::endl;
		eigen2ConnectedMatlab("nrmerr", nrm_err);
	}
	

	// fine-coarse topology information
	{
		PeriodicGridIndex idx(Eigen::Vector3d(-1, -1, -1), Eigen::Vector3d(2, 2, 2), 1e-3);
		for (auto vh : mlist.rbegin()->vertices()) {
			auto p = mlist.rbegin()->point(vh);
			idx.insert(p);
		}
		for (auto vh : mlist.begin()->vertices()) {
			auto p = mlist.begin()->point(vh);
			int qid = idx.query(p);
			if (qid == -1) { std::cout << "cannot find " << p << std::endl; }
			OM::SmartVertexHandle vnew(qid, &mlist[1]);
			v2fine[vh] = vnew;
		}
		for (auto eh : mlist.begin()->edges()) {
			auto he = eh.h0();
			auto v01 = make_period(mlist[0].calc_edge_vector(he), 2, 1);;
			auto p0 = mlist[0].point(he.from());
			auto pc = p0 + v01 / 2;
			for (int k = 0; k < 3; k++) { if (pc[k] < -1) pc[k] += 2; }
			int qid = idx.query(pc);
			if (qid == -1) { std::cout << "cannot find " << pc << std::endl; }
			OM::SmartVertexHandle vnew(qid, &mlist[1]);
			e2fine[eh] = vnew;
		}
		idx.clear();
		for (auto fh : mlist[1].faces()) {
			auto tri = mlist[1].getFacePeriodVertex(fh, 1);
			Eigen::Vector3d fc = tri.rowwise().mean();
			for (int k = 0; k < 3; k++) { if (fc[k] < -1) fc[k] += 2; }
			idx.insert(fc);
		}
		for (auto fh : mlist[0].faces()) {
			auto tri = mlist[0].getFacePeriodVertex(fh, 1);
			Eigen::Vector3d fc = tri.rowwise().mean();
			for (int k = 0; k < 3; k++) { if (fc[k] < -1) tri.row(k).array() += 2; }
			std::array<OM::SmartFaceHandle, 4> f4h;
			for (int k = 0; k < 3; k++) {
				Eigen::Vector3d c(1. / 6, 1. / 6, 1. / 6);
				c[k] = 2. / 3;
				Eigen::Vector3d p = tri * c;
				int qid = idx.query(p);
				if (qid == -1) { std::cout << "cannot find " << p.transpose() << std::endl; }
				f4h[k] = OM::SmartFaceHandle(qid, &mlist[1]);
			}
			Eigen::Vector3d p = tri.rowwise().mean(); int qid = idx.query(p);
			if (qid == -1) { std::cout << "cannot find " << p.transpose() << std::endl; }
			f4h[3] = OM::SmartFaceHandle(qid, &mlist[1]);
			f2fine[fh] = f4h;
		}
	}

	std::vector<std::vector<Compile1ring>> vringslist;
	std::vector<std::vector<Eigen::Matrix3d>> frlistlist;
	for (int mid = 0; mid < mlist.size(); mid++) {
		auto& m = mlist[mid];
		int ndof = m.n_vertices() * 3;
		std::vector<OM::SmartFaceHandle> fhlist(m.n_faces());
		auto [frlist, vrings] = generate_face_frame_and_vring(m);
		vringslist.push_back(vrings);
		for (auto fh : m.faces()) { fhlist[fh.idx()] = fh; }

#if 1
		auto [fm, um, Enmem, As] = asym_stif_fem(m, frlist, vrings, lam0, mu);
#else
		auto [fm, um, Enmem, As] = asym_stif_fem(m, frlist, vrings, lam0, mu, bvlist[mid]);
#endif

		frlistlist.push_back(frlist);
		ulist.emplace_back(um.col(0));

		Eigen::Matrix<double, 6, 6> CA;
		for (int i = 0; i < 6; i++) {
			for (int j = i; j < 6; j++) {
				CA(i, j) = (Enmem(i, j) - um.col(i).dot(fm.col(j))) / As;
				CA(j, i) = CA(i, j);
			}
		}
		std::cout << "CA = \n" << CA << std::endl;

		{
			std::ofstream ofs("Vp", std::ios::binary);
			for (auto vh : m.vertices()) {
				auto p = m.point(vh);
				ofs.write((const char*)p.data(), sizeof(p));
			}
			ofs.close(); ofs.open("um", std::ios::binary);
			ofs.write((const char*)um.col(0).data(), sizeof(double) * um.rows());
		}
	}

	

//	// interpolate coarse to fine
//	Eigen::VectorXd unew(mlist[1].n_vertices() * 3);
//	for (auto vh : mlist[0].vertices()) {
//		auto vh1 = v2fine[vh];
//		unew.block<3, 1>(vh1.idx() * 3, 0) = ulist[0].block<3, 1>(vh.idx() * 3, 0);
//	}
//	for (auto eh : mlist[0].edges()) {
//		auto he = eh.h0();
//		auto v0 = he.from();
//		auto v1 = he.to();
//		Eigen::Vector3d un = (ulist[0].block<3, 1>(v0.idx() * 3, 0) + ulist[0].block<3, 1>(v1.idx() * 3, 0)) / 2;
//		auto vh1 = e2fine[eh];
//		unew.block<3, 1>(vh1.idx() * 3, 0) = un;
//	}
//	{
//		std::ofstream ofs("um1", std::ios::binary);
//		ofs.write((const char*)unew.col(0).data(), sizeof(double)* unew.rows());
//	}
//	double Ecoarse2fine = 0;
//	for (auto fh : mlist[1].faces()) {
//		auto& vrings = vringslist[1];
//		auto& frlist = frlistlist[1];
//		auto fvh = mlist[1].getFaceVertexHandle(fh);
//		auto tri = mlist[1].getFacePeriodVertex(fh, 1);
//		Compile1ring v[3];
//		for (int k = 0; k < 3; k++) v[k] = vrings[fvh[k].idx()];
//		Eigen::Vector<double, 9> ue;
//		for (int k = 0; k < 3; k++) {
//			ue.block<3, 1>(k * 3, 0) = unew.block<3, 1>(fvh[k].idx() * 3, 0);
//		}
//		double ee = membrane_strain_energy(lam0, mu, tri, frlist[fh.idx()], v, ue, eps_M);
//		Ecoarse2fine += ee;
//	}
//
//	double As = 0;
//	for (auto vr : vringslist[1]) As += vr.As;
//	Ecoarse2fine /= As;
//	std::cout << "Ecoarse2fine = " << Ecoarse2fine << std::endl;
//
//	double Efine2coarse = 0;
//	unew.resize(mlist[0].n_vertices() * 3);
//	for (auto vh : mlist[0].vertices()) {
//		auto vh1 = v2fine[vh];
//		unew.middleRows(vh.idx() * 3, 3) = ulist[1].middleRows(vh1.idx() * 3, 3);
//	}
//	for (auto fh : mlist[0].faces()) {
//		auto& vrings = vringslist[0];
//		auto& frlist = frlistlist[0];
//		auto fvh = mlist[0].getFaceVertexHandle(fh);
//		auto tri = mlist[0].getFacePeriodVertex(fh, 1);
//		Compile1ring v[3];
//		for (int k = 0; k < 3; k++) v[k] = vrings[fvh[k].idx()];
//		Eigen::Vector<double, 9> ue;
//		for (int k = 0; k < 3; k++) {
//			ue.block<3, 1>(k * 3, 0) = unew.block<3, 1>(fvh[k].idx() * 3, 0);
//		}
//		double ee = membrane_strain_energy(lam0, mu, tri, frlist[fh.idx()], v, ue, eps_M);
//		Efine2coarse += ee;
//	}
//	Efine2coarse /= As;
//	std::cout << "Efine2coarse = " << Efine2coarse << std::endl;
//
//	std::vector<std::tuple<OM::SmartFaceHandle, double>> ferr;
//	for (auto fh : mlist[0].faces()) {
//		auto f4h = f2fine[fh];
//		double ee = calc_element_membrane_energy(lam0, mu, mlist[0], fh, frlistlist[0][fh.idx()], vringslist[0], ulist[0], eps_M);
//		double dncoarse = ee / frlistlist[0][fh.idx()].col(2).norm();
//		double e4 = 0, A4 = 0;
//		double dnfine = 0;
//		for (int k = 0; k < 4; k++) {
//			double ek = calc_element_membrane_energy(lam0, mu, mlist[1], f4h[k], frlistlist[1][f4h[k].idx()], vringslist[1], ulist[1], eps_M);
//			e4 += ek;
//			A4 += frlistlist[1][f4h[k].idx()].col(2).norm();
//		}
//		dnfine = e4 / A4;
//		ferr.emplace_back(fh, std::abs(dnfine - dncoarse));
//	}
//	std::sort(ferr.begin(), ferr.end(), [](const auto& f1, const auto& f2) {return std::get<1>(f1) < std::get<1>(f2); });
//	auto [fworst, err] = *ferr.rbegin();
//	{
//		auto f4h = f2fine[fworst];
//		double ee = calc_element_membrane_energy(lam0, mu, mlist[0], fworst, frlistlist[0][fworst.idx()], vringslist[0], ulist[0], eps_M);
//		double dncoarse = ee /*/ frlistlist[0][fworst.idx()].col(2).norm()*/;
//		std::cout << "dncoarse = " << dncoarse / frlistlist[0][fworst.idx()].col(2).norm() << ", err = " << err << std::endl;
//		std::cout << "bcoarse  = \n" << retrieve_element_bform(mlist[0], fworst, frlistlist[0], vringslist[0]) << std::endl;
//		auto tri_c = mlist[0].getFacePeriodVertex(fworst, 1);
//		auto u_c = retrieve_element_vector(mlist[0], fworst, ulist[0]);
//		std::ofstream ofs("badstrip");
//#if 0
//		ofs << tri_c.reshaped().transpose() << std::endl;
//		ofs << u_c.reshaped().transpose() << std::endl;
//		for (auto h4 : f4h) {
//			double e4 = calc_element_membrane_energy(lam0, mu, mlist[1], h4, frlistlist[1][h4.idx()], vringslist[1], ulist[1], eps_M);
//			std::cout << "bfine   = \n" << retrieve_element_bform(mlist[1], h4, frlistlist[1], vringslist[1]) << std::endl;
//			std::cout << "e4 fine = " << e4 / frlistlist[1][h4.idx()].col(2).norm() << std::endl;
//			auto tri_f = mlist[1].getFacePeriodVertex(h4, 1);
//			ofs << tri_f.reshaped().transpose() << std::endl;
//			auto u_f = retrieve_element_vector(mlist[1], h4, ulist[1]);
//			ofs << u_f.reshaped().transpose() << std::endl;
//#else
//		for (int i = 0; i < 100; i++) {
//			auto [fh, err] = ferr[ferr.size() - 1 - i];
//			auto tri = mlist[0].getFacePeriodVertex(fh, 1);
//			ofs << tri.reshaped().transpose() << std::endl;
//		}
//#endif
//	}
}

void test_bform_perturbation(void) {
	double Y = 1, nu = 0.3;
	auto [lam, mu] = mtl::lameCoeff(Y, nu);
	double lam0 = 2 * lam * mu / (lam + 2 * mu);
	std::cout << "lam = " << lam << ", lam0 = " << lam0 << ", mu = " << mu << std::endl;
	std::vector<std::string> mfilelist{
		"D:/projects/minisurf/image/siggraph/refine-mesh/cur1tube03.stl",
		"D:/projects/minisurf/image/siggraph/refine-mesh/cur1tube02.stl",
		"D:/projects/minisurf/image/siggraph/refine-mesh/cur1tube01.stl"
	};
	Eigen::Matrix3d eps_M; eps_M.setZero(); eps_M(0, 0) = 1;

	std::vector<PeriodSurfaceMesh> mlist;
	std::vector<std::vector<Eigen::Matrix3d>> bvlist;
	std::vector<PeriodicGridIndex> bvIndx;
	std::map<OM::SmartVertexHandle, OM::SmartVertexHandle> v2fine;
	std::map<OM::SmartEdgeHandle, OM::SmartVertexHandle> e2fine;
	std::map<OM::SmartFaceHandle, std::array<OM::SmartFaceHandle, 4>> f2fine;
	std::vector<Eigen::VectorXd> ulist;

	for (auto mfile : mfilelist) {
		PeriodSurfaceMesh m;
		m.readMergePeriodBoundary(mfile, false);
		mlist.emplace_back(m);
	}


	std::vector<std::vector<Compile1ring>> vringslist;
	std::vector<std::vector<Eigen::Matrix3d>> frlistlist;
	for (int mid = 0; mid < mlist.size(); mid++) {
		auto& m = mlist[mid];
		int ndof = m.n_vertices() * 3;
		std::vector<OM::SmartFaceHandle> fhlist(m.n_faces());
		auto [frlist, vrings] = generate_face_frame_and_vring(m);
		vringslist.push_back(vrings);
		for (auto fh : m.faces()) { fhlist[fh.idx()] = fh; }

#if 1
		auto [fm, um, Enmem, As] = asym_stif_fem(m, frlist, vrings, lam0, mu);
#else
		auto [fm, um, Enmem, As] = asym_stif_fem(m, frlist, vrings, lam0, mu, bvlist[mid]);
#endif

		frlistlist.push_back(frlist);
		ulist.emplace_back(um.col(0));

		Eigen::Matrix<double, 6, 6> CA;
		for (int i = 0; i < 6; i++) {
			for (int j = i; j < 6; j++) {
				CA(i, j) = (Enmem(i, j) - um.col(i).dot(fm.col(j))) / As;
				CA(j, i) = CA(i, j);
			}
		}
		std::cout << "CA = \n" << CA << std::endl;

		{
			std::ofstream ofs("Vp", std::ios::binary);
			for (auto vh : m.vertices()) {
				auto p = m.point(vh);
				ofs.write((const char*)p.data(), sizeof(p));
			}
			ofs.close(); ofs.open("um", std::ios::binary);
			ofs.write((const char*)um.col(0).data(), sizeof(double) * um.rows());
		}
	}
	
}

