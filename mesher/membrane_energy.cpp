#include <fstream>
#include "PeriodicMesher.h"
#include "material/materail.h"
#include "Eigen/PardisoSupport"
#include "matlab/matlab_utils.h"

using namespace msf;


extern void savePeriodicVector(std::string filename, msf::PeriodSurfaceMesh& mesh, Eigen::VectorXd& v);

Eigen::Matrix3d readStrain(std::string sfile) {
	std::ifstream ifs(sfile);
	if (!ifs) return Eigen::Matrix3d::Identity();
	Eigen::Matrix3d e;
	for (int i = 0; i < 9; i++) {
		ifs >> e.data()[i];
	}
	return e;
}

void solve_membrane_energy(PeriodSurfaceMesh& m, Real E, Real nu) {
	Eigen::Matrix3d eps_m = readStrain("epsm");
	eps_m = (eps_m + eps_m.transpose()).eval() / 2;
	std::cout << "eps_m = \n" << eps_m << std::endl;

	auto ff = m.getFaceFrame(2, 0.7);
	auto vnlist = m.getVertexNormal(2, 0.7);

	int ne = m.n_edges();
	int nv = m.n_vertices();
	//Eigen::SparseMatrix<double> K(ne + nv, ne + nv);
	Eigen::SparseMatrix<double> K(nv * 3, nv * 3);
	Eigen::VectorXd b(nv * 3); b.setZero();
	std::vector<Eigen::Triplet<double>> Kwise;
	auto D = mtl::planeElasticMatrix(E, nu);
	double Em = 0;
	double Am = 0;
	for (auto fh : m.faces()) {
		auto vt = m.getFacePeriodVertex(fh, 0.7);
		auto vh = m.getFaceVertexHandle(fh);
		auto h01 = m.findHe(vh[0], vh[1]);
		auto h12 = h01.next();
		auto h20 = h12.next();
		Eigen::Matrix3d fr = ff.block<3, 3>(fh.idx() * 3, 0).transpose();
		Eigen::Matrix<double, 3, 2> guv = fr.leftCols(2);
		Eigen::Vector3d v01_3 = vt.col(1) - vt.col(0);
		Eigen::Vector3d v02_3 = vt.col(2) - vt.col(0);
		Eigen::Vector3d v12_3 = v02_3 - v01_3;
		Eigen::Vector2d v01 = fr.leftCols(2).transpose() * v01_3; // e0
		Eigen::Vector2d v02 = fr.leftCols(2).transpose() * v02_3; // e1
		Eigen::Vector2d v12 = v02 - v01; // e2

		/****************************************************************/
		Eigen::Matrix3d Bt;
		Bt <<
			v01[0] * v01[0], v01[1] * v01[1], v01[0] * v01[1],
			v02[0] * v02[0], v02[1] * v02[1], v02[0] * v02[1],
			v12[0] * v12[0], v12[1] * v12[1], v12[0] * v12[1]; // no need to multiply 2 for 01-pair
		Bt = Bt.inverse().eval();
		/****************************************************************/
		//Eigen::Matrix3d vn;
		//vn << vnlist.row(vh[0].idx()), vnlist.row(vh[1].idx()), vnlist.row(vh[2].idx());
		//vn.transposeInPlace();
		//Eigen::Matrix3d Bn;
		//Bn << vn.col(0).dot(-v01), vn.col(0).dot(-v02), 0,
		//	vn.col(1).dot(v01), 0, vn.col(1).dot(-v12),
		//	0, vn.col(2).dot(v02), vn.col(2).dot(v12);
		//Bn = (Bt * Bn.transpose()).eval();
		/****************************************************************/
		Eigen::Matrix<double, 3, 9> Bv; Bv.setZero();
		auto ek = Eigen::Matrix3d::Identity();
		Bv.block<1, 3>(0, 0) = -v01_3.transpose();
		Bv.block<1, 3>(1,0) = -v02_3.transpose();
		Bv.block<1, 3>(0, 3) = v01_3.transpose();
		Bv.block<1, 3>(2, 3) = -v12_3.transpose();
		Bv.block<1, 3>(1, 6) = v02_3.transpose();
		Bv.block<1, 3>(2, 6) = v12_3.transpose();
		Bv = (Bt * Bv).eval();
		/****************************************************************/
#if 0
		//Eigen::Matrix<double, 6, 6> Ke; 
		//Eigen::Matrix<double, 3, 6> B; B << Bt, Bn;
		//Ke = B.transpose() * D * B;
#else
		Eigen::Matrix<double, 9, 9> Ke;
		Ke = Bv.transpose() * D * Bv;
		double A = (vt.col(1) - vt.col(0)).cross(vt.col(2) - vt.col(0)).norm() / 2;
		Ke *= A;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				for (int row = 0; row < 3; row++) {
					for (int col = 0; col < 3; col++) {
						Kwise.emplace_back(vh[i].idx() * 3 + row, vh[j].idx() * 3 + col, Ke(i * 3 + row, j * 3 + col));
					}
				}
			}
		}
#endif
		/****************************************************************/
		Eigen::Matrix2d ee = guv.transpose() * eps_m * guv;
		Eigen::Vector3d ee_vgt = { ee(0, 0), ee(1, 1), 2 * ee(0, 1) };
		Eigen::Vector3d sigm = D * ee_vgt;
		Eigen::Vector<double, 9> be = (sigm.transpose() * Bv).transpose();
		for (int i = 0; i < 3; i++) {
			for (int row = 0; row < 3; row++) {
				b[vh[i].idx() * 3 + row] -= be[i * 3 + row] * A;
			}
		}
		/****************************************************************/
		Em += sigm.dot(ee_vgt) * A;
		/****************************************************************/
		Am += A;
	}
	savePeriodicVector("bm", m, b);
	K.setFromTriplets(Kwise.begin(), Kwise.end());
	eigen2ConnectedMatlab("K", K);
	eigen2ConnectedMatlab("b", b);
	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> Kx(K);
	Eigen::VectorXd u = Kx.solve(b);
	eigen2ConnectedMatlab("u", u);
	savePeriodicVector("um", m, u);
	double Ea = Em - b.dot(u);
	std::cout << "Total area           :" << Am << std::endl;
	std::cout << "Pure Membrane Energy : " << Em << ", avg = " << Em / Am << std::endl;
	std::cout << "Total energy         :" << Ea <<  ", avg = " << Ea / Am << std::endl;
}

void test_membraine_energy(std::string meshfile) {
	PeriodSurfaceMesh m;
	m.readMergePeriodBoundary(meshfile);
	solve_membrane_energy(m, 1, 0.3);
}