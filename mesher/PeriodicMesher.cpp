#define _USE_MATH_DEFINES
#define EIGEN_USE_MKL_ALL
#define FMT_HEADER_ONLY
#include "PeriodicMesher.h"
//#include "nanoflann.hpp"
#include "PeriodMeshConvertion.h"
#include "grid/PeriodicGrid.h"
#include "unsupported/Eigen/AutoDiff"
#include "matlab/matlab_utils.h"
#include <bitset>
#include "Eigen/PardisoSupport"
//#include "AABB.hpp"
#include "cgal/homotopy_utils.h"
#include "igl/writeOBJ.h"
#include "igl/cotmatrix.h"
#include "fmt/core.h"
#include "igl/per_vertex_normals.h"
#include "igl/remove_duplicate_vertices.h"
#include "igl/readOBJ.h"
#include "igl/write_triangle_mesh.h"

std::string getPath(std::string s);

void msf::PeriodSurfaceMesh::setPeriod(void)
{
	// found boundary vertex
	std::map<OpenMesh::Vec3d, OpenMesh::VertexHandle> vblist;
	updateBb();
	auto bbscale = _bbox;
	bbscale.scale(1 - 1e-6);
	//std::vector<Point> plist, blist;
	for (auto v : vertices()) {
		if (is_boundary(v) && bbscale.isOut(toEigen(point(v)))) {
			vblist[point(v)] = v;
			//plist.push_back(point(v));
		}
	}
	//{ std::ofstream ofs(getPath("vblist")); for (int k = 0; k < plist.size(); k++) { ofs << plist[k] << std::endl; }}

	// map to 4D torus, and found overlapped point
	//auto p4list = nanoflann::embedTorusPoints(vbpos.data()->data(), vbpos.size());

	//using num_t = double;
	//using kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor <
	//	nanoflann::L2_Adaptor<num_t, nanoflann::PointCloud<num_t, 4>>,
	//	nanoflann::PointCloud<num_t, 4>, 4, size_t>;

	//kd_tree_t tree(4, p4list, nanoflann::KDTreeSingleIndexAdaptorParams());
	//tree.buildIndex();

	// rasterization
	std::unordered_map<Eigen::Vector3<size_t>, std::vector<OM::Vec3d>> buckets;
	double h = 1e-5;
	size_t N = 2 / h + 0.5;
	for (auto vb : vblist) {
		auto p = vb.first;
		Eigen::Vector3<size_t> pi;
		for (int k = 0; k < 3; k++) pi[k] = size_t((p[k] + 1) / h + 0.01) % N;
		buckets[pi].push_back(p);
	}
	
	// cluster overlapped point
	for (auto& pset : buckets) {
		auto& ps = pset.second;
		auto pi = pset.first;
		if (ps.empty()) continue;
		// check neighbor 
		Eigen::Vector3<size_t> pni;
		for (int ix = -1; ix <= 1; ix++) {
			pni[0] = (pi[0] + ix + N) % N;
			for (int iy = -1; iy <= 1; iy++) {
				pni[1] = (pi[1] + iy + N) % N;;
				for (int iz = -1; iz <= 1; iz++) {
					pni[2] = (pi[2] + iz + N) % N;;
					if (ix == 0 && iy == 0 && iz == 0) continue;
					if (buckets.count(pni)) {
						pset.second.insert(pset.second.end(), buckets[pni].begin(), buckets[pni].end());
						buckets[pni].clear();
					}
				}
			}
		}
	}

	_vflags.clear();
	_vpid.clear();
	_vdofid.clear();
	_vpid.clear();
	_vPeriodGroup.clear();
	_egperface.clear();
	_heTransfer.clear();
	_vflags.resize(n_vertices());
	// checkout cluster
	int periodVid = 0;
	for (auto pset : buckets) {
		const auto& ps = pset.second;
		auto pi = pset.first;
		if (ps.size() == 1 || ps.size() == 0) {
			continue;
		}
		else {
			OM::Vec3d mid = std::accumulate(ps.begin(), ps.end(), OM::Vec3d(0, 0, 0));
			mid /= ps.size();
			//for (int i = 0; i < ps.size(); i++) {
			//	std::cout << ps[i] << std::endl;
			//}
			//std::cout << mid << std::endl << std::endl;

			if (ps.size() == 2 || ps.size() == 4 || ps.size() == 8) { }
			// note that this situation potentialy occurred in normal case
			else {
				printf("Warning : unexpected number of overlaped points:\n");
				for (int k = 0; k < ps.size(); k++) { std::cout << ps[k] << std::endl; }
				//continue;
			}
			{
				int minid = -1;
				std::vector<OM::VertexHandle> vhper;
				for (int k = 0; k < ps.size(); k++) {
					auto omvid = vblist[ps[k]];
					vhper.emplace_back(omvid);
					int maxflg[3] = { ps[k][0] > mid[0] + 1e-6, ps[k][1] > mid[1] + 1e-6, ps[k][2] > mid[2] + 1e-6 };
					int minflg[3] = { ps[k][0] < mid[0] - 1e-6, ps[k][1] < mid[1] - 1e-6, ps[k][2] < mid[2] - 1e-6 };
					int flg[3] = { (minflg[0] ^ maxflg[0]) << maxflg[0], (minflg[1] ^ maxflg[1]) << maxflg[1],  (minflg[2] ^ maxflg[2]) << maxflg[2] };
					_vflags[omvid.idx()].set_period_boundary(flg[0], flg[1], flg[2]);
					_vpid[ps[k]] = periodVid;
					_vdofid[omvid] = periodVid;
					if (!maxflg[0] && !maxflg[1] && !maxflg[2]) {
						minid = omvid.idx();
					}
				}
				if (minid != -1) { }
				// this case potentialy occurrs
				else {
					printf("Warning : No representative (minimum) period vertex found\n");
				}
				periodVid++;
				_vPeriodGroup.emplace_back(vhper);
			}
		}
	}

	// output period vertex
	if (0) {
		std::ofstream ofs(getPath("vper"));
		for (int i = 0; i < _vPeriodGroup.size(); i++) {
			for (int j = 0; j < _vPeriodGroup[i].size(); j++) {
				ofs << point(_vPeriodGroup[i][j]) << " ";
			}
			ofs << std::endl;
		}
		ofs.close();

		ofs.open(getPath("doflist"));
		std::vector<std::vector<OM::Vec3d>> doflist(periodVid);
		for (auto v : vertices()){
			if (_vflags[v.idx()].is_period_boundary()) {
				doflist[_vdofid[v]].push_back(point(v));
			}
		}
		for (int i = 0; i < doflist.size(); i++) {
			for (int j = 0; j < doflist[i].size(); j++) {
				ofs << doflist[i][j] << " ";
			}
			ofs << std::endl;
		}
		ofs.close();
	}

	for (auto fh : faces()) {
		for (auto feh : fe_range(fh)) {
			auto v0h = feh.vertex(0);
			auto v1h = feh.vertex(1);
			if (_vflags[v0h.idx()].is_period_boundary() && _vflags[v1h.idx()].is_period_boundary()) {
				if (_egperface.count(feh)) {
					_egperface[feh].second = fh;
				}
				else {
					_egperface[feh].first = fh;
				}
			}
		}
	}

	// add other non period vertex dof
	for (auto v : vertices()) {
		if (_vflags[v.idx()].is_period_boundary()) continue;
		auto p = point(v);
		if (_vpid.count(p)) {
			_vdofid[v] = _vpid[p];
		} else {
			_vdofid[v] = periodVid;
			_vpid[p] = periodVid++;
		}
	}
	_nvdof = periodVid;

	// set periodic halfedge transfer
	for (auto vgroup : _vPeriodGroup) {
		for (auto vi : vgroup) {
			for (auto vj : vgroup) {
				if (vi == vj) continue;
				for (auto vioh : voh_ccw_range(vi)) {
					if (!vioh.is_boundary()) continue;
					auto vihe = this->calc_edge_vector(vioh);
					// find transfered halfedge
					for (auto vjoh : voh_ccw_range(vj)) {
#if 0
						if (vjoh.is_boundary()) continue;
						auto vjhe = this->calc_edge_vector(vjoh);
						bool paral = vihe.cross(vjhe).sqrnorm() < 1e-12;
						if (paral) _heTransfer[vioh] = vjoh;
#else
						if (vjoh.is_boundary()) continue; 
						auto vjhe = this->calc_edge_vector(vjoh);
						bool paral = (vihe - vjhe).sqrnorm() < 1e-12;
						if (paral) _heTransfer[vioh] = vjoh;
#endif
					}

				}
			}
		}
	}
}	

void msf::PeriodSurfaceMesh::normalizeBb(void)
{
	updateBb();
	Eigen::Vector3<Real> T = _bbox.getCenter();
	Eigen::Vector3<Real> S = Eigen::Vector3<Real>(2, 2, 2).cwiseQuotient(_bbox.diagonal());
	//std::cout << "scale = " << S.transpose() << std::endl;
	for (auto v : vertices()) {
		auto vp = point(v);
		for (int k = 0; k < 3; k++) {
			vp[k] = S[k] * (vp[k] - T[k]);
		// snap to boundary
			if (abs(vp[k]) < 1e-5) { vp[k] = 0; }
			if (abs(abs(vp[k]) - 1) < 1e-5) { vp[k] = vp[k] > 0 ? 1 : -1; }
		}
		point(v) = vp;
	}

	//OM::IO::write_mesh(*this, "unitmesh.obj");
}

bool msf::PeriodSurfaceMesh::read(std::string filepath, bool normalize, bool set_period, bool removeDup)
{
	try {
		OM::IO::read_mesh(*this, filepath);
	} catch (...) {
		printf("read mesh failed\n");
		return false;
	}
	//std::cout << "Read " << n_vertices() << " vertices, " << n_faces() << " faces" << std::endl;
	if (removeDup) remove_dup_vertices();
	//std::cout << "unduplicated " << n_vertices() << " vertices, " << n_faces() << " faces" << std::endl;
	//{ OM::IO::write_mesh(*this, getPath("removedup.obj")); }
	updateBb();
	if (normalize) { normalizeBb(); }
	if (set_period) setPeriod();

	request_normal_status();

	return true;
}

void msf::PeriodSurfaceMesh::faceVdof(OM::FaceHandle fh, int gid[3]) const
{
	int counter = 0;
	for (auto fv : fv_range(fh)) {
		gid[counter++] = vDofid(fv);
	}
}

bool msf::PeriodSurfaceMesh::save(std::string filepath)
{
	try {
		OM::IO::write_mesh(*this, filepath);
	}
	catch (...) {
		printf("OM write mesh failed\n");
		return false;
	}
	return true;
}

msf::OM::FaceHandle msf::PeriodSurfaceMesh::oppositePerFaces(OM::EdgeHandle eh, OM::FaceHandle fh) const
{
	if (!_egperface.count(eh))  return OM::FaceHandle();
	try {
		auto f2 = _egperface.at(eh);
		return f2.first == fh ? f2.second : f2.first;
	}
	catch (...) {
		return OM::FaceHandle();
	}
}

void msf::PeriodSurfaceMesh::saveVdofVector(std::string filename, const Eigen::VectorXd& u, Eigen::Vector3i beg, Eigen::Vector3i stride)
{
	std::ofstream ofs(filename, std::ios::binary);
	std::vector<Eigen::Vector3d> vdata(n_vertices());
	for (auto vh : vertices()) {
		int vdofid = vDofid(vh);
		Eigen::Vector3d u_v;
		for (int axis = 0; axis < 3; axis++) {
			u_v[axis] = u[beg[axis] + stride[axis] * vdofid];
		}
		vdata[vh.idx()] = u_v;
	}
	ofs.write((const char*)vdata.data(), sizeof(Eigen::Vector3d) * n_vertices());
	ofs.close();
}

void msf::PeriodSurfaceMesh::makeMinimal(void)
{
	Eigen::SparseMatrix<Real> gH(n_vdof() * 3, n_vdof() * 3);
	Eigen::VectorX<Real> H(n_vdof() * 3);
	H.setZero();
	constexpr int MaxAdj = 20;
	// newton's method
	for (int iter = 0; iter < 100; iter++) {
		if (iter % 10 == 0) {
			save("meshiter.obj");
		}
		//save("meshiter.obj");
		std::vector<Real> wcot(n_edges());
		std::vector<std::array<int, 4>> vdoflist(n_edges());
		std::vector<std::array<Eigen::Vector3<Real>, 4>> dofGradientlist(n_edges());
		// compute edge cot weight
		for (auto eh : edges()) {
			std::array<int, 4> vdof;
			std::array<Eigen::Vector3<Real>, 4> dofGradient;
			auto he = eh.halfedge(0);
			wcot[eh.idx()] = edgeCot(he, vdof.data(), dofGradient.data());
			//std::cout << "From " << point(he.from()) << ", vec = " << calc_edge_vector(he) << ", cot = " << wcot[eh.idx()] << std::endl;
			vdoflist[eh.idx()] = vdof;
			dofGradientlist[eh.idx()] = dofGradient;
		}
		std::vector<Eigen::Triplet<Real>> triplist;
		gH.setZero();
		for (auto vh : vertices()) {
			if (_vflags[vh.idx()].is_max_period()) continue;
			auto voh = *voh_begin(vh);
			if (voh.is_boundary()) voh = _heTransfer[voh];
			Real w_sum = 0;
			Eigen::Vector3<Real> w_sum_grad[MaxAdj];
			for (int k = 0; k < MaxAdj; k++) w_sum_grad[k].setZero();
			int nv_counter = 0;
			OM::SmartHalfedgeHandle vohlist[MaxAdj];
			int localvdoflist[MaxAdj];
			// compute weight sum
			do {
				vohlist[nv_counter] = voh;
				localvdoflist[nv_counter] = vDofid(voh.to());
				auto eh = voh.edge();
				Real w = wcot[eh.idx()];
				w_sum += w;
				if (voh.is_boundary()) { voh = _heTransfer[voh]; }
				// next halfedge
				voh = voh.prev().opp();
				nv_counter++;
				if (nv_counter > MaxAdj) throw std::runtime_error("adjacent buffer overflow");
			} while (voh != *voh_begin(vh));
			localvdoflist[MaxAdj - 1] = vDofid(vh);
			if (voh.is_boundary()) voh = _heTransfer[voh];
			// compute gradient
			Eigen::Vector3<Real> Hv(0, 0, 0);
			Eigen::Vector3<Real> hveclist[MaxAdj];
			for (int i = 0; i < nv_counter; i++) {
				voh = vohlist[i];
				auto hvec = this->calc_edge_vector(voh);
				hveclist[i] = Eigen::Vector3<Real>(hvec[0], hvec[1], hvec[2]);
				auto eh = voh.edge();
				Real w = wcot[eh.idx()];
				auto vdof = vdoflist[eh.idx()];
				auto dofGradient = dofGradientlist[eh.idx()];
				if (eh.halfedge(0) != voh) {
					std::swap(vdof[0], vdof[1]); std::swap(vdof[2], vdof[3]);
					std::swap(dofGradient[0], dofGradient[1]); std::swap(dofGradient[2], dofGradient[3]);
				}
				w_sum_grad[(i - 1 + nv_counter) % nv_counter] += dofGradient[3];
				w_sum_grad[i] += dofGradient[1];
				w_sum_grad[(i + 1) % nv_counter] += dofGradient[2];
				w_sum_grad[MaxAdj - 1] += dofGradient[0];
				// compute mean curvature 
				Hv += w / w_sum * hveclist[i];
				// next halfedge
				if (voh.is_boundary()) { voh = _heTransfer[voh]; }
				voh = voh.prev().opp();
			}
			H.block<3, 1>(vDofid(vh) * 3, 0) = Hv;
			// compute gradient
			Eigen::Matrix<Real, 3, MaxAdj * 3> localGrad;
			localGrad.setZero();
			for (int i = 0; i < nv_counter; i++) {
				voh = vohlist[i];
				auto hvec = hveclist[i];
				auto eh = voh.edge();
				auto vdof = vdoflist[eh.idx()];
				auto dofGradient = dofGradientlist[eh.idx()];
				Real w = wcot[eh.idx()];
				if (eh.halfedge(0) != voh) {
					std::swap(vdof[0], vdof[1]); std::swap(vdof[2], vdof[3]);
					std::swap(dofGradient[0], dofGradient[1]); std::swap(dofGradient[2], dofGradient[3]);
				}
				// weight gradient
				localGrad.block<3, 3>(0, ((i - 1 + nv_counter) % nv_counter) * 3) += hvec * (dofGradient[3].transpose() / w_sum);
				localGrad.block<3, 3>(0, i * 3) += hvec * (dofGradient[1].transpose() / w_sum);
				localGrad.block<3, 3>(0, ((i + 1) % nv_counter) * 3) += hvec * (dofGradient[2].transpose() / w_sum);
				localGrad.rightCols(3) += hvec * (dofGradient[0].transpose() / w_sum);
				// weight sum gradient
				Real wh = w / (w_sum * w_sum);
				for (int j = 0; j < nv_counter; j++) {
					localGrad.block<3, 3>(0, j * 3) -= (wh * hvec) * w_sum_grad[j].transpose();
				}
				//localGrad.block<3, 3>(0, i * 3) -= (wh * hvec) * w_sum_grad[i].transpose();
				localGrad.rightCols(3) -= (wh * hvec) * w_sum_grad[MaxAdj - 1].transpose();
				// out vector gradient 
				//for (int j = 0; j < nv_counter; j++) {
				//	localGrad.block<3, 3>(0, j * 3) += w / w_sum * Eigen::Matrix3<Real>::Identity();
				//}
				localGrad.block<3, 3>(0, i * 3) += w / w_sum * Eigen::Matrix3<Real>::Identity();
				localGrad.rightCols(3) -= w / w_sum * Eigen::Matrix3<Real>::Identity();
			}
			int vidof = localvdoflist[MaxAdj - 1];
			for (int i = 0; i < nv_counter; i++) {
				int vjdof = localvdoflist[i];
				for (int r = 0; r < 3; r++) {
					for (int c = 0; c < 3; c++) {
						triplist.emplace_back(vidof * 3 + r, vjdof * 3 + c, localGrad(r, c + i * 3));
					}
				}
			}
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) {
					triplist.emplace_back(vidof * 3 + r, vidof * 3 + c, localGrad(r, c + (MaxAdj - 1) * 3));
				}
			}
		}


		// setup Hessian matrix
		gH.setFromTriplets(triplist.begin(), triplist.end());

		//eigen2ConnectedMatlab("gH", gH);

		// solve equation
		//Eigen::SimplicialLDLT<Eigen::SparseMatrix<Real>> solver(gH);
#if 0
		//Eigen::SimplicialLDLT<Eigen::SparseMatrix<Real>> solver(gH.transpose()* gH);
		//Eigen::VectorX<Real> d = solver.solve(gH.transpose() * H);
		//Eigen::SparseLU<Eigen::SparseMatrix<Real>> solver(gH);
		Eigen::PardisoLU<Eigen::SparseMatrix<Real>> solver(gH);
		Eigen::VectorX<Real> d = solver.solve(H);
		if (d.dot(-H) < 0) {
			std::cout << "negative direction" << std::endl;
			d = -H;
		}
#else
		Eigen::VectorX<Real> d = -H;
#endif

		// project direction
		projectDirection(d);
		//eigen2ConnectedMatlab("d", d);

		//eigen2ConnectedMatlab("Hv", H);
		std::cout << "|H| = " << H.norm() << std::endl;
		if (iter % 10 == 0) {
			saveVdofVector("Hv", H, Eigen::Vector3i(0, 1, 2), Eigen::Vector3i(3, 3, 3));
		}
		//saveVdofVector("Hd", d, Eigen::Vector3i(0, 1, 2), Eigen::Vector3i(3, 3, 3));

		// move to new place
		auto Vnew = getVertexDofPosition();

		// test Gradient
		//auto H0 = meanCurvature(Vnew);
		//eigen2ConnectedMatlab("g0", H0);
		//auto H1 = meanCurvature(Vnew + 0.01 * Eigen::VectorX<Real>::Unit(Vnew.rows(), 80000));
		//eigen2ConnectedMatlab("g1", H1);

		Real alpha = searchStepMeanCurvature(Vnew, d);
		std::cout << "Newton step = " << alpha << std::endl;
		if (alpha < 1e-10) break;

		Vnew -=  alpha * d;
		setVertexDofPosition(Vnew);
	}
}

double msf::PeriodSurfaceMesh::edgeCot(OM::SmartHalfedgeHandle he, int vdof[4], Eigen::Vector3<Real> dofGradient[4])
{
	auto ophe = he.opp();
	OM::HalfedgeHandle heh0 = he, heh1 = ophe;
	if (he.is_boundary()) {
		heh0 = he.opp();
		heh1 = _heTransfer[he];
	} else if(ophe.is_boundary()){
		heh0 = he;
		heh1 = _heTransfer[ophe];
	}

	OM::VertexHandle vi0 = to_vertex_handle(heh0);
	OM::VertexHandle vj0 = from_vertex_handle(heh0);
	OM::VertexHandle vk0 = to_vertex_handle(next_halfedge_handle(heh0));

	OM::VertexHandle vi1 = to_vertex_handle(heh1);
	OM::VertexHandle vj1 = from_vertex_handle(heh1);
	OM::VertexHandle vk1 = to_vertex_handle(next_halfedge_handle(heh1));

	// vj0 vi0 vk0 vk1 
	vdof[0] = vDofid(vj0);
	vdof[1] = vDofid(vi0);
	vdof[2] = vDofid(vk0);
	vdof[3] = vDofid(vk1);

	OM::Vec3d pi0 = point(vi0);
	OM::Vec3d pj0 = point(vj0);
	OM::Vec3d pk0 = point(vk0);

	OM::Vec3d pi1 = point(vi1);
	OM::Vec3d pj1 = point(vj1);
	OM::Vec3d pk1 = point(vk1);

	OM::Vec3d vjk0 = pk0 - pj0;
	OM::Vec3d vik0 = pk0 - pi0;

	OM::Vec3d vjk1 = pk1 - pj1;
	OM::Vec3d vik1 = pk1 - pi1;

	Eigen::Vector3<Real> vjk0_vec(vjk0[0], vjk0[1], vjk0[2]), vjk0_grad;
	Eigen::Vector3<Real> vik0_vec(vik0[0], vik0[1], vik0[2]), vik0_grad;
	Eigen::Vector3<Real> vjk1_vec(vjk1[0], vjk1[1], vjk1[2]), vjk1_grad;
	Eigen::Vector3<Real> vik1_vec(vik1[0], vik1[1], vik1[2]), vik1_grad;

	//std::cout << "pi0 = " << pi0 << ", pj0 = " << pj0 << ", pk0 = " << pk0 << std::endl;
	//std::cout << "pi1 = " << pi1 << ", pj1 = " << pj1 << ", pk1 = " << pk1 << std::endl;
	//std::cout << "vjk0_vec = " << vjk0_vec.transpose() << ", vik0_vec = " << vik0_vec.transpose() << std::endl;
	//std::cout << "vjk1_vec = " << vjk1_vec.transpose() << ", vik1_vec = " << vik1_vec.transpose() << std::endl;

	double cot_alpha_i0 = _cot_gradient(vjk0_vec, vik0_vec, vjk0_grad, vik0_grad);
	double cot_alpha_i1 = _cot_gradient(vjk1_vec, vik1_vec, vjk1_grad, vik1_grad);

	dofGradient[0] = (-vjk0_grad - vik1_grad) / 2;
	dofGradient[1] = (-vik0_grad - vjk1_grad) / 2;
	dofGradient[2] = (vjk0_grad + vik0_grad) / 2;
	dofGradient[3] = (vik1_grad + vjk1_grad) / 2;

	return 0.5 * (cot_alpha_i0 + cot_alpha_i1);
}

msf::Real msf::PeriodSurfaceMesh::edgeCot(OM::SmartHalfedgeHandle he, Real period, Real deduce)
{
	auto vhe_sq = make_period(point(he.to()) - point(he.from()), period, deduce).sqrnorm();
	auto fv1 = he.next().to();
	auto fe1a_sq = (make_period(point(he.to()) - point(fv1))).sqrnorm();
	auto fe1b_sq = (make_period(point(he.from()) - point(fv1))).sqrnorm();
	auto fv2 = he.opp().next().to();
	auto fe2a_sq = (make_period(point(he.from()) - point(fv2))).sqrnorm();
	auto fe2b_sq = (make_period(point(he.to()) - point(fv2))).sqrnorm();
	Real cos1 = (fe1a_sq + fe1b_sq - vhe_sq) / (2 * std::sqrt(fe1a_sq * fe1b_sq));
	Real cos2 = (fe2a_sq + fe2b_sq - vhe_sq) / (2 * std::sqrt(fe2a_sq * fe2b_sq));
	Real cot1 = cos1 / std::sqrt(1 - cos1 * cos1);
	Real cot2 = cos2 / std::sqrt(1 - cos2 * cos2);
	return (cot1 + cot2) / 2;
}

msf::OM::HalfedgeHandle msf::PeriodSurfaceMesh::closedOpposite(OM::SmartHalfedgeHandle he)
{
	auto ophe = he.opp();
	if (!ophe.is_boundary()) {
		return ophe;
	} else {
		return _heTransfer[ophe];
	}
}

msf::Real msf::PeriodSurfaceMesh::_cot_gradient(const Eigen::Vector3<Real>& a, const Eigen::Vector3<Real>& b, Eigen::Vector3<Real>& agrad, Eigen::Vector3<Real>& bgrad) const
{
	using ADScalar = Eigen::AutoDiffScalar<Eigen::Vector<Real, 6>>;
	Eigen::Vector3<ADScalar> a_AD(a);
	Eigen::Vector3<ADScalar> b_AD(b);
	for (int i = 0; i < 3; i++) {
		a_AD[i].derivatives() = Eigen::Vector<Real, 6>::Unit(i);
		b_AD[i].derivatives() = Eigen::Vector<Real, 6>::Unit(i + 3);
	}
	ADScalar fcot = a_AD.dot(b_AD) / a_AD.cross(b_AD).norm();
	agrad = Eigen::Vector3<Real>(fcot.derivatives().topRows(3));
	bgrad = Eigen::Vector3<Real>(fcot.derivatives().bottomRows(3));
	return fcot.value();
}

Eigen::VectorX<msf::Real> msf::PeriodSurfaceMesh::getVertexDofPosition(void) const
{
	Eigen::VectorX<Real> V(n_vdof() * 3, 1);
	V.setZero();
	for (auto vh : vertices()) {
		if (_vflags[vh.idx()].is_max_period()) continue;
		int vdof = vDofid(vh);
		auto vp = point(vh);
		V.block<3, 1>(vdof * 3, 0) = Eigen::Vector3<Real>(vp[0], vp[1], vp[2]);
	}
	return V;
}

void msf::PeriodSurfaceMesh::setVertexDofPosition(const Eigen::VectorX<Real>& newV) {
	for (auto vh : vertices()) {
		if (_vflags[vh.idx()].is_max_period()) continue;
		int vdof = vDofid(vh);
		OM::Vec3d new_p(newV[vdof * 3], newV[vdof * 3 + 1], newV[vdof * 3 + 2]);
		this->set_point(vh, new_p);
	}
	for (auto vgroup : _vPeriodGroup) {
		OM::Vec3d vorbit;
		for (auto vh : vgroup) {
			if (!_vflags[vh.idx()].is_max_period()) {
				vorbit = point(vh);
				break;
			}
		}
		//std::cout << "vorbit = " << vorbit << std::endl;
		for (int i = 0; i < vgroup.size(); i++) {
			//if (!_vflags[vgroup[i].idx()].is_max_period()) continue;
			OM::Vec3d newp = vorbit;
			for (int axis = 0; axis < 3; axis++) {
				if (_vflags[vgroup[i].idx()].is_max_period(axis)) {
					newp[axis] = 1;
				}
				if (_vflags[vgroup[i].idx()].is_min_period(axis)) {
					newp[axis] = -1;
				}
			}
			//std::cout << "oldp = " << point(vgroup[i]) << ", newp = " << newp << ", flag = " << std::bitset<32>(_vflags[vgroup[i].idx()]._flag) << std::endl;
			this->set_point(vgroup[i], newp);
		}
	}
}

// only consider 3dof per-vertex
void msf::PeriodSurfaceMesh::projectDirection(Eigen::VectorX<Real>& d)
{
	Eigen::Matrix<Real, 3, -1> dmat = d.reshaped(3, d.size() / 3);
	Eigen::Vector3<Real> dmean = dmat.rowwise().sum();
	dmean /= dmat.cols();
	dmat.colwise() -= dmean;
	d = dmat.reshaped();
	saveVdofVector("dbalence", d, Eigen::Vector3i(0, 1, 2), Eigen::Vector3i(3, 3, 3));
	for (auto vgroup : _vPeriodGroup) {
		auto v0 = vgroup[0];
		int vdof = vDofid(v0);
		// remove x component
		for (int axis = 0; axis < 3; axis++) {
			if (_vflags[v0.idx()].is_period_boundary(axis)) {
				d[vdof * 3 + axis] = 0;
			}
		}
	}
}

Eigen::VectorX<msf::Real> msf::PeriodSurfaceMesh::meanCurvature(const Eigen::VectorX<Real>& v)
{
	Eigen::VectorX<Real> ubak = getVertexDofPosition();

	// update vertex position
	setVertexDofPosition(v);

	Eigen::VectorX<Real> H(n_vdof() * 3, 1);

	for (auto vh : vertices()) {
		if (_vflags[vh.idx()].is_max_period()) continue;
		auto voh = *voh_begin(vh);
		if (voh.is_boundary()) voh = _heTransfer[voh];
		Real wsum = 0;
		Eigen::Vector3<Real> Hv(0, 0, 0);
		do {
			int vdof[4];
			Eigen::Vector3<Real> dofgrad[4];
			auto hvec = calc_edge_vector(voh);
			Real w = edgeCot(voh, vdof, dofgrad);
			wsum += w;
			Hv += w * Eigen::Vector3<Real>(hvec[0], hvec[1], hvec[2]);
			if (voh.is_boundary()) voh = _heTransfer[voh];
			voh = voh.prev().opp();
		} while (voh != *voh_begin(vh));
		Hv /= wsum;
		int vidof = vDofid(vh);
		H.block<3, 1>(vidof * 3, 0) = Hv;
	}

	setVertexDofPosition(ubak);

	return H;
}

msf::Real msf::PeriodSurfaceMesh::searchStepMeanCurvature(const Eigen::VectorX<Real>& u, const Eigen::VectorX<Real>& d)
{
	Real step = 1;
	auto H0 = meanCurvature(u);
	Real H0norm = H0.norm();
	do {
		Eigen::VectorX<Real> u_step = u - step * d;
		auto H = meanCurvature(u_step);
		Real Hnorm = H.norm();
		if (Hnorm < H0norm) { break; }
		step /= 2;
	} while (step > 1e-12);
	return step;
}

msf::VertexFlag& msf::PeriodSurfaceMesh::getVFlag(OM::SmartVertexHandle vh)
{
	return _vflags[vh.idx()];
}

void msf::PeriodSurfaceMesh::extractAsymptoticLines(const std::vector<Eigen::Vector3d>& seeds, std::vector<std::vector<Eigen::Vector3d>>& plylines, std::vector<int>& faceid)
{
	
}

void msf::PeriodSurfaceMesh::initAABBTree(void)
{

}

void msf::PeriodSurfaceMesh::extractHarmonicField(std::vector<Eigen::VectorX<Real>>& val_vertex, std::vector<Eigen::VectorX<Real>>& vec_face)
{
	val_vertex.clear();
	vec_face.clear();
	auto L = getLaplacian();
	eigen2ConnectedMatlab("L", L);
	Eigen::SparseQR<Eigen::SparseMatrix<Real>, Eigen::COLAMDOrdering<int>> qr(L);
	std::vector<Eigen::VectorX<Real>> null_space;
	for (int i = 0; i < L.cols() - qr.rank(); i++) {
		Eigen::VectorXd e = Eigen::VectorXd::Unit(L.cols(), L.cols() - i);
		null_space.emplace_back(qr.matrixQ() * e);
	}
	for (int i = 0; i < null_space.size(); i++) {
		auto u_vertex = assembleVertexFromDof(null_space[i]);
		val_vertex.push_back(u_vertex);
		Eigen::VectorX<Real> vec_f(n_faces() * 3, 1);
		for (auto fh : faces()) {
			auto fv = getFaceVertexHandle(fh);
			auto fvpos = getFaceVertex(fh);
			Real valist[3] = { u_vertex[fv[0].idx()],u_vertex[fv[1].idx()],u_vertex[fv[2].idx()] };
			auto n_f = normal(fh);
			Eigen::Vector3<Real> n(n_f[0], n_f[1], n_f[2]);
			Eigen::Matrix3<Real> A;
			A << (fvpos.col(1) - fvpos.col(0)).transpose(),
				(fvpos.col(2) - fvpos.col(0)).transpose(),
				n.transpose();
			Eigen::Vector3<Real> g = A.lu().solve(Eigen::Vector3<Real>(valist[1] - valist[0], valist[2] - valist[0], 0));
			vec_f.block<3, 1>(fh.idx() * 3, 0) = g;
		}
		vec_face.emplace_back(vec_f);
	}
}

Eigen::SparseMatrix<msf::Real> msf::PeriodSurfaceMesh::getLaplacian(void)
{
	Eigen::SparseMatrix<Real> L(n_vdof(), n_vdof());

	std::vector<Eigen::Triplet<Real>> triplist;

	for (auto vh : vertices()) {
		if (_vflags[vh.idx()].is_max_period())  continue;
		int vdofid = vDofid(vh);
		auto voh = *voh_begin(vh);
		Real w_sum = 0;
		if (voh.is_boundary()) voh = _heTransfer[voh];
		do {
			auto v2 = voh.to();
			int vdof[4];
			Eigen::Vector3d dofgrad[4];
			auto w = edgeCot(voh, vdof, dofgrad);
			w_sum += w;
			int v2id = vDofid(v2);
			triplist.emplace_back(vdofid, v2id, w);
			if (voh.is_boundary()) { voh = _heTransfer[voh]; }
			// next halfedge
			voh = voh.prev().opp();
		} while (voh != *voh_begin(vh));
		triplist.emplace_back(vdofid, vdofid, -w_sum);
	}
	L.setFromTriplets(triplist.begin(), triplist.end());
	return L;
}

Eigen::VectorX<msf::Real> msf::PeriodSurfaceMesh::assembleVertexFromDof(const Eigen::VectorX<Real>& u)
{
	Eigen::VectorX<Real> uVertex(n_vertices() * 6, 1);
	uVertex.setZero();
	for (auto fh : faces()) {
		bool contains_max_period = false;
		VertexFlag vflag[3];
		int gid[3];
		int vid[3];
		int counter = 0;
		for (auto fv : fv_range(fh)) {
			vid[counter] = fv.idx();
			vflag[counter] = getVFlag(fv);
			gid[counter] = vDofid(fv);
			auto vp = point(fv);
			if (vflag[counter].is_max_period()) {
				contains_max_period = true;
			}
			counter++;
		}
		for (int i = 0; i < 3; i++) {
			uVertex.block<6, 1>(vid[i] * 6, 0) = u.block<6, 1>(gid[i] * 6, 0);
		}
	}
	return uVertex;
}

std::vector<msf::Real> msf::PeriodSurfaceMesh::edgeCotList(void)
{
	std::vector<Real> cotlist(n_edges());
	for (auto eh : edges()) {
		auto he = eh.h0();
		int vdof[4];
		Eigen::Vector3d grad[4];
		cotlist[eh.idx()] = edgeCot(he, vdof, grad);
	}
	return cotlist;
}


std::vector<Eigen::VectorXd> msf::PeriodSurfaceMesh::inextensionalDisplacement(void)
{
	
}

void msf::PeriodSurfaceMesh::extendPeriod(double ratio)
{
	std::vector<std::array<OM::Vec3d, 3>> ext;
	for (auto f : faces()) {
		auto fc = MeshSurface::calc_face_centroid(f);
		auto param = _bbox.localCoords(fc[0], fc[1], fc[2]);
		double off[3] = { 0,0,0 };
		for (int i = 0; i < 3; i++) {
			if (param[i] < ratio) {
				off[i] = 1;
			} else if (param[i] > 1 - ratio) {
				off[i] = -1;
			}
		}
		auto d = _bbox.diagonal();
		for (int z = 0; z < (off[2] != 0) + 1; z++) {
			for (int y = 0; y < (off[1] != 0) + 1; y++) {
				for (int x = 0; x < (off[0] != 0) + 1; x++) {
					OM::Vec3d poff(off[0] * d[0], off[1] * d[1], off[2] * d[2]);
					if (z == 0) poff[2] = 0;
					if (y == 0) poff[1] = 0;
					if (x == 0) poff[0] = 0;
					if (poff[0] == 0 && poff[1] == 0 && poff[2] == 0) { continue; }
					int counter = 0;
					OM::Vec3d fp[3];
					for (auto fvit : f.vertices()) {
						if (counter < 3) {
							fp[counter++] = MeshSurface::point(fvit) + poff;
						}
					}
					ext.push_back(std::array<OM::Vec3d, 3>{fp[0], fp[1], fp[2]});
				}
			}
		}
	}
	std::cout << "Found " << ext.size() << " faces on borderline" << std::endl;
	for (int k = 0; k < ext.size(); k++) {
		OM::VertexHandle vh[3];
		for (int j = 0; j < 3; j++) { vh[j] = MeshSurface::add_vertex(ext[k][j]); }
		MeshSurface::add_face(std::vector<OM::VertexHandle>{vh[0], vh[1], vh[2]});
	}
}

void msf::PeriodSurfaceMesh::updateBb(void)
{
	std::vector<msf::Point> bp;
	for (auto v : vertices()) {
		auto vp = point(v);
		bp.emplace_back(msf::Point(vp[0], vp[1], vp[2]));
	}
	_bbox = BBox(bp);
}

void msf::PeriodSurfaceMesh::embedMesh(std::string filename, std::vector<int> basepoints)
{
	OM::IO::read_mesh(*this, filename);
	auto [v, f] = getVFlist();
	{
		std::ofstream ofs(getPath("vlist")); ofs << v; ofs.close();
		ofs.open(getPath("flist")); ofs << f; ofs.close();
	}
	auto cyclelist = extract_shortest_homotopy_path(v, f, basepoints);
	{
		std::ofstream ofs(getPath("cyclelist"));
		for (int i = 0; i < cyclelist.size(); i++) {
			ofs << "cycle ";
			for (int k = 0; k < cyclelist[i].size(); k++) {
				ofs << cyclelist[i][k] << " ";
			}
			ofs << std::endl;
		}
	}
	std::vector<int> handle_tag(cyclelist.size(), 0);
	//handle_tag[0] = 1; handle_tag[1] = 2; handle_tag[2] = 3;
	handle_tag[0] = 1;
	relaxInTorus(cyclelist, handle_tag, 3);
}

//std::pair<Eigen::Matrix<msf::Real, -1, 3>, Eigen::Matrix<int, -1, 3>> msf::PeriodSurfaceMesh::getVFlist(void)
//{
//	Eigen::Matrix<Real, -1, 3> vlist(MeshSurface::n_vertices(), 3);
//	Eigen::Matrix<int, -1, 3> flist(n_faces(), 3);
//	int counter = 0;
//	for (auto v : MeshSurface::vertices()) {
//		auto pos = point(v);
//		for (int k = 0; k < 3; k++) { vlist(counter, k) = pos[k]; }
//		counter++;
//	}
//	counter = 0;
//	for (auto f : MeshSurface::faces()) {
//		if (status(f).deleted()) {
//			printf("\033[31mWarning ! Face has been deleted\033[0m\n");
//		}
//		//if (!f.is_valid()) continue;
//		auto fvlist = getFaceVertexHandle(f);
//		if (fvlist[0].is_valid() && fvlist[1].is_valid() && fvlist[2].is_valid()) {
//			for (int k = 0; k < 3; k++) flist(counter, k) = fvlist[k].idx();
//			counter++;
//		}
//	}
//	flist.conservativeResize(counter, 3);
//	return { vlist,flist };
//}

void msf::PeriodSurfaceMesh::relaxInTorus(const std::vector<std::vector<size_t>>& cutlocus, const std::vector<int>& handle_tag, int type)
{
	if (type == 2) {
		relaxInTorusH2(cutlocus, handle_tag);
	}
	else if (type == -1) {
		relaxInTorusLaplacian(cutlocus, handle_tag);
	}
	else if (type == 3) {
		relaxInTorusARAP(cutlocus, handle_tag);
	}
}

std::vector<Eigen::Vector3<msf::Real>> msf::PeriodSurfaceMesh::scaleTo(const BBox& new_bb) {
	std::vector<Eigen::Vector3<Real>> new_vlist;
	scaleTo(new_bb, new_vlist);
	for (int i = 0; i < n_vertices(); i++) {
		set_point(OM::VertexHandle(i), Point(new_vlist[i][0], new_vlist[i][1], new_vlist[i][2]));
	}
	return new_vlist;
}

void msf::PeriodSurfaceMesh::scaleTo(const BBox& new_bb, std::vector<Eigen::Vector3<Real>>& new_vlist)
{
	updateBb();
	Real s[3], b[3];
	auto dnew = new_bb.diagonal();
	auto bnew = new_bb.getCorner(0);
	auto dold = _bbox.diagonal();
	auto bold = _bbox.getCorner(0);
	for (int i = 0; i < 3; i++) {
		s[i] = dnew[i] / dold[i];
	}
	Real smin = *std::min_element(s, s + 3);
	s[0] = s[1] = s[2] = smin;
	for (int i = 0; i < 3; i++) {
		b[i] = (bnew[i] + dnew[i] / 2) - smin * (bold[i] + dold[i] / 2);
	}
	new_vlist.clear();
	for (auto v : vertices()) {
		auto p = point(v);
		for (int k = 0; k < 3; k++) { p[k] = s[k] * p[k] + b[k]; }
		new_vlist.push_back(Eigen::Vector3<Real>(p[0], p[1], p[2]));
	}
	return;
}

void msf::PeriodSurfaceMesh::scaleTo(const BBox& new_bb, Eigen::MatrixX3<Real>& new_vlist)
{
	std::vector<Eigen::Vector3<Real>> newvec(n_vertices());
	scaleTo(new_bb, newvec);
	new_vlist.resize(n_vertices(), 3);
	for (int i = 0; i < new_vlist.rows(); i++) {
		new_vlist.row(i) = newvec[i].transpose();
	}
}

#define USE_ANGLE_H2

std::tuple<msf::Real, msf::Real> msf::PeriodSurfaceMesh::meanH2(Eigen::Vector3<msf::Real>& o, std::vector<Eigen::Vector3<msf::Real>>& ring, bool area_averaged)
{
	using ADScalar = Eigen::AutoDiffScalar<Eigen::VectorX<Real>>;
	std::vector<Eigen::Vector3<ADScalar>> ring_ad(ring.size());
	Eigen::Vector3<ADScalar> o_ad;
	for (int i = 0; i < ring.size(); i++) {
		for (int k = 0; k < 3; k++) {
			ring_ad[i][k].derivatives() = Eigen::VectorX<Real>::Unit((ring.size() + 1) * 3, i * 3 + k);
			ring_ad[i][k].value() = ring[i][k];
		}
	}
	for (int k = 0; k < 3; k++) {
		o_ad[k].derivatives() = Eigen::VectorX<Real>::Unit((ring.size() + 1) * 3, ring.size() * 3 + k);
		o_ad[k].value() = o[k];
	}
#ifndef USE_ANGLE_H2
	std::vector<ADScalar> oh2, eh2;
	ADScalar A = 0;
	A.derivatives() = Eigen::VectorX<Real>::Zero((ring.size() + 1) * 3);
	for (int i = 0; i < ring.size(); i++) {
		oh2.push_back((ring_ad[i] - o_ad).squaredNorm());
		eh2.push_back((ring_ad[(i + 1) % ring.size()] - ring_ad[i]).squaredNorm());
		if (area_averaged) {
			A += (o_ad - ring_ad[i]).cross(ring_ad[(i + 1) % ring.size()] - ring_ad[i]).norm() / 6;
		} else {
			A.value() += (o - ring[i]).cross(ring[(i + 1) % ring.size()] - ring[i]).norm() / 6;
		}
	}
	Eigen::Vector3<ADScalar> Hv;
	for (int i = 0; i < 3; i++) {
		Hv[i].value() = 0;
		Hv[i].derivatives() = Eigen::VectorX<Real>::Zero((ring.size() + 1) * 3);
	}
	for (int i = 0; i < ring.size(); i++) {
		ADScalar cos_a = (oh2[i] + eh2[i] - oh2[(i + 1) % ring.size()]) / (2 * Eigen::sqrt(oh2[i] * eh2[i]));
		ADScalar cos_c = (oh2[(i + 2) % ring.size()] + eh2[(i + 1) % ring.size()] - oh2[(i + 1) % ring.size()]) / (2 * Eigen::sqrt(oh2[(i + 2) % ring.size()] * eh2[(i + 1) % ring.size()]));
		if (std::isnan(cos_a.value()) || std::isnan(cos_c.value())) {
			std::cout << i << "-th cos is Nan" << std::endl;
		}
		ADScalar cot_a = cos_a / Eigen::sqrt(1 - cos_a * cos_a);
		ADScalar cot_c = cos_c / Eigen::sqrt(1 - cos_c * cos_c);
		if (std::isnan(cot_a.value()) || std::isnan(cot_c.value())) {
			std::cout << i << "-th cot is Nan, " << "cos = " << cos_a.value() << ", " << cos_c.value() << std::endl;
		}
		Hv += (cot_a + cot_c) / 2 * (ring_ad[(i + 1) % ring.size()] - o_ad); // Note that this should be i+1 not i !!!
	}
	ADScalar H2 = Hv.squaredNorm();
	//A = A * A;
	//if (area_averaged) { H2 /= A; }
	for (int i = 0; i < ring.size(); i++) { ring[i] = H2.derivatives().block<3, 1>(i * 3, 0); }
	o = H2.derivatives().bottomRows(3);

	return { H2.value(), A.value() };
#else
	std::vector<Eigen::Vector3<ADScalar>> v_oh(ring.size()), v_eh(ring.size());
	for (int i = 0; i < ring.size(); i++) {
		v_oh[i] = (ring_ad[i] - o_ad).normalized();
		v_eh[i] = (ring_ad[(i + 1) % ring.size()] - ring_ad[i]).normalized();
	}

	ADScalar H2 = 0;
	H2.derivatives() = Eigen::VectorX<Real>::Zero((ring.size() + 1) * 3);
	int N = ring.size();
	for (int i = 0; i < ring.size(); i++) {
		const auto& a = v_eh[i];
		const auto& b = v_eh[(i + 1) % N];
		const auto& c = -v_oh[(i + 2) % N];
		const auto& d = v_oh[i];
		ADScalar cos_beta = a.dot(c) * b.dot(d) - a.dot(b) * (c).dot(d) - b.dot(c) * d.dot(a);
		if (1 - std::pow(cos_beta.value(), 2) < 1e-6 * 1e-6) {
			H2.value() += std::acos(std::clamp(cos_beta.value(), -1 + 1e-6, 1 - 1e-6));
		} else {
			H2 += Eigen::acos(cos_beta);
		}
		if (std::isnan(H2.value()) || H2.derivatives().hasNaN()) {
			std::cout << "\033[31m" << "H2 has NaN" << "\033[0m" << std::endl;
			std::cout << "dH2 = " << H2.derivatives().transpose() << std::endl;
			std::cout << "cos_beta = " << cos_beta.value() << ", " << "acos = " << std::acos(cos_beta.value()) << ", dcosbeta = " << cos_beta.derivatives().transpose() << std::endl;
		}
	}
	if (std::abs(H2.derivatives().maxCoeff()) > 200) {
		//for (int i = 0; i < ring.size(); i++) {
		//	std::cout << ring[i].transpose() << std::endl;
		//}
		//std::cout << o.transpose() << std::endl << std::endl;
		//std::cout << H2.derivatives().transpose() << std::endl;
	}
	for (int i = 0; i < ring.size(); i++) {
		ring[i] = H2.derivatives().block<3, 1>(i * 3, 0);
		//ring[i] = H2.derivatives().block<3, 1>(((i + 1) % ring.size()) * 3, 0);
	}
	o = H2.derivatives().bottomRows(3);
	return { H2.value(),1 };
#endif
}

std::tuple<msf::Real, msf::Real> msf::PeriodSurfaceMesh::meanH2_nograd(Eigen::Vector3<Real>& o, std::vector<Eigen::Vector3<Real>>& ring, bool area_averaged)
{
#ifndef USE_ANGLE_H2
	std::vector<Real> oh2, eh2;
	Real A = 0;
	for (int i = 0; i < ring.size(); i++) {
		oh2.push_back((ring[i] - o).squaredNorm());
		eh2.push_back((ring[(i + 1) % ring.size()] - ring[i]).squaredNorm());
		A += (o - ring[i]).cross(ring[(i + 1) % ring.size()] - ring[i]).norm() / 6;
	}
	Eigen::Vector3<Real> Hv;
	Hv.setZero();
	for (int i = 0; i < ring.size(); i++) {
		Real cos_a = (oh2[i] + eh2[i] - oh2[(i + 1) % ring.size()]) / (2 * std::sqrt(oh2[i] * eh2[i]));
		Real cos_c = (oh2[(i + 2) % ring.size()] + eh2[(i + 1) % ring.size()] - oh2[(i + 1) % ring.size()]) / (2 * std::sqrt(oh2[(i + 2) % ring.size()] * eh2[(i + 1) % ring.size()]));
		if (std::isnan(cos_a) || std::isnan(cos_c)) {
			std::cout << i << "-th cos is Nan" << std::endl;
		}
		Real cot_a = cos_a / std::sqrt(1 - cos_a * cos_a);
		Real cot_c = cos_c / std::sqrt(1 - cos_c * cos_c);
		if (std::isnan(cot_a) || std::isnan(cot_c)) {
			std::cout << i << "-th cot is Nan, " << "cos = " << cos_a << ", " << cos_c << std::endl;
		}
		Hv += (cot_a + cot_c) / 2 * (ring[(i + 1) % ring.size()] - o); // Note that here should be i+1 !!!
	}
	Real H2 = Hv.squaredNorm();
	//A = A * A;
	//if (area_averaged) { H2 /= A; }

	return { H2, A };
#else
	std::vector<Eigen::Vector3<Real>> v_oh(ring.size()), v_eh(ring.size());
	for (int i = 0; i < ring.size(); i++) {
		v_oh[i] = (ring[i] - o).normalized();
		v_eh[i] = (ring[(i + 1) % ring.size()] - ring[i]).normalized();
	}
	Real H2 = 0;
	int N = ring.size();
	for (int i = 0; i < ring.size(); i++) {
		const auto& a = v_eh[i];
		const auto& b = v_eh[(i + 1) % N];
		const auto& c = -v_oh[(i + 2) % N];
		const auto& d = v_oh[i];
		Real cos_beta = a.dot(c) * b.dot(d) - a.dot(b) * (c).dot(d) - b.dot(c) * d.dot(a);
		if (1 - std::pow(cos_beta, 2) < 1e-6 * 1e-6) {
			H2 += std::acos(std::clamp(cos_beta, -1 + 1e-6, 1 - 1e-6));
		} else {
			H2 += std::acos(cos_beta);
		}
	}
	return { H2, 1 };
#endif
}

std::tuple<Eigen::Vector3<msf::Real>, std::vector<Eigen::Vector3<msf::Real>>> msf::PeriodSurfaceMesh::find1ring(const Eigen::VectorX<Real>& vlist, OM::SmartVertexHandle vh, const std::vector<HETag>& het)
{
	std::vector<Eigen::Vector3<Real>> ring;
	for (auto voh : voh_ccw_range(vh)) {
		Eigen::Vector3<Real> p = vlist.block<3, 1>(voh.to().idx() * 3, 0);
		for (int k = 0; k < 3; k++) {
			if (het[voh.idx()].pass_boundary[k]) {
				//std::cout << p.transpose() << std::endl;
				p[k] -= het[voh.idx()].torus_trans[k];
				//std::cout << p.transpose() << std::endl;
			}
		}
		ring.push_back(Eigen::Vector3<Real>(p[0], p[1], p[2]));
		//std::cout << ring.rbegin()->transpose() << std::endl;
	}
	Eigen::Vector3<Real> op = vlist.block<3, 1>(vh.idx() * 3, 0);
	return { op,ring };
}

std::tuple<Eigen::Vector3<msf::Real>, std::vector<Eigen::Vector3<msf::Real>>> msf::PeriodSurfaceMesh::find1ring(const Eigen::MatrixX3<Real>& vlist, OM::SmartVertexHandle vh, const std::vector<HETag>& het)
{
	std::vector<Eigen::Vector3<Real>> ring;
	for (auto voh : voh_ccw_range(vh)) {
		Eigen::Vector3<Real> p = vlist.row(voh.to().idx()).transpose();
		for (int k = 0; k < 3; k++) {
			if (het[voh.idx()].pass_boundary[k]) {
				p[k] -= het[voh.idx()].torus_trans[k];
			}
		}
		ring.push_back(Eigen::Vector3<Real>(p[0], p[1], p[2]));
	}
	Eigen::Vector3<Real> op = vlist.row(vh.idx()).transpose();
	return { op,ring };
}

std::tuple<std::vector<OpenMesh::SmartVertexHandle>, std::vector<OpenMesh::SmartHalfedgeHandle>> msf::PeriodSurfaceMesh::find1ring(OM::SmartVertexHandle vh)
{
	std::vector<OM::SmartVertexHandle> vvlist;
	std::vector<OM::SmartHalfedgeHandle> vhelist;
	for (auto voh : voh_ccw_range(vh)) {
		vhelist.push_back(voh);
		vvlist.push_back(voh.to());
	}
	return { vvlist,vhelist };
}

std::tuple<Eigen::Vector3<msf::Real>, std::vector<Eigen::Vector3<msf::Real>>> msf::PeriodSurfaceMesh::find1ring(const Eigen::MatrixX3<Real>& vlist, OM::SmartVertexHandle vh)const
{
	std::vector<Eigen::Vector3<Real>> ring;
	for (auto voh : voh_ccw_range(vh)) {
		Eigen::Vector3<Real> p = vlist.row(voh.to().idx()).transpose();
		ring.push_back(Eigen::Vector3<Real>(p[0], p[1], p[2]));
	}
	Eigen::Vector3<Real> op = vlist.row(vh.idx()).transpose();
	return { op,ring };
}

std::vector<msf::Real> msf::PeriodSurfaceMesh::find1ring(OM::SmartVertexHandle vh, const Eigen::SparseMatrix<Real>& L)
{
	std::vector<Real> wlist;
	for (auto voh : voh_ccw_range(vh)) { wlist.push_back(L.coeff(vh.idx(), voh.to().idx())); }
	return wlist;
}

std::tuple<Eigen::Vector3<msf::Real>, std::vector<Eigen::Vector3<msf::Real>>> msf::PeriodSurfaceMesh::find1ring(OM::SmartVertexHandle vh, const std::vector<HETag>& het)
{
	std::vector<Eigen::Vector3<Real>> ring;
	for (auto voh : voh_ccw_range(vh)) {
		Eigen::Vector3<Real> p = toEigen(point(voh.to()));
		for (int k = 0; k < 3; k++) {
			if (het[voh.idx()].pass_boundary[k]) {
				p[k] -= het[voh.idx()].torus_trans[k];
			}
		}
		ring.push_back(Eigen::Vector3<Real>(p[0], p[1], p[2]));
	}
	Eigen::Vector3<Real> op = toEigen(point(vh));
	return { op,ring };

}

std::tuple<Eigen::Vector3<msf::Real>, std::vector<Eigen::Vector3<msf::Real>>> msf::PeriodSurfaceMesh::find1ring(OM::SmartVertexHandle vh, Real deduce /*= 0.5*/) const
{
	Eigen::Vector3<Real> op = toEigen(point(vh));
	std::vector<Eigen::Vector3<Real>> ring;
	for (auto voh : voh_ccw_range(vh)) {
		if (status(voh).deleted() || !voh.is_valid()) continue;
		if (status(voh.to()).deleted() || !voh.to().is_valid()) continue;
		Eigen::Vector3<Real> p = toEigen(point(voh.to()));
		p = op + make_period<Eigen::Vector3<Real>>(p - op, 2, deduce);
		ring.push_back(Eigen::Vector3<Real>(p[0], p[1], p[2]));
	}
	return { op,ring };
}

std::tuple<Eigen::Vector3<msf::Real>, Eigen::Vector3<msf::Real>, Eigen::Vector3<msf::Real>, Eigen::Vector3<msf::Real>> 
msf::PeriodSurfaceMesh::find1diag(OM::SmartHalfedgeHandle he, const std::vector<HETag>& hetags)
{
	auto d0 = toEigen(point(he.from()));
	auto d1 = toEigen(point(he.to()));
	d1 -= hetags[he.idx()].torus_trans;
	auto lef = toEigen(point(he.prev().from()));
	auto rig = toEigen(point(he.opp().next().to()));
	lef -= hetags[he.prev().opp().idx()].torus_trans;
	rig -= hetags[he.opp().next().idx()].torus_trans;
	return { d0,d1,lef,rig };
}

//  https://s2.loli.net/2024/07/11/8JdWhSCwj4eIL5M.png
std::tuple<Eigen::Vector3<msf::Real>, Eigen::Vector3<msf::Real>, Eigen::Vector3<msf::Real>, Eigen::Vector3<msf::Real>> msf::PeriodSurfaceMesh::find1diag(
	OM::SmartHalfedgeHandle he, Real deduce)
{
	auto d0 = toEigen(point(he.from()));
	auto d1 = toEigen(point(he.to()));
	d1 = d0 + make_period<Eigen::Vector3<Real>>(d1 - d0, 2, deduce);
	auto lef = toEigen(point(he.prev().from()));
	auto rig = toEigen(point(he.opp().next().to()));
	lef = d0 + make_period<Eigen::Vector3<Real>>(lef - d0, 2, deduce);
	rig = d0 + make_period<Eigen::Vector3<Real>>(rig - d0, 2, deduce);
	return { d0,d1,lef,rig };
}

std::tuple<msf::Real, msf::Real> msf::PeriodSurfaceMesh::meanH2_test(Eigen::Vector3<Real>& o, std::vector<Eigen::Vector3<Real>>& ring, bool area_averaged)
{
	auto o0 = o;
	auto ring0 = ring;
	auto [h0, a0] = meanH2(o0, ring0, area_averaged);
	std::cout << "Analytical gradient\n";
	for (int i = 0; i < ring.size(); i++) {
		std::cout << ring0[i].transpose() << std::endl;
	}
	std::cout << o0.transpose() << std::endl;
	
	Real step = 1e-6;
	auto ring_grad = ring;
	std::cout << "Numerical gradient\n";
	for (int i = 0; i < ring.size(); i++) {
		auto ring1 = ring;
		for (int j = 0; j < 3; j++) {
			auto o1 = o;
			ring1[i] = ring[i];
			ring1[i][j] += step;
			auto [h1, a1] = meanH2_nograd(o1, ring1, area_averaged);
			ring_grad[i][j] = (h1 - h0) / step;
		}
		std::cout << ring_grad[i].transpose() << std::endl;
	}
	auto o_grad = o;
	for (int j = 0; j < 3; j++) {
		auto o1 = o;
		o1[j] += step;
		auto [h1, a1] = meanH2_nograd(o1, ring, area_averaged);
		o_grad[j] = (h1 - h0) / step;
	}
	std::cout << o_grad.transpose() << std::endl;
	return { h0,a0 };
}

void msf::PeriodSurfaceMesh::relaxInTorusH2(const std::vector<std::vector<size_t>>& cutlocus, const std::vector<int>& handle_tag)
{
	auto [V, F] = getVFlist();
	std::vector<BBox::Vector> pdiag;
	pdiag.push_back(Eigen::Vector3<Real>(0.1, 0.1, 0.1));
	pdiag.push_back(Eigen::Vector3<Real>(0.9, 0.9, 0.9));
	BBox bb(pdiag);
	///  - - - - -
	std::vector<Eigen::Vector3<Real>> pt(n_vertices());
	std::vector<HETag> het(n_halfedges());
	///  - - - - -
	scaleTo(bb, pt);
	{
		std::ofstream ofs(getPath("ptinit"));
		for (int k = 0; k < pt.size(); k++) {
			ofs << pt[k] << std::endl;
		}
		ofs.close();
	}
	using VH = OM::VertexHandle;
#if 0
	for (int n = 0; n < cutlocus.size(); n++) {
		int handleTag = handle_tag[n];
		if (handleTag <= 0 || handleTag > 3) continue;
		for (int i = 0; i < cutlocus[n].size(); i++) {
			auto voh = findHe(VH(cutlocus[n][i]), VH(cutlocus[n][(i + 1) % cutlocus[n].size()]));
			if (!voh.is_valid()) { printf("\033[31mInvalid loop in cut locus"); break; }
			auto th = voh.next();
			het[th.idx()].pass_boundary[handleTag - 1] = true;
			het[th.idx()].torus_trans[handleTag - 1] = 1;
			het[th.opp().idx()].pass_boundary[handleTag - 1] = true;
			het[th.opp().idx()].torus_trans[handleTag - 1] = -1;
			if (pt[th.to().idx()][handleTag - 1] < 1) { pt[th.to().idx()][handleTag - 1] += 1; }

			th = voh.prev().opp();
			het[th.idx()].pass_boundary[handleTag - 1] = true;
			het[th.idx()].torus_trans[handleTag - 1] = 1;
			het[th.opp().idx()].pass_boundary[handleTag - 1] = true;
			het[th.opp().idx()].torus_trans[handleTag - 1] = -1;
			if (pt[th.to().idx()][handleTag - 1] < 1) { pt[th.to().idx()][handleTag - 1] += 1; }
		}
	}
#endif
	{
		std::ofstream ofs(getPath("cuthe"));
		for (auto he : halfedges()) {
			Eigen::Vector3<Real> p1 = pt[he.from().idx()];
			Eigen::Vector3<Real> p2 = pt[he.to().idx()];
			bool has_trans = false;
			for (int k = 0; k < 3; k++) {
				if (p1[k] > 1) { has_trans = true; p1[k] -= 1; } 
				if (p2[k] > 1) { has_trans = true; p2[k] -= 1; };
			}
			if (has_trans) {
				ofs << p1.transpose() << "; " << p2.transpose() << std::endl;
			}
		}
	}
	/// optimize Willmore energy
	std::cout << "\033[0m Begin iteration" << std::endl;
	Eigen::VectorX<Real> vlist(n_vertices() * 3);
	for (int i = 0; i < pt.size(); i++) { vlist.block<3, 1>(i * 3, 0) = pt[i]; }
	int max_iter = 1000;
	for (int iter = 0; iter < max_iter; iter++) {
		// compute gradient
		bool area_averg = true;
		Eigen::VectorX<Real> v_grad(n_vertices() * 3), v_newton(n_vertices() * 3);
		v_grad.setZero();
		v_newton.setZero();
		Real H2_total = 0;
		for (int vid = 0; vid < n_vertices(); vid++) {
			// find 1 ring
			std::vector<Eigen::Vector3<Real>> ring, old_ring;
			OM::SmartVertexHandle vh(vid, this);
			for (auto voh : voh_ccw_range(vh)) {
				Eigen::Vector3<Real> p = vlist.block<3, 1>(voh.to().idx() * 3, 0);
				for (int k = 0; k < 3; k++) {
					if (het[voh.idx()].pass_boundary[k]) {
						//std::cout << p.transpose() << std::endl;
						p[k] -= het[voh.idx()].torus_trans[k];
						//std::cout << p.transpose() << std::endl;
					}
				}
				ring.push_back(Eigen::Vector3<Real>(p[0], p[1], p[2]));
				//std::cout << ring.rbegin()->transpose() << std::endl;
			}
			Eigen::Vector3<Real> op = vlist.block<3, 1>(vh.idx() * 3, 0);
			Eigen::Vector3<Real> old_op;
			//std::cout << op.transpose() << std::endl;
			// compute gradient
			old_ring = ring;
			old_op = op;
			//auto [hi2, Ai] = meanH2_test(op, ring, area_averg);
			auto [hi2, Ai] = meanH2(op, ring, area_averg);
			bool hasNan = false;
			for (int i = 0; i < ring.size(); i++) {
				if (ring[i].hasNaN()) {
					if (!hasNan) {
						std::cout << "\033[31m" << "Ring gradient has NaN" << "\033[0m\n";
						hasNan = true; break;
					}
				}
			}
			if (op.hasNaN()) { hasNan = true; std::cout << "\033[31m" << "center gradient has NaN" << "\033[0m\n"; }
			if (hasNan) {
				std::cout << "ring : \n";
				for (int i = 0; i < ring.size(); i++) {
					std::cout << old_ring[i].transpose() << "; " << ring[i].transpose() << std::endl;
				}
				std::cout << "center \n" << old_op.transpose() << "; " << op.transpose() << std::endl;
			}
			Real J2 = op.squaredNorm(); if (J2 < 1e-12) J2 = 1e-12;
			for (int i = 0; i < ring.size(); i++) { J2 += ring[i].squaredNorm(); }
			// only counter gradient on center
#if 1
			{
				int counter = 0;
				for (auto voh : voh_ccw_range(vh)) {
					Eigen::Vector3<Real> g = ring[counter] / J2 * hi2 * 2;
					if (g.hasNaN()) {
						std::cout << "NaN,  J2 = " << J2 << ", hi2 = " << hi2 << ", v = " << ring[counter].transpose() << std::endl;
					}
					v_grad.block<3, 1>(voh.to().idx() * 3, 0) += area_averg ? ring[counter] : (ring[counter] * Ai);
					v_newton.block<3, 1>(voh.to().idx() * 3, 0) += area_averg ? g : (g * Ai);
					counter++;
				}
			}
#endif
			Eigen::Vector3<Real> g = op / J2 * hi2 * 2;
			v_grad.block<3, 1>(vh.idx() * 3, 0) += area_averg ? op : (op * Ai);
			v_newton.block<3, 1>(vh.idx() * 3, 0) += area_averg ? g : (g * Ai);
			H2_total += area_averg ? hi2 : (hi2 * Ai);
		}
		printf("Hsum = %.2e, g'd = %.2e, W = %f\n", H2_total, v_newton.dot(v_grad), H2_total - 2 * M_PI * n_vertices());
		if (v_newton.hasNaN()) {
			std::cout << "\033[31m" << "Energy gradient has NaN!" << "\033[0m" << std::endl;
		}
		if (iter % 20 == 0) {
			Eigen::MatrixX3d vigl(n_vertices(), 3);
			for (int i = 0; i < n_vertices(); i++) {
				vigl.row(i) = vlist.block<3, 1>(i * 3, 0).transpose();
			}
			for (int i = 0; i < vigl.size(); i++) {
				if (vigl.data()[i] > 1) { vigl.data()[i] -= 1; }
			}
			igl::writeOBJ(getPath("meshiter.obj"), vigl, F);
			std::ofstream ofs(getPath("Hgrad"));
			ofs << v_grad;
			eigen2ConnectedMatlab("Hgrad", v_grad);
		}
		// update vertex position
		//eigen2ConnectedMatlab("vlist", vlist);
		//eigen2ConnectedMatlab("vgrad", v_grad);
		//eigen2ConnectedMatlab("vnewt", v_newton);
		v_newton = v_grad;
		// line search
		Real step = 1;
		Real H2_new = 0;
		while (step > 1e-20) {
			Eigen::VectorX<Real> vlist_new = vlist - v_newton * step;
			Real d = v_grad.dot(v_newton);
			Real H2 = 0;
			for (auto vh : vertices()) {
				auto [o, ring] = find1ring(vlist_new, vh, het);
				auto [hi2, Ai] = meanH2_nograd(o, ring, true);
				H2 += area_averg ? hi2 : hi2 * Ai;
			}
			if (H2 < H2_total - 0.6 * step * d) { H2_new = H2; break; }
			//std::cout << "Searching step " << step << ", H2 = " << H2 << ", k = " << (H2 - H2_total) / step << ", d = " << d << std::endl;
			step *= 0.5;
		}
#if 1
		{
			std::cout << "g = " << v_newton.dot(v_grad) << std::endl;
			std::vector<Real> steplist, vallist;
			for (Real s = - step * 100; s < step * 100; s += step / 10) {
				Eigen::VectorX<Real> vlist_new = vlist - v_newton * s;
				Real H2 = 0;
				for (auto vh : vertices()) {
					auto [o, ring] = find1ring(vlist_new, vh, het);
					auto [hi2, Ai] = meanH2_nograd(o, ring, true);
					H2 += area_averg ? hi2 : hi2 * Ai;
				}
				steplist.push_back(s);
				vallist.push_back(H2);
			}
			eigen2ConnectedMatlab("steplist", Eigen::VectorX<Real>::Map(steplist.data(), steplist.size()));
			eigen2ConnectedMatlab("vallist", Eigen::VectorX<Real>::Map(vallist.data(), vallist.size()));
		}
#endif
		std::cout << "Searching step " << step << ", H2 = " << H2_new << std::endl;
		vlist -= v_newton * step;
	}

}

void msf::PeriodSurfaceMesh::relaxInTorusH1(const std::vector<std::vector<size_t>>& cutlocus, const std::vector<int>& handle_tag)
{
}

std::tuple<msf::Real, msf::Real> msf::PeriodSurfaceMesh::meanHs(Eigen::Vector3<Real>& o, std::vector<Eigen::Vector3<Real>>& ring, Real s, bool area_averaged)
{

}

void msf::PeriodSurfaceMesh::meanH(const Eigen::Vector3<Real>& o, const std::vector<Eigen::Vector3<Real>>& ring, Eigen::Vector3<Real>& Ho)
{
	auto [Hv, A] = meanH(o, ring);
	Ho = Hv;
}

std::tuple<Eigen::Vector3d, msf::Real> msf::PeriodSurfaceMesh::meanH(const Eigen::Vector3<Real>& o, const std::vector<Eigen::Vector3<Real>>& ring)
{
	std::vector<Real> oh2, eh2;
	Real A = 0;
	for (int i = 0; i < ring.size(); i++) {
		oh2.push_back((ring[i] - o).squaredNorm());
		eh2.push_back((ring[(i + 1) % ring.size()] - ring[i]).squaredNorm());
		A += (o - ring[i]).cross(ring[(i + 1) % ring.size()] - ring[i]).norm() / 6;
	}
	Eigen::Vector3<Real> Hv;
	Hv.setZero();
	for (int i = 0; i < ring.size(); i++) {
		Real cos_a = (oh2[i] + eh2[i] - oh2[(i + 1) % ring.size()]) / (2 * std::sqrt(oh2[i] * eh2[i]));
		Real cos_c = (oh2[(i + 2) % ring.size()] + eh2[(i + 1) % ring.size()] - oh2[(i + 1) % ring.size()]) / (2 * std::sqrt(oh2[(i + 2) % ring.size()] * eh2[(i + 1) % ring.size()]));
		if (std::isnan(cos_a) || std::isnan(cos_c)) {
			std::cout << i << "-th cos is Nan" << std::endl;
			std::cout << "o = " << o.transpose() << std::endl;
			for (int k = 0; k < ring.size(); k++) {
				std::cout << fmt::format("ring[{}] = ", k) << ring[k].transpose() << std::endl;
			}
			throw std::runtime_error("numerical issuse computing cos");
		}
		Real cot_a = cos_a / std::sqrt(1 - cos_a * cos_a);
		Real cot_c = cos_c / std::sqrt(1 - cos_c * cos_c);
		if (std::isnan(cot_a) || std::isnan(cot_c)) {
			std::cout << i << "-th cot is Nan, " << "cos = " << cos_a << ", " << cos_c << std::endl;
			continue;
			//throw std::runtime_error("numerical issuse computing cos");
		}
		Hv += (cot_a + cot_c) / 2 * (ring[(i + 1) % ring.size()] - o);
	}
	return std::make_tuple(Hv, A);
}

void msf::PeriodSurfaceMesh::relaxInTorusLaplacian(const std::vector<std::vector<size_t>>& cutlocus, const std::vector<int>& handle_tag)
{
	auto [V, F] = getVFlist();
	std::vector<BBox::Vector> pdiag;
	pdiag.push_back(Eigen::Vector3<Real>(0.1, 0.1, 0.1));
	pdiag.push_back(Eigen::Vector3<Real>(0.9, 0.9, 0.9));
	BBox bb(pdiag);
	///  - - - - -
	std::vector<Eigen::Vector3<Real>> pt(n_vertices());
	std::vector<HETag> het(n_halfedges());
	///  - - - - -
	scaleTo(pdiag, pt);

	{
		Eigen::MatrixX<Real> ptmat(pt.size(), 3);
		for (int i = 0; i < pt.size(); i++) { ptmat.row(i) = pt[i].transpose(); }
		igl::writeOBJ(getPath("rescale.obj"), ptmat, F);
	}
#if 0
	{
		std::ofstream ofs(getPath("ptinit"));
		for (int k = 0; k < pt.size(); k++) { ofs << pt[k] << std::endl; }
		ofs.close();
	}
#endif
	using VH = OM::SmartVertexHandle;
	std::vector<VH> vTlist;
	for (int i = 0; i < cutlocus.size(); i++) {
		for (int j = 0; j < cutlocus[i].size(); j++) {
			vTlist.push_back(VH(cutlocus[i][j], this));
		}
	}
	std::set<OM::SmartHalfedgeHandle> hetranslist;
#if 1
	for (int n = 0; n < cutlocus.size(); n++) {
		int handleTag = handle_tag[n];
		if (handleTag <= 0 || handleTag > 3) continue;
		for (int i = 0; i < cutlocus[n].size(); i++) {
			int NC = cutlocus[n].size();
			auto voh = findHe(VH(cutlocus[n][i]), VH(cutlocus[n][(i + 1) % NC], this));
			auto vohend = findHe(VH(cutlocus[n][i]), VH(cutlocus[n][(i + NC - 1) % NC], this));
			if (!voh.is_valid()|| !vohend.is_valid()) { printf("\033[31mInvalid loop in cut locus"); break; }
			voh = voh.prev().opp();
			while (voh != vohend) {
				hetranslist.insert(voh);
				het[voh.idx()].pass_boundary[handleTag - 1] = true;
				het[voh.idx()].torus_trans[handleTag - 1] = 1;
				het[voh.opp().idx()].pass_boundary[handleTag - 1] = true;
				het[voh.opp().idx()].torus_trans[handleTag - 1] = -1;
				if (pt[voh.to().idx()][handleTag - 1] < 1) { pt[voh.to().idx()][handleTag - 1] += 1; }
				voh = voh.prev().opp();
			}
		}
	}
#endif
#if 0
	{
		std::ofstream ofs(getPath("cuthe"));
		for (auto he : halfedges()) {
			Eigen::Vector3<Real> p1 = pt[he.from().idx()];
			Eigen::Vector3<Real> p2 = pt[he.to().idx()];
			bool has_trans = false;
			for (int k = 0; k < 3; k++) {
				if (p1[k] > 1) { has_trans = true; p1[k] -= 1; } 
				if (p2[k] > 1) { has_trans = true; p2[k] -= 1; };
			}
			if (has_trans) {
				ofs << p1.transpose() << "; " << p2.transpose() << std::endl;
			}
		}
	}
#endif
	/// optimize Willmore energy
	std::cout << "\033[0m Begin iteration" << std::endl;
	Eigen::MatrixX<Real> vlist(n_vertices(), 3);
	vlist.setZero();
	for (int i = 0; i < pt.size(); i++) { vlist.row(i) = pt[i].transpose(); }

	Eigen::SparseMatrix<Real> Linput;
	igl::cotmatrix(V, F, Linput);
	eigen2ConnectedMatlab("Li", Linput);
	
	Eigen::MatrixX<Real> x_sol(n_vertices(), 3);
	for (int iter = 0; iter < 200; iter++) {
		Eigen::SparseMatrix<Real> L;
		std::vector<OM::SmartHalfedgeHandle> transhe(hetranslist.begin(), hetranslist.end());
		cotweight(vlist, L, transhe, het);
		eigen2ConnectedMatlab("L", L);
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<Real>> Lsolver(L);
		for (int axis = 0; axis < 3; axis++) {
			Eigen::SparseMatrix<Real> Lp = L;
			Eigen::VectorX<Real> bp(Lp.rows(), 1);
			bp.setZero();
			Eigen::SparseMatrix<Real> A(hetranslist.size(), Lp.rows());
			for (auto he : hetranslist) {
				auto ht = het[he.idx()];
				if (ht.pass_boundary[axis]) {
					bp[he.from().idx()] -= Lp.coeffRef(he.from().idx(), he.to().idx()) * ht.torus_trans[axis];
					bp[he.to().idx()] += Lp.coeffRef(he.to().idx(), he.from().idx()) * ht.torus_trans[axis];
				}
			}
			x_sol.col(axis) = Lsolver.solve(bp);
		}
		Real blending = 3e-2;
		vlist = x_sol * blending + vlist * (1 - blending);
		{ std::ofstream ofs(getPath("x_sol")); ofs << x_sol; }
		//igl::writeOBJ(getPath("meshiter.obj"), vlist, F);
		savePeriodicMesh(getPath("meshiter.obj"), vlist, vTlist, het);
	}

}

void msf::PeriodSurfaceMesh::cotweight(const Eigen::MatrixX3<Real>& V, Eigen::SparseMatrix<Real>& L, const std::vector<OM::SmartHalfedgeHandle>& hetrans, const std::vector<HETag>& hetags)
{
	auto F = getFlist();
	L.setZero();
	igl::cotmatrix(V, F, L);
	for (int i = 0; i < hetrans.size(); i++) {
		auto v0h = hetrans[i].from();
		auto v1h = hetrans[i].to();
		//std::cout << " -  -  -  -  -  -  -  -  -  -  -  - " << std::endl;
		auto [v0, ring0] = find1ring(V, v0h, hetags);
		auto w0 = cot1ring(v0, ring0);
		auto [vhring, hering] = find1ring(v0h);
		//std::cout << "[" << v0h.idx() << "] " << v0.transpose() << std::endl;
		Real w01 = 0;
		for (int k = 0; k < vhring.size(); k++) {
			L.coeffRef(v0h.idx(), vhring[k].idx()) = w0[k]; 
			//printf("(%d, %d) : %f\n", v0h.idx(), vhring[k].idx(), w0[k]);
			if (vhring[k] == v1h) {
				//std::cout << "[" << vhring[k].idx() << "] " << ring0[k].transpose() << " " << w0[k] << std::endl;
				w01 = w0[k]; 
			}
		}
		L.coeffRef(v0h.idx(), v0h.idx()) = w0[vhring.size()];

		//std::cout << " -  -  -  -  -  -  -  -  -  -  -  - " << std::endl;
		auto [v1, ring1] = find1ring(V, v1h, hetags);
		auto w1 = cot1ring(v1, ring1);
		std::tie(vhring, hering) = find1ring(v1h);
		//std::cout << "[" << v1h.idx() << "] " << v1.transpose() << std::endl;
		Real w10 = 0;
		for (int k = 0; k < vhring.size(); k++) {
			L.coeffRef(v1h.idx(), vhring[k].idx()) = w1[k]; 
			//printf("(%d, %d) : %f\n", v1h.idx(), vhring[k].idx(), w1[k]);
			if (vhring[k] == v0h) {
				//std::cout << "[" << vhring[k].idx() << "] " << ring1[k].transpose() << " " << w1[k] << std::endl;
				w10 = w1[k]; 
			}
		}
		if (std::abs(w01 - w10) > 1e-1) {
			std::cout << "\033[31m" << "Unsymmetric pair" << "\033[0m" << std::endl;
			std::cout << "ring0" << std::endl;
			std::cout << v0.transpose() << std::endl;
			for (int k = 0; k < ring0.size(); k++) {
				std::cout << ring0[k].transpose() << std::endl;
			}
			std::cout << "ring1" << std::endl;
			std::cout << v1.transpose() << std::endl;
			for (int k = 0; k < ring1.size(); k++) {
				std::cout << ring1[k].transpose() << std::endl;
			}
		}
		L.coeffRef(v1h.idx(), v1h.idx()) = w1[vhring.size()];
	}
}

Eigen::Matrix<int, -1, 3> msf::PeriodSurfaceMesh::getFlist(void)
{
	Eigen::Matrix<int, -1, 3> flist(n_faces(), 3);
	int counter = 0;
	for (auto f : MeshSurface::faces()) {
		if (!f.is_valid()) continue;
		auto fvlist = getFaceVertexHandle(f);
		if (fvlist[0].is_valid() && fvlist[1].is_valid() && fvlist[2].is_valid()) {
			for (int k = 0; k < 3; k++) flist(counter, k) = fvlist[k].idx();
			counter++;
		}
	}
	flist.conservativeResize(counter, 3);
	return flist;
}

std::vector<msf::Real> msf::PeriodSurfaceMesh::cot1ring(const Eigen::Vector3<Real>& v0, const std::vector<Eigen::Vector3<Real>>& ring)
{
	std::vector<Real> w_cot(ring.size());
	int N = ring.size();
	for (int i = 0; i < ring.size(); i++) {
		const auto& va = ring[i] - v0;
		const auto& vb = ring[(i + 1) % N] - v0;
		const auto& vc = ring[(i + 2) % N] - v0;
		const auto& ab = ring[(i + 1) % N] - ring[i];
		const auto& cb = ring[(i + 1) % N] - ring[(i + 2) % N];
		Real cosa = (va.dot(va) + ab.dot(ab) - vb.dot(vb)) / (2 * va.norm() * ab.norm());
		Real cosc = (vc.dot(vc) + cb.dot(cb) - vb.dot(vb)) / (2 * vc.norm() * cb.norm());
		Real cota = cosa / std::sqrt(1 - cosa * cosa);
		Real cotc = cosc / std::sqrt(1 - cosc * cosc);
		w_cot[(i + 1) % N] = (cota + cotc) / 2;
	}
	w_cot.push_back(-std::accumulate(w_cot.begin(), w_cot.end(), 0));
	return w_cot;
}

void msf::PeriodSurfaceMesh::savePeriodicMesh(std::string filename, const Eigen::MatrixX<Real>& V, const std::vector<OM::SmartVertexHandle>& vcut, const std::vector<HETag>& hetags)
{
	MeshSurface mesh = *this;
	for (auto vh : mesh.vertices()) {
		OM::Vec3d op(V(vh.idx(), 0), V(vh.idx(), 1), V(vh.idx(), 2));
		mesh.point(vh) = op;
	}

	std::map<OM::VertexHandle, Eigen::Vector3<Real>> v2p;
	std::set<OM::SmartFaceHandle> f2v;
	for (int i = 0; i < vcut.size(); i++) {
		v2p[vcut[i]] = Eigen::Vector3<Real>();
	}
	for (int i = 0; i < hetags.size(); i++) {
		OM::SmartHalfedgeHandle he(i, &mesh);
		auto het = hetags[i];
		if (het.pass_boundary[0] + het.pass_boundary[1] + het.pass_boundary[2] == 0) {
			continue;
		}
		if (v2p.count(he.from())) {
			auto v0p = toEigen(mesh.point(he.from()));
			for (int k = 0; k < 3; k++) { v0p[k] += het.torus_trans[k]; }
			v2p[he.from()] = v0p;
			f2v.insert(he.face());
		}
		else if (v2p.count(he.to())) {
			auto v1p = toEigen(mesh.point(he.to()));
			for (int k = 0; k < 3; k++) { v1p[k] -= het.torus_trans[k]; }
			v2p[he.to()] = v1p;
			f2v.insert(he.face());
		} else {
			std::cout << "\033[31m" << "Unexpected translation edge" << "\033[0m" << std::endl;
		}
	}
	std::map<OM::SmartVertexHandle, OM::SmartVertexHandle> vold2new;
	for (auto fh : f2v) {
		std::vector<OM::VertexHandle> vlist;
		for (auto vh : fh.vertices()) {
			OM::VertexHandle newh = vh;
			if (v2p.count(vh)) {
				if (!vold2new.count(vh)) {
					OM::Vec3d p(v2p[vh][0], v2p[vh][1], v2p[vh][2]);
					vold2new[vh] = mesh.add_vertex(p);
				}
				newh = vold2new[vh];
			}
			vlist.push_back(newh);
		}
		mesh.delete_face(fh, true);
		mesh.add_face(vlist);
	}
	mesh.garbage_collection();
	OM::IO::write_mesh(mesh, filename);
}

void msf::PeriodSurfaceMesh::savePeriodicMesh(std::string filename, const Eigen::MatrixX<Real>& V, const std::set<OM::SmartVertexHandle>& vcut, const std::vector<HETag>& hetags)
{
	std::vector<OM::SmartVertexHandle> vcuth(vcut.begin(), vcut.end());
	savePeriodicMesh(filename, V, vcuth, hetags);
}

// [newv, newf, new2old]
std::tuple<Eigen::MatrixX3<msf::Real>, Eigen::MatrixX3i, Eigen::VectorXi> msf::PeriodSurfaceMesh::dupPeriodFaces(const Eigen::MatrixX3<Real>& V, const Eigen::MatrixX3i& F, Real detach /*= 0.5*/)
{
	auto indexer = generateIndexer(V);
	auto vlist = V;
	auto flist = F;
	std::unordered_map<Eigen::Vector3<Real>, int, Vector3Hash<Real>> poldset;
	{
		int counter = 0;
		for (int i = 0; i < vlist.rows(); i++) {
			Eigen::Vector3<Real> p = vlist.row(i).transpose();
			if(poldset.count(p)) continue ;
			poldset[p] = counter++;
		}
	}
	std::vector<Eigen::Vector3<Real>> vappend;
	for (int i = 0; i < flist.rows(); i++) {
		Eigen::Vector3i fv = flist.row(i).transpose();

		Eigen::Vector3<Real> e[3];
		Eigen::Vector3<Real> vnew[3];
		bool has_period_face = false;
		for (int j = 0; j < 3; j++) {
			vnew[j] = vlist.row(fv[j]); 
		}
		for (int j = 0; j < 3; j++) {
			e[j] = (vnew[(j + 1) % 3] - vnew[j]).transpose();
			for (int k = 0; k < 3; k++) {
				if (e[j][k] < -detach) {
					vnew[(j + 1) % 3][k] += 2;
					has_period_face = true;
				}
				else if (e[j][k] > detach) {
					vnew[j][k] += 2;
					has_period_face = true;
				}
			}
		}
		Eigen::Vector3d cmax = vnew[0].cwiseMax(vnew[1]).cwiseMax(vnew[2]);
		Eigen::Vector3d cmin = vnew[0].cwiseMin(vnew[1]).cwiseMin(vnew[2]);
		for (int k = 0; k < 3; k++) {
			if (cmax[k] > 1 + 1e-5) {
				vnew[0][k] -= 2; vnew[1][k] -= 2; vnew[2][k] -= 2;
			}
			else if (cmin[k] < -1 - 1e-5) {
				vnew[0][k] += 2; vnew[1][k] += 2; vnew[2][k] += 2;
			}
		}
		{
			//Eigen::Vector3d ccmax = vnew[0].cwiseMax(vnew[1]).cwiseMax(vnew[2]);
			//Eigen::Vector3d ccmin = vnew[0].cwiseMin(vnew[1]).cwiseMin(vnew[2]);
			//for (int k = 0; k < 3; k++) {
			//	if (ccmax[k] > 1 + 1e-5) {
			//		std::cout << "Out range " << vnew[0].transpose() << ", " << vnew[1].transpose() << ", " << vnew[2].transpose() << std::endl; 
			//		for (int j = 0; j < 3; j++) { vnew[j] = vlist.row(fv[j]); }
			//		std::cout << "vnew =  " << vnew[0].transpose() << ", " << vnew[1].transpose() << ", " << vnew[2].transpose() << std::endl;
			//		std::cout << "cmax = " << cmax.transpose() << ", cmin = " << cmin.transpose() << std::endl;
			//	}
			//	if (ccmin[k] < -1 - 1e-5) {
			//		std::cout << "Out range " << vnew[0].transpose() << ", " << vnew[1].transpose() << ", " << vnew[2].transpose() << std::endl; 
			//		for (int j = 0; j < 3; j++) { vnew[j] = vlist.row(fv[j]); }
			//		std::cout << "vnew =  " << vnew[0].transpose() << ", " << vnew[1].transpose() << ", " << vnew[2].transpose() << std::endl;
			//		std::cout << "cmax = " << cmax.transpose() << ", cmin = " << cmin.transpose() << std::endl;
			//	}
			//}
		}
		//. replace old face
		if (has_period_face) {
			for (int j = 0; j < 3; j++) {
				if (!poldset.count(vnew[j])) {
					int voldid = fv[j];
					fv[j] = poldset.size() + vappend.size();
					vappend.push_back(vnew[j]);
				}
				else {
					fv[j] = poldset[vnew[j]];
				}
			}
			flist.row(i) = fv.transpose();
		}
	}
	int n_old = vlist.rows();
	vlist.conservativeResize(n_old + vappend.size(), 3);
	for (int i = 0; i < vappend.size(); i++) {
		vlist.row(n_old + i) = vappend[i].transpose();
	}
	auto [new_v, new_f] = remove_dup_vertices(vlist, flist, 1e-5);

	Eigen::VectorXi new2old(new_v.rows()); new2old.setConstant(-1);
	for (int i = 0; i < new_v.rows(); i++) {
		Eigen::Vector3<Real> p = new_v.row(i).transpose();
		int old_id = indexer.query(p);
		new2old[i] = old_id;
	}
	return { new_v,new_f, new2old };
}

std::tuple<Eigen::MatrixX3<msf::Real>, Eigen::MatrixX3i> msf::PeriodSurfaceMesh::getPeriodicMesh(void)
{
	auto [vlist, flist] = getVFlist();

	auto [newvlist, newflist, new2old] = dupPeriodFaces(vlist, flist, 0.5);
	return{ newvlist, newflist };
}

void msf::PeriodSurfaceMesh::savePeriodicMesh(std::string filename, const std::vector<OM::SmartVertexHandle>& vcut, Real detach /*= 0.5*/)
{
	auto [vlist, flist] = getVFlist();
	Eigen::VectorXi new2old;
	std::tie(vlist, flist, new2old) = dupPeriodFaces(vlist, flist, detach);
	//igl::writeOBJ(filename, vlist, flist);
	igl::write_triangle_mesh(filename, vlist, flist);
}

void msf::PeriodSurfaceMesh::savePeriodicMesh(std::string filename, const std::set<OM::SmartVertexHandle>& vcut, Real detach /*= 0.5*/)
{
	std::vector<OM::SmartVertexHandle> vcuth(vcut.begin(), vcut.end());
	savePeriodicMesh(filename, vcuth, detach);
}

void msf::PeriodSurfaceMesh::relaxInTorusARAP(const std::vector<std::vector<size_t>>& cutlocus, const std::vector<int>& handle_tag)
{
	auto [V, F] = getVFlist();
	std::vector<BBox::Vector> pdiag{ Eigen::Vector3<Real>(0.1, 0.1, 0.1),Eigen::Vector3<Real>(0.9, 0.9, 0.9) };
	BBox bb(pdiag);
	Eigen::MatrixX3<Real> vlist;
	scaleTo(pdiag, vlist);
	Eigen::MatrixX3<Real> vlist0(vlist);
	auto [hetrans, het] = toTransHETag(cutlocus, handle_tag, vlist);
	{
		std::ofstream ofs(getPath("hetrans"));
		for (int i = 0; i < hetrans.size(); i++) {
			ofs << point(hetrans[i].from()) << " " << point(hetrans[i].to()) << std::endl;
		}
	}
	std::cout << "hetrans size = " << hetrans.size() << std::endl;
	std::vector<OM::SmartVertexHandle> vTlist;
	for (int i = 0; i < cutlocus.size(); i++) {
		for (int j = 0; j < cutlocus[i].size(); j++) {
			vTlist.push_back(OM::SmartVertexHandle(cutlocus[i][j], this));
		}
	}

	auto Flist = getFlist();

	for (int iter = 0; iter < 100; iter++) {
		std::vector<Eigen::Matrix3<Real>> Rv(n_vertices());
		Eigen::SparseMatrix<Real> L;
		//cotweight(vlist0, L, hetrans, het);
		igl::cotmatrix(vlist0, Flist, L);
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<Real>> Lsolver(L);
		// compute R
		for (auto vh : vertices()) {
			auto [o, ring] = find1ring(vlist0, vh);
			auto [o1, ring1] = find1ring(vlist, vh, het);
			auto wij = find1ring(vh, L);
			Rv[vh.idx()] = matchRot(wij, o, ring, o1, ring1);
			if ((Rv[vh.idx()] - Eigen::Matrix3<Real>::Identity()).norm() > 1e-3) {
				//std::cout << Rv[vh.idx()] << std::endl;
				//std::cout << "ring\n";
				//for (int i = 0; i < ring1.size(); i++) {
				//	std::cout << ring1[i].transpose() << std::endl;
				//}
			}
		}
		// ARAP update
		Eigen::MatrixX3<Real> b(n_vertices(), 3);
		b.setZero();
		for (auto vh : vertices()) {
			Eigen::Vector3<Real> bv(0, 0, 0);
			auto [o, ring] = find1ring(vlist0, vh);
			auto [vnei, henei] = find1ring(vh);
			for (int i = 0; i < vnei.size(); i++) {
				auto vv = vnei[i];
				Real w = L.coeff(vh.idx(), vv.idx());
				bv += w / 2 * (Rv[vh.idx()] + Rv[vv.idx()]) * (ring[i] - o);
			}
			b.row(vh.idx()) = bv.transpose();
		}
		// get periodic translation force
		auto bperioid = laplacianTranslationForce(L, hetrans, het);
#if 0
		Eigen::MatrixX3<Real> x = Lsolver.solve(b + bperioid);
#elif 0
		Eigen::MatrixX3<Real> x = Lsolver.solve(b);
		Eigen::MatrixX3<Real> xtrans = Lsolver.solve(bperioid);
		vlist = x + xtrans;
		{std::ofstream ofs(getPath("x_sol")); ofs << xtrans; }
#elif 1
		Eigen::MatrixX3<Real> x = Lsolver.solve(b + bperioid);
		vlist = x;
#endif
		if (iter % 100 == 0) {
			savePeriodicMesh(getPath("meshiter.obj"), vlist, vTlist, het);
		}
	}
	savePeriodicMesh(getPath("meshiter.obj"), vlist, vTlist, het);

	// then use willmore flow
	//willmore_flow(vlist, vTlist, hetrans, het);

	// use mean curvature flow
	meanCurvatureFlow(vlist, vTlist, hetrans, het);
}

std::tuple<std::vector<OpenMesh::SmartHalfedgeHandle>, std::vector<msf::PeriodSurfaceMesh::HETag>> 
msf::PeriodSurfaceMesh::toTransHETag(const std::vector<std::vector<size_t>>& cutlocus, const std::vector<int>& handle_tag, Eigen::MatrixX3<Real>& vtrans)
{
	std::vector<HETag> het(n_halfedges());
	using VH = OM::VertexHandle;
	std::vector<OM::SmartHalfedgeHandle> hetrans;
	for (int n = 0; n < cutlocus.size(); n++) {
		int handleTag = handle_tag[n];
		if (handleTag <= 0 || handleTag > 3) continue;
		for (int i = 0; i < cutlocus[n].size(); i++) {
			int NC = cutlocus[n].size();
			auto voh = findHe(VH(cutlocus[n][i]), VH(cutlocus[n][(i + 1) % NC]));
			auto vohend = findHe(VH(cutlocus[n][i]), VH(cutlocus[n][(i + NC - 1) % NC]));
			if (!voh.is_valid() || !vohend.is_valid()) { printf("\033[31mInvalid loop in cut locus"); break; }
			voh = voh.prev().opp();
			while (voh != vohend) {
				hetrans.push_back(voh);
				het[voh.idx()].pass_boundary[handleTag - 1] = true;
				het[voh.idx()].torus_trans[handleTag - 1] = 1;
				het[voh.opp().idx()].pass_boundary[handleTag - 1] = true;
				het[voh.opp().idx()].torus_trans[handleTag - 1] = -1;
				if (vtrans(voh.to().idx(), handleTag - 1) < 1) { vtrans(voh.to().idx(), handleTag - 1) += 1; }
				voh = voh.prev().opp();
			}
		}
	}
	std::sort(hetrans.begin(), hetrans.end());
	auto uit = std::unique(hetrans.begin(), hetrans.end());
	hetrans.erase(uit, hetrans.end());
	return { hetrans, het };
}

Eigen::Matrix3<msf::Real> msf::PeriodSurfaceMesh::matchRot(
	const Eigen::Vector3<Real>& o, const std::vector<Eigen::Vector3<Real>>& ring,
	const Eigen::Vector3<Real>& o1, const std::vector<Eigen::Vector3<Real>>& ring1
)
{
	auto wij = cot1ring(o, ring);
	matchRot(wij, o, ring, o1, ring1);
}

Eigen::Matrix3<msf::Real> msf::PeriodSurfaceMesh::matchRot(const std::vector<Real>& wij, const Eigen::Vector3<Real>& o, const std::vector<Eigen::Vector3<Real>>& ring, const Eigen::Vector3<Real>& o1, const std::vector<Eigen::Vector3<Real>>& ring1)
{
	Eigen::Matrix3<Real> S;
	S.setZero();
	for (int i = 0; i < ring.size(); i++) {
		S += wij[i] * (ring[i] - o) * (ring1[i] - o1).transpose();
	}
	Eigen::BDCSVD<Eigen::Matrix3<Real>> svd(S, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix3<Real> R = svd.matrixV() * svd.matrixU().transpose();
	return R;
}

Eigen::MatrixX3<msf::Real> msf::PeriodSurfaceMesh::laplacianTranslationForce(
const Eigen::SparseMatrix<Real>& L, const std::vector<OM::SmartHalfedgeHandle>& hetrans, const std::vector<HETag>& hetags)
{
	Eigen::MatrixX3<Real> b(n_vertices(), 3);
	for (int axis = 0; axis < 3; axis++) {
		Eigen::VectorX<Real> bp(L.rows(), 1);
		bp.setZero();
		for (auto he : hetrans) {
			auto ht = hetags[he.idx()];
			if (ht.pass_boundary[axis]) {
				bp[he.from().idx()] += L.coeff(he.from().idx(), he.to().idx()) * ht.torus_trans[axis];
				bp[he.to().idx()] -= L.coeff(he.to().idx(), he.from().idx()) * ht.torus_trans[axis];
			}
		}
		b.col(axis) = bp;
	}
	return b;
}

std::vector<int> msf::PeriodSurfaceMesh::getVertexValence(void)
{
	std::vector<int> valencelist(n_vertices());
	for (auto vh : vertices()) { valencelist[vh.idx()] = valence(vh); }
	return valencelist;
}

std::vector<msf::Real> msf::PeriodSurfaceMesh::getFaceArea(void)
{
	std::vector<Real> facearea(n_faces());
	for (auto fh : faces()) {
		facearea[fh.idx()] = calc_face_area(fh);
	}
	return facearea;
}

std::tuple<std::array<msf::OM::SmartVertexHandle, 3>, std::array<Eigen::Vector3<msf::Real>, 3>> msf::PeriodSurfaceMesh::findFaceVInc(
	const Eigen::MatrixX3<Real>& vlist, OM::SmartFaceHandle fh, const std::vector<HETag>& hetags)
{
	auto he = fh.halfedge();
	OM::SmartHalfedgeHandle vhe[2] = { he, he.prev().opp() };
	std::array<OM::SmartVertexHandle, 3> v = { he.from(), he.next().from(), he.next().next().from() };

	std::array<Eigen::Vector3<Real>, 3> vp;
	for (int i = 0; i < 3; i++) { vp[i] = vlist.row(v[i].idx()).transpose(); }

	if (hetags[vhe[0].idx()].has_trans()) { vp[1] -= hetags[vhe[0].idx()].torus_trans; }

	if (hetags[vhe[1].idx()].has_trans()) { vp[2] -= hetags[vhe[1].idx()].torus_trans; }

	return { v, vp };
}


Eigen::VectorX<msf::Real> msf::PeriodSurfaceMesh::getFaceArea(const Eigen::MatrixX3<Real>& vlist, const std::vector<HETag>& het)
{
	Eigen::VectorX<msf::Real> facearea(n_faces());
	for (auto fh : faces()) {
		auto [vh, vp] = findFaceVInc(vlist, fh, het);
		facearea[fh.idx()] = (vp[1] - vp[0]).cross(vp[2] - vp[1]).norm() / 2;
	}
	for (int i = 0; i < het.size(); i++) {
		if (het[i].has_trans()) {
			auto he0 = OM::SmartHalfedgeHandle(i, this);
			Eigen::Vector3<Real> v_he0 = (vlist.row(he0.to().idx()) - vlist.row(he0.from().idx())).transpose();
			v_he0 -= het[i].torus_trans;
			auto he1 = he0.next();
			Eigen::Vector3<Real> v_he1 = (vlist.row(he1.to().idx()) - vlist.row(he1.from().idx())).transpose();
			v_he1 -= het[he1.idx()].torus_trans;
			auto fh = he0.face();
			facearea[fh.idx()] = v_he0.cross(v_he1).norm() / 2;
		}
	}
	return facearea;
}

Eigen::Vector3<msf::Real> msf::PeriodSurfaceMesh::getFaceArea(OM::FaceHandle fh, Real deduce)
{
	auto fv = getFaceVertex(fh);
	fv.col(1) = make_period<Eigen::Vector3<Real>>(fv.col(1) - fv.col(0), 2, deduce) + fv.col(0);
	fv.col(2) = make_period<Eigen::Vector3<Real>>(fv.col(2) - fv.col(0), 2, deduce) + fv.col(0);
	return (fv.col(1) - fv.col(0)).cross(fv.col(2) - fv.col(0)) / 2;
}

Eigen::Matrix3<msf::Real> msf::PeriodSurfaceMesh::getFacePeriodVertex(OM::FaceHandle fh, Real deduce) const
{
	auto fv = getFaceVertex(fh);
	fv.col(1) = make_period<Eigen::Vector3<Real>>(fv.col(1) - fv.col(0), 2, deduce) + fv.col(0);
	fv.col(2) = make_period<Eigen::Vector3<Real>>(fv.col(2) - fv.col(0), 2, deduce) + fv.col(0);
	return fv;
}

Eigen::Matrix3<msf::Real> msf::PeriodSurfaceMesh::getFacePeriodVertex(OM::FaceHandle fh, OM::VertexHandle vhsync, Real deduce) const
{
	Eigen::Matrix3<msf::Real> tri = getFacePeriodVertex(fh, 1);
	auto fvh = getFaceVertexHandle(fh);
	for (int i = 0; i < 3; i++) {
		if (fvh[i] == vhsync) {
			auto p = toEigen(point(vhsync));
			if ((p - tri.col(i)).isZero()) break;
			tri.col((i + 1) % 3) = tri.col(i) + make_period((tri.col((i + 1) % 3) - tri.col(i)).eval(), 2, 1);
			tri.col((i + 2) % 3) = tri.col(i) + make_period((tri.col((i + 2) % 3) - tri.col(i)).eval(), 2, 1);
			break;
		}
	}
	return tri;
}

Eigen::MatrixX3<msf::Real> msf::PeriodSurfaceMesh::getVertexNormals(int type,
	const Eigen::MatrixX3<Real>& vlist, const std::vector<OM::SmartVertexHandle>& vcut, const std::vector<OM::SmartHalfedgeHandle>& hetrans, const std::vector<HETag>& hetags)
{
	auto flist = getFlist();
	Eigen::MatrixX3<Real> vn;
	if (type == 1) {
		igl::per_vertex_normals(vlist, flist, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA, vn);
	} else if (type == 2) {
		igl::per_vertex_normals(vlist, flist, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE, vn);
	}
	// fix it on transition region
	for (auto vh : vertices()) {
		bool has_trans = false;
		for (auto voh : voh_ccw_range(vh)) {
			if (hetags[voh.idx()].pass_boundary[0] + hetags[voh.idx()].pass_boundary[1] +
				hetags[voh.idx()].pass_boundary[2] > 0) {
				has_trans = true; break;
			}
		}
		if (has_trans) {
			auto [o, ring] = find1ring(vlist, vh, hetags);
			int N = ring.size();
			Eigen::Vector3<Real> n(0, 0, 0);
			Real As = 0;
			for (int i = 0; i < ring.size(); i++) {
				if (type == 1) {
					Eigen::Vector3<Real> An = (ring[i] - o).cross(ring[(i + 1) % N] - o) / 2;
					Real A = An.norm(); //  
					As += A;
					n += An;
				} else if (type == 2) {
					Eigen::Vector3<Real> e0 = (ring[i] - o).normalized(), e1 = (ring[(i + 1) % N] - o).normalized();
					Real A = std::acos(e0.dot(e1)); //  
					n += A * e0.cross(e1).normalized();
				}
			}
			// n.normalize();
			n /= (As / 3); // normal density
			vn.row(vh.idx()) = n.transpose();
		}
	}
	return vn;
}

std::tuple<Eigen::Vector3<msf::Real>, Eigen::Vector3<msf::Real>> msf::PeriodSurfaceMesh::getVertexNormal(const Eigen::Vector3<Real>& o, const std::vector<Eigen::Vector3<Real>>& ring, int type /*= 1*/)
{
	Eigen::Vector3<Real> vn; vn.setZero();
	Real As = 1e-30;
	Eigen::Vector3<Real> barycenter(0, 0, 0);;
	for (int i = 0; i < ring.size(); i++) {
		Eigen::Vector3<Real> fn = (ring[(i + 1) % ring.size()] - ring[i]).cross(ring[i] - o) / 2;
		Real A = fn.norm();
		As += A;
		barycenter += A * (ring[(i + 1) % ring.size()] + ring[i] + o) / 3;
		vn += fn;
	}
	barycenter /= As;
	CHECK_NAN(vn);
	return { vn.normalized(),barycenter };
}

Eigen::MatrixX3<msf::Real> msf::PeriodSurfaceMesh::getVertexNormal(Real period, Real deduce)
{
	Eigen::MatrixX3<Real> vnlist(n_vertices(), 3);
	for (auto vh : vertices()) {
		auto [o, ring] = find1ring(vh, deduce);
		auto [N, c] = getVertexNormal(o, ring);
		vnlist.row(vh.idx()) = N.transpose();
	}
	return vnlist;
}

Eigen::MatrixX<msf::Real> msf::PeriodSurfaceMesh::vecs2quats(const Eigen::MatrixX<Real>& v3list)
{
	int n = v3list.rows() / 3;
	Eigen::MatrixX<Real> q4list(n * 4, v3list.cols());
	q4list.setZero();
	for (int i = 0; i < v3list.cols(); i++) {
		for (int j = 0; j < n; j++) { q4list.block<3, 1>(j * 4 + 1, i) = v3list.block<3, 1>(j * 3, i); }
	}
	return q4list;
}

Eigen::MatrixX<msf::Real> msf::PeriodSurfaceMesh::quats2vecs(const Eigen::MatrixX<Real>& q4list)
{
	int n = q4list.rows() / 4;
	Eigen::MatrixX<Real> v3list(n * 3, q4list.cols());
	v3list.setZero();
	for (int i = 0; i < q4list.cols(); i++) {
		for (int j = 0; j < n; j++) { v3list.block<3, 1>(j * 3, i) = q4list.block<3, 1>(j * 4 + 1, i); }
	}
	return v3list;
}

msf::Real msf::PeriodSurfaceMesh::angle(const Eigen::Vector3<Real>& va, const Eigen::Vector3<Real>& vo, const Eigen::Vector3<Real>& vb)
{
	auto a = (va - vo).normalized();
	auto b = (vb - vo).normalized();
	return std::acos(a.dot(b));
}

//Eigen::Vector3<msf::Real> msf::PeriodSurfaceMesh::make_period(const Eigen::Vector3<Real>& vec, Real period /*= 2*/, Real deduce /*= 0.5*/)

msf::Real msf::PeriodSurfaceMesh::cosangle(const Eigen::Vector3<Real>& vo, OM::SmartHalfedgeHandle he)
{
	OM::Vec3d o(vo[0], vo[1], vo[2]);
	auto o1 = point(he.from()) - o;
	auto o2 = point(he.to()) - o;
	return o1.normalized().dot(o2.normalized());
}

std::vector<msf::OM::SmartVertexHandle> msf::PeriodSurfaceMesh::mergePeriodBoundary(void)
{
	mergePeriodEdges();
	return mergePeriodVertices();
}

Eigen::VectorXi msf::PeriodSurfaceMesh::findUniqueCorrespond(const Eigen::MatrixX3<Real>& p_uniq, const Eigen::MatrixX3<Real>& pdup)
{
	int n_old = p_uniq.rows();
	std::vector<int> newid2old(n_old);
	torus_kd_tree_t tree(2, 2, 2, 1e-7);
	tree.build(p_uniq);
	Eigen::VectorXi idlist(pdup.rows());
	for (int i = 0; i < pdup.rows(); i++) {
		Eigen::Vector3<Real> p = pdup.row(i).transpose();
		int oldid = tree.query(p);
		if (oldid == -1) {
			std::cout << "not point found at " << p.transpose() << std::endl;
			throw std::runtime_error("torus kd tree search failed"); 
		}
		idlist[i] = oldid;
	}
	return idlist;
}

msf::PeriodicGridIndex msf::PeriodSurfaceMesh::generateIndexer(const Eigen::MatrixX3<Real> vlist)
{
	PeriodicGridIndex indexer(Eigen::Vector3<Real>(-1, -1, -1), Eigen::Vector3<Real>(2, 2, 2));
	for (int i = 0; i < vlist.rows(); i++) {
		Eigen::Vector3<Real> p = vlist.row(i).transpose();
		indexer.insert(p);
	}
	return indexer;
}

Eigen::SparseMatrix<msf::Real> msf::PeriodSurfaceMesh::getPeriodicLaplacian(Real period, Real deduce)
{
	Eigen::SparseMatrix<Real> L(n_vertices(), n_vertices());
	std::vector<Eigen::Triplet<Real>> triplist;
	for (auto vh : vertices()) {
		Real w_sum = 0;
		auto [vv, voh] = find1ring(vh);
		auto [v0, ring] = find1ring(vh, 0.5);
		auto wlist = cot1ring(v0, ring);
		for (int i = 0; i < vv.size(); i++) {
			w_sum += wlist[i];
			triplist.emplace_back(vh.idx(), vv[i].idx(), wlist[i]);
		}
		triplist.emplace_back(vh.idx(), vh.idx(), -w_sum);
	}
	L.setFromTriplets(triplist.begin(), triplist.end());
	return L;
}

std::tuple<Eigen::MatrixX3<msf::Real>, Eigen::MatrixX3<msf::Real>> msf::PeriodSurfaceMesh::getFaceNormal(Real period, Real deduce)
{
	Eigen::MatrixX3<Real> fn(n_faces(), 3);
	Eigen::MatrixX3<Real> fc(n_faces(), 3);
	for (auto fh : faces()) {
		Eigen::Vector3<Real> hv[3], v[3];
		int counter = 0;
		for (auto he : fh_range(fh)) {
			hv[counter] = make_period(toEigen(calc_edge_vector(he)));
			v[counter] = toEigen(point(he.from()));
			counter++;
		}
		v[1] = v[0] + hv[0];
		v[2] = v[1] + hv[1];
		Eigen::Vector3<Real> c = (v[0] + v[1] + v[2]) / 3;
		for (int i = 0; i < 3; i++) {
			if (c[i] < -1) { c[i] += 2; }
		}
		fc.row(fh.idx()) = c.transpose();
		fn.row(fh.idx()) = hv[0].cross(hv[1]).normalized().transpose();
	}
	return { fc,fn };
}

bool msf::PeriodSurfaceMesh::has_boundary(void)
{
	for (auto vh : vertices()) {
		if (status(vh).deleted() || !vh.is_valid()) continue;
		if (vh.is_boundary())return true;
	}
	return false;
}

std::tuple<msf::Real, msf::Real> msf::PeriodSurfaceMesh::normalDistribution(const Eigen::Vector3<Real>& q)
{
	Real qs = 0;
	Real s = 0;
	//std::cout << "q = " << q.transpose() << std::endl;
	for (auto vh : vertices()) {
		auto [o, ring] = find1ring(vh, 0.5);
		auto [n, c] = getVertexNormal(o, ring);
		//std::cout << "n = " << n.transpose() << std::endl;
		auto [A, N] = area1ring(o, ring, 2);
		qs += (std::pow)(q.dot(n), 2) * A;
		s += A;
	}
	return { qs, s };
}

msf::Real msf::PeriodSurfaceMesh::period_cote(OM::EdgeHandle e)
{
	OM::SmartEdgeHandle eh(e.idx(), this);
	if (eh.is_boundary()) {
		throw std::runtime_error("evaluating cot on boundary edge");
	}
	auto he = eh.halfedge(0);

	auto vh = he.from();

	auto p = point(vh);

	auto q = p + make_period(calc_edge_vector(he));

	auto vlef = p - make_period(calc_edge_vector(he.prev()));
	auto vrig = p + make_period(calc_edge_vector(he.opp().next()));

	Real pq = (q - p).norm();
	Real lef_p = (p - vlef).norm();
	Real lef_q = (q - vlef).norm();
	Real rig_p = (p - vrig).norm();
	Real rig_q = (q - vrig).norm();
	Real cos_lef = (std::pow(lef_p, 2) + std::pow(lef_q, 2) - std::pow(pq, 2)) / (2 * lef_p * lef_q);
	Real cos_rig = (std::pow(rig_p, 2) + std::pow(rig_q, 2) - std::pow(pq, 2)) / (2 * rig_p * rig_q);

	Real cot_lef = cos_lef / (std::sqrt)(1 - cos_lef * cos_lef);
	Real cot_rig = cos_rig / (std::sqrt)(1 - cos_rig * cos_rig);
	return (cot_lef + cot_rig) / 2;
}

msf::Real msf::PeriodSurfaceMesh::periodicArea(msf::Real deduce)
{
	Real As = 0;
	for (auto fh : faces()) { As += getFaceArea(fh, deduce).norm(); }
	return As;
}

Eigen::MatrixX3<msf::Real> msf::PeriodSurfaceMesh::getFaceFrame(Real period, Real deduce)
{
	Eigen::MatrixX3<msf::Real> ff(n_faces() * 3, 3);
	for (auto fh : faces()) {
		auto vh = getFaceVertexHandle(fh);
		auto vt = getFacePeriodVertex(fh, deduce);
		Eigen::Vector3d n = (vt.col(1) - vt.col(0)).cross(vt.col(2) - vt.col(0)).normalized();
		Eigen::Vector3d t1 = (vt.col(1) - vt.col(0)).normalized();
		ff.row(fh.idx() * 3) = t1;
		ff.row(fh.idx() * 3 + 1) = n.cross(t1);
		ff.row(fh.idx() * 3 + 2) = n;
	}
	return ff;
}

Eigen::Vector3<msf::Real> msf::PeriodSurfaceMesh::eval_period_edge(OM::SmartHalfedgeHandle he, double t) const
{
	auto vec = calc_edge_vector(he);
	vec = make_period(vec, 2, 1);
	auto p0 = point(he.from());
	return toEigen(p0 + t * vec);
}

double msf::PeriodSurfaceMesh::period_sector_angle(OM::HalfedgeHandle h) const
{
	auto he = OM::SmartHalfedgeHandle(h.idx(), this);
	auto a = point(he.from());
	auto b = point(he.to());
	auto c = point(he.next().to());
	b = a + make_period(b - a);
	c = a + make_period(c - a);
	double ct = (c - b).dot(a - b) / (c - b).norm() / (a - b).norm();
	return std::acos(ct);
}

double msf::PeriodSurfaceMesh::period_face_area(OM::FaceHandle fh) const
{
	auto tri = getFacePeriodVertex(fh, 1);
	return (tri.col(1) - tri.col(0)).cross(tri.col(2) - tri.col(0)).norm() / 2;
}

double msf::PeriodSurfaceMesh::period_vertex_dual_area(OM::VertexHandle vh0) const
{
	OM::SmartVertexHandle vh(vh0.idx(), this);
	auto [o, ring] = find1ring(vh, 1.);

	double As = 0;
	int N = ring.size();
	for (int k = 0; k < ring.size(); k++) {
		As += (ring[(k + 1) % N] - o).cross(ring[k] - o).norm();
	}
	return As / 6;
}

void msf::PeriodSurfaceMesh::saveUnitCell(std::string filename)
{
	auto m_new = *this;
	m_new.split_unit_cell();
	m_new.savePeriodicMesh(filename, std::set<OM::SmartVertexHandle>{}, 1);
}

void savePeriodicVector(std::string filename, msf::PeriodSurfaceMesh& mesh, Eigen::VectorXd& v) {
	std::ofstream ofs(filename, std::ios::binary);
	for (auto vh : mesh.vertices()) {
		auto p = mesh.point(vh);
		Eigen::Vector3d pp;
		Eigen::Vector3d up = v.block<3, 1>(vh.idx() * 3, 0);
		for (int i = 0; i < (p[0] == -1) + 1; i++) {
			pp[0] = p[0] + i * 2;
			for (int j = 0; j < (p[1] == -1) + 1; j++) {
				pp[1] = p[1] + j * 2;
				for (int k = 0; k < (p[2] == -1) + 1; k++) {
					pp[2] = p[2] + k * 2;
					ofs.write((const char*)pp.data(), sizeof(pp));
					ofs.write((const char*)up.data(), sizeof(up));
				}
			}
		}
	}
}


