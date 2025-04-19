#define _USE_MATH_DEFINES
#include "PeriodicMesher.h"
#include "igl/read_triangle_mesh.h"
#include "queue"
#include <Eigen/Geometry>
#include "dir_utils.h"

using namespace msf;

std::string getPath(std::string s);

void test_conjugate_surf(std::string meshfile) {
	//Eigen::MatrixX3<Real> V; Eigen::MatrixX3i F;
	//igl::read_triangle_mesh(meshfile, V, F);
	
	MeshSurface mesh;
	mesh.read(meshfile, false);

	{
		std::ofstream ofs("vbd");
		for (auto vp : mesh.vertices()) {
			if (vp.is_boundary()) ofs << mesh.point(vp) << std::endl;
		}
	}

	std::vector<OM::Vec3d> Vdual(mesh.n_faces());

	std::queue<OM::SmartHalfedgeHandle> front;

	std::set<OM::EdgeHandle> ehset;

	for (auto eh : mesh.edges()) { ehset.emplace(eh); }

	auto f0 = *mesh.faces_begin();
	OM::SmartHalfedgeHandle he0;
	for (auto he : f0.halfedges()) {
		if (!he.opp().is_boundary() && !he.is_boundary()) { he0 = he; break; }
	}

	Vdual[f0.idx()] = OM::Vec3d(0, 0, 0);
	front.emplace(he0);
	ehset.erase(he0.edge());
	while (!front.empty()) {
		auto fr = front.front();
		front.pop();
		if (fr.is_boundary() || fr.opp().is_boundary()) continue;
		OM::Vec3d ev = mesh.calc_edge_vector(fr) * mesh.cote(fr.edge());
		Vdual[fr.opp().face().idx()] = Vdual[fr.face().idx()] + ev;
		for (auto fh : fr.opp().face().halfedges()) {
			if (fh == fr.opp()) { continue; }
			if (ehset.count(fh.edge())) {
				front.push(fh); 
				ehset.erase(fh.edge());
			}
		}
	}

	//{ std::ofstream ofs(getPath("Vdual"), std::ios::binary); ofs.write((const char*)Vdual.data(), Vdual.size() * sizeof(Vdual[0])); }

	MeshSurface conj = mesh;
	
	for (auto vh : conj.vertices()) {
		Eigen::MatrixX3<Real> A;
		Eigen::VectorX<Real> b;
		OM::Vec3d cdual(0, 0, 0);
		if (vh.is_boundary()) {
			//if (vh.valence() == 3) { std::cout << "Found corner" << std::endl; }
			int valen = 0;
			Eigen::Matrix3<Real> translayer = Eigen::Matrix3<Real>::Identity();
			OM::Vec3d res(0, 0, 0);
			bool ccw = true;
			auto h0 = *vh.outgoing_halfedges_ccw().begin();
			Eigen::Vector3<Real> p0;
			if (h0.is_boundary()) {
				res = mesh.cothe(h0.opp()) * mesh.calc_edge_vector(h0);
				p0 = toEigen(Vdual[h0.opp().face().idx()]);
			}
			else if (h0.opp().is_boundary()) {
				res = mesh.cothe(h0) * mesh.calc_edge_vector(h0);
				// plus res because we need to propagate dual points based on inner points, instead of boundary
				p0 = toEigen(Vdual[h0.face().idx()] + res);
			}
			else {
				res = mesh.cote(h0.edge()) * mesh.calc_edge_vector(h0);
				p0 = toEigen(Vdual[h0.opp().face().idx()]);
			}
			valen++;
			auto h = h0;
			std::vector<Eigen::Vector3<Real>> plist{ p0 };
			cdual += toOM(p0);
			std::vector<OM::Vec3d> res_his;
			//std::cout << "V = " << mesh.point(vh) << " p0 = " << p0.transpose() << std::endl;
			do {
				res_his.push_back(res);
				if (valen > 15) {
					std::cout << "vp = " << conj.point(vh) << std::endl;
					std::cout << "plist = \n";
					for (auto pk : plist) { std::cout << pk.transpose() << std::endl; }
					std::cout << "res = \n";
					for (auto rp : res_his) { std::cout << rp << std::endl; }
					throw std::runtime_error("boundary loop overflow, possibly imperfect minimal surface");
				}
				// minus res because we need to propagate dual points based on inner dual points, instead of on boundary
				plist.emplace_back(p0 - toEigen(res));
				//std::cout << "p = " << plist.rbegin()->transpose() << std::endl;
				cdual += toOM(*plist.rbegin());
				if ((h.is_boundary() && ccw) || (h.opp().is_boundary() && !ccw)) {
					ccw = !ccw; 
					translayer = (Eigen::AngleAxis<Real>(M_PI, (translayer * toEigen(mesh.calc_edge_vector(h))).normalized()) * translayer).eval();
				}
				if (ccw) {
					h = h.prev().opp();
				} else {
					h = h.opp().next();
				}
				Eigen::Vector3<Real> hv;
				if (h.is_boundary()) {
					hv = toEigen(mesh.cothe(h.opp()) * mesh.calc_edge_vector(h));
				} else if (h.opp().is_boundary()) {
					hv = toEigen(mesh.cothe(h) * mesh.calc_edge_vector(h));
				} else {
					hv = toEigen(mesh.cote(h.edge()) * mesh.calc_edge_vector(h));
				}
				res += toOM((translayer * hv).eval());
				//std::cout << "h = " << (translayer * hv).transpose() << std::endl;
				valen++;
			} while (res.sqrnorm() > 4e-6);
			cdual /= plist.size();
			A.resize(plist.size(), 3); b.resize(plist.size());
			for (int i = 0; i < plist.size(); i++) {
				Eigen::Vector3<Real> di = plist[(i + 1) % plist.size()] - plist[i];
				Eigen::Vector3<Real> ci = (plist[(i + 1) % plist.size()] + plist[i]) / 2;
				A.row(i) = di.transpose();
				b[i] = di.dot(ci);
			}
		}
		else {
			int counter = 0;
			A.resize(vh.valence(), 3); b.resize(vh.valence());
			for (auto voh : vh.outgoing_halfedges_ccw()) {
				OM::Vec3d vi_1 = Vdual[voh.opp().face().idx()];
				OM::Vec3d vi = Vdual[voh.face().idx()];
				cdual += vi;
				A.row(counter) = toEigen(vi - vi_1).transpose();
				b[counter] = (vi - vi_1).dot(vi + vi_1) / 2;
				counter++;
			}
			cdual /= counter;
		}
		conj.set_point(vh, cdual);
		//Eigen::Vector3<Real> x = (A.transpose() * A).fullPivLu().solve(A.transpose() * b);
		//conj.set_point(vh, toOM(x));
	}
	auto fn = dir_utils::path2filename(meshfile);
	conj.write(getPath(fn + ".conj.obj"));
}
