#include "PeriodMeshConvertion.h"
#include <iostream>
#include <fstream>
//#include <utils/HashTypes.h>
#include <algorithm>
#include <numeric>
//#include "cmdline.hpp"

using namespace msf;

std::pair<std::vector<int>, std::vector<VertexFlag>>
msf::removeDupDof(
	std::vector<Eigen::Vector3d>& vlist, std::vector<std::array<int, 4>>& elist,
	std::unordered_set<std::array<int, 3>>& boundaryFaces, std::map<int, int>& vMuxMap,
	std::map<int, int>& vnew2old
) {
	// checkout boundary faces and vertices
	std::unordered_map<Eigen::Vector3d, Eigen::Vector3d> periodVertices;
	boundaryFaces.clear();
	for (int i = 0; i < elist.size(); i++) {
		std::array<int, 3> f_arr;
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 3; k++) {
				f_arr[k] = elist[i][(j + k) % 4];
			}
			std::sort(f_arr.begin(), f_arr.end());
			if (boundaryFaces.count(f_arr)) boundaryFaces.erase(f_arr);
			else boundaryFaces.insert(f_arr);
		}
	}
	for (auto f : boundaryFaces) { for (auto v : f)  periodVertices[vlist[v]] = vlist[v]; }
	// DEBUG
	if(0) {
		std::ofstream ofs("periodVertices.txt");
		for (auto it = periodVertices.begin(); it != periodVertices.end(); it++) {
			ofs << it->first.transpose() << std::endl;
		}
		ofs.close();
	}
	// mark periodic points
	std::unordered_map<Eigen::Vector3d, VertexFlag> vperiodFlags;
	{
		using num_t = double;
#if 0
		nanoflann::PointCloud<num_t> cloud;
		cloud.pts.resize(periodVertices.size());
		{
			int counter = 0;
			for (auto it = periodVertices.begin(); it != periodVertices.end(); it++) {
				// only find points near boundary
				/*if (it->first.norm() > 0.6)*/ {
					auto p = it->first;
					cloud.pts[counter++] = { p[0],p[1],p[2] };
				}
			}
		}
		using kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
			nanoflann::T3_Adaptor<num_t, nanoflann::PointCloud<num_t>>,
			nanoflann::PointCloud<num_t>, 3 /* dim */ >;
#else
		using kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor <
			nanoflann::L2_Adaptor<num_t, nanoflann::PointCloud<num_t, 4>>,
			nanoflann::PointCloud<num_t, 4>, 4, size_t>;
		std::vector<Eigen::Vector3d> torusParamlist;
		{
			for (auto it = periodVertices.begin(); it != periodVertices.end(); it++) {
				torusParamlist.push_back(it->first);
			}
		}
		auto cloud = nanoflann::embedTorusPoints((double*)torusParamlist.data(), periodVertices.size());
#endif
		kd_tree_t tree(4, cloud);
		tree.buildIndex();
		// checkout period vertices on boundary
		for (auto v : periodVertices) {
			const size_t num_results = 2;
			std::vector<size_t> ret_indices(num_results);
			std::vector<num_t> out_dist_sqr(num_results);
#if 0
			nanoflann::KNNResultSet<num_t> resultSet(num_results);
			resultSet.init(ret_indices.data(), out_dist_sqr.data());
			tree.findNeighbors(resultSet, &p[0], nanoflann::SearchParams());
			if (out_dist_sqr[0] < 1e-10 && out_dist_sqr[1] < 1e-10) {
				// period vertex
				Eigen::Vector3d q[2] = {
					{cloud.pts[ret_indices[0]].x, cloud.pts[ret_indices[0]].y, cloud.pts[ret_indices[0]].z},
					{cloud.pts[ret_indices[1]].x, cloud.pts[ret_indices[1]].y, cloud.pts[ret_indices[1]].z}
				};
				if (q[0] == q[1]) { printf("\033[31mWarning duplicated points in Kd-tree\033[0m\n"); }
				Eigen::Vector3d qq = q[1] - q[0];
				if (q[0] < q[1]) {
					periodVertices[q[1]] = q[0];
				}
			}
			num_t p[3] = { v.first[0], v.first[1], v.first[2] };
#else
			Eigen::Vector4d p = nanoflann::embedTorusPoint(v.first.data());
			std::vector<nanoflann::ResultItem<size_t, num_t>> indices_dists;
			//nanoflann::RadiusResultSet<num_t> resultSet(1e-4, indices_dists);
			tree.radiusSearch(&p[0], 1e-12, indices_dists, nanoflann::SearchParameters());
			std::vector<Eigen::Vector3d> p_dup;
			Eigen::Vector3d p_min; p_min.setOnes(); p_min *= 1e30;
			// non period vertex
			if (indices_dists.size() == 1) {}
			// face,edge, corner period vertex
			else if (indices_dists.size() == 2 || indices_dists.size() == 4 || indices_dists.size() == 8) {
				std::vector<Eigen::Vector3d> plist(indices_dists.size());
				for (int i = 0; i < indices_dists.size(); i++) {
					//auto p = cloud.pts[indices_dists[i].first]; plist[i] = Eigen::Vector3d(p.x, p.y, p.z);
					auto p = torusParamlist[indices_dists[i].first]; plist[i] = p;
					if (!(p_min < plist[i]))  p_min = plist[i];
				}
				for (int i = 0; i < indices_dists.size(); i++) { periodVertices[plist[i]] = p_min; }
				bool p_xyz_min[3] = { false,false,false };
				for (int i = 0; i < plist.size(); i++) {
					Eigen::Vector3d pp = plist[i] - p_min;
					if (pp[0] < 0 || pp[1] < 0 || pp[2] < 0) { printf("\033[31mUnexpected error\033[0m\n"); }
					bool q_xyz_max[3] = { pp[0] > 1e-4, pp[1] > 1e-4, pp[2] > 1e-4 };
					for (int j = 0; j < 3; j++) p_xyz_min[j] = p_xyz_min[j] || q_xyz_max[j];
				}
				vperiodFlags[p_min].set_period_boundary(p_xyz_min[0], p_xyz_min[1], p_xyz_min[2]);
				for (int i = 0; i < plist.size(); i++) {
					Eigen::Vector3d pp = plist[i] - p_min;
					bool q_xyz_max[3] = { pp[0] > 1e-4, pp[1] > 1e-4, pp[2] > 1e-4 };
					bool q_xyz_min[3] = { abs(pp[0]) <= 1e-4 && p_xyz_min[0],
						abs(pp[1]) <= 1e-4 && p_xyz_min[1], abs(pp[2]) <= 1e-4 && p_xyz_min[2] };
					int flg[3];
					for (int j = 0; j < 3; j++) flg[j] = int(q_xyz_min[j]) | (int(q_xyz_max[j]) << 1);
					vperiodFlags[plist[i]].set_period_boundary(flg[0], flg[1], flg[2]);
				}
			}
			else {
				printf("\033[31mUnexpected period duplicate number %d, on Point : (%e, %e, %e)\033[0m",
					int(indices_dists.size()), p[0], p[1], p[2]);
				printf("[ ");
				for (int kk = 0; kk < indices_dists.size(); kk++) {
					printf("\033[33m%d \033[0m", indices_dists[kk].first);
				}
				printf("]\n");
			}
#endif
		}
	}
	// DEBUG
	if(0) {
		std::ofstream ofs("vperiod.txt");
		for (auto vf : vperiodFlags) {
			if (vf.second.is_min_period()) {
				ofs << "min " << vf.first.transpose() << std::endl;
			}
			if (vf.second.is_max_period()) {
				ofs << "max " << vf.first.transpose() << std::endl;
			}
		}
		ofs.close();
	}
	std::unordered_map<Eigen::Vector3d, int> vindex;
	std::unordered_map<Eigen::Vector3d, int> vdofs;
	// remove duplicated vertices and dof
	for (int i = 0; i < elist.size(); i++) {
		for (int j = 0; j < 4; j++) {
			if (!vindex.count(vlist[elist[i][j]])) {
				auto siz = vindex.size();
				// add vertex index
				vindex[vlist[elist[i][j]]] = siz;
				// add vertex dof
				if (!vdofs.count(vlist[elist[i][j]])) {
					auto dofsiz = vdofs.size();
					auto it = vperiodFlags.find(vlist[elist[i][j]]);
					if (it != vperiodFlags.end()) {
						if ((it->second.is_min_period() && !it->second.is_max_period()) || !it->second.is_period_boundary()) {
							vdofs[it->first] = dofsiz;
						}
					} else {
						vdofs[vlist[elist[i][j]]] = dofsiz;
					}
				}
			}
		}
	}
	// add period vertices dof to its equivalent
	for (auto it = periodVertices.begin(); it != periodVertices.end(); it++) {
		if (it->first != it->second) {
			if (vdofs.find(it->second) == vdofs.end()) {
				printf("\033[31mPeriod source is not set yet!\033[0m\n");
			}
			vdofs[it->first] = vdofs[it->second];
		}
	}
	// remap vertex index in cells
	for (int i = 0; i < elist.size(); i++) {
		for (auto vit = elist[i].begin(); vit != elist[i].end(); vit++) {
			*vit = vindex[vlist[*vit]];
		}
	}
	// remove duplicated cells
	for (int i = 0; i < elist.size(); i++) {
		std::sort(elist[i].begin(), elist[i].end());
	}
	{ auto uit = std::unique(elist.begin(), elist.end()); elist.erase(uit, elist.end()); }
	// new boundary vertices
	std::unordered_set<std::array<int, 3>> new_boundaryFaces;
	for (auto it = boundaryFaces.begin(); it != boundaryFaces.end(); it++) {
		auto arr = *it;
		for (int i = 0; i < 3; i++) arr[i] = vindex[vlist[arr[i]]];
		std::sort(arr.begin(),arr.end());
		new_boundaryFaces.insert(arr);
	}
	std::swap(boundaryFaces, new_boundaryFaces);
	// reorder vertices
	std::vector<Eigen::Vector3d> new_vlist;
	vnew2old.clear();
	for (auto vpos2id : vindex) {
		if (new_vlist.size() <= vpos2id.second) {
			new_vlist.resize(vpos2id.second + 1);
		}
		new_vlist[vpos2id.second] = vpos2id.first;
	}
	for (int i = 0; i < vlist.size(); i++) {
		auto it = vindex.find(vlist[i]);
		if (it != vindex.end()) { vnew2old[it->second] = i; }
	}
	// update vertices list
	std::swap(vlist, new_vlist);
	// get dof vector
	std::vector<int> vdoflist(vlist.size());
	std::vector<VertexFlag> vflags(vlist.size());
	for (int i = 0; i < vdoflist.size(); i++) {
		auto it = vdofs.find(vlist[i]);
		if (it == vdofs.end()) {
			printf("\033[31mWarning: unrecorded vertex [%d](%e, %e, %e)\033[0m\n", i, vlist[i][0], vlist[i][1], vlist[i][2]);
		}
		vdoflist[i] = vdofs[vlist[i]];
	}
	for (int i = 0; i < vlist.size(); i++) {
		if (vperiodFlags.count(vlist[i])) {
			vflags[i] = vperiodFlags[vlist[i]];
		}
	}
	for (auto it = periodVertices.begin(); it != periodVertices.end(); it++) {
		if (it->first != it->second) {
			vMuxMap[vindex[it->first]] = vindex[it->second];
		}
	}
	return { vdoflist, vflags };
}

struct Ratio {
	union {
		int Numerator = 1;
		int b;
	};
	union {
		int denominator = 1;
		int a;
	};
	operator double() const{
		return double(Numerator) / denominator;
	}
	bool nonNegative(void) const{
		return !(Numerator >= 0 ^ denominator >= 0);
	}
	void reduce(void) {
		int g = std::gcd(Numerator, denominator);
		Numerator /= g;
		denominator /= g;
	}
	bool isLarge(void) const{
		return Numerator > 1e8 || denominator > 1e8 || Numerator < -1e8 || denominator < -1e8;
	}
	Ratio operator*(const Ratio& z) const{
		Ratio res;
		res.Numerator = Numerator * z.Numerator;
		res.denominator = denominator * z.denominator;
		return res;
	}
	bool operator==(const Ratio& z) const {
		return Numerator * z.denominator == denominator * z.Numerator;
	}
	bool operator<(const Ratio& z) const {
		return (b * z.a - a * z.b) * (a * z.a) < 0;
	}
	bool operator>(const Ratio& z) const {
		return (b * z.a - a * z.b) * (a * z.a) > 0;
	}
	bool operator>=(const Ratio& z) const {
		return (b * z.a - a * z.b) * (a * z.a) >= 0;
	}
	bool operator<=(const Ratio& z) const {
		return (b * z.a - a * z.b) * (a * z.a) <= 0;
	}
	Ratio operator/(const Ratio& z) const{
		Ratio res;
		res.Numerator = Numerator * z.denominator;
		res.denominator = denominator * z.Numerator;
	}
	Ratio operator-() const noexcept {
		Ratio z{-Numerator, -denominator};
		return z;
	}
	Ratio& operator/=(const Ratio& z) {
		Numerator *= z.denominator;
		denominator *= z.Numerator;
		return *this;
	}
	Ratio& operator*=(const Ratio& z) {
		Numerator *= z.Numerator;
		denominator *= z.denominator;
		return *this;
	}
	Ratio& operator+=(const Ratio& z) {
		// ToDo : add overflow test
		Numerator = Numerator * z.denominator + denominator * z.Numerator;
		denominator *= z.denominator;
		reduce();
		return *this;
	}
	Ratio& operator-=(const Ratio& z) {
		Ratio a = -z;
		return (*this) += a;
	}
	Ratio operator-(const Ratio &z) const {
		Ratio a = *this;
		a -= z;
		return a;
	}
	Ratio operator+(const Ratio &z) const {
		Ratio a = *this;
		a += z;
		return a;
	}
};

template<typename T>
struct Box {
	std::array<Eigen::Vector3<T>, 8> coorners;
	Box(Eigen::Vector3<T> mi, Eigen::Vector3<T> ma) {
		coorners[0] = Eigen::Vector3<T>(mi[0], mi[1], mi[2]);
		coorners[1] = Eigen::Vector3<T>(ma[0], mi[1], mi[2]);
		coorners[2] = Eigen::Vector3<T>(ma[0], ma[1], mi[2]);
		coorners[3] = Eigen::Vector3<T>(mi[0], ma[1], mi[2]);
		coorners[4] = Eigen::Vector3<T>(mi[0], mi[1], ma[2]);
		coorners[5] = Eigen::Vector3<T>(ma[0], mi[1], ma[2]);
		coorners[6] = Eigen::Vector3<T>(ma[0], ma[1], ma[2]);
		coorners[7] = Eigen::Vector3<T>(mi[0], ma[1], ma[2]);
	}
	template<typename Map>
	std::array<int, 8> getIndex(Map& map) {
		std::array<int, 8> ind;
		for (int i = 0; i < 8; i++)
			ind[i] = map[coorners[i]];
		return ind;
	}
};

std::pair<std::vector<Eigen::Vector3d>,std::vector<std::array<int,4>>> sampleTetmesh(int id)
{
	std::vector<Eigen::Vector3d> vertices;
	std::vector<std::array<int, 4>> tets;
#if 0
	if (id == 0) {
		// 3D cross
		std::vector<Box<int>> boxes;
		boxes.emplace_back(Eigen::Vector3<int>{1, 1, 0}, Eigen::Vector3<int>{2, 2, 1});
		boxes.emplace_back(Eigen::Vector3<int>{1, 1, 1}, Eigen::Vector3<int>{2, 2, 2});
		boxes.emplace_back(Eigen::Vector3<int>{1, 1, 2}, Eigen::Vector3<int>{2, 2, 3});
		boxes.emplace_back(Eigen::Vector3<int>{1, 0, 1}, Eigen::Vector3<int>{2, 1, 2});
		boxes.emplace_back(Eigen::Vector3<int>{1, 2, 1}, Eigen::Vector3<int>{2, 3, 2});
		boxes.emplace_back(Eigen::Vector3<int>{0, 1, 1}, Eigen::Vector3<int>{1, 2, 2});
		boxes.emplace_back(Eigen::Vector3<int>{2, 1, 1}, Eigen::Vector3<int>{3, 2, 2});
		std::unordered_map<Eigen::Vector3i, int> coords;
		for (int i = 0; i < boxes.size(); i++) {
			for (int j = 0; j < 8; j++) {
				auto c = boxes[i].coorners[j]; 
				if(!coords.count(c)) { coords[c] = coords.size(); }
			}
		}
		vertices.resize(coords.size());
		for(auto v2id : coords) {
			Eigen::Vector3d p = v2id.first.cast<double>() * 2 / 3 - Eigen::Vector3d(1., 1., 1.);
			vertices[v2id.second] = p;
		}
		for (int i = 0; i < boxes.size(); i++) {
			auto btets = splitBox2Tets(boxes[i].getIndex(coords));
			tets.insert(tets.end(), btets.begin(), btets.end());
		}
	}

	for (int i = 0; i < vertices.size(); i++) {
		if (getConfig().usr_number == 1) {
			std::swap(vertices[i][0], vertices[i][1]);
		} else if (getConfig().usr_number == 2) {
			std::swap(vertices[i][1], vertices[i][2]);
		} else if(getConfig().usr_number == 3) {
			std::swap(vertices[i][0], vertices[i][2]);
		}
		else if (getConfig().usr_number == 4)
		{
			double t = vertices[i][1];
			vertices[i][1] = vertices[i][0];
			vertices[i][0] = - t;
		}
		else if (getConfig().usr_number == 5)
		{
			double t = vertices[i][2];
			vertices[i][2] = vertices[i][1];
			vertices[i][1] = - t;
		}
		else if (getConfig().usr_number == 6)
		{
			double t = vertices[i][2];
			vertices[i][2] = vertices[i][0];
			vertices[i][0] = - t;
		}
		else if (getConfig().usr_number == 7) {
			double t = vertices[i][2];
			vertices[i][2] = vertices[i][0];
			vertices[i][0] = - t;
			t = vertices[i][2];
			vertices[i][2] = vertices[i][1];
			vertices[i][1] = - t;
		}
	}

#endif
	return {vertices, tets};
}

void torus_kd_tree_t::build(const Eigen::MatrixX3<Real>& pmat)
{
	Eigen::Matrix3X<Real> pcol(3, pmat.rows());
	for (int k = 0; k < 3; k++)
		pcol.row(k) = pmat.col(k).transpose() / _period[k] * 2;
	cloud = nanoflann::embedTorusPoints(pcol.data(), pmat.rows());
	_tree = std::make_unique<kd_tree_t>(4, cloud);
	_tree->buildIndex();
}

int torus_kd_tree_t::query(const Eigen::Vector3<Real>& p_torus)
{
	Eigen::Vector3<Real> p_parm = p_torus.cwiseQuotient(_period / 2);
	Eigen::Vector4<Real> p = nanoflann::embedTorusPoint(p_parm.data());
	std::vector<nanoflann::ResultItem<size_t, Real>> indices_dists;
	//nanoflann::RadiusResultSet<num_t> resultSet(1e-4, indices_dists);
	_tree->radiusSearch(&p[0], _tol, indices_dists, nanoflann::SearchParameters());
	if (indices_dists.empty())return -1;
	else { return indices_dists[0].first; }
}
