#include "reorient_mesh.h"
#include "igl/triangle_triangle_adjacency.h"
#include "igl/read_triangle_mesh.h"
#include "igl/writeSTL.h"
#include "dir_utils.h"
#include <queue>
#include <set>

std::string getPath(std::string);

void msf::reorient_mesh(const Eigen::MatrixX3<Real>& V, Eigen::MatrixX3i& F)
{
	//
	{
		//igl::writeSTL(getPath("befori.stl"), V, F, igl::FileEncoding::Binary);
	}

	// find triangle with correct orientation
	Real xmax = -1e30;
	int v_xmax = -1;
	for (int i = 0; i < V.rows(); i++) {
		if (V(i, 0) > xmax) {
			v_xmax = i; xmax = V(i, 0);
		}
	}
	std::set<int> f_xmax;
	for (int j = 0; j < 3; j++) {
		for (int i = 0; i < F.rows(); i++) {
			if (F(i, j) == v_xmax) {
				f_xmax.insert(i);
			}
		}
	}
	Real cmax = -1e30;
	int cmax_id = -1;
	bool flip_max = false;
	for (auto f : f_xmax) {
		Eigen::Matrix3<Real> fv;
		for (int i = 0; i < 3; i++) { fv.col(i) = V.row(F(f, i)).transpose(); }
		Eigen::Vector3<Real> c = fv.rowwise().sum();
		if (c[0] > cmax) {
			cmax = c[0]; cmax_id = f; 
			flip_max = (fv.col(1) - fv.col(0)).cross(fv.col(2) - fv.col(0)).dot(Eigen::Vector3<Real>::UnitX()) < 0;
		}
	}
	if (flip_max) { std::swap(F(cmax_id, 0), F(cmax_id, 1)); }

	// unify orientation
	Eigen::MatrixX3i FF;
	igl::triangle_triangle_adjacency(F, FF);
	std::vector<bool> passed(F.rows(), false);
	std::queue<std::pair<int, int>> frontier;
	frontier.emplace(cmax_id, cmax_id);
	passed[cmax_id] = true;
	do {
		auto fr = frontier.front().second;
		auto fr_src = frontier.front().first;
		frontier.pop();
		if (fr != fr_src) {
			Eigen::Vector3i src = F.row(fr_src).transpose();
			Eigen::Vector3i dst = F.row(fr).transpose();
			int common_in_src[2], common_in_dst[2];
			int counter = 0;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					if (src[i] == dst[j]) {
						common_in_src[counter] = i; common_in_dst[counter] = j;
						counter++; break;
					}
				}
			}
			//common_in_dst[1] = (common_in_dst[1] + 3 + (common_in_dst[0] - common_in_src[0])) % 3;
			if (((common_in_dst[1] - common_in_dst[0]) - (common_in_src[1] - common_in_src[0])) % 3 == 0) {
				// this means wrong orientation
				std::swap(F(fr, common_in_dst[0]), F(fr, common_in_dst[1]));
			}
		}
		// add neighbor frontier 
		for (int i = 0; i < 3; i++) {
			try {
				if (passed.at(FF(fr, i))) continue;
				passed[FF(fr, i)] = true;
				frontier.emplace(fr, FF(fr, i));
			}
			catch (...) {
				printf("\033[31mOut of Range : FF(%d, %d) = %d\033[0m\n", fr, i, (int)FF(fr, i));
				std::cout << "F =\n" << V.row(F(fr, 0)) << "\n" <<
					V.row(F(fr, 1)) << "\n" <<
					V.row(F(fr, 2)) << "\n";
			}
		}
	} while (!frontier.empty());
}

extern std::string meshfile;

std::string getPath(std::string);

std::tuple<Eigen::MatrixX3<msf::Real>, Eigen::MatrixX3i> removeDupVertices(const Eigen::MatrixX3<msf::Real>& v, const Eigen::MatrixX3i& f, msf::Real eps /*= 1e-5*/);

void test_orient_mesh(void) {
	Eigen::MatrixX3<msf::Real> V;
	Eigen::MatrixX3i F;
	igl::read_triangle_mesh(meshfile, V, F);
	std::tie(V, F) = removeDupVertices(V, F, 1e-5);
	msf::reorient_mesh(V, F);
	auto fn = dir_utils::path2filename(meshfile);
	igl::writeSTL(getPath(fn) + "_ro.stl", V, F, igl::FileEncoding::Binary);
}
