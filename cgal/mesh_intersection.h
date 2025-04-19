#pragma once

#include "Config.h"
#include "Eigen/Eigen"

namespace msf {

	struct UDEdge : std::pair<int, int> {
		UDEdge(int v0, int v1) {
			if (v0 > v1) {
				first = v1; second = v0;
			}
			else { first = v0; second = v1; }
		}
		bool connect(const UDEdge& e2) {
			return first == e2.first || first == e2.second || second == e2.first || second == e2.second;
		}
	};


	std::map<UDEdge, char> extractNonManifoldUDEdges(const Eigen::MatrixX3i& Flist);

	std::vector<std::vector<int>> connected_components(const Eigen::MatrixX3<Real>& Vlist, const Eigen::MatrixX3i& Flist);


	void mesh_split(Eigen::MatrixX3<Real>& V, Eigen::MatrixX3i& F, const  Eigen::MatrixX3<Real>& spliter_V2, const Eigen::MatrixX3i& spliter_F2);
	void mesh_split_non_manifold_vertices(Eigen::MatrixX3<Real>& V, Eigen::MatrixX3i& F);

};

