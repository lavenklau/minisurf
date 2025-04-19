#pragma once
#include "tuple"
#include "Eigen/Eigen"
#include <set>

template<typename T>
std::tuple<Eigen::MatrixX3<T>, Eigen::MatrixX3i> cgal_remesh(const Eigen::MatrixX3<T>& V, const Eigen::MatrixX3i& F, T len = 0.03);

template<typename T>
std::vector<int> cgal_kd_tree_search(const Eigen::MatrixX3<T>& Vcloud, const Eigen::MatrixX3<T>& toSearch);

std::tuple<Eigen::MatrixX3d, Eigen::MatrixX3i> patch_holes(const Eigen::MatrixX3d& V, const Eigen::MatrixX3i& F, const std::set<int>& cand_hole_edges);
