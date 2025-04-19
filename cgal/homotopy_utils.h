#pragma once

#include <vector>
#include <Eigen/Eigen>
#include <functional>

std::vector<std::vector<size_t>> extract_shortest_homotopy_path(const Eigen::Matrix<double, -1, 3>& V, const Eigen::Matrix<int, -1, 3>& F, const std::vector<int>& baselist);
std::vector<size_t> extract_shortest_homotopy_path(const Eigen::Matrix<double, -1, 3>& V, const Eigen::Matrix<int, -1, 3>& F, double period, double deduce);



