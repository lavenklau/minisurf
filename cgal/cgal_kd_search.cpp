#include "cgal_utils.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <list>
#include <cmath>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_d;
typedef CGAL::Search_traits_3<K> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree Tree;

template<typename T>
std::vector<int> cgal_kd_tree_search(const Eigen::MatrixX3<T>& Vcloud, const Eigen::MatrixX3<T>& toSearch) {
	const unsigned int N = 1;
	std::list<Point_d> points;
	std::map<Point_d, int> p2id;
	for (int i = 0; i < Vcloud.rows(); i++) {
		points.emplace_back(Vcloud(i, 0), Vcloud(i, 1), Vcloud(i, 2));
		p2id[*points.rbegin()] = i;
	}
	Tree tree(points.begin(), points.end());
	// Initialize the search structure, and search all N points
	std::vector<int> pid;
	for (int i = 0; i < toSearch.rows(); i++) {
		Point_d query(toSearch(i, 0), toSearch(i, 1), toSearch(i, 2));
		Neighbor_search search(tree, query, N);
		// report the N nearest neighbors and their distance
	   // This should sort all N points by increasing distance from origin
		if (search.begin()->second > 1e-5) {
			std::cout << "Warning : possible unmatch point " << query << std::endl;
		}
		pid.push_back(p2id.at(search.begin()->first));
	}
	return pid;
}

template std::vector<int> cgal_kd_tree_search<double>(const Eigen::MatrixX3<double>& Vcloud, const Eigen::MatrixX3<double>& toSearch);
template std::vector<int> cgal_kd_tree_search<float>(const Eigen::MatrixX3<float>& Vcloud, const Eigen::MatrixX3<float>& toSearch);