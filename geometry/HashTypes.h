#pragma once

#include <Eigen/Eigen>
#include "geometry/surface/mesh.h"

namespace std {
	template <class T>
	inline void std_hash_combine(size_t& seed, const T& v) {
		hash<T> hasher;
		seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
	}

	template<typename T, int N>
	struct hash<Eigen::Matrix<T, N, 1>> {
		size_t operator()(Eigen::Matrix<T, N, 1> const& arr) const noexcept {
			size_t hs = 0;
			for (auto iter = arr.begin(); iter != arr.end(); iter++) {
				std_hash_combine(hs, *iter);
			}
			return hs;
		}
	};

	template<typename T, size_t N>
	struct hash<array<T, N>> {
		size_t operator()(array<T, N> const& arr) const noexcept {
			size_t hs = 0;
			for (auto iter = arr.begin(); iter != arr.end(); iter++) {
				std_hash_combine(hs, *iter);
			}
			return hs;
		}
	};

	template<typename T1, typename T2>
	struct hash<pair<T1, T2>> {
		size_t operator()(pair<T1,T2> const& arr) const noexcept {
			size_t hs = 0;
			std_hash_combine(hs, arr.first);
			std_hash_combine(hs, arr.second);
			return hs;
		}
	};



	// Custom hash function for an array of 3 double values
	template<>
	struct hash<OpenMesh::Vec3d> {
		std::size_t operator()(const OpenMesh::Vec3d& arr) const {
			std::size_t seed = 0;
			// Combine the hash values of the individual elements
			for (const double& value : arr) {
				// Use the hash function from std::hash for doubles
				seed ^= std::hash<double>{}(value)+0x9e3779b9 + (seed << 6) + (seed >> 2);
			}
			return seed;
		}
	};

	template<>
	struct hash<OpenMesh::SmartHalfedgeHandle> : public hash<OpenMesh::HalfedgeHandle> {};
}


