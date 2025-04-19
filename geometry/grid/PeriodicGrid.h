#pragma once

#include "Config.h"
#include <vector>
#include <Eigen/Eigen>
#include <unordered_map>
#include "HashTypes.h"

namespace msf {
	class PeriodicGridIndex {
		Eigen::Vector3d zero, diag;
		double h;
		std::unordered_map<Eigen::Vector3i, int> activeLattice;
		std::vector<Eigen::Vector3d> plist;
		Eigen::Vector3i Ncell;
		int insert_counter = 0;
	protected:
		template<typename P>
		void periodize(P& p) const {
			for (int i = 0; i < 3; i++) {
				p[i] = (p[i] % Ncell[i] + Ncell[i]) % Ncell[i];
			}
		}
		template<typename P>
		auto raster(const P& p) const {
			Eigen::Vector3i p_ind;
			for (int i = 0; i < 3; i++) {
				p_ind[i] = (p[i] - zero[i]) / h + 0.5;
				p_ind[i] = (p_ind[i] % Ncell[i] + Ncell[i]) % Ncell[i];
			}
			return p_ind;
		}
		Eigen::Vector3d deraster(const Eigen::Vector3i& pi) const {
			Eigen::Vector3d p;
			for (int i = 0; i < 3; i++) p[i] = Real(pi[i]) / Ncell[i] * diag[i] + zero[i];
			return p;
		}
	public:
		PeriodicGridIndex(const Eigen::Vector3d& o, const Eigen::Vector3d& d, double eps_ = 0) {
			zero = o; diag = d;
			if (eps_ == 0) {
				h = d.minCoeff() / 1e6;
			} else {
				h = eps_;
			}
			for (int i = 0; i < 3; i++) { Ncell[i] = d[i] / h; }
		}
		Eigen::MatrixX3<Real> dumpPoints(void) const {
			Eigen::Matrix3X<Real> vlist(3, plist.size());
			//for (auto p_id : activeLattice) {
			//	vlist.col(p_id.second) = deraster(p_id.first);
			//}
			for (int i = 0; i < plist.size(); i++) vlist.col(i) = plist[i];
			return vlist.transpose();
		}
		template<typename P>
		Eigen::Vector3i insert(const P& p) {
			Eigen::Vector3i p_ind = raster(p);
			if (!activeLattice.count(p_ind)) {
				plist.push_back(deraster(p_ind));
				for (int i = 0; i < 27; i++) {
					Eigen::Vector3i pnear(p_ind[0] + i % 3 - 1, p_ind[1] + i / 3 % 3 - 1, p_ind[2] + i / 9 - 1);
					periodize(pnear);
					activeLattice[pnear] = insert_counter;
				}
				insert_counter++;
			}
			return p_ind;
		}
		void clear(void) {
			insert_counter = 0;
			activeLattice.clear();
			plist.clear();
		}
		template<typename P>
		int query(const P& p) const {
			auto p_ind = raster(p);
			if (activeLattice.count(p_ind)) { return activeLattice.at(p_ind); }
			return -1;
		}
		
		template<typename P> Eigen::Vector3i i3d(const P& p) const { return raster(p); }
	};

	template<typename P>
	static P make_period(const P& vec, Real period = 2, Real deduce = 0.5) {
		P T01 = vec;
		for (int k = 0; k < 3; k++) { if (T01[k] < -deduce) { T01[k] += period; } else if (T01[k] > deduce) { T01[k] -= period; } }
		return T01;
	}

	template<typename P>
	static P sync_period(const P& fr, const P& ach, Real period = 2, Real deduce = 0.5) {
		return ach + make_period(P(fr - ach), period, deduce);
	}
}

