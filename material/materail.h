#pragma once

#include <tuple>
#include <Eigen/Eigen>

namespace msf {

	namespace mtl {

		inline auto lameCoeff(double E, double nu) { return std::make_tuple(E * nu / (1 + nu) / (1 - 2 * nu), E / 2 / (1 + nu)); }

		inline auto planeLameCoeff(double E, double nu) {
			auto [lam, mu] = lameCoeff(E, nu);
			return std::make_tuple(2 * lam * mu / (lam + 2 * mu), mu);
		}

		inline auto elasticMatrix(double E, double nu) {
			auto [lam, mu] = lameCoeff(E, nu);
			Eigen::Matrix<double, 6, 6> D;
			D <<
				lam + 2 * mu, lam, lam, 0, 0, 0,
				lam, lam + 2 * mu, lam, 0, 0, 0,
				lam, lam, lam + 2 * mu, 0, 0, 0,
				0, 0, 0, mu, 0, 0,
				0, 0, 0, 0, mu, 0,
				0, 0, 0, 0, 0, mu;
			return D;
		}

		inline auto planeElasticMatrix(double E, double nu) {
			auto [lam, mu] = planeLameCoeff(E, nu);
			Eigen::Matrix3d D;
			D <<
				lam + 2 * mu, lam, 0,
				lam, lam + 2 * mu, 0,
				0, 0, mu;
			return D;
		}

		Eigen::Matrix<double, 6, 6> isotropic_projection_s(const Eigen::Matrix<double, 6, 6>& S, double& lam, double& mu);
	}
}
