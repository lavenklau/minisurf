#include "Config.h"
#include <igl/copyleft/marching_cubes.h>
#include "tbb/tbb.h"
#include "Eigen/Eigen"


using namespace msf;

std::pair<Eigen::MatrixX3<Real>, Eigen::MatrixX3i> isosurf_mesh(
	int reso, const std::pair<Eigen::Vector3<Real>, Eigen::Vector3<Real>>& bbox,
	std::function<Real(Real, Real, Real)> sdf, const double isolevel)
{
	size_t gs = reso + 1;

	const size_t nsamples = gs * gs * gs;

	Eigen::MatrixXd sampleLocations(nsamples, 3);
	{
		size_t i = 0;
		for (size_t zi = 0; zi < gs; ++zi) {
			for (size_t yi = 0; yi < gs; ++yi) {
				for (size_t xi = 0; xi < gs; ++xi) {
					Eigen::Vector3<Real> p =
						bbox.first + (bbox.second - bbox.first).cwiseProduct(Eigen::Vector3<Real>(xi, yi, zi) / reso);
					sampleLocations.row(i) = p.transpose();
					++i;
				}
			}
		}
	}

	// Evaluate signed distances at each grid point
	Eigen::VectorXd signedDistances(nsamples);

	tbb::parallel_for(tbb::blocked_range<size_t>(0, nsamples),
		[&](const tbb::blocked_range<size_t>& r) {
			for (size_t i = r.begin(); i < r.end(); ++i) {
				const auto& p = sampleLocations.row(i);
				signedDistances(i) = sdf(p[0], p[1], p[2]) - isolevel;
			}
		}
	);

	Eigen::MatrixX3<Real> V;
	Eigen::MatrixX3i F;
	igl::copyleft::marching_cubes(signedDistances, sampleLocations,
		gs, gs, gs, V, F);

	return { V,F };
}

std::pair<Eigen::MatrixX3<msf::Real>, Eigen::MatrixX3i> isosurf_mesh(
	int reso, const Eigen::MatrixX3<msf::Real>& sampleLocations,
	const Eigen::VectorXd& signedDistances)
{
	size_t gs = reso + 1;
	const size_t nsamples = gs * gs * gs;

	Eigen::MatrixX3<Real> V;
	Eigen::MatrixX3i F;
	igl::copyleft::marching_cubes(signedDistances, sampleLocations,
		gs, gs, gs, V, F);

	return { V,F };
}
