#pragma once
#include "Geometry.h"
#include <vector>
#include "Eigen/Eigen"
#include "bbox/BBox.h"
#include <memory>
#include "facet/Rect.h"
#include <functional>
//#include "fft/fft_util.h"

BEGIN_MINSURF_NAMESPACE

class Grid;

class Form;

class VectorField;

class Grid {
private:
	Eigen::SparseMatrix<Real> _D;
	int _D_major;
	// 1. position major  2. direction component major
	Eigen::SparseMatrix<Real> differenceMatrix(int major) const;
public:
	bool periodic = true;
	BBox bb;
	Real hx, hy, hz;
	int nx, ny, nz;
	Grid(Real ox, Real oy, Real oz, int nx, int ny, int nz, Real hx, Real hy, Real hz) : nx(nx), ny(ny), nz(nz), hx(hx), hy(hy), hz(hz) {
		bb = BBox(Point(ox, oy, oz), Point(nx * hx, ny * hy, nz * hz));
		_D = differenceMatrix(1);
		_D_major = 1;
	}
	Grid(int nx, int ny, int nz) : nx(nx), ny(ny), nz(nz) {
		Real ox = 0, oy = 0, oz = 0;
		hx = 1 / nx; hy = 1 / ny; hz = 1 / nz;
		bb = BBox(Point(ox, oy, oz), Point(nx * hx, ny * hy, nz * hz));
		_D = differenceMatrix(1);
		_D_major = 1;
	}
	void saveField(std::string filename, const Eigen::Matrix<Real, -1, 1>* pNodeData, int n);
	void saveField(std::string filename, const Eigen::Matrix<Real, -1, 3>& nodedata);
	const Eigen::SparseMatrix<Real>& getDifferenceMatrix(int major) const { if (_D_major == major) return _D; else { return differenceMatrix(major); } };
	Real cellVolume(void) const { return hx * hy * hz; }
	size_t n_nodes(void) const { return nx * ny * nz; }
	bool isUniform(void) const { return std::abs(hx - hy) / std::abs(hx) < 1e-6 && std::abs(hy - hz) / std::abs(hy) < 1e-6; }
	auto createScalarField(Real initVal = 0) {
		return std::make_unique<Form>(*this, 0, initVal);
	}
	auto getRange(void) { return Eigen::Vector3i(nx, ny, nz); }
	Rect getFacet(int axis, int i, int y, int z);
	Point getNode(int x, int y, int z) { return Point(x * hx + bb.getCorner(0)[0], y * hy + bb.getCorner(0)[1], z * hz + bb.getCorner(0)[2]); }
	Real access(const Eigen::Matrix<Real, -1, 1>& valnode, int i, int j, int k) const {
		size_t memloc = ((i + nx) % nx) + ((j + ny) % ny) * nx + ((k + nz) % nz) * nx * ny;
		return valnode[memloc];
	}
	Real& access(Eigen::Matrix<Real, -1, 1>& valnode, int i, int j, int k) const {
		size_t memloc = ((i + nx) % nx) + ((j + ny) % ny) * nx + ((k + nz) % nz) * nx * ny;
		return valnode[memloc];
	}
	void accessRange(Eigen::Vector3i mi, Eigen::Vector3i ma, std::function<void(Eigen::Vector3i)>);

	int fft_freqsize(void) { return nx * ny * (nz / 2 + 1); }
	// the output layout is [Nx][Ny][Nz/2+1]
	template<typename Scalar>
	Eigen::VectorX<std::complex<Scalar>> fft(const Eigen::VectorX<Scalar>& src) const {
		Eigen::VectorX<std::complex<Scalar>> fftres(nx * ny * (nz / 2 + 1));
		std::vector<Scalar> val(src.begin(), src.end());
		//printf("valmin = %f, valmax = %f\n", *std::min_element(val.begin(), val.end()), *std::max_element(val.begin(), val.end()));
		std::vector<std::complex<Scalar>> result;
		cuda_fft(val, result, nx, ny, nz);
		std::copy(result.begin(), result.end(), fftres.data());
		return fftres;
	}
	// only first [Nx][Ny][Nz/2+1] is used
	template<typename Scalar>
	Eigen::VectorX<Scalar> ifft(const Eigen::VectorX<std::complex<Scalar>>& src) const {
		Eigen::VectorX<Scalar> ifftres(n_nodes());
		std::vector<std::complex<Scalar>> val(src.size());
		memcpy(val.data(), src.data(), sizeof(std::complex<Scalar>) * src.size());
		std::vector<Scalar> result;
		cuda_ifft(val, result, nx, ny, nz);
		std::copy(result.begin(), result.end(), ifftres.data());
		return ifftres / (n_nodes());
	}
	Eigen::VectorX<Real> poissonSolve(const Eigen::VectorX<Real>& src) const;
	void checkFFT(const Eigen::VectorX<Real>& src) const;
	Eigen::Vector3i locate(Point p) {
		Real L[3] = { hx * nx, hy * ny, hz * nz };
		for (int i = 0; i < 3; i++) {
			while (p[i] < bb.getCorner(0)[i]) p[i] += L[i];
		}
		Eigen::Vector3i ipos;
		ipos[0] = (p[0] - bb.getCorner(0)[0]) / hx; ipos[0] %= nx;
		ipos[1] = (p[1] - bb.getCorner(0)[1]) / hy; ipos[1] %= ny;
		ipos[2] = (p[2] - bb.getCorner(0)[2]) / hz; ipos[2] %= nz;
		return ipos;
	}
};

class Form {
	int order = 0;
	const Grid& owner;
	Eigen::Matrix<Real, -1, 1> nodeData[3];
public:
	Form(const Grid& g, int ord, Real defaultVal);
	const Eigen::VectorX<Real>* getNodeValue(void) const { return nodeData; }
	Eigen::VectorX<Real>* getNodeValue(void) { return nodeData; }
	Eigen::Matrix<Real, -1, 3> getNodeValueMatrix(void) {
		Eigen::Matrix<Real, -1, 3> mat(owner.n_nodes(), 3);
		for (int i = 0; i < 3; i++) {
			mat.col(i) = nodeData[i];
		}
		return mat;
	}
	Real innerProd(const Form& g2) const {
		if (order != g2.order) throw std::runtime_error("unmatched form order for inner product");
		if (order == 0 || order == 3) {
			return nodeData[0].dot(g2.nodeData[0]) * owner.cellVolume();
		} else if (order < 3) {
			Real s = 0;
			for (int k = 0; k < 3; k++) {
				s += nodeData[k].dot(g2.nodeData[k]);
			}
			return s * owner.cellVolume();
		}
		return 0;
	}
	// Eigen::Matrix<std::complex<Real>, -1, 3> fft(void) const;
	enum DiffOption {
		FORWARD,
		CENTRAL
	};
	Eigen::Vector3<Real> at(int i, int j, int k) const {
		Eigen::Vector3<Real> v;
		v[0]  = owner.access(nodeData[0], i, j, k);
		v[1]  = owner.access(nodeData[1], i, j, k);
		v[2]  = owner.access(nodeData[2], i, j, k);
		return v;
	}
	void set(int i, int j, int k, const Eigen::Vector3<Real>& v) {
		owner.access(nodeData[0], i, j, k) = v[0];
		owner.access(nodeData[1], i, j, k) = v[1];
		owner.access(nodeData[2], i, j, k) = v[2];
	}
	Form d(void) const;
	Real norm(void) const;
	VectorField difference(void) const;
};

class VectorField 
	:public Eigen::Matrix<Real, -1, 3>
{
	const Grid& owner;
	//Eigen::Matrix<Real, -1, 3>& nodeData = E;
	friend class Form;
public:
	VectorField(const Grid& g) : owner(g), Eigen::Matrix<Real, -1, 3>(g.n_nodes(), 3) {
		fill(0);
	}
	Real innerProd(const VectorField& g2) const {
		Real s = 0;
		return cwiseProduct(g2).sum() / owner.cellVolume();
	}
};


END_MINSURF_NAMESPACE
