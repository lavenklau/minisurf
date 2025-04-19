#pragma once

#include "Config.h"
#include <Eigen/Eigen>

BEGIN_MINSURF_NAMESPACE

//template<typename Scalar>
//struct Quaternion 
//	: public Eigen::Quaternion<Scalar>
//{
//	using Quat = Eigen::Quaternion<Scalar>;
//	Eigen::Matrix4<Scalar> toMatrix(void) const {
//		Eigen::Matrix4<Scalar> m;
//		const auto a = w(), b = x(), c = y(), d = z();
//		m << a, -b, -c, -d, b, a, -d, c, c, d, a, -b, d, -c, b, a;
//		return m;
//	}
//	Quaternion(const Eigen::Vector3<Scalar>& v) {
//		Quat::w() = 0; Quat::x() = v[0]; Quat::y() = v[1]; Quat::z() = v[2];
//	}
//	Quaternion(void) = default;
//	Quaternion(const Scalar wval) {
//		Quat::w() = wval; Quat::x() = 0; Quat::y() = 0; Quat::z() = 0;
//	}
//	Scalar real(void) const { return w(); }
//	Scalar& real(void) { return w(); }
//	Eigen::Vector3<Scalar>& imag(void) const { return Quat::vec(); }
//	Eigen::Vector3<Scalar>& imag(void) { return Quat::vec(); }
//};

template<typename Scalar>
Eigen::Matrix4<Scalar> q2mat(const Eigen::Quaternion<Scalar>& q) {
	Eigen::Matrix4<Scalar> m;
	const auto a = q.w(), b = q.x(), c = q.y(), d = q.z();
	m << 
		a, -b, -c, -d,
		b, a, -d, c, 
		c, d, a, -b,
		d, -c, b, a;
	return m;
}

template<typename Scalar>
Eigen::Quaternion<Scalar> vec2q(const Eigen::Vector4<Scalar>& x) {
	Eigen::Quaternion<Scalar> q;
	q.w() = x[0]; q.vec() = x.template tail<3>();
	return q;
}

template<typename Scalar>
Eigen::Vector4<Scalar> q2vec(const Eigen::Quaternion<Scalar>& x) {
	Eigen::Vector4<Scalar> q;
	q[0] = x.w(); q[1] = x.x(); q[2] = x.y(); q[3] = x.z();
	return q;
}

template<typename Scalar>
Eigen::Quaternion<Scalar> operator*(Scalar s, const Eigen::Quaternion<Scalar>& q) {
	Eigen::Quaternion<Scalar> sq = q; sq.w() *= s; sq.vec() *= s;
	return sq;
}

template<typename Scalar>
Eigen::Quaternion<Scalar> operator-(const Eigen::Quaternion<Scalar>& p, const Eigen::Quaternion<Scalar>& q) {
	Eigen::Quaternion<Scalar> pq;
	pq.w() = p.w() - q.w(); pq.vec() = p.vec() - q.vec();
	return pq;
}

template<typename Scalar>
Eigen::Quaternion<Scalar> operator+(const Eigen::Quaternion<Scalar>& p, const Eigen::Quaternion<Scalar>& q) {
	Eigen::Quaternion<Scalar> pq;
	pq.w() = q.w() + p.w(); pq.vec() = q.vec() + p.vec();
	return pq;
}


END_MINSURF_NAMESPACE
 
