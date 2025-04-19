#pragma once
#include "Config.h"
#include "Eigen/Eigen"

BEGIN_MINSURF_NAMESPACE

//class Point
//	: public Eigen::Matrix<Real, 3, 1>
//{
//	using Base = Eigen::Matrix<Real, 3, 1>;
//public:
//	template<typename... Args>
//	Point(Args&&... args) : Eigen::Vector3d(std::forward<Args>(args)...) { }
//
//};
using Point = Eigen::Vector3<Real>;

END_MINSURF_NAMESPACE
