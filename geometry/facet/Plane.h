#pragma once

#include "Config.h"
#include "Point/Point.h"
#include "Geometry.h"

BEGIN_MINSURF_NAMESPACE

class Plane 
{
public:
	using Vector = Point;
private:
	Point origin;
	Vector x,y;
public:
	static Plane XYPlane(const Point& p) {
		return Plane(p, Vector(1, 0, 0), Vector(0, 1, 0));
	}
	static Plane YZPlane(const Point& p) {
		return Plane(p, Vector(0, 1, 0), Vector(0, 0, 1));
	}
	static Plane ZXPlane(const Point& p) {
		return Plane(p, Vector(0, 0, 1), Vector(1, 0, 0));
	}
	static Plane AxisPlane(int axis, const Point& p) {
		switch (axis) {
		case 0:
			return YZPlane(p);
		case 1:
			return ZXPlane(p);
		case 2:
			return XYPlane(p);
		default:
			return XYPlane(p);
		}
	}
	const auto& Origin(void) const { return origin; }
	auto X(void) const { return x; }
	auto Y(void) const { return y; }
	Plane(const Point& p, const Vector& x, const Vector& y) : origin(p), x(x), y(y) {}
	Point eval(const Real coord[2]) const {
		Point p = coord[0] * x + coord[1] * y + origin;
		return p;
	}
};

END_MINSURF_NAMESPACE
