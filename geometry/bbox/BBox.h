#pragma once

#include "Point/Point.h"


BEGIN_MINSURF_NAMESPACE

class Geometry;

class BBox
{
	Point corner[2];
public:
	using Vector = Point;
	BBox(const Point& mi, const Point& ma);
	BBox(std::vector<Point> pointlist);
	BBox(void) { corner[0].setConstant(max_real); corner[1].setConstant(min_real); }
	bool intersect(const Geometry& geo);
	bool intersect(const BBox& bb);
	BBox unionBox(const BBox& b2);
	Real volume(void) const { return (corner[1] - corner[0]).prod(); }
	bool isValid(void) const { return corner[0][0] < corner[1][1] + eps_real; }
	BBox dilate2Cube(void);
	Point getCenter(void) { return evaluate(0.5, 0.5, 0.5); }
	Point localCoords(Real px, Real py, Real pz);
	void scale(const Point& center, Real s);
	void scale(Real s);
	bool isOut(const Point& p);
	Vector diagonal(void) const;
	Point evaluate(Real u, Real v, Real w);
	const auto& getCorner(int i) const { return corner[i]; }
};

END_MINSURF_NAMESPACE
