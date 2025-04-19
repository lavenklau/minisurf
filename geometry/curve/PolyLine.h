#pragma once
#include "Curve.h"

BEGIN_MINSURF_NAMESPACE

class Line : public Curve {
	Point _p[2];
public:
	Line(Point p0, Point p1) { _p[0] = p0; _p[1] = p1; }
	const auto& start(void) const { return _p[0]; }
	const auto& end(void) const { return _p[1]; }
	Eigen::Vector3<Real> vec(void) { return _p[1] - _p[0]; }
	BBox getBBox(void) const override {
		Point mi, ma;
		mi.setConstant(max_real);
		ma.setConstant(min_real);
		mi = mi.cwiseMin(_p[0]).cwiseMin(_p[1]);
		ma = ma.cwiseMax(_p[0]).cwiseMax(_p[1]);
		return BBox(mi, ma);
	}
	Real signedIntersection(const Rect& rect);
};

class PolyLine :
	public Curve
{
	std::vector<Point> pointSeq;
	bool isClosed;
public:
	PolyLine(const std::vector<Point>& pointlist);
	std::vector<Real> signedIntersection(const std::vector<Point>& pointlist);
	Real signedIntersection(const Rect& rect);
	BBox getBBox(void) const override;
	std::vector<Line> explode(void) const;
};


END_MINSURF_NAMESPACE
