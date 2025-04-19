#include "PolyLine.h"
#include "facet/Rect.h"

using namespace msf;

msf::PolyLine::PolyLine(const std::vector<msf::Point>& pointlist)
	: pointSeq(pointlist)
{
	if (*pointSeq.rbegin() == *pointSeq.begin()) {
		isClosed = true;
	}
}

std::vector<Real> PolyLine::signedIntersection(const std::vector<Point>& pointlist)
{
	
	return std::vector<Real>();
}

msf::Real PolyLine::signedIntersection(const Rect& rect)
{
	auto segs = explode();
	Real s = 0;
	for (int i = 0; i < segs.size(); i++) {
		s += segs[i].signedIntersection(rect);
	}
	return s;
}

msf::BBox PolyLine::getBBox(void) const
{
	Point mi(max_real, max_real, max_real), ma(min_real, min_real, min_real);
	for (auto p : pointSeq) {
		mi = mi.cwiseMin(p);
		ma = ma.cwiseMax(p);
	}
	return BBox(mi, ma);
}

std::vector<msf::Line> PolyLine::explode(void) const
{
	std::vector<Line> linelist;
	for (int i = 0; i < pointSeq.size() - 1; i++) {
		linelist.emplace_back(pointSeq[i], pointSeq[i + 1]);
	}
	return linelist;
}

msf::Real Line::signedIntersection(const Rect& rect)
{
	Eigen::Vector3<Real> param;
	Eigen::Matrix3<Real> A;
	A << rect.getXdir(), rect.getYdir(), _p[1] - _p[0];
	Eigen::Vector3<Real> b = _p[1] - rect.getPlane().Origin();
	Real Adet = A.determinant();
	if (std::abs(Adet) < 1e-30) return 0;
	param = A.lu().solve(b);
	if (rect.contains(param.topRows(2)) && param[2] >= 0 && param[2] <= 1) {
		return Adet < 0 ? -1 : 1;
	}
	else {
		return 0;
	}
}
