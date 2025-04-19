#include "BBox.h"
#include "Geometry.h"



bool msf::BBox::intersect(const Geometry& geo)
{
	if (!isValid()) return false;
	return geo.intersect(*this);
}

bool msf::BBox::intersect(const BBox& bb)
{
	if (!isValid() || !bb.isValid()) return false;
	for (int i = 0; i < 3; i++) {
		if (corner[0][i] <= bb.corner[0][i] && corner[1][i] <= bb.corner[0][i]) {
			return true;
		}
		else if (bb.corner[0][i] <= corner[0][i] && bb.corner[1][i] <= corner[0][i]) {
			return true;
		}
	}
	return false;
}

msf::BBox::BBox(const Point& mi, const Point& ma)
{
	corner[0] = mi;
	corner[1] = ma;
}

msf::BBox::BBox(std::vector<Point> pointlist)
{
	corner[0].setConstant(max_real); corner[1].setConstant(min_real);
	for (int i = 0; i < pointlist.size(); i++) {
		corner[0] = corner[0].cwiseMin(pointlist[i]);
		corner[1] = corner[1].cwiseMax(pointlist[i]);
	}
}

msf::BBox msf::BBox::unionBox(const BBox& b2)
{
	BBox b(*this);
	b.corner[0] = b.corner[0].cwiseMin(b2.corner[0]).cwiseMin(b2.corner[1]);
	b.corner[1] = b.corner[1].cwiseMax(b2.corner[0]).cwiseMax(b2.corner[1]);
	return b;
}

msf::BBox msf::BBox::dilate2Cube(void)
{
	auto diag = diagonal();
	Real Mlen = diag.maxCoeff();
	for (int i = 0; i < 3; i++) {
		corner[0][i] -= (Mlen - diag[i]) / 2;
		corner[1][i] += (Mlen - diag[i]) / 2;
	}
	return *this;
}

msf::Point msf::BBox::evaluate(Real u, Real v, Real w)
{
	Point p;
	Point L(u, v, w);
	p = corner[0] + L.cwiseProduct(corner[1] - corner[0]);
	return p;
}

msf::BBox::Vector msf::BBox::diagonal(void) const
{
	Vector diag = corner[1] - corner[0];
	return diag;
}

void msf::BBox::scale(const Point& center, Real s)
{
	corner[0] = s * (corner[0] - center) + center;
	corner[1] = s * (corner[1] - center) + center;
}

void msf::BBox::scale(Real s)
{
	Point center = (corner[0] + corner[1]) / 2;
	scale(center, s);
}

msf::Point msf::BBox::localCoords(Real px, Real py, Real pz)
{
	Point pp;
	pp[0] = (px - corner[0][0]) / (corner[1][0] - corner[0][0]);
	pp[1] = (py - corner[0][1]) / (corner[1][1] - corner[0][1]);
	pp[2] = (pz - corner[0][2]) / (corner[1][2] - corner[0][2]);
	return pp;
}

bool msf::BBox::isOut(const Point& p)
{
	bool out = false;
	for (int k = 0; k < 3; k++) {
		out = out || p[k]<corner[0][k] || p[k]>corner[1][k];
	}
	return out;
}
