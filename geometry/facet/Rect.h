#pragma once
#include "Facet.h"
#include "Plane.h"

BEGIN_MINSURF_NAMESPACE

class Rect
	: public Facet
{
	Plane pl;
	Real mi[2], ma[2];
public:
	using Vector = Eigen::Vector3<Real>;
	const Plane& getPlane(void) const { return pl; }
	Vector getXdir(void) const { return pl.X(); }
	Vector getYdir(void) const { return pl.Y(); }
	Rect(Plane pl, Eigen::Vector<Real, 2> vmi, Eigen::Vector<Real, 2> vma) : pl(pl) {
		mi[0] = vmi[0]; mi[1] = vmi[1];
		ma[0] = vma[0]; ma[1] = vma[1];
	};
	bool contains(Eigen::Vector2<Real> p) const {
		return p[0] >= mi[0] && p[0] <= ma[0] && p[1] >= mi[1] && p[1] <= ma[1];
	}
	// bool intersect(const BBox& bb) const override;
	BBox getBBox() const override {
		Point p = pl.eval(mi);
		Point q = pl.eval(ma);
		return BBox(p, q);
	}
	
};

END_MINSURF_NAMESPACE
