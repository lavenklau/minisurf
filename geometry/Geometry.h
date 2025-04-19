#pragma once

#include "Config.h"
#include "bbox/BBox.h"

BEGIN_MINSURF_NAMESPACE


class BBox;
class Curve;
class PolyLine;
class Facet;
class Plane;
class Rect;
//class Point;
class Grid;
class Facet;

class Geometry
{
public:
	virtual bool intersect(const BBox& bb) const { return false; };
	virtual BBox getBBox(void) const { return BBox(); };
	virtual operator BBox()  const { return getBBox(); }
	virtual bool intersect(const Rect& rect) { return false; }
};


END_MINSURF_NAMESPACE
