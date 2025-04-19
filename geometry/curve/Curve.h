#pragma once
#include "Config.h"
#include "Geometry.h"
#include "Point/Point.h"
#include <vector>

BEGIN_MINSURF_NAMESPACE

class Curve
	: public Geometry
{
public:
	virtual std::vector<Real> signedIntersection(const std::vector<Point>& pointlist) { return {}; };
};

END_MINSURF_NAMESPACE
