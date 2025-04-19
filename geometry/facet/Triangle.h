#pragma once

#include "Geometry.h"
#include "Facet.h"

BEGIN_MINSURF_NAMESPACE

class Triangle
	: public Facet
{
	Point v[3];
public:
	BBox getBBox() const override;
};

END_MINSURF_NAMESPACE
