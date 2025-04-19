#include "Triangle.h"

msf::BBox msf::Triangle::getBBox() const
{
	return BBox({ v[0],v[1],v[2] });
}

