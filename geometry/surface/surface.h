#pragma once
#include "Geometry.h"
#include <stdint.h>

BEGIN_MINSURF_NAMESPACE

enum class SurfaceType : std::int16_t {
	Mesh,
	Nurbs,
	Manifold
};

class Surface {
	SurfaceType _type;
	
};

END_MINSURF_NAMESPACE
