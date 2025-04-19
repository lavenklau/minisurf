#pragma once
#include "Config.h"
#include <tuple>

namespace msf {
	//extern Real STRANG3_weight[4];
	//extern Real STRANG3_uv[4][2];

	//extern Real STRANG7_weight[7];
	//extern Real STRANG7_uv[7][2];

	//extern Real STRANG10_uv[13][2];
	//extern Real STRANG10_w[13];

	//extern Real TOMS584_uv[19][2];
	//extern Real TOMS584_w[19];

	//extern Real TOMS612_19uv[19][2];
	//extern Real TOMS612_19w[19];

	//extern Real TOMS612_28uv[28][2];
	//extern Real TOMS612_28w[28];

	//extern Real TOMS706_uv[37][2];
	//extern Real TOMS706_w[37];

	std::tuple<Real(*)[2], Real*, int> getTriQuadrature(int order);
}

