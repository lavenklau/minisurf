#pragma once

#include <vector>

namespace qd {

	namespace tri {
		struct QuadraturePoint {
			double c[3];
			double w;
		};
		namespace symq {
			std::vector<QuadraturePoint> get(int precision);
			int getSize(int precision);
		}
	}


}
