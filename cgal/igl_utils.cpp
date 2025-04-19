#include "igl_utils.h"
#include "igl/readOBJ.h"
#include "igl/readOFF.h"
#include "igl/readSTL.h"
#include "igl/readPLY.h"
#include "igl/read_triangle_mesh.h"

using namespace msf;

std::tuple<Eigen::MatrixX3<msf::Real>, Eigen::MatrixX3i> readMesh(const std::string& filename)
{
	Eigen::MatrixX3<Real> V;
	Eigen::MatrixX3i F;
	igl::read_triangle_mesh(filename, V, F);
}

