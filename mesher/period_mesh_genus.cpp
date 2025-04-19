#include "PeriodicMesher.h"

using namespace msf;

void test_period_mesh_genus(std::string mfile) {
	PeriodSurfaceMesh mesher;
	mesher.readMergePeriodBoundary(mfile);
	int genus = 1 - (mesher.n_vertices() + mesher.n_faces() - mesher.n_edges()) / 2;
	std::cout << "genus = " << genus << std::endl;
}