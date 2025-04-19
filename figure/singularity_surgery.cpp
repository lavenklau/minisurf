#include "mesher/asymptotic_analysis.h"

using namespace msf;



void test_surgery(void) {
	//std::string mfile = "D:/projects/minisurf/image/siggraph/surgery/befsur.obj";
	std::string mfile = "D:/projects/minisurf/image/siggraph/surgery/test4.obj";
	PeriodSurfaceMesh mesher;
	mesher.readMergePeriodBoundary(mfile);
	mesher.delaunayRemesh(20, 0.2, 0.02, 1);
	mesher.saveUnitCell(getPath("surinput.obj"));

	mesher.surgery();
	mesher.saveUnitCell(getPath("suroutput.obj"));
	mesher.delaunayRemesh(5, 0.2, 0.02, 1);
	mesher.saveUnitCell(getPath("surremsh.obj"));
}



