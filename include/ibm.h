#ifndef IBM_H
	#define IBM_H

	#include "mesh.h"
	#include "fluid.h"

	using namespace std;

	void interpolation(fluid fluido, mesh membrana, int x, int y, int z);
	void spread(fluid fluido, mesh membrana, int x, int y, int z);
	double dirac_2(double *x);
	double dirac_3(double *x);
	double dirac_4(double *x);

#endif
