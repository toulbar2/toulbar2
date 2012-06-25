
/**
 * Test toulbar2 API
 */

#include "toulbar2.hpp"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char * argv[])
{
	Cost ub = 1000; // MAX_COST; // BUG!
	ToulBar2::verbose = 1;
	ToulBar2::showSolutions = true;

	/// no preprocessing
	ToulBar2::elimDegree = -1;
	ToulBar2::elimDegree_preprocessing = -1;
	ToulBar2::preprocessTernaryRPC = 0;
	ToulBar2::preprocessFunctional  = 0;
	ToulBar2::preprocessNary  = 0;
	ToulBar2::costfuncSeparate = false;

	WeightedCSPSolver *solver = WeightedCSPSolver::makeWeightedCSPSolver(STORE_SIZE, ub);
	int x = solver->getWCSP()->makeEnumeratedVariable("x",0,4);
	int y = solver->getWCSP()->makeEnumeratedVariable("y",0,4);
	int z = solver->getWCSP()->makeEnumeratedVariable("z",0,4);
	int w = solver->getWCSP()->makeEnumeratedVariable("w",0,4);
	solver->getWCSP()->postSupxyc(x,y,1,0);
	solver->getWCSP()->postSupxyc(y,z,1,0);
	solver->getWCSP()->postSupxyc(z,w,1,0);
	int scope[4] = {x,y,z,w};
	string alldiff = "salldiff";
	string parameters_("var " + to_string(ub));
	istringstream parameters(parameters_);
	solver->getWCSP()->postGlobalConstraint(scope, 4, alldiff, parameters);
	solver->solve();
	cout << "end." << endl;
	return 0;
}
