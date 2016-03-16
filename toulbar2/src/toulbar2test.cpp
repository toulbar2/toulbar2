
/**
 * Test toulbar2 API
 */

#include "toulbar2lib.hpp"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// INCOP default command line option
const string Incop_cmd = "0 1 3 idwa 100000 cv v 0 200 1 0 0";

int main(int argc, char * argv[])
{
    mysrand(getpid());

    tb2init(); // must be call before setting specific ToulBar2 options and creating a model

    ToulBar2::verbose = -1;  // change to 0 to see more information

    // uncomment if Virtual Arc Consistency (equivalent to Augmented DAG algorithm) enable
    //	ToulBar2::vac = 1; // option -A
    //	ToulBar2::vacValueHeuristic = true; // option -V
    // uncomment if partial Limited Discrepancy Search enable
    //	ToulBar2::lds = 1;  // option -l=1
    // uncomment if INCOP local search enable
    //	ToulBar2::incop_cmd = Incop_cmd; // option -i

    // create a problem with three 0/1 variables
    WeightedCSPSolver *solver = WeightedCSPSolver::makeWeightedCSPSolver(STORE_SIZE, MAX_COST);
    int x = solver->getWCSP()->makeEnumeratedVariable("x",0,1); // note that for efficiency issue, I assume domain values start at zero (otherwise remove flag -DWCSPFORMATONLY in Makefile)
    int y = solver->getWCSP()->makeEnumeratedVariable("y",0,1);
    int z = solver->getWCSP()->makeEnumeratedVariable("z",0,1);

    // add random unary cost functions on each variable
    {
        vector<Cost> costs(2, 0);
        costs[0] = randomCost(0,100);
        costs[1] = randomCost(0,100);
        solver->getWCSP()->postUnary(x, costs);
        costs[0] = randomCost(0,100);
        costs[1] = randomCost(0,100);
        solver->getWCSP()->postUnary(y, costs);
        costs[0] = randomCost(0,100);
        costs[1] = randomCost(0,100);
        solver->getWCSP()->postUnary(z, costs);
    }

    // add binary cost functions (Ising) on each pair of variables
    {
        vector<Cost> costs;
        for (unsigned int i = 0 ; i < 2; i++) {
            for (unsigned int j = 0 ; j < 2; j++) {
                costs.push_back((solver->getWCSP()->toValue(x, i) == solver->getWCSP()->toValue(y, j)) ? 0 : 30);  // penalizes by a cost=30 if variables are assigned to different values
            }
        }
        solver->getWCSP()->postBinaryConstraint(x,y,costs);
        solver->getWCSP()->postBinaryConstraint(x,z,costs);
        solver->getWCSP()->postBinaryConstraint(y,z,costs);
    }

    // add a ternary hard constraint (x+y=z)
    {
        vector<Cost> costs;
        for (unsigned int i = 0 ; i < 2; i++) {
            for (unsigned int j = 0 ; j < 2; j++) {
                for (unsigned int k = 0 ; k < 2; k++) {
                    costs.push_back((solver->getWCSP()->toValue(x, i) + solver->getWCSP()->toValue(y, j) == solver->getWCSP()->toValue(z, k)) ? 0 : MAX_COST);
                }
            }
        }
        solver->getWCSP()->postTernaryConstraint(x,y,z,costs);
    }

    solver->getWCSP()->sortConstraints(); // to be done before search

    //	int verbose = ToulBar2::verbose;
    //	ToulBar2::verbose = 5;  // high verbosity to see the cost functions
    //	solver->getWCSP()->print(cout);
    //	ToulBar2::verbose = verbose;

    if (solver->solve()) {
        // show optimal solution
        vector<Value> sol;
        Cost optimum = solver->getSolution(sol);
        cout << "Optimum=" << optimum << endl;
        cout << "Solution: x=" << sol[x] << " ,y=" << sol[y] << " ,z=" << sol[z] << endl;
    } else {
        cout << "No solution found!" << endl;
    }
    // cout << "Problem lower bound: " << solver->getWCSP()->getLb() << endl; // initial problem lower bound possibly enhanced by value removals at the root during search

    return 0;
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

