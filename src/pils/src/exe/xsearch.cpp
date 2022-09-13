#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <math.h>
#include <fstream>

#include <core/tb2wcsp.hpp>
#include <search/tb2solver.hpp>
#include <core/tb2binconstr.hpp>

// protect name conflicts
namespace PILS {
#include <base/solution.h>
#include <base/costFunction.h>
#include <init/init.h>
#include <base/fullNeighborEval.h>
#include <base/incrNeighborEval.h>
#include <algo/doubleIncr_biHC.h>
#include <algo/ils.h>
#include <algo/staticPerturb.h>
#include <algo/adaptivePertub.h>
#include <algo/rndPerturb.h>
#include <algo/xsearch.h>
}

using namespace std;
using namespace PILS;

/*
  PILS

*/

Cost Solver::pils(string cmd, vector<Value>& solution)
{
    int nbruns = 3;
    int perturb_id = 0;
    double perturb_s = 0.333;
    Long flatMax = 100;
    Long nEvalHC = 500;
    Long nEvalMax = 10000;
    double strengthMin = 0.1;
    double strengthMax = 0.5;
    double incrFactor = 0.1;
    double decrFactor = 0.1;
    istringstream ss_cmd(cmd);
    if (ss_cmd) ss_cmd >> nbruns;
    if (ss_cmd) ss_cmd >> perturb_id;
    if (ss_cmd) ss_cmd >> perturb_s;
    if (ss_cmd) ss_cmd >> flatMax;
    if (ss_cmd) ss_cmd >> nEvalHC;
    if (ss_cmd) ss_cmd >> nEvalMax;
    if (ss_cmd) ss_cmd >> strengthMin;
    if (ss_cmd) ss_cmd >> strengthMax;
    if (ss_cmd) ss_cmd >> incrFactor;
    if (ss_cmd) ss_cmd >> decrFactor;
    return pils(solution, nbruns, perturb_id, perturb_s, flatMax, nEvalHC, nEvalMax, strengthMin, strengthMax, incrFactor, decrFactor);
}

Cost Solver::pils(vector<Value>& solution, int nbruns, int perturb_id, double perturb_s, unsigned long long flatMax, unsigned long long nEvalHC, unsigned long long nEvalMax, double strengthMin, double strengthMax, double incrFactor, double decrFactor)
{
    if (ToulBar2::verbose >= 1) {
        cout << "PILS parameters: " << nbruns;
        cout << " " <<  perturb_id;
        cout << " " <<  perturb_s;
        cout << " " <<  flatMax;
        cout << " " <<  nEvalHC;
        cout << " " <<  nEvalMax;
        cout << " " <<  strengthMin;
        cout << " " <<  strengthMax;
        cout << " " <<  incrFactor;
        cout << " " <<  decrFactor << endl;
    }

    Cost initialLowerBound = wcsp->getLb();
    Cost initialUpperBound = wcsp->getUb();

    CostFunction eval(wcsp); // read a binary WCSP without its lower bound

    // random initialization
    Init init(myrandom_generator, eval);
    IncrNeighborEval neighborEval(eval);

    // local search
    DoubleIncr_biHC hc(myrandom_generator, eval, neighborEval, nEvalHC, flatMax, std::cout);


    // perturbation
    Perturbation * perturbation;

    switch (perturb_id) {
    case 0: // static perturbation
        unsigned strength;

        if (perturb_s < 1) // proportion of the total length
            strength = round(perturb_s * eval.n_variables);
        else
            strength = (unsigned) perturb_s;
        perturbation = new StaticPerturb(myrandom_generator, eval, neighborEval, strength);
        break;

    case 1: // random perturbation
        unsigned strengthmin;
        if (strengthMin < 1) // proportion of the total length
            strengthmin = round(strengthMin * eval.n_variables);
        else
            strengthmin = (unsigned) strengthMin;

        unsigned strengthmax;
        if (strengthMax < 1) // proportion of the total length
            strengthmax = round(strengthMax * eval.n_variables);
        else
            strengthmax = (unsigned) strengthMax;

        perturbation = new RandomPerturb(myrandom_generator, eval, neighborEval, strengthmin, strengthmax);
        break;

    case 2: // adaptive perturbation
        perturbation = new AdaptivePerturb(myrandom_generator, eval, neighborEval, strengthMin, strengthMax, incrFactor, decrFactor);
        break;

    default: // wrong parameter value
        cerr << "Sorry! Wrong perturb mode (should be 0, 1 or 2): " << perturb_id << endl;
        throw BadConfiguration();
    }


    // CrossOVer algorithm search initialization
    Xsearch xs(eval, neighborEval, hc, *perturbation, nEvalMax, flatMax);


    //------------------------------------------------------
    // computation and output

    Solution x;

    for (int nessai = 0; nessai < nbruns && wcsp->getLb() < wcsp->getUb(); nessai++) {
        // random initialization except first run
        init(x);
        if (nessai==0) {
            for(size_t i = 0; i < eval.size(); i++) {
                int val = eval.getPILSValueIndex(i, solution[eval.getWCSPIndex(i)]);
                if (val >= 0 && val < (int)eval.variable_size(i)) {
                    x[i] = val;
                }
            }
        }
        eval(x);
        xs.nEval = 1;
        if (ToulBar2::verbose >= 1) {
            x.print();
            cout << endl;
        }

        // do it
        xs(x);

        // verify result
        eval(x);

        if (initialLowerBound + x.fitness()  < wcsp->getUb()) {
            vector<Value> bestsolution(wcsp->numberOfVariables(), 0);
            int depth = Store::getDepth();
            try {
                Store::store();
                vector<int> tabvars(eval.size());
                vector<Value> solution(eval.size());
                for (unsigned i = 0; i < eval.size(); i++) {
                    tabvars[i] = eval.getWCSPIndex(i);
                    solution[i] = eval.getValue(i, x[i]);
                }
                wcsp->assignLS(tabvars, solution);
                newSolution();
                for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
                    bestsolution[i] = wcsp->getValue(i);
                    ((WCSP *)wcsp)->setBestValue(i, bestsolution[i]);
                }
            } catch (const Contradiction&) {
                wcsp->whenContradiction();
            }
            Store::restore(depth);
        }
    }

    delete perturbation;

    if (wcsp->getUb() < initialUpperBound) {
        wcsp->enforceUb();
        wcsp->propagate();
    }

    return (wcsp->getUb() < initialUpperBound)?wcsp->getSolutionCost():MAX_COST;
}
