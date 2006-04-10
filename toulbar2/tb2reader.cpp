/*
 * **************** Read wcsp format files **************************
 * 
 */

#include "tb2system.hpp"
#include "tb2wcsp.hpp"
#include "tb2enumvar.hpp"

void WCSP::read_wcsp(const char *fileName)
{
    string pbname;
    int nbvar,nbval,nbconstr;
    Cost top;
    int i,j,c,k;
    string varname;
    int domsize;
    unsigned int a;
    unsigned int b;
    Cost defval;
    Cost cost;
    int ntuples;
    int arity;
    string funcname;
    Value funcparam1;
            
    // open the file
    ifstream file(fileName);
    if (!file) {
        cerr << "Could not open file " << fileName << endl;
        exit(EXIT_FAILURE);
    }
    
    // read problem name and sizes
    file >> pbname;
    file >> nbvar;
    file >> nbval;
    file >> nbconstr;
    file >> top;
    if (ToulBar2::verbose) cout << "Read problem: " << pbname << endl;
    
    assert(vars.empty());
    assert(constrs.empty());
    updateUb(top);
    NCBucketSize = cost2log2(getUb()) + 1;
    
    // read variable domain sizes
    for (i = 0; i < nbvar; i++) {
        char varname_[128];
        sprintf(varname_, "%d", i);
        string varname(varname_);
        file >> domsize;
        if (ToulBar2::verbose >= 3) cout << "read variable " << i << " of size " << domsize << endl;
        int theindex = -1;
        if (domsize >= 0) theindex = makeEnumeratedVariable(varname,0,domsize-1);
        else theindex = makeIntervalVariable(varname,0,-domsize-1);
        assert(theindex == i);
    }
    
    // read each constraint
    for (c = 0; c < nbconstr; c++) {
        file >> arity;
        if (arity == 2) {
            file >> i;
            file >> j;
        	if (i == j) {
    	       cerr << "Error: binary constraint with only one variable in its scope!" << endl;
               exit(EXIT_FAILURE);
            }
            assert(vars[i]->enumerated());
            assert(vars[j]->enumerated());
            EnumeratedVariable *x = (EnumeratedVariable *) vars[i];
            EnumeratedVariable *y = (EnumeratedVariable *) vars[j];
            file >> defval;
            if (defval >= 0) {
                if (ToulBar2::verbose >= 3) cout << "read binary constraint " << c << " on " << i << "," << j << endl;
                file >> ntuples;
                vector<Cost> costs;
                for (a = 0; a < x->getDomainInitSize(); a++) {
                    for (b = 0; b < y->getDomainInitSize(); b++) {
                        costs.push_back(defval);
                    }
                }
            	for (k = 0; k < ntuples; k++) {
                    file >> a;
                    file >> b;
                    file >> cost;
                    costs[a * y->getDomainInitSize() + b] = cost;
                }
                postBinaryConstraint(i,j,costs);
            } else {
                file >> funcname;
                if (funcname == ">=") {
                    file >> funcparam1;
                    postSupxyc(i,j,funcparam1);
                } else if (funcname == ">") {
                    file >> funcparam1;
                    postSupxyc(i,j,funcparam1 + 1);
                } else if (funcname == "<=") {
                    file >> funcparam1;
                    postSupxyc(j,i, -funcparam1);
                } else if (funcname == "<") {
                    file >> funcparam1;
                    postSupxyc(j,i, -funcparam1 + 1);
                } else if (funcname == "=") {
                    file >> funcparam1;
                    postSupxyc(i,j,funcparam1);
                    postSupxyc(j,i,-funcparam1);
                } else {
                    cerr << "Error: function " << funcname << " not implemented!" << endl;
                    exit(EXIT_FAILURE);
                }
            }
        } else if (arity == 1) {
            file >> i;
            if (ToulBar2::verbose >= 3) cout << "read unary constraint " << c << " on " << i << endl;
            assert(vars[i]->enumerated());
            EnumeratedVariable *x = (EnumeratedVariable *) vars[i];
            file >> defval;
            file >> ntuples;
            for (a = 0; a < x->getDomainInitSize(); a++) {
                x->project(a, defval);
            }
            for (k = 0; k < ntuples; k++) {
                file >> a;
                file >> cost;
                x->extend(a, defval);
                x->project(a, cost);
            }
            x->findSupport();
//            x->propagateNC();       // Let the initial propagation be done only once in solver.cpp
//            propagate();
            x->queueNC();
        } else if (arity == 0) {
            file >> defval;
            file >> ntuples;
            if (ToulBar2::verbose >= 3) cout << "read global lower bound contribution " << c << " of " << defval << endl;
        	if (ntuples != 0) {
                cerr << "Error: global lower bound contribution with several tuples!" << endl;
                exit(EXIT_FAILURE);
            }
            increaseLb(getLb() + defval);
//            propagate();       // Let the initial propagation be done only once in solver.cpp
        } else {
            cerr << "Error: not implemented for this constraint arity " << arity << "!" << endl;
            exit(EXIT_FAILURE);
        }
    }
    sortConstraints();
    if (ToulBar2::verbose >= 0) {
        cout << "Read " << nbvar << " variables, with " << nbval << " values at most, and " << nbconstr << " constraints." << endl;
    }
}
