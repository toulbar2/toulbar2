/*
 * **************** Read wcsp format files **************************
 * 
 */

#include "tb2system.hpp"
#include "tb2wcsp.hpp"

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
            
    // open the file
    ifstream file(fileName);
    
    // read problem name and sizes
    file >> pbname;
    file >> nbvar;
    file >> nbval;
    file >> nbconstr;
    file >> top;
    if (ToulBar2::verbose) cout << "Read problem: " << pbname << endl;
    
    assert(vars.empty());
    assert(constrs.empty());
    objective->decrease(top - 1);
    NCBucketSize = cost2log2(objective->getSup()) + 1;
    
    // read variable domain sizes
    for (i = 0; i < nbvar; i++) {
        char varname_[128];
        sprintf(varname_, "%d", i);
        string varname(varname_);
        file >> domsize;
        if (ToulBar2::verbose >= 3) cout << "read variable " << i << " of size " << domsize << endl;
        Variable *x = new Variable(varname,0,domsize-1,storeData,true);
        readVars.push_back(x);
        link(x);
    }
    
    // read each constraint
    for (c = 0; c < nbconstr; c++) {
        file >> arity;
        if (arity == 2) {
            file >> i;
            file >> j;
            if (ToulBar2::verbose >= 3) cout << "read binary constraint " << c << " on " << i << "," << j << endl;
        	if (i == j) {
    	       cerr << "Error: binary constraint with only one variable in its scope!" << endl;
               exit(EXIT_FAILURE);
            }
            file >> defval;
            file >> ntuples;
            vector<Cost> costs;
            for (a = 0; a < readVars[i]->getDomainInitSize(); a++) {
                for (b = 0; b < readVars[j]->getDomainInitSize(); b++) {
                    costs.push_back(defval);
                }
            }
        	for (k = 0; k < ntuples; k++) {
                file >> a;
                file >> b;
                file >> cost;
                costs[a * readVars[j]->getDomainInitSize() + b] = cost;
            }
            addBinaryConstraint(readVars[i],readVars[j],costs);
        } else if (arity == 1) {
            file >> i;
            if (ToulBar2::verbose >= 3) cout << "read unary constraint " << c << " on " << i << endl;
            file >> defval;
            file >> ntuples;
            for (a = 0; a < readVars[i]->getDomainInitSize(); a++) {
                vars[i]->project(a, defval);
            }
            for (k = 0; k < ntuples; k++) {
                file >> a;
                file >> cost;
                vars[i]->extend(a, defval);
                vars[i]->project(a, cost);
            }
            vars[i]->findSupport();
            vars[i]->propagateNC();
            propagate();
        } else if (arity == 0) {
            file >> defval;
            file >> ntuples;
            if (ToulBar2::verbose >= 3) cout << "read global lower bound contribution " << c << " of " << defval << endl;
        	if (ntuples != 0) {
                cerr << "Error: global lower bound contribution with several tuples!" << endl;
                exit(EXIT_FAILURE);
            }
            objective->increase(getLb() + defval);
            propagate();
        } else {
            cerr << "Error: not implemented for this constraint arity " << arity << endl;
            exit(EXIT_FAILURE);
        }
    }
    if (ToulBar2::verbose) {
        cout << "Read " << nbvar << " variables, with " << nbval << " values at most, and " << nbconstr << " constraints." << endl;
    }
}
