/*
 * **************** Read wcsp format files **************************
 * 
 */

#include "tb2system.hpp"
#include "tb2wcsp.hpp"
#include "tb2enumvar.hpp"
#include "tb2pedigree.hpp"
#include "tb2naryconstr.hpp"

typedef struct {
    EnumeratedVariable *var;
    vector<Cost> costs;
} TemporaryUnaryConstraint;



#define MAX_ARITY 50

void WCSP::read_wcsp(const char *fileName)
{
    if (ToulBar2::pedigree) {
      if (!ToulBar2::bayesian) ToulBar2::pedigree->read(fileName, this);
      else ToulBar2::pedigree->read_bayesian(fileName, this);
      return;
    }
    string pbname;
    int nbvar,nbval,nbconstr;
    Cost top;
    int i,j,k,t, ic;
    string varname;
    int domsize;
    unsigned int a;
    unsigned int b;
    unsigned int c;
    Cost defval;
    Cost cost;
    int ntuples;
    int arity;
    string funcname;
    Value funcparam1;
    vector<TemporaryUnaryConstraint> unaryconstrs;
    Cost inclowerbound = 0;
    
    
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

    // read variable domain sizes
    for (i = 0; i < nbvar; i++) {
        string varname;
        varname = to_string(i);
        file >> domsize;
        if (ToulBar2::verbose >= 3) cout << "read variable " << i << " of size " << domsize << endl;
        int theindex = -1;
         
        if (domsize >= 0) theindex = makeEnumeratedVariable(varname,0,domsize-1);
        else theindex = makeIntervalVariable(varname,0,-domsize-1);
        assert(theindex == i);   
    }
    
    // read each constraint
    for (ic = 0; ic < nbconstr; ic++) {
        file >> arity;
        if (!file) {
            cerr << "Warning: EOF reached before reading all the constraints (initial number of constraints too large?)" << endl;
            break;
        }
        if (arity > 3) {
        	EnumeratedVariable* scope[MAX_ARITY];
			for(i=0;i<arity;i++) {
	            file >> j;
	            scope[i] = (EnumeratedVariable*) vars[j];
			}     	
            file >> defval;
		    file >> ntuples;
    
		    if((defval != 0) || (ntuples > 0))           
		    { 
	            NaryConstraint* nary = postNaryConstraint(scope,arity,defval);
                    
	            char buf[MAX_ARITY];
	            for (t = 0; t < ntuples; t++) {
					for(i=0;i<arity;i++) {
			            file >> j;
			            buf[i] = j + CHAR_FIRST;
					}
					buf[i] = '\0';
				    file >> cost;
				
					string tup = buf;
					nary->setTuple(tup, cost, NULL);
	            } 	           
		    }
        } else if (arity == 3) {
            file >> i;
            file >> j;
            file >> k;
        	if ((i == j) || (i == k) || (k == j)) {
    	       cerr << "Error: ternary constraint!" << endl;
               exit(EXIT_FAILURE);
            }
            file >> defval;
            if (defval >= 0) {
                assert(vars[i]->enumerated());
                assert(vars[j]->enumerated());
                assert(vars[k]->enumerated());
                EnumeratedVariable *x = (EnumeratedVariable *) vars[i];
                EnumeratedVariable *y = (EnumeratedVariable *) vars[j];
                EnumeratedVariable *z = (EnumeratedVariable *) vars[k];
                if (ToulBar2::verbose >= 3) cout << "read ternary constraint " << ic << " on " << i << "," << j << "," << k << endl;
                file >> ntuples;
                vector<Cost> costs;
                for (a = 0; a < x->getDomainInitSize(); a++) {
                    for (b = 0; b < y->getDomainInitSize(); b++) {
	                    for (c = 0; c < z->getDomainInitSize(); c++) {
	                        costs.push_back(defval);
						}
                    }
                }
            	for (t = 0; t < ntuples; t++) {
                    file >> a;
                    file >> b;
                    file >> c;
                    file >> cost;
                    costs[a * y->getDomainInitSize() * z->getDomainInitSize() + b * z->getDomainInitSize() + c] = cost;
                }
                if((defval != 0) || (ntuples > 0)) {
                		TernaryConstraint* ctr = postTernaryConstraint(i,j,k,costs); 
                }
            }
		} else if (arity == 2) {
            file >> i;
            file >> j;
        	if (i == j) {
    	       cerr << "Error: binary constraint with only one variable in its scope!" << endl;
               exit(EXIT_FAILURE);
            }
            file >> defval;
            if (defval >= 0) {
                assert(vars[i]->enumerated());
                assert(vars[j]->enumerated());
                EnumeratedVariable *x = (EnumeratedVariable *) vars[i];
                EnumeratedVariable *y = (EnumeratedVariable *) vars[j];
                if (ToulBar2::verbose >= 3) cout << "read binary constraint " << ic << " on " << i << "," << j << endl;
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
                if((defval != 0) || (ntuples > 0)) {
                	BinaryConstraint* ctr =  postBinaryConstraint(i,j,costs);
                }
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
            TemporaryUnaryConstraint unaryconstr;
            unaryconstr.var = x;
            for (a = 0; a < x->getDomainInitSize(); a++) {
                unaryconstr.costs.push_back(defval);
            }
            for (k = 0; k < ntuples; k++) {
                file >> a;
                file >> cost;
                unaryconstr.costs[a] = cost;
            }
            unaryconstrs.push_back(unaryconstr);
            x->queueNC();
        } else if (arity == 0) {
            file >> defval;
            file >> ntuples;
            if (ToulBar2::verbose >= 3) cout << "read global lower bound contribution " << c << " of " << defval << endl;
        	if (ntuples != 0) {
                cerr << "Error: global lower bound contribution with several tuples!" << endl;
                exit(EXIT_FAILURE);
            }
            inclowerbound += defval;
        } 
    }

    sortConstraints();
    // apply basic initial propagation AFTER complete network loading
    increaseLb(getLb() + inclowerbound);
    for (unsigned int u=0; u<unaryconstrs.size(); u++) {
        for (a = 0; a < unaryconstrs[u].var->getDomainInitSize(); a++) {
            if (unaryconstrs[u].costs[a] > 0) unaryconstrs[u].var->project(a, unaryconstrs[u].costs[a]);
        }
        unaryconstrs[u].var->findSupport();
    }
    if (ToulBar2::verbose >= 0) {
        cout << "Read " << nbvar << " variables, with " << nbval << " values at most, and " << nbconstr << " constraints." << endl;
    }
}

