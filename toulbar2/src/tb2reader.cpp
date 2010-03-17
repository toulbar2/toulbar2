/*
 * **************** Read wcsp format files **************************
 * 
 */

#include "tb2wcsp.hpp"
#include "tb2enumvar.hpp"
#include "tb2pedigree.hpp"
#include "tb2bep.hpp"
#include "tb2naryconstr.hpp"
#include "tb2randomgen.hpp"
#include <list>


typedef struct {
    EnumeratedVariable *var;
    vector<Cost> costs;
} TemporaryUnaryConstraint;


void WCSP::read_wcsp(const char *fileName)
{
    if (ToulBar2::pedigree) {
      if (!ToulBar2::bayesian) ToulBar2::pedigree->read(fileName, this);
      else ToulBar2::pedigree->read_bayesian(fileName, this);
      return;
    }
    else if (ToulBar2::uai) {
	    read_uai2008(fileName);
	    return;
    }
    else if (ToulBar2::xmlflag) {
	    read_XML(fileName);
	    return;
    }
    else if (ToulBar2::bep) {
	  ToulBar2::bep->read(fileName, this);
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
	Value funcparam2;
    vector<TemporaryUnaryConstraint> unaryconstrs;
    Cost inclowerbound = MIN_COST;

   
    
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
    
	Cost K = ToulBar2::costMultiplier;    
	if(top < MAX_COST / K)	top = top * K;
	else top = MAX_COST;
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
        if(arity > MAX_ARITY)  { cerr << "Nary constraints of arity > " << MAX_ARITY << " not supported" << endl; exit(EXIT_FAILURE); }       
        if (!file) {
            cerr << "Warning: EOF reached before reading all the constraints (initial number of constraints too large?)" << endl;
            break;
        }
        if (arity > 3) {
		    if (ToulBar2::verbose >= 3) cout << "read " << arity << "-ary constraint " << ic << " on";
        	int scopeIndex[MAX_ARITY];
			for(i=0;i<arity;i++) {
	            file >> j;
				if (ToulBar2::verbose >= 3) cout << " " << j;
	            scopeIndex[i] = j;
			}
			if (ToulBar2::verbose >= 3) cout << endl;
            file >> defval;
		    file >> ntuples;
    
		    if((defval != MIN_COST) || (ntuples > 0))           
		    { 
			    Cost tmpcost = defval*K;
				if(CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST*getUb())) tmpcost *= MEDIUM_COST;
	            int naryIndex = postNaryConstraint(scopeIndex,arity,tmpcost);
                NaryConstraint *nary = (NaryConstraint *) constrs[naryIndex];

	            Char buf[MAX_ARITY];
	            for (t = 0; t < ntuples; t++) {
					for(i=0;i<arity;i++) {
			            file >> j;			            
			            buf[i] = j + CHAR_FIRST;
					}
					buf[i] = '\0';
				    file >> cost;
				    Cost tmpcost = cost * K;
					if(CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST*getUb())) tmpcost *= MEDIUM_COST;
					String tup = buf;
					nary->setTuple(tup, tmpcost, NULL);
	            }
				
				Cost minc = nary->getMinCost();
				if (minc > MIN_COST) {
				  Cost defcost = nary->getDefCost();
				  if (CUT(defcost, minc)) nary->setDefCost(defcost - minc);
				  String tuple;
				  Cost cost;
				  nary->first();
				  while (nary->next(tuple,cost)) {
					nary->setTuple(tuple, cost-minc, NULL);
				  }
				  if (ToulBar2::verbose >= 2) cout << "IC0 performed for constraint " << nary << " with initial minimum cost " << minc << endl;
				  inclowerbound += minc;
				}

	            //((NaryConstraintMap*) nary)->changeDefCost( top );
	            //((NaryConstraintMap*) nary)->preprojectall2();
				nary->propagate();

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
            if (defval >= MIN_COST) {
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
						    Cost tmpcost = defval*K;
							if(CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST*getUb())) tmpcost *= MEDIUM_COST;
	                        costs.push_back(tmpcost);
						}
                    }
                }
            	for (t = 0; t < ntuples; t++) {
                    file >> a;
                    file >> b;
                    file >> c;
                    file >> cost;
					Cost tmpcost = cost*K;
					if(CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST*getUb())) tmpcost *= MEDIUM_COST;
                    costs[a * y->getDomainInitSize() * z->getDomainInitSize() + b * z->getDomainInitSize() + c] = tmpcost;                    
                }
                if(ToulBar2::vac) {
	                for (a = 0; a < x->getDomainInitSize(); a++) {
	                    for (b = 0; b < y->getDomainInitSize(); b++) {
		                    for (c = 0; c < z->getDomainInitSize(); c++) {
		             			Cost co = costs[a * y->getDomainInitSize() * z->getDomainInitSize() + b * z->getDomainInitSize() + c];
		                        histogram(co);
		                    }
	                    }
	                }               	
                }                
                if((defval != MIN_COST) || (ntuples > 0)) postTernaryConstraint(i,j,k,costs);
            }
		} else if (arity == 2) {
            file >> i;
            file >> j;
			if (ToulBar2::verbose >= 3) cout << "read binary constraint " << ic << " on " << i << "," << j << endl;
        	if (i == j) {
    	       cerr << "Error: binary constraint with only one variable in its scope!" << endl;
               exit(EXIT_FAILURE);
            }
            file >> defval;
            if (defval >= MIN_COST) {
                assert(vars[i]->enumerated());
                assert(vars[j]->enumerated());
                EnumeratedVariable *x = (EnumeratedVariable *) vars[i];
                EnumeratedVariable *y = (EnumeratedVariable *) vars[j];
                file >> ntuples;
                vector<Cost> costs;
                for (a = 0; a < x->getDomainInitSize(); a++) {
                    for (b = 0; b < y->getDomainInitSize(); b++) {
					    Cost tmpcost = defval*K;
						if(CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST*getUb())) tmpcost *= MEDIUM_COST;
                        costs.push_back(tmpcost);
                    }
                }
            	for (k = 0; k < ntuples; k++) {
                    file >> a;
                    file >> b;
                    file >> cost;
					Cost tmpcost = cost*K;
					if(CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST*getUb())) tmpcost *= MEDIUM_COST;
                    costs[a * y->getDomainInitSize() + b] = tmpcost;
                }
                
                if(ToulBar2::vac) {
	                for (a = 0; a < x->getDomainInitSize(); a++) {
	                    for (b = 0; b < y->getDomainInitSize(); b++) {
	             			Cost c = costs[a * y->getDomainInitSize() + b];
	                        histogram(c);
	                    }
	                }               	
                }
                if((defval != MIN_COST) || (ntuples > 0)) postBinaryConstraint(i,j,costs);
            } else {
                file >> funcname;
                if (funcname == ">=") {
                    file >> funcparam1;
                    file >> funcparam2;
                    postSupxyc(i,j,funcparam1,funcparam2);
                } else if (funcname == ">") {
                    file >> funcparam1;
                    file >> funcparam2;
                    postSupxyc(i,j,funcparam1 + 1,funcparam2);
                } else if (funcname == "<=") {
                    file >> funcparam1;
                    file >> funcparam2;
                    postSupxyc(j,i, -funcparam1,funcparam2);
                } else if (funcname == "<") {
                    file >> funcparam1;
                    file >> funcparam2;
                    postSupxyc(j,i, -funcparam1 + 1,funcparam2);
                } else if (funcname == "=") {
                    file >> funcparam1;
                    file >> funcparam2;
                    postSupxyc(i,j,funcparam1,funcparam2);
                    postSupxyc(j,i,-funcparam1,funcparam2);
                } else if (funcname == "disj") {
				  Cost funcparam3;
				  file >> funcparam1;
				  file >> funcparam2;
				  file >> funcparam3;
				  postDisjunction(i,j,funcparam1,funcparam2,funcparam3);
                } else if (funcname == "sdisj") {
				  Value funcparam3;
				  Value funcparam4;
				  Cost funcparam5;
				  Cost funcparam6;
				  file >> funcparam1;
				  file >> funcparam2;
				  file >> funcparam3;
				  file >> funcparam4;
				  file >> funcparam5;
				  file >> funcparam6;
				  postSpecialDisjunction(i,j,funcparam1,funcparam2,funcparam3,funcparam4,funcparam5,funcparam6);
                } else {
                    cerr << "Error: function " << funcname << " not implemented!" << endl;
                    exit(EXIT_FAILURE);
                }
            }
        } else if (arity == 1) {
            file >> i;
            if (ToulBar2::verbose >= 3) cout << "read unary constraint " << ic << " on " << i << endl;
			if (vars[i]->enumerated()) {
			  EnumeratedVariable *x = (EnumeratedVariable *) vars[i];
			  file >> defval;
			  file >> ntuples;
			  TemporaryUnaryConstraint unaryconstr;
			  unaryconstr.var = x;
			  for (a = 0; a < x->getDomainInitSize(); a++) {
				Cost tmpcost = defval*K;
     			if(CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST*getUb())) tmpcost *= MEDIUM_COST;
                unaryconstr.costs.push_back(tmpcost);
			  }
			  for (k = 0; k < ntuples; k++) {
                file >> a;
                file >> cost;
				Cost tmpcost = cost*K;
 	 	    	if(CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST*getUb())) tmpcost *= MEDIUM_COST;
                unaryconstr.costs[a] = tmpcost;
			  }
			  if(ToulBar2::vac) {
                for (a = 0; a < x->getDomainInitSize(); a++) {
				  Cost c = unaryconstr.costs[a];
				  histogram(c);
                }               	
			  }
			  unaryconstrs.push_back(unaryconstr);
			  x->queueNC();
			} else {
			  file >> defval;
			  if (defval == MIN_COST) {
				cerr << "Error: unary cost function with zero penalty cost!" << endl;
				exit(EXIT_FAILURE);
			  }
			  file >> ntuples;
			  Value *dom = new Value[ntuples];
			  for (k = 0; k < ntuples; k++) {
                file >> dom[k];
                file >> cost;
				if (cost != MIN_COST) {
				  cerr << "Error: unary cost function with non-zero cost tuple!" << endl;
				  exit(EXIT_FAILURE);
				}
			  }			  
			  postUnary(i,dom,ntuples,defval);
			  delete [] dom;
			}
        } else if (arity == 0) {
            file >> defval;
            file >> ntuples;
            if (ToulBar2::verbose >= 3) cout << "read global lower bound contribution " << ic << " of " << defval << endl;
        	if (ntuples > 1) {
                cerr << "Error: global lower bound contribution with several tuples!" << endl;
                exit(EXIT_FAILURE);
            }
            if (ntuples == 1) file >> cost;
              else cost = defval;
            inclowerbound += cost*K;
        } 
    }
    
	file >> funcname;
	if (file) {
	  cerr << "Warning: EOF not reached after reading all the constraints (initial number of constraints too small?)" << endl;
	}
    sortVariables();
    sortConstraints();
    // apply basic initial propagation AFTER complete network loading
    increaseLb(inclowerbound);
    
    for (unsigned int u=0; u<unaryconstrs.size(); u++) {
        for (a = 0; a < unaryconstrs[u].var->getDomainInitSize(); a++) {
            if (unaryconstrs[u].costs[a] > MIN_COST) unaryconstrs[u].var->project(a, unaryconstrs[u].costs[a]);
        }
        unaryconstrs[u].var->findSupport();
    }
    if (ToulBar2::verbose >= 0) {
        cout << "Read " << nbvar << " variables, with " << nbval << " values at most, and " << nbconstr << " constraints." << endl;
    }   
    histogram();
}


void WCSP::read_random(int n, int m, vector<int>& p, int seed, bool forceSubModular ) 
{
	naryRandom randwcsp(this,seed);
    randwcsp.Input(n,m,p,forceSubModular);	    
 
 	unsigned int nbconstr = numberOfConstraints();
    sortVariables();
    sortConstraints();
    
    if (ToulBar2::verbose >= 0) {
        cout << "Generated random problem " << n << " variables, with " << m << " values, and " << nbconstr << " constraints." << endl;
    }  
}




void WCSP::read_uai2008(const char *fileName)
{
    ToulBar2::NormFactor = (-Log( (TProb)10.)/Log1p( - Pow( (TProb)10., -(TProb)ToulBar2::resolution)));
    if (ToulBar2::NormFactor > (Pow( (TProb)2., (TProb)INTEGERBITS)-1)/-Log10(Pow( (TProb)10., -(TProb)ToulBar2::resolution))) {
	   cerr << "This resolution cannot be ensured on the data type used to represent costs." << endl;
	   exit(EXIT_FAILURE);
    }
	
	//    Cost inclowerbound = MIN_COST;
	string uaitype;
	ifstream file(fileName);
  	if (!file) { cerr << "Could not open file " << fileName << endl; exit(EXIT_FAILURE); }

	Cost inclowerbound = MIN_COST;
    updateUb( (MAX_COST-UNIT_COST)/MEDIUM_COST/MEDIUM_COST );

    int nbval = 0;
    int nbvar,nbconstr;
    int i,j,k,ic;
    string varname;
    int domsize;
    EnumeratedVariable *x;
    EnumeratedVariable *y;
    EnumeratedVariable *z;
    unsigned int a;
    unsigned int b;
    unsigned int c;
    Cost cost;
    int ntuples;
    int arity;
    vector<TemporaryUnaryConstraint> unaryconstrs;
	                
	list<int> lctrs;

	file >> uaitype;
	
	if (ToulBar2::verbose >= 3) cout << "Reading " << uaitype << "  file." << endl;
	
	
	bool markov = uaitype == string("MARKOV");
	//bool bayes = uaitype == string("BAYES");

    

    file >> nbvar;
    // read variable domain sizes
    for (i = 0; i < nbvar; i++) {
        string varname;
        varname = to_string(i);
        file >> domsize;
        if (ToulBar2::verbose >= 3) cout << "read variable " << i << " of size " << domsize << endl;
	    if(domsize > nbval) nbval = domsize; 
        int theindex = -1;
        if (domsize >= 0) theindex = makeEnumeratedVariable(varname,0,domsize-1);
        else theindex = makeIntervalVariable(varname,0,-domsize-1);
        assert(theindex == i);   
    }


    file >> nbconstr;
    // read each constraint
    for (ic = 0; ic < nbconstr; ic++) {
        file >> arity;

        if(arity > MAX_ARITY)  { cerr << "Nary constraints of arity > " << MAX_ARITY << " not supported" << endl; exit(EXIT_FAILURE); }       
        if (!file) {
            cerr << "Warning: EOF reached before reading all the constraints (initial number of constraints too large?)" << endl;
            break;
        }
        if (arity > 3) {
        	int scopeIndex[MAX_ARITY];        	
            if (ToulBar2::verbose >= 3) cout << "read nary constraint on ";

			for(i=0;i<arity;i++) {
	            file >> j;
	            scopeIndex[i] = j;
	            if (ToulBar2::verbose >= 3) cout << j << " ";
			}     	
			if (ToulBar2::verbose >= 3) cout << endl;
            lctrs.push_back( postNaryConstraint(scopeIndex,arity,MAX_COST) );
        } 
        else if (arity == 3) {
            file >> i;
            file >> j;
            file >> k;
        	if ((i == j) || (i == k) || (k == j)) {
    	       cerr << "Error: ternary constraint!" << endl;
               exit(EXIT_FAILURE);
            }
            x = (EnumeratedVariable *) vars[i];
            y = (EnumeratedVariable *) vars[j];
            z = (EnumeratedVariable *) vars[k];
            if (ToulBar2::verbose >= 3) cout << "read ternary constraint " << ic << " on " << i << "," << j << "," << k << endl;
            vector<Cost> costs;
            for (a = 0; a < x->getDomainInitSize(); a++) {
                for (b = 0; b < y->getDomainInitSize(); b++) {
                    for (c = 0; c < z->getDomainInitSize(); c++) {
                        costs.push_back(MIN_COST);
					}
                }
            }
			lctrs.push_back( postTernaryConstraint(i,j,k,costs) );
		} 
		else if (arity == 2) {
            file >> i;
            file >> j;
			if (ToulBar2::verbose >= 3) cout << "read binary constraint " << ic << " on " << i << "," << j << endl;
        	if (i == j) {
    	       cerr << "Error: binary constraint with only one variable in its scope!" << endl;
               exit(EXIT_FAILURE);
            }
            x = (EnumeratedVariable *) vars[i];
            y = (EnumeratedVariable *) vars[j];
            vector<Cost> costs;
            for (a = 0; a < x->getDomainInitSize(); a++) {
                for (b = 0; b < y->getDomainInitSize(); b++) {
                    costs.push_back(MIN_COST);
                }
            }
            lctrs.push_back( postBinaryConstraint(i,j,costs) );          
        } 
        else if (arity == 1) {
            file >> i;
            if (ToulBar2::verbose >= 3) cout << "read unary constraint " << ic << " on " << i << endl;
		    x = (EnumeratedVariable *) vars[i];
			TemporaryUnaryConstraint unaryconstr;
			unaryconstr.var = x;
	    	unaryconstrs.push_back(unaryconstr);
		    x->queueNC();
            lctrs.push_back(-1);            
        } else if(arity == 0) {
            lctrs.push_back(-2);            
        }
        
    }
    
	int iunaryctr = 0;
	int ictr = 0;
	Constraint*	ctr = NULL;
	TernaryConstraint* tctr = NULL;
	BinaryConstraint* bctr = NULL;
	NaryConstraint* nctr = NULL;
	String s;

	ToulBar2::markov_log = 0;   // for the MARKOV Case   
					
	list<int>::iterator it = lctrs.begin();
	while(it !=  lctrs.end()) {
	
		int arity;
		if(*it == -1) { ctr = NULL; arity = 1; }
		else if(*it == -2) { ctr = NULL; arity = 0; }
		else { ctr = getCtr(*it); arity = ctr->arity(); }
		
		file >> ntuples;

		TProb p;
		vector<Cost>  costs;
		vector<TProb> costsProb;
	
		TProb maxp = 0.;	
		for (k = 0; k < ntuples; k++) {
	        file >> p;
	        costsProb.push_back( p );
	        if(p > maxp) maxp = p;
	    }

		Cost minc = MAX_COST;	
		for (k = 0; k < ntuples; k++) {
			p = costsProb[k];
			Cost cost;
	        if(markov) cost = Prob2Cost(p / maxp);
	        else 	   cost = Prob2Cost(p);
			costs.push_back(cost);
			if(cost < minc) minc = cost;
	        }

		if(minc > MIN_COST) {	    
			for (k = 0; k < ntuples; k++) {
				costs[k] -= minc;
			}
		    inclowerbound += minc;
		}
		
		if(markov) ToulBar2::markov_log += log10( maxp );	
			
		switch(arity) {
			case 0:     inclowerbound += costs[0];
						break;
							
			case 1: 	unaryconstrs[iunaryctr].costs.clear();
						for (a = 0; a < unaryconstrs[iunaryctr].var->getDomainInitSize(); a++) {
						      unaryconstrs[iunaryctr].costs.push_back(costs[a]);
						}			
						iunaryctr++; 
						if (ToulBar2::verbose >= 3) cout << "read unary costs."  << endl;							
						break;	
	
			case 2: bctr = (BinaryConstraint*) ctr;
					x = (EnumeratedVariable*) bctr->getVar(0);
            		y = (EnumeratedVariable*) bctr->getVar(1);
            		bctr->addCosts( x,y, costs );
					ictr++;
					if (ToulBar2::verbose >= 3) cout << "read binary costs."  << endl;							
					break;

			case 3: tctr = (TernaryConstraint*) ctr;
					x = (EnumeratedVariable*) tctr->getVar(0);
            		y = (EnumeratedVariable*) tctr->getVar(1);
            		z = (EnumeratedVariable*) tctr->getVar(2);
		    		tctr->addCosts( x,y,z, costs );
        			ictr++;
					if (ToulBar2::verbose >= 3) cout << "read ternary costs." << endl;							
					break;
			 
			default: nctr = (NaryConstraint*) ctr;
					j = 0;
					nctr->firstlex();
					while(nctr->nextlex(s,cost)) {
						nctr->setTuple(s, costs[j]);
						j++;
					}
					ictr++; 
		            //((NaryConstraintMap*) nctr)->preprojectall2();
		            //((NaryConstraintMap*) nctr)->preproject3();
				    if (ToulBar2::verbose >= 3) cout << "read arity " << arity << " table costs."  << endl;						nctr->propagate();
					break;
			
		}
		++it;
	}

    sortVariables();
    sortConstraints();
    // apply basic initial propagation AFTER complete network loading
	//    increaseLb(inclowerbound);

    for (unsigned int u=0; u<unaryconstrs.size(); u++) {
        for (a = 0; a < unaryconstrs[u].var->getDomainInitSize(); a++) {
            if (unaryconstrs[u].costs[a] > MIN_COST) unaryconstrs[u].var->project(a, unaryconstrs[u].costs[a]);
        }
        unaryconstrs[u].var->findSupport();
    }
    if (ToulBar2::verbose >= 0) {
        cout << "Read " << nbvar << " variables, with " << nbval << " values at most, and " << nbconstr << " constraints." << endl;
    }   
 
 	int nevi = 0;	
	ifstream fevid(ToulBar2::evidence_file.c_str());
  	if (!fevid) 
  	{ 
  		string strevid(string(fileName) + string(".evid"));
  		fevid.open(strevid.c_str());
  		cerr << "No evidence file specified. Trying " << strevid << endl;
		if(!fevid) cerr << "No evidence file. " << endl;			 
  	}
  	if(fevid) {
		fevid >> nevi;
	 	while(nevi) {
	 		if(!fevid) {
				cerr << "Error: incorrect number of evidences." << endl;
	            exit(EXIT_FAILURE);
	 		}
	 		fevid >> i;
	 		fevid >> j;
	 		getVar(i)->assign(j);
	 		nevi--;
	 	}
  	}
 	
    increaseLb(inclowerbound);
 
    histogram();
}




void WCSP::solution_UAI(Cost res)
{
 	if (!ToulBar2::uai) return;
    cout << "t " << cpuTime() - ToulBar2::startCpuTime << endl;
	cout << "s " << Cost2LogLike(res) + ToulBar2::markov_log << " ";
	ifstream sol;
	sol.open("sol");	
	if(sol) { 	
		cout << vars.size() << " ";
	    for (unsigned int i=0; i<vars.size(); i++) {
			int value;
	    	sol >> value;
	    	cout << value << " ";
	    }
	}
	cout << endl;
	sol.close();
}





#ifdef XMLFLAG
#include "./xmlcsp/xmlcsp.h"
#endif


void WCSP::read_XML(const char *fileName)
{
	 #ifdef XMLFLAG
		 MyCallback xmlCallBack; 
		 xmlCallBack.wcsp = this;	
	 	 xmlCallBack.fname = string(fileName);
	  	 xmlCallBack.convertWCSP = true;
		 try {
		    XMLParser_libxml2<> parser( xmlCallBack );
		    parser.setPreferredExpressionRepresentation(INFIX_C);
		    parser.parse(fileName); 
		  } catch (exception &e) {
		    cout.flush();
		    cerr << "\n\tUnexpected exception in XML parsing\n";
		    cerr << "\t" << e.what() << endl;
		    exit(1);
		  }
	#else	
		cerr << "\nXML format without including in Makefile flag XMLFLAG and files ./xmlcsp\n" << endl;			   
	    exit(1);
	#endif
}


void WCSP::solution_XML(bool opt)
  {
	 #ifdef XMLFLAG
	 	if (!ToulBar2::xmlflag) return;
	
		if(opt)  cout << "s OPTIMUM FOUND" << endl;	

		//ofstream fsol;
		ifstream sol;
		sol.open("sol");	
		//if(!sol) { cout << "cannot open solution file to translate" << endl; exit(1); }
		//fsol.open("solution");	
		//fsol << "SOL ";


		cout << "v "; 
	    for (unsigned int i=0; i<vars.size(); i++) {
			int value;
	    	sol >> value;
			int index = ((EnumeratedVariable*) getVar(i))->toIndex(value);	    	
		  	cout << Doms[varsDom[i]][ index ] << " ";
	    }
		cout << endl;
			    
		//fsol << endl;	
		//fsol.close();
		sol.close();
	#endif
  }





