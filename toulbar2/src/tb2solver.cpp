/*
 * **************** Generic solver *******************
 *
 */

#include "tb2solver.hpp"
#include "tb2domain.hpp"
#include "tb2pedigree.hpp"
#include "tb2haplotype.hpp"
#include "tb2bep.hpp"
#include "tb2clusters.hpp"
#include <strings.h>

extern void setvalue(int wcspId, int varIndex, Value value);


/*
 * Solver constructors
 *
 */

Solver *Solver::currentSolver = NULL;

Solver::Solver(int storeSize, Cost initUpperBound) : store(NULL), nbNodes(0), nbBacktracks(0), nbBacktracksLimit(LONGLONG_MAX), wcsp(NULL),
                                                     allVars(NULL), unassignedVars(NULL), lastConflictVar(-1), nbSol(0.), nbSGoods(0), nbSGoodsUse(0)
{
    store = new Store(storeSize);
    wcsp = WeightedCSP::makeWeightedCSP(store, initUpperBound);
}

Solver::~Solver()
{
    delete store;
    delete wcsp;
    delete unassignedVars;
	delete[] allVars;
}

void Solver::initVarHeuristic()
{
    unassignedVars = new BTList<Value>(&store->storeDomain);
    allVars = new DLink<Value>[wcsp->numberOfVariables()];
    for (unsigned int j=0; j<wcsp->numberOfVariables(); j++) {
	  unsigned int i = wcsp->getDACOrder(j);
      allVars[i].content = j;
	}
    for (unsigned int i=0; i<wcsp->numberOfVariables(); i++) {
        unassignedVars->push_back(&allVars[i], false);
		if (wcsp->assigned(allVars[i].content) || (ToulBar2::nbDecisionVars > 0 && allVars[i].content >=  ToulBar2::nbDecisionVars)) unassignedVars->erase(&allVars[i], false);
    }
    // Now function setvalue can be called safely!
	ToulBar2::setvalue = setvalue;
}

void Solver::read_wcsp(const char *fileName)
{
    ToulBar2::setvalue = NULL;
    wcsp->read_wcsp(fileName);
}

void Solver::read_random(int n, int m, vector<int>& p, int seed, bool forceSubModular)
{
    ToulBar2::setvalue = NULL;
    wcsp->read_random(n,m,p,seed, forceSubModular);
}

void Solver::read_solution(const char *filename)
{
    currentSolver = this;

	if (ToulBar2::btdMode>=2) wcsp->propagate();

	int depth = store->getDepth();
    store->store();

    // open the file
    ifstream file(filename);
    if (!file) {
        cerr << "Solution file " << filename << " not found!" << endl;
        exit(EXIT_FAILURE);
    }

	vector<int> variables;
	vector<Value> values;
    int i = 0;
    while (!file.eof()) {
        if ((unsigned int) i >= wcsp->numberOfVariables()) break;
        Value value = 0;
        file >> value;
        variables.push_back(i);
        values.push_back(value);
        // side-effect: remember last solution
        wcsp->setBestValue(i, value);
//        if (wcsp->unassigned(i)) {
//		  assign(i, value);
//		  // side-effect: remember last solution
//		  wcsp->setBestValue(i, value);
//        } else {
//		    if (wcsp->getValue(i) != value) {
//			  THROWCONTRADICTION;
//			} else {
//			  wcsp->setBestValue(i, value); // side-effect: remember last solution
//			}
//        }
        i++;
    }
    wcsp->assignLS(variables, values);
    wcsp->propagate();
    cout << " Solution cost: [" << wcsp->getLb() << "," << wcsp->getUb() << "] (nb. of unassigned variables: " << wcsp->numberOfUnassignedVariables() << ")" << endl;
	assert(wcsp->numberOfUnassignedVariables() == 0);
	wcsp->updateUb(wcsp->getLb()+UNIT_COST);
    store->restore(depth);
}

void Solver::parse_solution(const char *certificate)
{
    currentSolver = this;

    if (ToulBar2::btdMode>=2) wcsp->propagate();

	//  int depth = store->getDepth();
	//    store->store();

    //certif2 = index(certif2,',');
   char *certif2;
   char sep[]=",";
    certif2 = strdup(certificate);
    certif2= strstr(certif2,sep);

    if (certif2) certif2++;
    
	vector<int> variables;
	vector<Value> values;
    int var;
	Value value;
	int items;
    while ((certif2 != NULL) && (certif2[0] != 0)) {
        items = sscanf(certif2,"%d=%d",&var,&value);
        if ((items != 2) || ((unsigned int)var >= wcsp->numberOfVariables())) {
             cerr << "Certificate " << certif2 << " incorrect!" << endl;
             exit(EXIT_FAILURE);
        }
        certif2 = strstr(certif2,sep);
        if (certif2) certif2++;

        variables.push_back(var);
        values.push_back(value);
        // side-effect: remember last solution
        wcsp->setBestValue(var, value);
//        if (wcsp->unassigned(var)) {
//          assign(var, value);
//          // side-effect: remember last solution
//          wcsp->setBestValue(var, value);
//        } else {
//		  if (wcsp->getValue(var) != value) {
//			THROWCONTRADICTION;
//		  } else {
//			wcsp->setBestValue(var, value); // side-effect: remember last solution
//		  }
//        }
    }
    wcsp->assignLS(variables, values);
    wcsp->propagate();

    cout << " Solution cost: [" << wcsp->getLb() << "," << wcsp->getUb() << "] (nb. of unassigned variables: " << wcsp->numberOfUnassignedVariables() << ")" << endl;
    
//    if (ToulBar2::btdMode>=2) wcsp->updateUb(wcsp->getLb()+UNIT_COST);
//    store->restore(depth);
}

void Solver::dump_wcsp(const char *fileName, bool original)
{
    ofstream pb(fileName);
    if (pb) wcsp->dump(pb, original);
}

/*
 * Link between solver and wcsp: maintain a backtrackable list of unassigned variable indexes
 *
 */

void setvalue(int wcspId, int varIndex, Value value)
{
//    assert(wcspId == 0); // WARNING! assert not compatible with sequential execution of solve() method
    assert(Solver::currentSolver != NULL);
	unsigned int i = Solver::currentSolver->getWCSP()->getDACOrder(varIndex);
    if(!Solver::currentSolver->allVars[i].removed) {
	  Solver::currentSolver->unassignedVars->erase(&Solver::currentSolver->allVars[i], true);
	}
}

/*
 * Variable ordering heuristics
 *
 */

/// \defgroup heuristics Variable and value search ordering heuristics
/// \see <em> Boosting Systematic Search by Weighting Constraints </em>. Frédéric Boussemart, Fred Hemery, Christophe Lecoutre, Lakhdar Sais. Proc. of ECAI 2004, pages 146-150. Valencia, Spain, 2004.
/// \see <em> Last Conflict Based Reasoning </em>. Christophe Lecoutre, Lakhdar Sais, Sébastien Tabary, Vincent Vidal. Proc. of ECAI 2006, pages 133-137. Trentino, Italy, 2006.

int Solver::getNextUnassignedVar()
{
  //    if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar)) return lastConflictVar;
    return (unassignedVars->empty())?-1:(*unassignedVars->begin());
}

int Solver::getVarMinDomainDivMaxDegree()
{
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;

    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
        double heuristic = (double) wcsp->getDomainSize(*iter) / (double) (wcsp->getDegree(*iter)+1);
        if (varIndex < 0 || heuristic < best - 1./100001.
            || (heuristic < best + 1./100001. && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
        }
    }
    return varIndex;
}

int Solver::getVarMinDomainDivMaxDegreeRandomized()
{
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;
	int ties[unassignedVars->getSize()];
	int nbties = 0;

    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
        double heuristic = (double) wcsp->getDomainSize(*iter) / (double) (wcsp->getDegree(*iter)+1);
        if (varIndex < 0 || heuristic < best - 1./1000001.
            || (heuristic < best + 1./1000001. && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
			nbties = 1;
			ties[0] = varIndex;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
        } else if (heuristic < best + 1./1000001. && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) {
		 ties[nbties] = *iter;
		 nbties++;
		}
    }
    if (nbties>1) return ties[myrand()%nbties];
    else return varIndex;
}

int Solver::getVarMinDomainDivMaxDegreeLastConflict()
{
    if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar)) return lastConflictVar;
	// int varIndexVAC = wcsp->getVACHeuristic();
	// if(varIndexVAC != -1) return varIndexVAC;
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;

    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
        // remove following "+1" when isolated variables are automatically assigned
        double heuristic = (double) wcsp->getDomainSize(*iter) / (double) (wcsp->getDegree(*iter) + 1);
        if (varIndex < 0 || heuristic < best - 1./100001.
            || (heuristic < best + 1./100001. && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
        }
    }
    return varIndex;
}

int Solver::getVarMinDomainDivMaxDegreeLastConflictRandomized()
{
    if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar)) return lastConflictVar;
	// int varIndexVAC = wcsp->getVACHeuristic();
	// if(varIndexVAC != -1) return varIndexVAC;
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;
	int ties[unassignedVars->getSize()];
	int nbties = 0;

    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
        // remove following "+1" when isolated variables are automatically assigned
        double heuristic = (double) wcsp->getDomainSize(*iter) / (double) (wcsp->getDegree(*iter) + 1);
        if (varIndex < 0 || heuristic < best - 1./1000001.
            || (heuristic < best + 1./1000001. && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
			nbties = 1;
			ties[0] = varIndex;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
			//        } else if ((heuristic < best + 1./1000001. && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) || ((myrand()%100)==0)) {
        } else if (heuristic < best + 1./1000001. && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) {
		 ties[nbties] = *iter;
		 nbties++;
	   }
    }
    if (nbties>1) {if (ToulBar2::debug>1) cout << "RAND VAR " << nbties << endl; return ties[myrand()%nbties];}
    else return varIndex;
}

int Solver::getVarMinDomainDivMaxWeightedDegree()
{
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;

    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
	  double heuristic = (double) wcsp->getDomainSize(*iter) / (double) (wcsp->getWeightedDegree(*iter)+1);
	  if (varIndex < 0 || heuristic < best - 1./100001.
		  || (heuristic < best + 1./100001. && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
		best = heuristic;
		varIndex = *iter;
		worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
	  }
    }
    return varIndex;
}

int Solver::getVarMinDomainDivMaxWeightedDegreeRandomized()
{
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;
	int ties[unassignedVars->getSize()];
	int nbties = 0;

    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
	  double heuristic = (double) wcsp->getDomainSize(*iter) / (double) (wcsp->getWeightedDegree(*iter)+1);
	  if (varIndex < 0 || heuristic < best - 1./1000001.
		  || (heuristic < best + 1./1000001. && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
		best = heuristic;
		varIndex = *iter;
		nbties = 1;
		ties[0] = varIndex;
		worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
	  } else if (heuristic < best + 1./1000001. && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) {
		 ties[nbties] = *iter;
		 nbties++;
	   }
    }
	if (nbties>1) return ties[myrand()%nbties];
    else return varIndex;
}

int Solver::getVarMinDomainDivMaxWeightedDegreeLastConflict()
{
   if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar)) return lastConflictVar;
   int varIndex = -1;
   Cost worstUnaryCost = MIN_COST;
   double best = MAX_VAL - MIN_VAL;

   for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
       // remove following "+1" when isolated variables are automatically assigned
	 double heuristic = (double) wcsp->getDomainSize(*iter) / (double) (wcsp->getWeightedDegree(*iter) + 1);
       if (varIndex < 0 || heuristic < best - 1./100001.
           || (heuristic < best + 1./100001. && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
           best = heuristic;
           varIndex = *iter;
           worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
       }
   }
   return varIndex;
}

int Solver::getVarMinDomainDivMaxWeightedDegreeLastConflictRandomized()
{
   if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar)) return lastConflictVar;
   int varIndex = -1;
   Cost worstUnaryCost = MIN_COST;
   double best = MAX_VAL - MIN_VAL;
   int ties[unassignedVars->getSize()];
   int nbties = 0;

   for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
       // remove following "+1" when isolated variables are automatically assigned
	 double heuristic = (double) wcsp->getDomainSize(*iter) / (double) (wcsp->getWeightedDegree(*iter) + 1);
       if (varIndex < 0 || heuristic < best - 1./1000001.
           || (heuristic < best + 1./1000001. && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
           best = heuristic;
           varIndex = *iter;
		   nbties = 1;
		   ties[0] = varIndex;
           worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
		   //       } else if ((heuristic < best + 1./1000001. && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) || ((myrand()%100)==0)) {
       } else if (heuristic < best + 1./1000001. && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) {
		 ties[nbties] = *iter;
		 nbties++;
	   }
   }
   if (nbties>1) {if (ToulBar2::debug>1) cout << "RAND VAR " << nbties << endl; return ties[myrand()%nbties];}
   else return varIndex;
}

int Solver::getMostUrgent()
{
    int varIndex = -1;
	Value best = MAX_VAL;
    Cost worstUnaryCost = MIN_COST;

    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
        if (varIndex < 0 || wcsp->getInf(*iter) < best ||
			(wcsp->getInf(*iter) == best && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
		  best = wcsp->getInf(*iter);
		  worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
		  varIndex = *iter;
        }
	}
    return varIndex;
}

/*
 * Choice points
 *
 */

/// \brief Enforce WCSP upper-bound and backtrack if ub <= lb or in the case of probabilistic inference if the contribution is too small
void Solver::enforceUb()
{
	wcsp->enforceUb();
	if (ToulBar2::isZ) {
		Cost newCost = wcsp->getLb() + wcsp->getNegativeLb() + wcsp->LogLike2Cost(unassignedVars->getSize() * Log10(wcsp->getMaxDomainSize()));
		Double newlogz = wcsp->SumLogLikeCost(ToulBar2::logZ, newCost);
		if (newlogz - ToulBar2::logZ < 0.001) {
			THROWCONTRADICTION;
		}
	}
}

void Solver::increase(int varIndex, Value value)
{
    enforceUb();
    nbNodes++;
    if (ToulBar2::verbose >= 1) {
        if (ToulBar2::verbose >= 2) cout << *wcsp;
		if (ToulBar2::debug >= 3) {
		  string pbname = "problem" + to_string(nbNodes) + ".wcsp";
		  ofstream pb(pbname.c_str());
		  wcsp->dump(pb);
		}
        cout << "[" << store->getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum();
		if (wcsp->getTreeDec()) cout << ",C" << wcsp->getTreeDec()->getCurrentCluster()->getId();
		cout << "] Try " << wcsp->getName(varIndex) << " >= " << value << " (s:" << wcsp->getSupport(varIndex) << ")" << endl;
    }
    wcsp->increase(varIndex, value);
    wcsp->propagate();
}

void Solver::decrease(int varIndex, Value value)
{
    enforceUb();
    nbNodes++;
    if (ToulBar2::verbose >= 1) {
        if (ToulBar2::verbose >= 2) cout << *wcsp;
		if (ToulBar2::debug >= 3) {
		  string pbname = "problem" + to_string(nbNodes) + ".wcsp";
		  ofstream pb(pbname.c_str());
		  wcsp->dump(pb);
		}
        cout << "[" << store->getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum();
		if (wcsp->getTreeDec()) cout << ",C" << wcsp->getTreeDec()->getCurrentCluster()->getId();
		cout << "] Try " << wcsp->getName(varIndex) << " <= " << value << " (s:" << wcsp->getSupport(varIndex) << ")" << endl;
    }
    wcsp->decrease(varIndex, value);
    wcsp->propagate();
}

void Solver::assign(int varIndex, Value value)
{
    enforceUb();
    nbNodes++;
	if (ToulBar2::debug && ((nbNodes % 128) == 0)) {
	  cout << "\r" << store->getDepth();
	  if (wcsp->getTreeDec()) cout << " C" << wcsp->getTreeDec()->getCurrentCluster()->getId();
	  cout << "   ";
	  cout.flush();
	}
    if (ToulBar2::verbose >= 1) {
        if (ToulBar2::verbose >= 2) cout << *wcsp;
		if (ToulBar2::debug >= 3) {
		  string pbname = "problem" + to_string(nbNodes) + ".wcsp";
		  ofstream pb(pbname.c_str());
		  wcsp->dump(pb);
		}
        cout << "[" << store->getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum();
		if (wcsp->getTreeDec()) cout << ",C" << wcsp->getTreeDec()->getCurrentCluster()->getId();
		cout << "] Try " << wcsp->getName(varIndex) << " == " << value << endl;
    }
    wcsp->assign(varIndex, value);
    wcsp->propagate();
}

void Solver::remove(int varIndex, Value value)
{
    enforceUb();
    nbNodes++;
    if (ToulBar2::verbose >= 1) {
        if (ToulBar2::verbose >= 2) cout << *wcsp;
		if (ToulBar2::debug >= 3) {
		  string pbname = "problem" + to_string(nbNodes) + ".wcsp";
		  ofstream pb(pbname.c_str());
		  wcsp->dump(pb);
		}
        cout << "[" << store->getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum();
		if (wcsp->getTreeDec()) cout << ",C" << wcsp->getTreeDec()->getCurrentCluster()->getId();
		cout << "] Try " << wcsp->getName(varIndex) << " != " << value << endl;
    }
    wcsp->remove(varIndex, value);
    wcsp->propagate();
}

void Solver::binaryChoicePoint(int varIndex, Value value)
{
    assert(wcsp->unassigned(varIndex));
    assert(wcsp->canbe(varIndex,value));
    bool dichotomic = (ToulBar2::dichotomicBranching && ToulBar2::dichotomicBranchingSize < wcsp->getDomainSize(varIndex));
    Value middle = value;
    bool increasing = true;
	//	bool reverse = ((ToulBar2::restart>0) && (myrand()%10==0));
    if (dichotomic) {
      middle = (wcsp->getInf(varIndex) + wcsp->getSup(varIndex)) / 2;
	  //      if (value <= middle || reverse) increasing = true;
      if (value <= middle) increasing = true;
      else increasing = false;
    }
    try {
        store->store();
        lastConflictVar = varIndex;
        if (dichotomic) {
    	  if (increasing) decrease(varIndex, middle); else increase(varIndex, middle+1);
    	// } else if (reverse) {
		//   remove(varIndex, value);
		} else assign(varIndex, value);
        lastConflictVar = -1;
        recursiveSolve();
    } catch (Contradiction) {
        wcsp->whenContradiction();
    }
    store->restore();
    enforceUb();
    nbBacktracks++;
	if (ToulBar2::restart>0 && nbBacktracks > nbBacktracksLimit) throw NbBacktracksOut();

    if (dichotomic) {
      if (increasing) increase(varIndex, middle+1); else decrease(varIndex, middle);
    // } else if (reverse) {
	//   assign(varIndex, value);
	} else remove(varIndex, value);
    recursiveSolve();

}

void Solver::binaryChoicePointLDS(int varIndex, Value value, int discrepancy)
{
    assert(wcsp->unassigned(varIndex));
    assert(wcsp->canbe(varIndex,value));
    bool dichotomic = (ToulBar2::dichotomicBranching && ToulBar2::dichotomicBranchingSize < wcsp->getDomainSize(varIndex));
    Value middle = value;
    bool increasing = true;
    if (dichotomic) {
      middle = (wcsp->getInf(varIndex) + wcsp->getSup(varIndex)) / 2;
      if (value <= middle) increasing = true;
      else increasing = false;
    }
    if (discrepancy > 0) {
        try {
            store->store();
            lastConflictVar = varIndex;
            if (dichotomic) {
              if (increasing) increase(varIndex, middle+1); else decrease(varIndex, middle);
            } else remove(varIndex, value);
            lastConflictVar = -1;
            recursiveSolveLDS(discrepancy - 1);
        } catch (Contradiction) {
            wcsp->whenContradiction();
        }
        store->restore();
        enforceUb();
        nbBacktracks++;
		if (ToulBar2::restart>0 && nbBacktracks > nbBacktracksLimit) throw NbBacktracksOut();
        if (dichotomic) {
          if (increasing) decrease(varIndex, middle); else increase(varIndex, middle+1);
        } else assign(varIndex, value);
        recursiveSolveLDS(discrepancy);
    } else {
        ToulBar2::limited = true;
        lastConflictVar = varIndex;
        if (dichotomic) {
          if (increasing) decrease(varIndex, middle); else increase(varIndex, middle+1);
        } else assign(varIndex, value);
        lastConflictVar = -1;
        recursiveSolveLDS(0);
    }
}

Value Solver::postponeRule(int varIndex)
{
  assert(ToulBar2::bep);
  Value best = ToulBar2::bep->latest[varIndex] + 1;

  for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
	if (*iter != varIndex) {
	  Value time = wcsp->getInf(*iter) + ToulBar2::bep->duration[*iter] + ToulBar2::bep->delay[*iter * ToulBar2::bep->size + varIndex];
	  if (time < best) {
		best = time;
	  }
	}
  }
  return best;
}

void Solver::scheduleOrPostpone(int varIndex)
{
    assert(wcsp->unassigned(varIndex));
	Value xinf = wcsp->getInf(varIndex);
	Value postponeValue = postponeRule(varIndex);
	postponeValue = max(postponeValue, xinf+1);
	assert(postponeValue <= ToulBar2::bep->latest[varIndex]+1);
	bool reverse = (wcsp->getUnaryCost(varIndex,xinf) > MIN_COST)?true:false;
    try {
        store->store();
        if (reverse) increase(varIndex, postponeValue);
		else assign(varIndex, xinf);
        recursiveSolve();
    } catch (Contradiction) {
        wcsp->whenContradiction();
    }
    store->restore();
    enforceUb();
    nbBacktracks++;
	if (ToulBar2::restart>0 && nbBacktracks > nbBacktracksLimit) throw NbBacktracksOut();
	if (reverse) assign(varIndex, xinf);
	else increase(varIndex, postponeValue);
    recursiveSolve();
}

int cmpValueCost(const void *p1, const void *p2)
{
    Cost c1 = ((ValueCost *) p1)->cost;
    Cost c2 = ((ValueCost *) p2)->cost;
    Value v1 = ((ValueCost *) p1)->value;
    Value v2 = ((ValueCost *) p2)->value;
    if (c1 < c2) return -1;
    else if (c1 > c2) return 1;
    else if (v1 < v2) return -1;
    else if (v1 > v2) return 1;
    else return 0;
}

void Solver::narySortedChoicePoint(int varIndex)
{
    assert(wcsp->enumerated(varIndex));
    int size = wcsp->getDomainSize(varIndex);
    //ValueCost sorted[size];
    ValueCost* sorted = new ValueCost [size];
    wcsp->getEnumDomainAndCost(varIndex, sorted);
    qsort(sorted, size, sizeof(ValueCost), cmpValueCost);
    for (int v = 0; wcsp->getLb() < wcsp->getUb() && v < size; v++) {
        try {
            store->store();
            assign(varIndex, sorted[v].value);
            recursiveSolve();
        } catch (Contradiction) {
            wcsp->whenContradiction();
        }
        store->restore();
    }
	delete [] sorted;
    enforceUb();
    nbBacktracks++;
	if (ToulBar2::restart>0 && nbBacktracks > nbBacktracksLimit) throw NbBacktracksOut();
}

/* returns true if a limit has occured during the search */
void Solver::narySortedChoicePointLDS(int varIndex, int discrepancy)
{
    assert(wcsp->enumerated(varIndex));
    int size = wcsp->getDomainSize(varIndex);
    ValueCost* sorted = new ValueCost [size];
	wcsp->getEnumDomainAndCost(varIndex, sorted);
    qsort(sorted, size, sizeof(ValueCost), cmpValueCost);
    if (discrepancy < size-1) ToulBar2::limited = true;
    for (int v = min(size-1, discrepancy); wcsp->getLb() < wcsp->getUb() && v >= 0; v--) {
        try {
            store->store();
            assign(varIndex, sorted[v].value);
            recursiveSolveLDS(discrepancy - v);
        } catch (Contradiction) {
            wcsp->whenContradiction();
        }
        store->restore();
    }
	delete [] sorted;
    enforceUb();
    nbBacktracks++;
	if (ToulBar2::restart>0 && nbBacktracks > nbBacktracksLimit) throw NbBacktracksOut();
}

void Solver::singletonConsistency()
{
    bool deadend;
    bool done = false;
    while(!done) {
    	done = true;
	    for (unsigned int varIndex = 0; varIndex < ((ToulBar2::nbDecisionVars>0)?ToulBar2::nbDecisionVars:wcsp->numberOfVariables()); varIndex++) {
			  int size = wcsp->getDomainSize(varIndex);
		      ValueCost* sorted = new ValueCost [size];
			  wcsp->iniSingleton();
			  wcsp->getEnumDomainAndCost(varIndex, sorted);
			  qsort(sorted, size, sizeof(ValueCost), cmpValueCost);
			  for (int a = 0; a < size; a++) {
					deadend = false;
			        try {
			            store->store();
			            assign(varIndex, sorted[a].value);
			        } catch (Contradiction) {
			            wcsp->whenContradiction();
			            deadend = true;
			            done = false;
			        }
			        store->restore();
					wcsp->updateSingleton();
					//cout << "(" << varIndex << "," << a <<  ")" << endl;
					if(deadend) {
					  remove(varIndex, sorted[a].value);
					  cout << "."; flush(cout);
// WARNING!!! can we stop if the variable is assigned, what about removeSingleton after???
					}
		      }
			  wcsp->removeSingleton();
			  delete [] sorted;
	    }
    }
    cout << "Done Singleton Consistency" << endl;
}

/*
 * Depth-First Branch and Bound
 *
 */

void Solver::newSolution()
{
    assert(unassignedVars->empty());
#ifndef NDEBUG
    bool allVarsAssigned = true;
    for (unsigned int i=0; i<wcsp->numberOfVariables(); i++) {
        if (wcsp->unassigned(i)) {
            allVarsAssigned = false;
            break;
        }
    }
    assert(allVarsAssigned);
#endif
    if (!ToulBar2::allSolutions && !ToulBar2::isZ) wcsp->updateUb(wcsp->getLb());
	else if (!ToulBar2::btdMode) nbSol += 1.;
	if (ToulBar2::isZ) {
	  ToulBar2::logZ = wcsp->SumLogLikeCost(ToulBar2::logZ, wcsp->getLb() + wcsp->getNegativeLb());
	  if (ToulBar2::debug && (nbBacktracks % 10000LL)==0) cout << (ToulBar2::logZ + ToulBar2::markov_log) << endl;
	}
  	if((!ToulBar2::allSolutions && !ToulBar2::isZ) || ToulBar2::debug>=2) {
  		if (ToulBar2::verbose>=0 || ToulBar2::showSolutions) {
  			if(ToulBar2::haplotype) cout <<  "***New solution: " <<  wcsp->getLb() << " log10like: " << ToulBar2::haplotype->Cost2LogLike(wcsp->getLb())<< " logProb: " << ToulBar2::haplotype->Cost2Prob( wcsp->getLb()) << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << store->getDepth() << ")" << endl;
  			else if(!ToulBar2::bayesian) cout << "New solution: " <<  wcsp->getLb() << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << store->getDepth() << ")" << endl;
  			else cout << "New solution: " <<  wcsp->getLb() << " log10like: " << wcsp->Cost2LogLike(wcsp->getLb() + wcsp->getNegativeLb()) + ToulBar2::markov_log << " prob: " << wcsp->Cost2Prob( wcsp->getLb() + wcsp->getNegativeLb() ) * Exp10(ToulBar2::markov_log) << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << store->getDepth() << ")" << endl;
  		}
  	}
  	else {
  		if(ToulBar2::xmlflag) {
		  cout << "o " << wcsp->getLb() << endl; //" ";
		  //	((WCSP*)wcsp)->solution_XML();
  		}
  	}

    wcsp->restoreSolution();

    if (ToulBar2::showSolutions) {

        if (ToulBar2::verbose >= 2) cout << *wcsp << endl;

        if(ToulBar2::allSolutions) {
        	cout << nbSol << " solution: ";
        }

        for (unsigned int i=0; i<wcsp->numberOfVariables(); i++) {
	        cout << " ";
            if (ToulBar2::pedigree) {
                cout <<  wcsp->getName(i) << ":";
                ToulBar2::pedigree->printGenotype(cout, wcsp->getValue(i));
            } else if (ToulBar2::haplotype) {
            	ToulBar2::haplotype->printHaplotype(cout,wcsp->getValue(i),i);
			} else {
                cout << wcsp->getValue(i);
            }
        }
        cout << endl;
		if (ToulBar2::bep) ToulBar2::bep->printSolution((WCSP *) wcsp);
    }
    if (ToulBar2::pedigree) {
      ToulBar2::pedigree->printCorrection((WCSP *) wcsp);
    }
    if (ToulBar2::writeSolution) {
        if (ToulBar2::pedigree) {
            ToulBar2::pedigree->save("pedigree_corrected.pre", (WCSP *) wcsp, true, false);
            ToulBar2::pedigree->printSol((WCSP*) wcsp);
            ToulBar2::pedigree->printCorrectSol((WCSP*) wcsp);
        } else if (ToulBar2::haplotype) {
		  ToulBar2::haplotype->printSol((WCSP*) wcsp);
		}
//        else {
	        ofstream file("sol");
	        if (!file) {
	          cerr << "Could not write file " << "solution" << endl;
	          exit(EXIT_FAILURE);
	        }
	        for (unsigned int i=0; i<wcsp->numberOfVariables(); i++) {
	            file << " " << wcsp->getValue(i);
	        }
	        file << endl;
//        }
    }
	if(ToulBar2::uai && !ToulBar2::isZ) {
	  ((WCSP*)wcsp)->solution_UAI(wcsp->getLb());
	}

	if (ToulBar2::newsolution) (*ToulBar2::newsolution)(wcsp->getIndex());

	if (ToulBar2::restart==0 && !ToulBar2::lds && !ToulBar2::isZ) throw NbBacktracksOut();
}

void Solver::recursiveSolve()
{
	int varIndex = -1;
	if (ToulBar2::bep) varIndex = getMostUrgent();
	else if (ToulBar2::Static_variable_ordering) varIndex = getNextUnassignedVar();
	else if(ToulBar2::weightedDegree && ToulBar2::lastConflict) varIndex = ((ToulBar2::restart>0)?getVarMinDomainDivMaxWeightedDegreeLastConflictRandomized():getVarMinDomainDivMaxWeightedDegreeLastConflict());
	else if(ToulBar2::lastConflict) varIndex = ((ToulBar2::restart>0)?getVarMinDomainDivMaxDegreeLastConflictRandomized():getVarMinDomainDivMaxDegreeLastConflict());
	else if(ToulBar2::weightedDegree) varIndex = ((ToulBar2::restart>0)?getVarMinDomainDivMaxWeightedDegreeRandomized():getVarMinDomainDivMaxWeightedDegree());
	else varIndex = ((ToulBar2::restart>0)?getVarMinDomainDivMaxDegreeRandomized():getVarMinDomainDivMaxDegree());
    if (varIndex >= 0) {
	    if (ToulBar2::bep) scheduleOrPostpone(varIndex);
        else if (wcsp->enumerated(varIndex)) {
            if (ToulBar2::binaryBranching) {
			  assert(wcsp->canbe(varIndex, wcsp->getSupport(varIndex)));
			  // Reuse last solution found if available
			  Value bestval = wcsp->getBestValue(varIndex);
			  binaryChoicePoint(varIndex, (wcsp->canbe(varIndex, bestval))?bestval:wcsp->getSupport(varIndex));
            } else narySortedChoicePoint(varIndex);
        } else {
            binaryChoicePoint(varIndex, wcsp->getInf(varIndex));
        }
    } else newSolution();

}

/* returns true if a limit has occured during the search */
void Solver::recursiveSolveLDS(int discrepancy)
{
	int varIndex = -1;
	if (ToulBar2::bep) varIndex = getMostUrgent();
	else if(ToulBar2::weightedDegree && ToulBar2::lastConflict) varIndex = ((ToulBar2::restart>0)?getVarMinDomainDivMaxWeightedDegreeLastConflictRandomized():getVarMinDomainDivMaxWeightedDegreeLastConflict());
	else if(ToulBar2::lastConflict) varIndex =  ((ToulBar2::restart>0)?getVarMinDomainDivMaxDegreeLastConflictRandomized():getVarMinDomainDivMaxDegreeLastConflict());
	else if(ToulBar2::weightedDegree) varIndex = ((ToulBar2::restart>0)?getVarMinDomainDivMaxWeightedDegreeRandomized():getVarMinDomainDivMaxWeightedDegree());
	else varIndex = ((ToulBar2::restart>0)?getVarMinDomainDivMaxDegreeRandomized():getVarMinDomainDivMaxDegree());
    if (varIndex >= 0) {
	    if (ToulBar2::bep) scheduleOrPostpone(varIndex);
        else if (wcsp->enumerated(varIndex)) {
            if (ToulBar2::binaryBranching) {
			  assert(wcsp->canbe(varIndex, wcsp->getSupport(varIndex)));
			  // Reuse last solution found if available
			  Value bestval = wcsp->getBestValue(varIndex);
			  binaryChoicePointLDS(varIndex, (wcsp->canbe(varIndex, bestval))?bestval:wcsp->getSupport(varIndex), discrepancy);
            } else {
                narySortedChoicePointLDS(varIndex, discrepancy);
            }
        } else {
            binaryChoicePointLDS(varIndex, wcsp->getInf(varIndex), discrepancy);
        }
    } else newSolution();
}

static Cost UpperBound = MAX_COST;
static WeightedCSP *CurrentWeightedCSP = NULL;
static bool IsASolution = false;
static int *CurrentSolution = NULL;

void solution_restart(int wcspIndex)
{
  assert(CurrentWeightedCSP->getIndex() == wcspIndex);
  IsASolution = true;
  if (CurrentWeightedCSP->getUb() < UpperBound) UpperBound = CurrentWeightedCSP->getUb();
}

Long luby(Long r) {
  int j = cost2log2(r+1);
  if (r+1 == (1L << j)) return (1L << (j-1));
  else return luby(r - (1L << j) + 1);
}

bool Solver::solve()
{
    currentSolver = this;
    Cost initialUpperBound = wcsp->getUb();
    nbBacktracks = 0;
    nbNodes = 0;
	lastConflictVar = -1;
	int tailleSep = 0;

	if (ToulBar2::isZ) {
	  ToulBar2::logZ = -numeric_limits<TProb>::infinity();
	}

    try {
//        store->store();       // if uncomment then solve() does not change the problem but all preprocessing operations will allocate in backtrackable memory
        wcsp->decreaseUb(initialUpperBound);

        wcsp->propagate();                // initial propagation
        wcsp->preprocessing();            // preprocessing after initial propagation
		cout << "Preprocessing Time               : " << cpuTime() - ToulBar2::startCpuTime << " seconds." << endl;

		cout << wcsp->numberOfUnassignedVariables() << " unassigned variables, " << wcsp->getDomainSizeSum() << " values in all current domains (med. size:" << wcsp->medianDomainSize() << ") and " << wcsp->numberOfConnectedConstraints() << " non-unary cost functions (med. degree:" << wcsp->medianDegree() << ")" << endl;
		cout << "Initial lower and upper bounds: [" << wcsp->getLb() << "," << wcsp->getUb() << "[" << endl;

        if (ToulBar2::singletonConsistency) {
		  singletonConsistency();
		  wcsp->propagate();
		}

		if (ToulBar2::isZ && ToulBar2::verbose >= 1) cout << "NegativeShiftingCost= " << wcsp->getNegativeLb() << endl;

	    if (ToulBar2::btdMode) {
		  if(wcsp->numberOfUnassignedVariables()==0 || wcsp->numberOfConnectedConstraints()==0)	ToulBar2::approximateCountingBTD = 0;
		  wcsp->buildTreeDecomposition();
	    } else if (ToulBar2::weightedDegree && (((Long) wcsp->numberOfConnectedConstraints()) >= ((Long) ToulBar2::weightedDegree))) {
		  cout << "Weighted degree heuristic disabled (#costfunctions=" << wcsp->numberOfConnectedConstraints() << " >= " << ToulBar2::weightedDegree << ")" << endl;
		  ToulBar2::weightedDegree = 0;
		}
		
		if (ToulBar2::dumpWCSP) {dump_wcsp("problem.wcsp",false); cout << "end." << endl; exit(0);}

        // special data structure to be initialized for variable ordering heuristics
	    initVarHeuristic();

		if (ToulBar2::restart>=0) {
		  if (ToulBar2::restart>0)nbBacktracksLimit = 1;
		  CurrentWeightedCSP = wcsp;
		  UpperBound = wcsp->getUb();
		  ToulBar2::newsolution = solution_restart;
		}
		bool nbbacktracksout = false;
		int nbrestart = 0;
		Long currentNbBacktracksLimit = 1;
		Long nbBacktracksLimitTop = 1;
		int storedepth = store->getDepth();
		do {
		  store->store();
		  if (ToulBar2::restart>=0) {
			nbbacktracksout = false;
			nbrestart++;
			// currentNbBacktracksLimit = max(currentNbBacktracksLimit + 1, (Long) (1.2 * (Double) currentNbBacktracksLimit + 0.5));
			// if (ToulBar2::lds) currentNbBacktracksLimit *= 4;
			currentNbBacktracksLimit = luby(nbrestart);
			if (currentNbBacktracksLimit > nbBacktracksLimitTop || IsASolution) {
			  nbBacktracksLimitTop = currentNbBacktracksLimit;
			  currentNbBacktracksLimit = 1;
			}
			//			if (!IsASolution && nbNodes >= ToulBar2::restart) {
			if (nbNodes >= ToulBar2::restart) {
			  nbBacktracksLimit = LONGLONG_MAX;
			  ToulBar2::restart = 0;
			  cout << "****** Restart " << nbrestart << " with no backtrack limit and UB=" << UpperBound << " ****** (" << nbNodes << " nodes)" << endl;
			} else {
			  nbBacktracksLimit = nbBacktracks + currentNbBacktracksLimit;
			  cout << "****** Restart " << nbrestart << " with " << currentNbBacktracksLimit << " backtracks max and UB=" << UpperBound << " ****** (" << nbNodes << " nodes)" << endl;
			}
			IsASolution = false;
			wcsp->setUb(UpperBound);
			enforceUb();
			wcsp->propagate();
			if (ToulBar2::isZ) ToulBar2::logZ = -numeric_limits<TProb>::infinity();
		  }
		  try {
			if (ToulBar2::restart <= 0 && ToulBar2::lds) {
			  int discrepancy = 0;
			  do {
				if (discrepancy > ToulBar2::lds) {cout << "--- [" << store->getDepth() << "] Search with no discrepancy limit --- (" << nbNodes << " nodes)" << endl;}
				else {cout << "--- [" << store->getDepth() << "] LDS " << discrepancy << " --- (" << nbNodes << " nodes)" << endl;}
				ToulBar2::limited = false;
				try {
				  store->store();
				  if (ToulBar2::isZ) ToulBar2::logZ = -numeric_limits<TProb>::infinity();
				  if (discrepancy > ToulBar2::lds) {ToulBar2::lds = false; recursiveSolve();} else recursiveSolveLDS(discrepancy);
				} catch (Contradiction) {
				  wcsp->whenContradiction();
				}
				store->restore();
				if (discrepancy > 0) discrepancy *= 2;
				else discrepancy++;
			  } while (ToulBar2::limited);
			} else {
			  TreeDecomposition* td = wcsp->getTreeDec();
			  if(td) {
				Cost ub = wcsp->getUb();
				Cluster* start = td->getRoot();
				assert(start->getLbRec() == MIN_COST); // local lower bounds (and delta costs) must be zero!
				if(ToulBar2::btdSubTree >= 0) start = td->getCluster(ToulBar2::btdSubTree);
				td->setCurrentCluster(start);
				if (start==td->getRoot()) start->setLb(wcsp->getLb()); // initial lower bound found by propagation is associated to tree decompostion root cluster
				switch(ToulBar2::btdMode) {
				case 0:case 1:
				  if(ToulBar2::allSolutions)
					{
					  timeDeconnect = 0.;
					  BigInteger cartesianProduct = 1;
					  nbSol=(wcsp->numberOfConnectedConstraints() == 0)?(wcsp->cartProd(cartesianProduct),cartesianProduct):sharpBTD(start);
					  if(ToulBar2::approximateCountingBTD && nbSol>0. && td->getRoot()->getNbVars()==0)
						{ //if there are several parts
						  approximate(nbSol,td);
						}
					  // computation of maximal separator size
					  for(int i=0;i<td->getNbOfClusters();i++)
						{
						  if(td->getCluster(i)->sepSize()>tailleSep)
							tailleSep=td->getCluster(i)->sepSize();
						}
					}
				  else
					ub = recursiveSolve(start, wcsp->getLb(), ub);
				  break;
				case 2:case 3:
				  russianDollSearch(start, ub);
				  ub = start->getLbRDS();
				  break;
				default:
				  cerr << "Unknown search method B" << ToulBar2::btdMode << endl;
				  exit(EXIT_FAILURE);
				}
				wcsp->setUb(ub);
				if(ToulBar2::debug) start->printStatsRec();
			  } else recursiveSolve();
			}
		  } catch (NbBacktracksOut) {
			nbbacktracksout = true;
		  }
		  store->restore(storedepth);
		} while (nbbacktracksout);
    } catch (Contradiction) {
        wcsp->whenContradiction();
    }
    if(ToulBar2::isZ) {
	  if (ToulBar2::verbose >= 1) cout << "NegativeShiftingCost= " << wcsp->getNegativeLb() << endl;
	  if (ToulBar2::uai) {
		if (ToulBar2::uai_firstoutput) ToulBar2::uai_firstoutput = false;
		else ToulBar2::solution_file << "-BEGIN-" << endl;
		ToulBar2::solution_file << "1" << endl;
		ToulBar2::solution_file << (ToulBar2::logZ + ToulBar2::markov_log) << endl;
		ToulBar2::solution_file.flush();
	  }
	  cout << "Log10(Z)= ";
	  cout <<  (ToulBar2::logZ + ToulBar2::markov_log) << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes and " << cpuTime() - ToulBar2::startCpuTime << " seconds" << endl;
	  return true;
    }
	if(ToulBar2::allSolutions) {
	  if(ToulBar2::approximateCountingBTD)
		cout << "Number of solutions    : ~= " << nbSol << endl;
	  else
		cout << "Number of solutions    : =  " << nbSol << endl;
	  cout << "Number of #goods       :    " << nbSGoods << endl;
	  cout << "Number of used #goods  :    " << nbSGoodsUse << endl;
	  cout << "Size of sep            :    " << tailleSep << endl;
	  cout << "Time                   :    " << cpuTime() - ToulBar2::startCpuTime << " seconds" << endl;
	  cout << "... in " <<nbBacktracks << " backtracks and " << nbNodes << " nodes" << endl;
	  return true;
	}
//  store->restore();         // see above for store->store()

   	if(ToulBar2::vac) wcsp->printVACStat();

    if (wcsp->getUb() < initialUpperBound) {
	  if(!ToulBar2::uai && !ToulBar2::xmlflag) {
	        if(ToulBar2::haplotype) cout <<  "\nOptimum: " <<  wcsp->getUb() << " log10like: " << ToulBar2::haplotype->Cost2LogLike(wcsp->getUb())<< " logProb: " << ToulBar2::haplotype->Cost2Prob( wcsp->getUb()) << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes, and " << cpuTime() - ToulBar2::startCpuTime << " seconds." << endl;
    		else if(!ToulBar2::bayesian) cout << "Optimum: " << wcsp->getUb() << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << " and " << cpuTime() - ToulBar2::startCpuTime << " seconds." << endl;
			else cout << "Optimum: " << wcsp->getUb() << " log10like: " << wcsp->Cost2LogLike(wcsp->getUb()) + ToulBar2::markov_log << " prob: " << wcsp->Cost2Prob( wcsp->getUb() ) * Exp10(ToulBar2::markov_log) << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << " and " << cpuTime() - ToulBar2::startCpuTime << " seconds." << endl;
    	} else {
    		if(ToulBar2::xmlflag) ((WCSP*)wcsp)->solution_XML(true);
    		else if(ToulBar2::uai && !ToulBar2::isZ) {
			  ((WCSP*)wcsp)->solution_UAI(wcsp->getUb(), NULL, true);
			  cout << "Optimum: " << wcsp->getUb() << " log10like: " << wcsp->Cost2LogLike(wcsp->getUb()) + ToulBar2::markov_log << " prob: " << wcsp->Cost2Prob( wcsp->getUb() ) * Exp10(ToulBar2::markov_log) << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << " and " << cpuTime() - ToulBar2::startCpuTime << " seconds." << endl;
			}
    	}
        return true;
    } else {
	    cout << "No solution in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << " and " << cpuTime() - ToulBar2::startCpuTime << " seconds." << endl;
        return false;
    }
}


void Solver::approximate(BigInteger &nbsol, TreeDecomposition* td)
{
	BigInteger cartesianProduct = 1;
	wcsp->cartProd(cartesianProduct);
	for(map<int, BigInteger>:: iterator it = ubSol.begin(); it != ubSol.end(); ++it){
		(it->second) *= cartesianProduct;
	}
	BigInteger nbSolInter = nbsol*cartesianProduct;
	BigInteger subCartesianProduct = 1.;
	for( int i = 0; i< td->getNbOfClusters();i++)
	{
		BigInteger ssCartProd = 1.;
		if((td->getCluster(i)->getParent()!=NULL) && (td->getCluster(i)->getParent()->getParent()==NULL))
		{
			/* on considere seulement les clusters fils de la racine */
			Cluster * c = td->getCluster(i);
			c->cartProduct(ssCartProd);
			subCartesianProduct *= ssCartProd;
			(ubSol.find(c->getPart())->second) /= ssCartProd;

		}
	}
	nbsol = (nbSolInter/subCartesianProduct);
	if(nbsol < 1)
		nbsol = 1;
	// the minimum upper bound of solutions number
	cout << "\nCartesian product \t\t   :    " << cartesianProduct << endl;
	BigInteger minUBsol = cartesianProduct;
	for(map<int, BigInteger> :: iterator it = ubSol.begin(); it != ubSol.end(); ++it)
	{
		if(it->second < minUBsol) minUBsol = it->second;
	}
	cout << "Upper bound of number of solutions : <= " << minUBsol << endl;

}

void solution_symmax2sat(int wcspIndex)
{
  IsASolution = true;
  for (unsigned int i=0; i<CurrentWeightedCSP->numberOfVariables(); i++) {
	assert(CurrentWeightedCSP->assigned(i));
	if (CurrentWeightedCSP->getValue(i) == 0) {
	  CurrentSolution[i] = 1;
	} else {
	  CurrentSolution[i] = -1;
	}
  }
}

// Maximize h' W h where W is expressed by all its
// non-zero half squared matrix costs (can be positive or negative, with posx <= posy)
// notice that costs for posx <> posy are multiplied by 2 by this method

// convention: h = 1 <=> x = 0 and h = -1 <=> x = 1

// warning! does not allow infinite costs (no forbidden assignments)

// returns true if at least one solution has been found (array sol being filled with the best solution)
bool Solver::solve_symmax2sat(int n, int m, int *posx, int *posy, double *cost, int *sol)
{
  if (n == 0 || m == 0) return true;
  ToulBar2::setvalue = NULL;

  // create Boolean variables
  for (int i=0; i<n; i++) {
	wcsp->makeEnumeratedVariable(to_string(i), 0, 1);
  }

  vector<Cost> unaryCosts0(n, 0);
  vector<Cost> unaryCosts1(n, 0);

  // find total cost
  Double sumcost = 0.;
  for (int e=0; e<m; e++) {
	sumcost += 2. * abs(cost[e]);
  }
  Double multiplier = ((Double) MAX_COST) / sumcost;
  multiplier /= MEDIUM_COST;

  // create weighted binary clauses
  for (int e=0; e<m; e++) {
	if (posx[e] != posy[e]) {
	  vector<Cost> costs(4, 0);
	  if (cost[e] > 0) {
		costs[1] = (Cost) (multiplier * 2. * cost[e]);
		costs[2] = costs[1];
	  } else {
		costs[0] = (Cost) (multiplier * -2. * cost[e]);
		costs[3] = costs[0];
	  }
	  wcsp->postBinaryConstraint(posx[e] - 1, posy[e] - 1, costs);
	} else {
	  if (cost[e] > 0) {
		unaryCosts1[posx[e] - 1] += (Cost) (multiplier * cost[e]);
	  } else {
		unaryCosts0[posx[e] - 1] += (Cost) (multiplier * -cost[e]);
	  }
	}
  }

  wcsp->sortConstraints();

  // create weighted unary clauses
  for (int i=0; i<n; i++) {
	if (unaryCosts0[i] > 0 || unaryCosts1[i] > 0) {
	  vector<Cost> costs(2, 0);
	  costs[0] = unaryCosts0[i];
	  costs[1] = unaryCosts1[i];
	  wcsp->postUnary(i, costs);
    }
  }

  wcsp->histogram();  

  cout << "Read " << n << " variables, with " << 2 << " values at most, and " << m << " cost functions." << endl;
  // dump_wcsp("mydebug.wcsp", true);

  // solve using BTD exploiting a lexicographic elimination order with a path decomposition

  ToulBar2::btdMode = 3;
  ToulBar2::minProperVarSize = 4;
  ToulBar2::elimDegree_preprocessing = 12; // Prefer variable elimination than search (do not impose a limit on maximum separator size)

  IsASolution = false;
  CurrentSolution = sol;
  CurrentWeightedCSP = wcsp;
  ToulBar2::newsolution = solution_symmax2sat;
  solve();
  return IsASolution;
}

/// \brief interface for Fortran call
/// \code
/// integer     ,dimension(sW),intent(out)        :: H
/// integer     ,dimension(:)    ,allocatable       :: posx,posy ! On exit dimension value is m
/// real(kind=dp),dimension(:)   ,allocatable   :: cost
/// logical :: ok
/// allocate (posx(sW*sW),posy(sW*sW),cost(sW*sW))
/// ret = solvesymmax2sat_(n,m,posx,posy,cost,H)
/// ok = ( ret /= 0 )
/// deallocate(posx,posy,cost)
/// \endcode
int solvesymmax2sat_(int *n, int *m, int *posx, int *posy, double *cost, int *sol)
{return solveSymMax2SAT(*n,*m,posx,posy,cost,sol);}

int solveSymMax2SAT(int n, int m, int *posx, int *posy, double *cost, int *sol)
{
  // select verbosity during search
  ToulBar2::verbose = -1;

  initCosts(MAX_COST);
  Solver solver(STORE_SIZE, MAX_COST);

  ToulBar2::startCpuTime = cpuTime();
  return solver.solve_symmax2sat(n , m, posx, posy, cost, sol);
}
