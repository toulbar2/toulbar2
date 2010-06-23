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

extern void setvalue(int wcspId, int varIndex, Value value);


/*
 * Solver constructors
 *
 */

Solver *Solver::currentSolver = NULL;

Solver::Solver(int storeSize, Cost initUpperBound) : store(NULL), nbNodes(0), nbBacktracks(0), nbBacktracksLimit(-1), wcsp(NULL),
                                                     unassignedVars(NULL), lastConflictVar(-1), nbSol(0.), nbSGoods(0), nbSGoodsUse(0)
{
    store = new Store(storeSize);
    wcsp = WeightedCSP::makeWeightedCSP(store, initUpperBound);
}

Solver::~Solver()
{
    delete store;
    delete wcsp;
    delete unassignedVars;
}

void Solver::initVarHeuristic()
{
    unassignedVars = new BTList<Value>(&store->storeDomain);
    allVars = new DLink<Value>[wcsp->numberOfVariables()];
    for (unsigned int i=0; i<wcsp->numberOfVariables(); i++) {
        allVars[i].content = i;
        unassignedVars->push_back(&allVars[i], false);
        if (wcsp->assigned(i)) unassignedVars->erase(&allVars[i], false);
    }
}

void Solver::read_wcsp(const char *fileName)
{
    ToulBar2::setvalue = NULL;
    wcsp->read_wcsp(fileName);
	initVarHeuristic();
    ToulBar2::setvalue = setvalue;
}

void Solver::read_random(int n, int m, vector<int>& p, int seed, bool forceSubModular)
{
    ToulBar2::setvalue = NULL;
    wcsp->read_random(n,m,p,seed, forceSubModular);
	initVarHeuristic();
    ToulBar2::setvalue = setvalue;
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

    int i = 0;
    while (!file.eof()) {
        if ((unsigned int) i >= wcsp->numberOfVariables()) break;
        int value = 0;
        file >> value;
        if (wcsp->unassigned(i)) {
		  assign(i, value);
		  // side-effect: remember last solution
		  wcsp->setBestValue(i, value);
        } else {
		    if (wcsp->getValue(i) != value) {
			  THROWCONTRADICTION;
			} else {
			  wcsp->setBestValue(i, value); // side-effect: remember last solution
			}
        }
        i++;
    }
    cout << " Solution cost: [" << wcsp->getLb() << "," << wcsp->getUb() << "] (nb. of unassigned variables: " << wcsp->numberOfUnassignedVariables() << ")" << endl;
	if (ToulBar2::btdMode>=2) wcsp->updateUb(wcsp->getLb()+UNIT_COST);
    store->restore(depth);
}

void Solver::parse_solution(const char *certificate)
{
    currentSolver = this;

    if (ToulBar2::btdMode>=2) wcsp->propagate();

	//  int depth = store->getDepth();
	//    store->store();

    certificate = index(certificate,',');
    if (certificate) certificate++;
    
    int var;
	Value value;
	int items;
    while ((certificate != NULL) && (certificate[0] != 0)) {

        items = sscanf(certificate,"%d=%d",&var,&value);
        
        if ((items != 2) || ((unsigned int)var >= wcsp->numberOfVariables())) {
             cerr << "Certificate " << certificate << " incorrect!" << endl;
             exit(EXIT_FAILURE);
        }
        certificate = index(certificate,',');
        if (certificate) certificate++;

        if (wcsp->unassigned(var)) {
          assign(var, value);
          // side-effect: remember last solution
          wcsp->setBestValue(var, value);
        } else {
		  if (wcsp->getValue(var) != value) {
			THROWCONTRADICTION;
		  } else {
			wcsp->setBestValue(var, value); // side-effect: remember last solution
		  }
        }
    }
    cout << " Solution cost: [" << wcsp->getLb() << "," << wcsp->getUb() << "] (nb. of unassigned variables: " << wcsp->numberOfUnassignedVariables() << ")" << endl;
    
    if (ToulBar2::btdMode>=2) wcsp->updateUb(wcsp->getLb()+UNIT_COST);
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
    assert(wcspId == 0);
    if(!Solver::currentSolver->allVars[varIndex].removed) {
	  Solver::currentSolver->unassignedVars->erase(&Solver::currentSolver->allVars[varIndex], true);
	}
}

/*
 * Variable ordering heuristics
 *
 */

int Solver::getNextUnassignedVar()
{
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
        } else if (heuristic < best + 1./1000001. && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) {
		 ties[nbties] = *iter;
		 nbties++;
	   }
    }
    if (nbties>1) return ties[myrand()%nbties];
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
       } else if (heuristic < best + 1./1000001. && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) {
		 ties[nbties] = *iter;
		 nbties++;
	   }
   }
   if (nbties>1) return ties[myrand()%nbties];
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

void Solver::increase(int varIndex, Value value)
{
    wcsp->enforceUb();
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
    wcsp->enforceUb();
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
    wcsp->enforceUb();
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
    wcsp->enforceUb();
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
    if (dichotomic) {
      middle = (wcsp->getInf(varIndex) + wcsp->getSup(varIndex)) / 2;
      if (value <= middle) increasing = true;
      else increasing = false;
    }
    try {
        store->store();
        lastConflictVar = varIndex;
        if (dichotomic) {
    	  if (increasing) decrease(varIndex, middle); else increase(varIndex, middle+1);
    	} else assign(varIndex, value);
        lastConflictVar = -1;
        recursiveSolve();
    } catch (Contradiction) {
        wcsp->whenContradiction();
    }
    store->restore();
    nbBacktracks++;
	if (ToulBar2::restart && nbBacktracks > nbBacktracksLimit) throw NbBacktracksOut();

    if (dichotomic) {
      if (increasing) increase(varIndex, middle+1); else decrease(varIndex, middle);
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
        nbBacktracks++;
		if (ToulBar2::restart && nbBacktracks > nbBacktracksLimit) throw NbBacktracksOut();
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
    nbBacktracks++;
	if (ToulBar2::restart && nbBacktracks > nbBacktracksLimit) throw NbBacktracksOut();
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
    nbBacktracks++;
	if (ToulBar2::restart && nbBacktracks > nbBacktracksLimit) throw NbBacktracksOut();
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
    nbBacktracks++;
	if (ToulBar2::restart && nbBacktracks > nbBacktracksLimit) throw NbBacktracksOut();
}

void Solver::singletonConsistency()
{
    bool deadend;
    bool done = false;
    while(!done) {
    	done = true;
	    for (unsigned int varIndex = 0; varIndex < wcsp->numberOfVariables(); varIndex++) {
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
    if (!ToulBar2::allSolutions) wcsp->updateUb(wcsp->getLb());
	else if (!ToulBar2::btdMode) nbSol += 1.;

  	if(!ToulBar2::xmlflag && !ToulBar2::uai && (!ToulBar2::allSolutions || ToulBar2::debug)) {
  		if(ToulBar2::haplotype) cout <<  "***New solution: " <<  wcsp->getLb() << " log10like: " << ToulBar2::haplotype->Cost2LogLike(wcsp->getLb())<< " logProb: " << ToulBar2::haplotype->Cost2Prob( wcsp->getLb()) << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << store->getDepth() << ")" << endl;
		else if(!ToulBar2::bayesian) cout << "New solution: " <<  wcsp->getLb() << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << store->getDepth() << ")" << endl;
		else cout << "New solution: " <<  wcsp->getLb() << " log10like: " << wcsp->Cost2LogLike(wcsp->getLb()) << " prob: " << wcsp->Cost2Prob( wcsp->getLb() ) << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << store->getDepth() << ")" << endl;
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
	if(ToulBar2::uai) {
	  ((WCSP*)wcsp)->solution_UAI(wcsp->getLb());
	}

	if (ToulBar2::newsolution) (*ToulBar2::newsolution)(wcsp->getIndex());
}

void Solver::recursiveSolve()
{
	int varIndex = -1;
	if (ToulBar2::bep) varIndex = getMostUrgent();
	else if(ToulBar2::weightedDegree && ToulBar2::lastConflict) varIndex = ((ToulBar2::restart)?getVarMinDomainDivMaxWeightedDegreeLastConflictRandomized():getVarMinDomainDivMaxWeightedDegreeLastConflict());
	else if(ToulBar2::lastConflict) varIndex = ((ToulBar2::restart)?getVarMinDomainDivMaxDegreeLastConflictRandomized():getVarMinDomainDivMaxDegreeLastConflict());
	else if(ToulBar2::weightedDegree) varIndex = ((ToulBar2::restart)?getVarMinDomainDivMaxWeightedDegreeRandomized():getVarMinDomainDivMaxWeightedDegree());
	else varIndex = ((ToulBar2::restart)?getVarMinDomainDivMaxDegreeRandomized():getVarMinDomainDivMaxDegree());
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
	else if(ToulBar2::weightedDegree && ToulBar2::lastConflict) varIndex = ((ToulBar2::restart)?getVarMinDomainDivMaxWeightedDegreeLastConflictRandomized():getVarMinDomainDivMaxWeightedDegreeLastConflict());
	else if(ToulBar2::lastConflict) varIndex =  ((ToulBar2::restart)?getVarMinDomainDivMaxDegreeLastConflictRandomized():getVarMinDomainDivMaxDegreeLastConflict());
	else if(ToulBar2::weightedDegree) varIndex = ((ToulBar2::restart)?getVarMinDomainDivMaxWeightedDegreeRandomized():getVarMinDomainDivMaxWeightedDegree());
	else varIndex = ((ToulBar2::restart)?getVarMinDomainDivMaxDegreeRandomized():getVarMinDomainDivMaxDegree());
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

    try {
//        store->store();       // if uncomment then solve() does not change the problem but all preprocessing operations will allocate in backtrackable memory
        wcsp->decreaseUb(initialUpperBound);

        wcsp->propagate();                // initial propagation
        wcsp->preprocessing();            // preprocessing after initial propagation
        if (ToulBar2::verbose >= 1||(!ToulBar2::xmlflag && !ToulBar2::uai)) cout << wcsp->numberOfUnassignedVariables() << " unassigned variables, " << wcsp->getDomainSizeSum() << " values in all current domains and " << wcsp->numberOfConnectedConstraints() << " constraints." << endl;

        if (ToulBar2::singletonConsistency) {
		  singletonConsistency();
		  wcsp->propagate();
		}

	    if (ToulBar2::btdMode) {
	    	if(wcsp->numberOfUnassignedVariables() == 0)
	    		ToulBar2::approximateCountingBTD = 0;
	    	wcsp->buildTreeDecomposition();
	    }

		if (ToulBar2::dumpWCSP) {dump_wcsp("problem.wcsp",false); cout << "end." << endl; exit(0);}

		if (ToulBar2::restart) {
		  nbBacktracksLimit = 1;
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
		  if (ToulBar2::restart) {
			nbbacktracksout = false;
			nbrestart++;
			//			currentNbBacktracksLimit = luby(nbrestart);
			currentNbBacktracksLimit *= 4;
			if (currentNbBacktracksLimit > nbBacktracksLimitTop || IsASolution) {
			  nbBacktracksLimitTop = currentNbBacktracksLimit;
			  currentNbBacktracksLimit = 1;
			}
			if (!IsASolution && nbNodes > ToulBar2::restart) {
			  nbBacktracksLimit = LONGLONG_MAX;
			  ToulBar2::restart = 0;
			  cout << "****** Restart " << nbrestart << " with no backtrack limit and UB=" << UpperBound << " ****** (" << nbNodes << " nodes)" << endl;
			} else {
			  nbBacktracksLimit = nbBacktracks + currentNbBacktracksLimit;
			  cout << "****** Restart " << nbrestart << " with " << currentNbBacktracksLimit << " backtracks max and UB=" << UpperBound << " ****** (" << nbNodes << " nodes)" << endl;
			  store->store();
			}
			IsASolution = false;
			wcsp->setUb(UpperBound);
			wcsp->enforceUb();
			wcsp->propagate();
		  }
		  try {
			if (ToulBar2::lds) {
			  int discrepancy = 0;
			  do {
				if (discrepancy > ToulBar2::lds) cout << "--- [" << store->getDepth() << "] Search with no discrepancy limit --- (" << nbNodes << " nodes)" << endl;
				else cout << "--- [" << store->getDepth() << "] LDS " << discrepancy << " --- (" << nbNodes << " nodes)" << endl;
				ToulBar2::limited = false;
				try {
				  store->store();
				  if (discrepancy > ToulBar2::lds) recursiveSolve(); else recursiveSolveLDS(discrepancy);
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
					  nbSol=sharpBTD(start);
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
			store->restore(storedepth);
		  }
		} while (nbbacktracksout);
    } catch (Contradiction) {
        wcsp->whenContradiction();
    }
    if(ToulBar2::allSolutions)
    {
    	if(ToulBar2::approximateCountingBTD)
    		cout << "Number of solutions    : ~= " << nbSol << endl;
    	else
    		cout << "Number of solutions    : =  " << nbSol << endl;
    	cout << "Number of #goods       :    " << nbSGoods << endl;
    	cout << "Number of used #goods  :    " << nbSGoodsUse << endl;
    	cout << "Size of sep            :    " << tailleSep << endl;
    	cout << "Time                   :    " << cpuTime() - ToulBar2::startCpuTime << " seconds" << endl;
    	cout << "... with " <<nbBacktracks << " backtracks and " << nbNodes << " nodes" << endl;
    	return true;
    }
//  store->restore();         // see above for store->store()

   	if(ToulBar2::vac) wcsp->printVACStat();

    if (wcsp->getUb() < initialUpperBound) {
    	if(!ToulBar2::uai && !ToulBar2::xmlflag) {
	        if(ToulBar2::haplotype) cout <<  "\nOptimum: " <<  wcsp->getUb() << " log10like: " << ToulBar2::haplotype->Cost2LogLike(wcsp->getUb())<< " logProb: " << ToulBar2::haplotype->Cost2Prob( wcsp->getLb()) << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes, and " << cpuTime() - ToulBar2::startCpuTime << " seconds." << endl;
    		else if(!ToulBar2::bayesian) cout << "Optimum: " << wcsp->getUb() << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << " and " << cpuTime() - ToulBar2::startCpuTime << " seconds." << endl;
			else cout << "Optimum: " << wcsp->getUb() << " log10like: " << wcsp->Cost2LogLike(wcsp->getUb()) << " prob: " << wcsp->Cost2Prob( wcsp->getUb() ) << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << " and " << cpuTime() - ToulBar2::startCpuTime << " seconds." << endl;
    	} else {
    		if(ToulBar2::xmlflag) ((WCSP*)wcsp)->solution_XML(true);
    		else if(ToulBar2::uai) ((WCSP*)wcsp)->solution_UAI(wcsp->getUb());
    	}
        return true;
    } else {
	    if(ToulBar2::debug || (!ToulBar2::uai && !ToulBar2::xmlflag)) {
 	       cout << "No solution in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << " and " << cpuTime() - ToulBar2::startCpuTime << " seconds." << endl;
    	}
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
  double sumcost = 0.;
  for (int e=0; e<m; e++) {
	sumcost += 2. * abs(cost[e]);
  }
  double multiplier = MAX_COST / sumcost;
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

  // create weighted unary clauses
  for (int i=0; i<n; i++) {
	if (unaryCosts0[i] > 0 || unaryCosts1[i] > 0) {
	  vector<Cost> costs(2, 0);
	  costs[0] = unaryCosts0[i];
	  costs[1] = unaryCosts1[i];
	  wcsp->postUnary(i, costs);
    }
  }

  cout << "Read " << n << " variables, with " << 2 << " values at most, and " << m << " constraints." << endl;

  // special data structure to be initialized for variable ordering heuristics
  initVarHeuristic();
  ToulBar2::setvalue = setvalue;

  // solve using BTD exploiting a lexicographic elimination order, a path decomposition, and only small separators of size ToulBar2::smallSeparatorSize

  ToulBar2::btdMode = 3;
  ToulBar2::maxSeparatorSize = ToulBar2::smallSeparatorSize;

  IsASolution = false;
  CurrentSolution = sol;
  CurrentWeightedCSP = wcsp;
  ToulBar2::newsolution = solution_symmax2sat;
  solve();
  return IsASolution;
}

int solveSymMax2SAT(int n, int m, int *posx, int *posy, double *cost, int *sol)
{
  // select verbosity during search
  ToulBar2::verbose = 0;
  //  ToulBar2::elimDegree = -1;

  initCosts(MAX_COST);
  Solver solver(STORE_SIZE, MAX_COST);

  ToulBar2::startCpuTime = cpuTime();
  return solver.solve_symmax2sat(n , m, posx, posy, cost, sol);
}
