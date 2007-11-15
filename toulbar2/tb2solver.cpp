/*
 * **************** Generic solver *******************
 * 
 */

#include "tb2solver.hpp"
#include "tb2domain.hpp"
#include "tb2pedigree.hpp"
#include "tb2bep.hpp"

extern void setvalue(int wcspId, int varIndex, Value value);

/*
 * Solver constructors
 * 
 */

Solver *Solver::currentSolver = NULL;

Solver::Solver(int storeSize, Cost initUpperBound) : store(NULL), nbNodes(0), nbBacktracks(0), wcsp(NULL), 
                                                     unassignedVars(NULL), lastConflictVar(-1)
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

void Solver::read_wcsp(const char *fileName)
{
    ToulBar2::setvalue = NULL;
    wcsp->read_wcsp(fileName);
    unassignedVars = new BTList<Value>(&store->storeDomain);
    allVars = new DLink<Value>[wcsp->numberOfVariables()];
    for (unsigned int i=0; i<wcsp->numberOfVariables(); i++) {
        allVars[i].content = i;
        unassignedVars->push_back(&allVars[i], false);
        if (wcsp->assigned(i)) unassignedVars->erase(&allVars[i], false);
    }
    ToulBar2::setvalue = setvalue;
}

void Solver::read_random(int n, int m, vector<int>& p, int seed, bool forceSubModular)
{
    ToulBar2::setvalue = NULL;
    wcsp->read_random(n,m,p,seed, forceSubModular);
    unassignedVars = new BTList<Value>(&store->storeDomain);
    allVars = new DLink<Value>[wcsp->numberOfVariables()];
    for (unsigned int i=0; i<wcsp->numberOfVariables(); i++) {
        allVars[i].content = i;
        unassignedVars->push_back(&allVars[i], false);
        if (wcsp->assigned(i)) unassignedVars->erase(&allVars[i], false);
    }
    ToulBar2::setvalue = setvalue;
}


void Solver::read_solution(const char *filename)
{
    currentSolver = this;
    store->store();

    // open the file
    ifstream file(filename);
    if (!file) {
        cerr << "Solution file " << filename << " not found!" << endl;
        exit(EXIT_FAILURE);
    }
    
    int i = 0;
    while (file) {
        if ((unsigned int) i >= wcsp->numberOfVariables()) break;
        int value = 0;
        file >> value;
        if (wcsp->unassigned(i)) {
		  assign(i, value);
        } else {
            if (wcsp->getValue(i) != value) THROWCONTRADICTION;
        }
        i++;
    }
    cout << " Solution cost: [" << wcsp->getLb() << "," << wcsp->getUb() << "] (nb. of unassigned variables: " << wcsp->numberOfUnassignedVariables() << ")" << endl;
}

void Solver::dump_wcsp(const char *fileName)
{
    ofstream pb(fileName);
    if (pb) wcsp->dump(pb);
}

/*
 * Link between solver and wcsp: maintain a backtrackable list of unassigned variable indexes
 * 
 */

void setvalue(int wcspId, int varIndex, Value value)
{
    assert(wcspId == 0);
    assert(!Solver::currentSolver->allVars[varIndex].removed);
    Solver::currentSolver->unassignedVars->erase(&Solver::currentSolver->allVars[varIndex], true);
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
        // remove following "+1" when isolated variables are automatically assigned
        double heuristic = (double) wcsp->getDomainSize(*iter) / (wcsp->getDegree(*iter) + 1);
        if (varIndex < 0 || heuristic < best - 1./100001.
            || (heuristic < best + 1./100001. && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
        }
    }
    return varIndex;   
}

int Solver::getVarMinDomainDivMaxDegreeLastConflict()
{
    if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar)) return lastConflictVar;
    int varIndex = -1;
    Cost worstUnaryCost = MIN_COST;
    double best = MAX_VAL - MIN_VAL;
    
    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
        // remove following "+1" when isolated variables are automatically assigned
        double heuristic = (double) wcsp->getDomainSize(*iter) / (wcsp->getDegree(*iter) + 1);
        if (varIndex < 0 || heuristic < best - 1./100001.
            || (heuristic < best + 1./100001. && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
        }
    }
    return varIndex;
}

int Solver::getVarMinDomainDivMaxWeightedDegree()
{
    int varIndex = -1;
    double best = MAX_VAL - MIN_VAL;
    
    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
	  double size = wcsp->getDomainSize(*iter);
	  double heuristic = size / (wcsp->getWeightedDegree(*iter) + to_double(wcsp->getMaxUnaryCost(*iter)));
	  if (varIndex  < 0 || heuristic < best) {
		best = heuristic;
		varIndex = *iter;
	  }
    }
    return varIndex;
}

int Solver::getVarMinDomainDivMaxWeightedDegreeLastConflict()
{
    if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar)) return lastConflictVar;
    int varIndex = -1;
    double best = MAX_VAL - MIN_VAL;
    
    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
	  double size = wcsp->getDomainSize(*iter);
	  double heuristic = size / (wcsp->getWeightedDegree(*iter) + to_double(wcsp->getMaxUnaryCost(*iter)));
	  if (varIndex  < 0 || heuristic < best) {
		best = heuristic;
		varIndex = *iter;
	  }
    }
    return varIndex;
}

// First experiments are not convincing. Wait to try on CSP benchmark of [Lecoutre ECAI-04]
//int Solver::getVarMinDomainDivWeightedDegreeLastConflict()
//{
//    if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar)) return lastConflictVar;
//    int varIndex = -1;
//    Cost worstUnaryCost = MIN_COST;
//    double best = MAX_VAL - MIN_VAL;
//    
//    for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
//        // remove following "+1" when isolated variables are automatically assigned
//        double heuristic = (double) wcsp->getDomainSize(*iter) / (wcsp->getWeightedDegree(*iter) + 1);
//        if (varIndex < 0 || heuristic < best - 1./100001.
//            || (heuristic < best + 1./100001. && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
//            best = heuristic;
//            varIndex = *iter;
//            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
//        }
//    }
//    return varIndex;
//}

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
        cout << "[" << store->getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum() << "] Try " << wcsp->getName(varIndex) << " >= " << value << " (s:" << wcsp->getSupport(varIndex) << ")" << endl;
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
        cout << "[" << store->getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum() << "] Try " << wcsp->getName(varIndex) << " <= " << value << " (s:" << wcsp->getSupport(varIndex) << ")" << endl;
    }
    wcsp->decrease(varIndex, value);
    wcsp->propagate();
}

void Solver::assign(int varIndex, Value value)
{
    wcsp->enforceUb();
    nbNodes++;
    if (ToulBar2::verbose >= 1) {
        if (ToulBar2::verbose >= 2) cout << *wcsp;
        cout << "[" << store->getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum() << "] Try " << wcsp->getName(varIndex) << " == " << value << endl;
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
        cout << "[" << store->getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum() << "] Try " << wcsp->getName(varIndex) << " != " << value << endl;
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
    
    if(!ToulBar2::bayesian) cout << "New solution: " <<  wcsp->getLb() << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << store->getDepth() << ")" << endl;
	else cout << "New solution: " <<  wcsp->getLb() << " log10like: " << wcsp->Cost2LogLike(wcsp->getLb()) << " prob: " << wcsp->Cost2Prob( wcsp->getLb() ) << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << store->getDepth() << ")" << endl;
    wcsp->restoreSolution();
    if (ToulBar2::showSolutions) {
        if (ToulBar2::verbose >= 2) cout << *wcsp << endl;
        for (unsigned int i=0; i<wcsp->numberOfVariables(); i++) {
            cout << " ";
            if (ToulBar2::pedigree) {
                cout <<  wcsp->getName(i) << ":";
                ToulBar2::pedigree->printGenotype(cout, wcsp->getValue(i));
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
}

void Solver::recursiveSolve()
{		
	int varIndex = -1;
	if (ToulBar2::bep) varIndex = getMostUrgent();
	else if(ToulBar2::weightedDegree && ToulBar2::lastConflict) varIndex = getVarMinDomainDivMaxWeightedDegreeLastConflict();
	else if(ToulBar2::lastConflict) varIndex = getVarMinDomainDivMaxDegreeLastConflict();
	else if(ToulBar2::weightedDegree) varIndex = getVarMinDomainDivMaxWeightedDegree();
	else varIndex = getVarMinDomainDivMaxDegree();
    if (varIndex >= 0) {
	    if (ToulBar2::bep) scheduleOrPostpone(varIndex);
        else if (wcsp->enumerated(varIndex)) {
            if (ToulBar2::binaryBranching) {
			  assert(wcsp->canbe(varIndex, wcsp->getSupport(varIndex)));
			  binaryChoicePoint(varIndex, wcsp->getSupport(varIndex));
			  //			    assert(wcsp->canbe(varIndex, wcsp->getMaxUnaryCostValue(varIndex)));
			  //			    binaryChoicePoint(varIndex, wcsp->getMaxUnaryCostValue(varIndex));
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
	else if(ToulBar2::weightedDegree && ToulBar2::lastConflict) varIndex = getVarMinDomainDivMaxWeightedDegreeLastConflict();
	else if(ToulBar2::lastConflict) varIndex = getVarMinDomainDivMaxDegreeLastConflict();
	else if(ToulBar2::weightedDegree) varIndex = getVarMinDomainDivMaxWeightedDegree();
	else varIndex = getVarMinDomainDivMaxDegree();
    if (varIndex >= 0) {
	    if (ToulBar2::bep) scheduleOrPostpone(varIndex);
        else if (wcsp->enumerated(varIndex)) {
            if (ToulBar2::binaryBranching) {
			  assert(wcsp->canbe(varIndex, wcsp->getSupport(varIndex)));
			  binaryChoicePointLDS(varIndex, wcsp->getSupport(varIndex), discrepancy);
			  //			    assert(wcsp->canbe(varIndex, wcsp->getMaxUnaryCostValue(varIndex)));
			  //			    binaryChoicePointLDS(varIndex, wcsp->getMaxUnaryCostValue(varIndex), discrepancy);
            } else {
                narySortedChoicePointLDS(varIndex, discrepancy);
            }
        } else {
            binaryChoicePointLDS(varIndex, wcsp->getInf(varIndex), discrepancy);
        }
    } else newSolution();
}

bool Solver::solve()
{
    currentSolver = this;
    Cost initialUpperBound = wcsp->getUb();
    nbBacktracks = 0;
    nbNodes = 0;
    try {
//        store->store();       // if uncomment then solve() does not change the problem but all preprocessing operations will allocate in backtrackable memory
        wcsp->decreaseUb(initialUpperBound);
           
        wcsp->propagate();                // initial propagation
        wcsp->preprocessing();            // preprocessing after initial propagation
        cout << wcsp->numberOfUnassignedVariables() << " unassigned variables, " << wcsp->getDomainSizeSum() << " values in all current domains and " << wcsp->numberOfConnectedConstraints() << " constraints." << endl;

        if (ToulBar2::singletonConsistency) singletonConsistency();

        if (ToulBar2::lds) {
            int discrepancy = 0;
            do {
                cout << "--- [" << store->getDepth() << "] LDS " << discrepancy << " ---" << endl;
                ToulBar2::limited = false;
                try {
                    store->store();
                    recursiveSolveLDS(discrepancy);
                } catch (Contradiction) {
                    wcsp->whenContradiction();
                }
                store->restore();
                if (discrepancy > 0) discrepancy *= 2;
                else discrepancy++;
            } while (ToulBar2::limited);
        } else {
            recursiveSolve();
        }
    } catch (Contradiction) {
        wcsp->whenContradiction();
    }
//    store->restore();         // see above for store->store()

   	if(ToulBar2::vac) wcsp->printVACStat();

    if (wcsp->getUb() < initialUpperBound) {
        if(!ToulBar2::bayesian) cout << "Optimum: " << wcsp->getUb() << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << endl;
		else cout << "Optimum: " << wcsp->getUb() << " log10like: " << wcsp->Cost2LogLike(wcsp->getUb()) << " prob: " << wcsp->Cost2Prob( wcsp->getUb() ) << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << endl; 
        return true;
    } else {
        cout << "No solution in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << endl;
        return false;
    }
    
}
