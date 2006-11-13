/*
 * **************** Generic solver *******************
 * 
 */

#include "tb2solver.hpp"
#include "tb2domain.hpp"
#include "tb2pedigree.hpp"

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

int Solver::getVarMinDomainDivMaxDegree()
{
    int varIndex = -1;
    double best = MAX_VAL - MIN_VAL;
    Cost worstUnaryCost = 0;
    
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
    Cost worstUnaryCost = 0;
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

// First experiments are not convincing. Wait to try on CSP benchmark of [Lecoutre ECAI-04]
//int Solver::getVarMinDomainDivWeightedDegreeLastConflict()
//{
//    if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar)) return lastConflictVar;
//    int varIndex = -1;
//    Cost worstUnaryCost = 0;
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

int Solver::getNextUnassignedVar()
{
    return (unassignedVars->empty())?-1:(*unassignedVars->begin());
}

/*
 * Choice points
 * 
 */

void Solver::binaryChoicePoint(int varIndex, Value value)
{
    assert(wcsp->unassigned(varIndex));
    assert(wcsp->canbe(varIndex,value));
    bool dichotomic = (ToulBar2::dichotomicBranching && ToulBar2::dichotomicBranchingSize < wcsp->getDomainSize(varIndex));
    Value middle;
    bool increasing;
    if (dichotomic) {
      middle = (wcsp->getInf(varIndex) + wcsp->getSup(varIndex)) / 2;
      if (value <= middle) increasing = true;
      else increasing = false;
    }
    try {
        store->store();
        if (ToulBar2::verbose >= 2) {
            cout << *wcsp;
        }
        nbNodes++;
        lastConflictVar = varIndex;
        if (dichotomic) {
	  if (increasing) {
	    if (ToulBar2::verbose >= 1) cout << "[" << store->getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum() << "] Try " << wcsp->getName(varIndex) << " <= " << middle << " (s:" << value << ")" << endl;
	    wcsp->decrease(varIndex, middle);
	  } else {
	    if (ToulBar2::verbose >= 1) cout << "[" << store->getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum() << "] Try " << wcsp->getName(varIndex) << " >= " << middle+1 << " (s:" << value << ")" << endl;
	    wcsp->increase(varIndex, middle+1);
	  }
	} else {
	  if (ToulBar2::verbose >= 1) cout << "[" << store->getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum() << "] Try " << wcsp->getName(varIndex) << " = " << value << endl;
	  wcsp->assign(varIndex, value);
	}
        wcsp->propagate();
        lastConflictVar = -1;
        recursiveSolve();
    } catch (Contradiction) {
        wcsp->whenContradiction();
    }
    store->restore();
    nbBacktracks++;
    wcsp->enforceUb();
    if (ToulBar2::verbose >= 2) {
        cout << *wcsp;
    }
    nbNodes++;
    if (dichotomic) {
      if (increasing) {
	if (ToulBar2::verbose >= 1) cout << "[" << store->getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum() << "] Refute " << wcsp->getName(varIndex) << " >= " << middle+1 << " (s:" << value << ")" << endl;
	wcsp->increase(varIndex, middle+1);
      } else {
	if (ToulBar2::verbose >= 1) cout << "[" << store->getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum() << "] Refute " << wcsp->getName(varIndex) << " <= " << middle << " (s:" << value << ")" << endl;
	wcsp->decrease(varIndex, middle);
      }
    } else {
      if (ToulBar2::verbose >= 1) cout << "[" << store->getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum() << "] Refute " << wcsp->getName(varIndex) << " != " << value << endl;
      wcsp->remove(varIndex, value);
    }
    wcsp->propagate();
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
    ValueCost sorted[size]; // replace size by MAX_DOMAIN_SIZE in case of compilation problem
    wcsp->getEnumDomainAndCost(varIndex, sorted);
    qsort(sorted, size, sizeof(ValueCost), cmpValueCost);
    for (int v = 0; wcsp->getLb() < wcsp->getUb() && v < size; v++) {
        try {
            store->store();
            nbNodes++;
            if (ToulBar2::verbose >= 2) {
                cout << *wcsp << endl;
            }
            if (ToulBar2::verbose >= 1) cout << "[" << store->getDepth() << "," << wcsp->getLb() << "," << wcsp->getUb() << "," << wcsp->getDomainSizeSum() << "] Try " << wcsp->getName(varIndex) << " = " << sorted[v].value << endl;
            wcsp->enforceUb();
            wcsp->assign(varIndex, sorted[v].value);
            wcsp->propagate();
            recursiveSolve();
        } catch (Contradiction) {
            wcsp->whenContradiction();
        }
        store->restore();
    }
    nbBacktracks++;
}

/*
 * Depth-First Branch and Bound
 * 
 */

void Solver::recursiveSolve()
{
//    int varIndex = getNextUnassignedVar();
    int varIndex = (ToulBar2::lastConflict)?getVarMinDomainDivMaxDegreeLastConflict():getVarMinDomainDivMaxDegree();
    if (varIndex >= 0) {
        if (wcsp->enumerated(varIndex)) {
            if (ToulBar2::binaryBranching) {
                assert(wcsp->canbe(varIndex, wcsp->getSupport(varIndex)));
                binaryChoicePoint(varIndex, wcsp->getSupport(varIndex));
            } else {
                narySortedChoicePoint(varIndex);
            }
        } else {
            binaryChoicePoint(varIndex, wcsp->getInf(varIndex));
        }
    } else {
        assert(unassignedVars->empty());
#ifndef NDEBUG
        bool allVarsAssigned = true;
        for (unsigned int i=0; i<wcsp->numberOfVariables(); i++) allVarsAssigned &= wcsp->assigned(i);
        assert(allVarsAssigned);
#endif
        wcsp->updateUb(wcsp->getLb());
        cout << "New solution: " <<  wcsp->getUb() << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << store->getDepth() << ")" << endl;
        wcsp->restoreSolution();
        if (ToulBar2::showSolutions) {
            if (ToulBar2::verbose >= 2) cout << *wcsp << endl;
            for (unsigned int i=0; i<wcsp->numberOfVariables(); i++) {
                cout << " ";
                if (ToulBar2::pedigree) {
                    ToulBar2::pedigree->printGenotype(cout, wcsp->getValue(i));
                } else {
                    cout << wcsp->getValue(i);
                }
            }
            cout << endl;
        }
    	if (ToulBar2::pedigree) {
    	  ToulBar2::pedigree->printCorrection((WCSP *) wcsp);
        }
        if (ToulBar2::writeSolution) {
            if (ToulBar2::pedigree) {
                ToulBar2::pedigree->save("pedigree_corrected.pre", (WCSP *) wcsp, true);
            } else {
                ofstream file("solution");
                if (!file) {
                  cerr << "Could not write file " << "solution" << endl;
                  exit(EXIT_FAILURE);
                }
                for (unsigned int i=0; i<wcsp->numberOfVariables(); i++) {
                    file << " " << wcsp->getValue(i);
                }
                file << endl;
            }
        }
    }
}

bool Solver::solve()
{
    currentSolver = this;
    Cost initialUpperBound = wcsp->getUb();
    nbBacktracks = 0;
    nbNodes = 0;
    try {
        store->store();
        wcsp->decreaseUb(initialUpperBound);
        wcsp->propagate();                // initial propagation
        wcsp->preprocessing();            // preprocessing after initial propagation
        cout << wcsp->numberOfUnassignedVariables() << " unassigned variables, " << wcsp->getDomainSizeSum() << " values in current domains and " << wcsp->numberOfConnectedConstraints() << " constraints." << endl;
        recursiveSolve();
    } catch (Contradiction) {
        wcsp->whenContradiction();
    }
    store->restore();
    if (wcsp->getUb() < initialUpperBound) {
        cout << "Optimum: " << wcsp->getUb() << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << endl;
        return true;
    } else {
        cout << "No solution in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << endl;
        return false;
    }
}
