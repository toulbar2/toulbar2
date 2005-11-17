/*
 * **************** Heuristics and search methods *******************
 * 
 */

#include "tb2wcsp.hpp"


/*
 * Variable ordering heuristics
 * 
 */

CostVariable *WCSP::getVarMinDomainDivMaxDegree()
{
    CostVariable *var = NULL;
    double best = MAX_VAL - MIN_VAL;

    for (unsigned int i=0; i<vars.size(); i++) {
        if (vars[i]->unassigned()) {
            double heuristic = (double) vars[i]->getDomainSize() / vars[i]->getDegree();
            if (var == NULL || heuristic < best) {
                best = heuristic;
                var = vars[i];
            }
        }
    }
    return var;
}

/*
 * Choice points
 * 
 */

void WCSP::store()
{
    storeData->store();
    NC.clear();
    AC.clear();
    objectiveChanged = false;
    nbNodes++;
}

void WCSP::restore()
{
    storeData->restore();
    NC.clear();
    AC.clear();
    objectiveChanged = false;
    nbNodes++;
}

void WCSP::binaryChoicePoint(CostVariable *x, Value value)
{
    assert(x->unassigned());
    assert(x->canbe(value));
    try {
        store();
        if (ToulBar2::verbose >= 2) {
            for (unsigned int i=0; i<vars.size(); i++) cout << *vars[i] << endl;
        }
        if (ToulBar2::verbose >= 1) cout << "[" << storeData->getDepth() << "] Try " << x->getName() << " = " << value << endl;
        x->assign(value);
        propagate();
        recursiveSolve();
    } catch (Contradiction) {
    }
    restore();
    getObjective()->decrease(upperBound - 1);
    nbBacktracks++;
    if (ToulBar2::verbose >= 2) {
        for (unsigned int i=0; i<vars.size(); i++) cout << *vars[i] << endl;
    }
    if (ToulBar2::verbose >= 1) cout << "[" << storeData->getDepth() << "] Try " << x->getName() << " != " << value << endl;
    x->remove(value);
    propagate();
    recursiveSolve();
}

/*
 * Depth-First Branch and Bound
 * 
 */
            
void WCSP::recursiveSolve()
{
    CostVariable *var = getVarMinDomainDivMaxDegree();
    if (var != NULL) {
        assert(var->canbe(var->getSupport()));
        binaryChoicePoint(var, var->getSupport());
    } else {
        upperBound = getLb();
        getObjective()->decrease(upperBound);
//        propagate();                  // not needed ???????????????????????????????????????????????????
        cout << "New solution found!" << endl << *this;
    }
}

bool WCSP::solve()
{
    Cost initialUpperBound = getUb() + 1;
    upperBound = initialUpperBound;
    nbBacktracks = 0;
    nbNodes = 0;
    try {
        assert(!objectiveChanged);
        assert(NC.empty());            // (initial) propagation is already done
        assert(AC.empty());
        store();
        recursiveSolve();
    } catch (Contradiction) {
    }
    restore();
    if (upperBound < initialUpperBound) {
        cout << "Optimun: " << upperBound << " in " << nbBacktracks << " backtracks" << endl;
        return true;
    } else {
        cout << "No solution in " << nbBacktracks << " backtracks" << endl;
        return false;
    }
}
