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

CostVariable *WCSP::getNextUnassignedVar()
{
    for (unsigned int i=0; i<vars.size(); i++) {
        if (vars[i]->unassigned()) return vars[i];
    }
    return NULL;
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
            // cout << *this;
            for (unsigned int i=0; i<vars.size(); i++) cout << *vars[i] << endl;
        }
        if (ToulBar2::verbose >= 1) cout << "[" << storeData->getDepth() << "," << getLb() << "," << upperBound << "," << getDomainSizeSum() << "] Try " << x->getName() << " = " << value << endl;
        x->assign(value);
        propagate();
        recursiveSolve();
    } catch (Contradiction) {
    }
    restore();
    nbBacktracks++;
    getObjective()->decrease(upperBound - 1);
    if (ToulBar2::verbose >= 2) {
        // cout << *this;
        for (unsigned int i=0; i<vars.size(); i++) cout << *vars[i] << endl;
    }
    if (ToulBar2::verbose >= 1) cout << "[" << storeData->getDepth() << "," << getLb() << "," << upperBound << "," << getDomainSizeSum() << "] Refute " << x->getName() << " != " << value << endl;
    x->remove(value);
    propagate();
    recursiveSolve();
}

void WCSP::naryChoicePoint(CostVariable *x)
{
    assert(x->unassigned());
    for (Variable::iterator iter = x->begin(); iter != x->end(); ++iter) {
        try {
            store();
            if (ToulBar2::verbose >= 2) {
                for (unsigned int i=0; i<vars.size(); i++) cout << *vars[i] << endl;
            }
            if (ToulBar2::verbose >= 1) cout << "[" << storeData->getDepth() << "," << getLb() << "," << upperBound << "," << getDomainSizeSum() << "] Try " << x->getName() << " = " << *iter << endl;
            getObjective()->decrease(upperBound - 1);
            x->assign(*iter);
            propagate();
            recursiveSolve();
        } catch (Contradiction) {
        }
        restore();
        nbBacktracks++;
    }
}

/*
 * Depth-First Branch and Bound
 * 
 */
            
void WCSP::recursiveSolve()
{
    CostVariable *var = getVarMinDomainDivMaxDegree();
//    CostVariable *var = getNextUnassignedVar();
    if (var != NULL) {
        assert(var->canbe(var->getSupport()));
        binaryChoicePoint(var, var->getSupport());
//        binaryChoicePoint(var, var->getInf());
//        naryChoicePoint(var);
    } else {
        upperBound = getLb();
        getObjective()->decrease(upperBound);
//        propagate();                  // not needed ???????????????????????????????????????????????????
        cout << "New solution: " <<  upperBound << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes)" << endl;
        if (ToulBar2::showSolutions) cout << *this;
    }
}

bool WCSP::solve()
{
    Cost initialUpperBound = getUb() + 1;
    upperBound = initialUpperBound;
    nbBacktracks = 0;
    nbNodes = 0;
    try {
        assert(!objectiveChanged);            // (initial) propagation is already done
        assert(NC.empty());
        assert(AC.empty());
        store();
        recursiveSolve();
    } catch (Contradiction) {
    }
    restore();
    if (upperBound < initialUpperBound) {
        cout << "Optimun: " << upperBound << " in " << nbBacktracks << " backtracks (and " << nbNodes << " nodes)" << endl;
        return true;
    } else {
        cout << "No solution in " << nbBacktracks << " backtracks (and " << nbNodes << " nodes)" << endl;
        return false;
    }
}
