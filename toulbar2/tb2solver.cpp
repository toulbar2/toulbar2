/*
 * **************** Generic solver *******************
 * 
 */

#include "tb2solver.hpp"


/*
 * Solver constructors
 * 
 */

Solver::Solver(int storeSize, Cost initUpperBound) : store(storeSize), nbNodes(0), nbBacktracks(0), 
        unassignedVars(&store.storeVariable),
        upperBound(initUpperBound), objective("obj",0,MAX_COST,this,OBJ_VAR), wcsp(&objective, &store)
{
}
        
Solver::~Solver()
{
    for (unsigned int i=0; i<vars.size(); i++) delete vars[i];
}

void Solver::read_wcsp(const char *fileName)
{
    wcsp.read_wcsp(fileName, this);
}

/*
 * Variable ordering heuristics
 * 
 */

Variable *Solver::getVarMinDomainDivMaxDegree()
{
    Variable *var = NULL;
    double best = MAX_VAL - MIN_VAL;

    for (VariableList::iterator iter = unassignedVars.begin(); iter != unassignedVars.end(); ++iter) {
        // remove following "+1" when isolated variables are automatically assigned
        double heuristic = (double) (*iter)->getDomainSize() / ((*iter)->getDegree() + 1);
        if (var == NULL || heuristic < best) { // - 1./100001.) {
            best = heuristic;
            var = *iter;
        }
    }
    return var;
}

Variable *Solver::getNextUnassignedVar()
{
    return (unassignedVars.empty())?NULL:*unassignedVars.begin();
}

/*
 * Choice points
 * 
 */

void Solver::binaryChoicePoint(Variable *x, Value value)
{
    assert(x->unassigned());
    assert(x->canbe(value));
    try {
        store.store();
        if (ToulBar2::verbose >= 2) {
            cout << wcsp;
        }
        if (ToulBar2::verbose >= 1) cout << "[" << store.getDepth() << "," << wcsp.getLb() << "," << upperBound << "," << wcsp.getDomainSizeSum() << "] Try " << x->getName() << " = " << value << endl;
        nbNodes++;
        x->assign(value);
        wcsp.propagate();
        recursiveSolve();
    } catch (Contradiction) {
        whenContradiction();
    }
    store.restore();
    nbBacktracks++;
    wcsp.decreaseUb(upperBound - 1);
    if (ToulBar2::verbose >= 2) {
        cout << wcsp;
    }
    if (ToulBar2::verbose >= 1) cout << "[" << store.getDepth() << "," << wcsp.getLb() << "," << upperBound << "," << wcsp.getDomainSizeSum() << "] Refute " << x->getName() << " != " << value << endl;
    nbNodes++;
    x->remove(value);
    wcsp.propagate();
    recursiveSolve();
}

struct CostValue {
    Value val;
    Cost cost;
};

int cmpCostValue(const void *p1, const void *p2)
{
    Cost c1 = ((CostValue *) p1)->cost;
    Cost c2 = ((CostValue *) p2)->cost;
    Value v1 = ((CostValue *) p1)->val;
    Value v2 = ((CostValue *) p2)->val;
    if (c1 < c2) return -1;
    else if (c1 > c2) return 1;
    else if (v1 < v2) return -1;
    else if (v1 > v2) return 1;
    else return 0;
}

void Solver::narySortedChoicePoint(Variable *x)
{
    assert(x->unassigned());
    assert(x->getEnumerated());
    int size = x->getDomainSize();
    CostValue sorted[size]; 
    int v = 0;
    int currentVariableIndex = x->findWCSPIndex(&wcsp);
    for (Variable::iterator iter = x->begin(); iter != x->end(); ++iter) {
        sorted[v].val = *iter;
        sorted[v].cost = wcsp.getUnaryCost(currentVariableIndex,*iter);
        v++;
    }
    qsort(sorted, size, sizeof(CostValue), cmpCostValue);
    for (v = 0; v < size; v++) {
        try {
            store.store();
            nbNodes++;
            if (ToulBar2::verbose >= 2) {
                cout << wcsp << endl;
            }
            if (ToulBar2::verbose >= 1) cout << "[" << store.getDepth() << "," << wcsp.getLb() << "," << upperBound << "," << wcsp.getDomainSizeSum() << "] Try " << x->getName() << " = " << sorted[v].val << endl;
            wcsp.decreaseUb(upperBound - 1);
            x->assign(sorted[v].val);
            wcsp.propagate();
            recursiveSolve();
        } catch (Contradiction) {
            whenContradiction();
        }
        store.restore();
    }
    nbBacktracks++;
}

void Solver::naryChoicePoint(Variable *x)
{
    assert(x->unassigned());
    for (Variable::iterator iter = x->begin(); iter != x->end(); ++iter) {
        try {
            store.store();
            nbNodes++;
            if (ToulBar2::verbose >= 2) {
                cout << wcsp;
            }
            if (ToulBar2::verbose >= 1) cout << "[" << store.getDepth() << "," << wcsp.getLb() << "," << upperBound << "," << wcsp.getDomainSizeSum() << "] Try " << x->getName() << " = " << *iter << endl;
            wcsp.decreaseUb(upperBound - 1);
            x->assign(*iter);
            wcsp.propagate();
            recursiveSolve();
        } catch (Contradiction) {
            whenContradiction();
        }
        store.restore();
    }
    nbBacktracks++;
}

/*
 * Depth-First Branch and Bound
 * 
 */

void Solver::whenContradiction()
{
    wcsp.whenContradiction();
}
            
void Solver::recursiveSolve()
{
//    Variable *var = getNextUnassignedVar();
    Variable *var = getVarMinDomainDivMaxDegree();
    if (var != NULL) {
        if (var->getEnumerated()) {
            if (ToulBar2::binaryBranching) {
//                binaryChoicePoint(var, var->getInf());
                assert(var->canbe(wcsp.getSupport(var->findWCSPIndex(&wcsp))));
                binaryChoicePoint(var, wcsp.getSupport(var->findWCSPIndex(&wcsp)));
            } else {
//                naryChoicePoint(var);
                narySortedChoicePoint(var);
            }
        } else {
            binaryChoicePoint(var, var->getInf());
        }
    } else {
        upperBound = wcsp.getLb();
        wcsp.decreaseUb(upperBound);        //  do not need to propagate after this ?
        cout << "New solution: " <<  upperBound << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes)" << endl;
        if (ToulBar2::showSolutions) cout << wcsp;
    }
}

bool Solver::solve()
{
    Cost initialUpperBound = min(wcsp.getUb() + 1, upperBound);
    upperBound = initialUpperBound;
    nbBacktracks = 0;
    nbNodes = 0;
    try {
        store.store();
        wcsp.decreaseUb(upperBound - 1);
        wcsp.propagate();                // initial propagation
        recursiveSolve();
    } catch (Contradiction) {
        whenContradiction();
    }
    store.restore();
    if (upperBound < initialUpperBound) {
        cout << "Optimun: " << upperBound << " in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << endl;
        return true;
    } else {
        cout << "No solution in " << nbBacktracks << " backtracks and " << nbNodes << " nodes" << endl;
        return false;
    }
}
