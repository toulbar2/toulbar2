/*
 * **************** Globals and WCSP methods ************************
 * 
 */

#include "tb2system.hpp"
#include "tb2wcsp.hpp"
#include "tb2binconstr.hpp"


/*
 * Global variables
 * 
 */
 
int ToulBar2::verbose  = 0;


/*
 * WCSP constructors
 * 
 */

WCSP::WCSP(Variable *obj, Store *s) : objective(obj), storeData(s), 
        NCBucketSize(cost2log2(objective->getSup()) + 1),
        NCBuckets(NCBucketSize, CostVariableList(&s->storeCostVariable)),
        upperBound(objective->getSup()), objectiveChanged(false),
        nbNodes(0), nbBacktracks(0)
{
    objective->addWCSP(this, -1);
}

// create a CostVariable object AND link it and the original variable to the wcsp
// return the index in wcsp.vars
int WCSP::link(Variable *x)
{
    if (x == objective) {
        cerr << "Cannot use the objective of a WCSP as a variable in a constraint of the same WCSP!" << endl;
        exit(EXIT_FAILURE);
    }
    int idx = x->findWCSPIndex(this);
    if (idx < 0) {
        idx = vars.size();
        x->addWCSP(this, idx);
        CostVariable *xx = new CostVariable(x,&storeData->storeCost,&storeData->storeConstraint,&storeData->storeValue,&storeData->storeValue);
        xx->wcsp = this;
        xx->wcspIndex = idx;
        vars.push_back(xx);        
    } else assert(vars[idx]->getVar() == x);
    return idx;
}

// link the constraint to the wcsp
// return the index in wcsp.constrs
int WCSP::link(Constraint *c)
{
    int idx = constrs.size();
    c->wcsp = this;
    c->wcspIndex = idx;
    constrs.push_back(c);
    return idx;
}

int WCSP::addBinaryConstraint(Variable *x, Variable *y, vector<Cost> &tab)
{
    CostVariable *xx = vars[link(x)];
    CostVariable *yy = vars[link(y)];
    Constraint *c = new BinaryConstraint(xx,yy,tab,&storeData->storeCost);
    int index = link(c);
   c->propagate();
    propagate();
    return index;
}

ostream& operator<<(ostream& os, WCSP &wcsp)
{
    os << "Objective: " << *wcsp.objective << endl;
    os << "Variables:" << endl;
    for (unsigned int i=0; i<wcsp.vars.size(); i++) os << *wcsp.vars[i] << endl;
    os << "Constraints:" << endl;
    for (unsigned int i=0; i<wcsp.constrs.size(); i++) if (wcsp.constrs[i]->connected()) os << *wcsp.constrs[i];
    return os;
}


/*
 * WCSP propagation methods
 * 
 */

void WCSP::printNCBuckets()
{
    for (int bucket = 0; bucket < NCBucketSize; bucket++) {
        cout << "NC " << bucket << ":";
        for (CostVariableList::iterator iter = NCBuckets[bucket].begin (); iter != NCBuckets[bucket].end(); ++iter) {
           cout << " " << (*iter)->getName() << "," << (*iter)->maxCostValue << "," << (*iter)->maxCost;
           assert((*iter)->canbe((*iter)->maxCostValue));
           assert((*iter)->getCost((*iter)->maxCostValue) == (*iter)->maxCost);
        }
        cout << endl;
    }
}

bool WCSP::verify()
{
    for (unsigned int i=0; i<vars.size(); i++) {
        if (vars[i]->unassigned() && !vars[i]->verifyNC()) return false;
    }
    for (unsigned int i=0; i<constrs.size(); i++) {
        if (constrs[i]->connected() &&  !constrs[i]->verify()) return false;
    }
    return true;
}

void WCSP::propagateNC()
{
    if (ToulBar2::verbose >= 2) cout << "NCQueue size: " << NC.getSize() << endl;
    while (!NC.empty()) {
        CostVariable *x = NC.pop();
        if (x->unassigned()) x->propagateNC();
    }
    if (ToulBar2::verbose >= 3) {
        for (unsigned int i=0; i<vars.size(); i++) cout << *vars[i] << endl;
    }
    if (ToulBar2::verbose >= 2) printNCBuckets();

    if (objectiveChanged) {
        objectiveChanged = false;
        for (int bucket = cost2log2(getUb() - getLb() + 1); bucket < NCBucketSize; bucket++) {
            for (CostVariableList::iterator iter = NCBuckets[bucket].begin(); iter != NCBuckets[bucket].end(); ++iter) {
                CostVariable *x = *iter;
                if (x->unassigned() && x->maxCost + getLb() > getUb()) x->propagateNC();
            }
        }
    }
}

void WCSP::propagateAC()
{
    if (ToulBar2::verbose >= 2) cout << "ACQueue size: " << AC.getSize() << endl;
    while (!AC.empty()) {
        CostVariable *x = AC.pop();
        if (x->unassigned()) x->propagateAC();
    }
}

void WCSP::propagate()
{
    while (!AC.empty() || objectiveChanged || !NC.empty()) {
        propagateAC();
        while (objectiveChanged || !NC.empty()) propagateNC();
    }
    if (ToulBar2::verbose >= 2) verify();
}
