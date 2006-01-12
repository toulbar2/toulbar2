/*
 * ****** Global soft constraint representing a weighted CSP ********
 * 
 * Contains also ToulBar2 global variable definitions
 */

#include "tb2system.hpp"
#include "tb2wcsp.hpp"
#include "tb2binconstr.hpp"
#include "tb2arithmetic.hpp"


/*
 * Global variables
 * 
 */
 
int ToulBar2::verbose  = 0;
bool ToulBar2::showSolutions  = false;
bool ToulBar2::binaryBranching = false;

/*
 * WCSP constructors
 * 
 */

WCSP::WCSP(Variable *obj, Store *s) : objective(obj), storeData(s), 
        NCBucketSize(cost2log2(objective->getSup()) + 1),
        NCBuckets(NCBucketSize, CostVariableList(&s->storeCostVariable)),
        objectiveChanged(false),
        nbNodes(0)
{
    objective->addWCSP(this, -1);
}

WCSP::~WCSP()
{
    for (unsigned int i=0; i<vars.size(); i++) delete vars[i];
    for (unsigned int i=0; i<constrs.size(); i++) delete constrs[i];
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

int WCSP::postBinaryConstraint(Variable *x, Variable *y, vector<Cost> &tab)
{
    CostVariable *xx = vars[link(x)];
    CostVariable *yy = vars[link(y)];
    Constraint *c = new BinaryConstraint(xx,yy,tab,&storeData->storeCost);
    int index = link(c);
//    c->propagate();       // Let the initial propagation be done only once in solver.cpp
//    propagate();
    xx->queueAC();
    xx->queueDAC();
    yy->queueAC();
    yy->queueDAC();
    return index;
}

int WCSP::postSupxyc(Variable *x, Variable *y, Value cst)
{
    CostVariable *xx = vars[link(x)];
    CostVariable *yy = vars[link(y)];
    Constraint *c = new Supxyc(xx,yy,cst,&storeData->storeCost,&storeData->storeValue);
    int index = link(c);
//    c->propagate();       // Let the initial propagation be done only once in solver.cpp
//    propagate();
    xx->queueInc();
    xx->queueDec();
    yy->queueInc();
    yy->queueDec();
    return index;
}

void WCSP::sortConstraints()
{
    for (unsigned int i=0; i<vars.size(); i++) {
        vars[i]->sortConstraints();
    }
}

Value WCSP::getDomainSizeSum()
{
    Value sum = 0;
    for (unsigned int i=0; i<vars.size(); i++) {
        if (vars[i]->unassigned()) sum += vars[i]->getDomainSize();
    }
    return sum;
}

void WCSP::printNCBuckets()
{
    for (int bucket = 0; bucket < NCBucketSize; bucket++) {
        cout << "NC " << bucket << ":";
        for (CostVariableList::iterator iter = NCBuckets[bucket].begin (); iter != NCBuckets[bucket].end(); ++iter) {
           cout << " " << (*iter)->getName() << "," << (*iter)->maxCostValue << "," << (*iter)->maxCost;
           assert((*iter)->canbe((*iter)->maxCostValue));
           assert((*iter)->getUnaryCost((*iter)->maxCostValue) == (*iter)->maxCost);
        }
        cout << endl;
    }
}

ostream& operator<<(ostream& os, WCSP &wcsp)
{
    os << "Objective: " << *wcsp.objective << endl;
    os << "Variables:" << endl;
    for (unsigned int i=0; i<wcsp.vars.size(); i++) os << *wcsp.vars[i] << endl;
    if (ToulBar2::verbose >= 4) {
        os << "Constraints:" << endl;
        for (unsigned int i=0; i<wcsp.constrs.size(); i++) if (wcsp.constrs[i]->connected()) os << *wcsp.constrs[i];
    }
    return os;
}


/*
 * WCSP propagation methods
 * 
 */

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

void WCSP::whenContradiction()
{
    NC.clear();
    IncDec.clear();
    AC.clear();
    DAC.clear();
    objectiveChanged = false;
    nbNodes++;
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
            for (CostVariableList::iterator iter = NCBuckets[bucket].begin(); iter != NCBuckets[bucket].end();) {
                CostVariable *x = *iter;
                ++iter; // Warning! the iterator could be moved to another place by propagateNC
                if (x->unassigned() && x->maxCost + getLb() > getUb()) x->propagateNC();
            }
        }
    }
}

void WCSP::propagateIncDec()
{
    if (ToulBar2::verbose >= 2) cout << "IncDecQueue size: " << IncDec.getSize() << endl;
    while (!IncDec.empty()) {
        int incdec;
        CostVariable *x = IncDec.pop(&incdec);
        if (x->unassigned()) x->propagateIncDec(incdec);
    }
}

void WCSP::propagateAC()
{
    if (ToulBar2::verbose >= 2) cout << "ACQueue size: " << AC.getSize() << endl;
    while (!AC.empty()) {
        CostVariable *x = AC.pop_min();
        if (x->unassigned()) x->propagateAC();
        // Warning! propagateIncDec() necessary to transform inc/dec event into remove event
        propagateIncDec();          // always examine inc/dec events before remove events
    }
}

void WCSP::propagateDAC()
{
    if (ToulBar2::verbose >= 2) cout << "DACQueue size: " << DAC.getSize() << endl;
    while (!DAC.empty()) {
        CostVariable *x = DAC.pop_max();
        if (x->unassigned()) x->propagateDAC();
        propagateIncDec();          // always examine inc/dec events before projectFromZero events
    }
}

void WCSP::propagate()
{
    while (!IncDec.empty() || !AC.empty() || !DAC.empty() || !NC.empty() || objectiveChanged) {
        propagateIncDec();
        propagateAC();
        assert(IncDec.empty());
        propagateDAC();
        assert(IncDec.empty());
        propagateNC();
    }
    assert( verify() );
    assert(!objectiveChanged);
    assert(NC.empty());
    assert(IncDec.empty());
    assert(AC.empty());
    assert(DAC.empty());
    nbNodes++;
}
