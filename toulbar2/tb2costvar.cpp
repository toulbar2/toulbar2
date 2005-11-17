/*
 * **************** Variable extended with unary costs **************
 */

#include "tb2system.hpp"
#include "tb2costvar.hpp"
#include "tb2wcsp.hpp"


/*
 * Constructors and misc.
 * 
 */

CostVariable::CostVariable(Variable *v, StoreStack<Cost,Cost> *storeCost, 
        StoreStack<ConstraintList, DLink<ConstraintLink> *> *storeConstraint,
        StoreStack<Value,Value> *storeValue, StoreStack<int,int> *storeInt) : 
        enumerated(v->getEnumerated()), var(v), constrs(storeConstraint),
        infCost(0, storeCost), supCost(0, storeCost), deltaCost(0, storeCost), support(v->getInf(), storeValue),
        maxCost(0, storeCost), maxCostValue(v->getInf(), storeValue), NCBucket(-1, storeInt)
{
    if (enumerated) {
        costs = vector<StoreCost>(v->getDomainInitSize(), StoreCost(0, storeCost));
    }
    linkNCBucket.content = this;
    linkNCQueue.content.var = this;
    linkNCQueue.content.timeStamp = -1;
    linkACQueue.content.var = this;
    linkACQueue.content.timeStamp = -1;
}

DLink<ConstraintLink> *CostVariable::addConstraint(Constraint *c, int index)
{
    ConstraintLink e;
    e.constr = c;
    e.scopeIndex = index;
    DLink<ConstraintLink> *elt = new DLink<ConstraintLink>;
    elt->content = e;
    constrs.push_back(elt,true);
    return elt;
}

ostream& operator<<(ostream& os, CostVariable &var)
{
    os << *var.var;
    if (var.unassigned()) {
        os << " <";
        if (var.enumerated) {
            for (CostVariable::iterator iter=var.begin(); iter != var.end(); ++iter) {
                os << " " << var.getCost(*iter);
            }
        } else {
            os << " " << var.getInfCost() << "," << var.getSupCost();
        }
        os << " > s:" << var.support;
        assert(var.canbe(var.support) && var.getCost(var.support)==0);
    }
    for (ConstraintList::iterator iter=var.constrs.begin(); iter != var.constrs.end(); ++iter) {
        os << " (" << (*iter).constr << "," << (*iter).scopeIndex << ")";
    }
    return os;
}


/*
 * Propagation methods
 * 
 */

void CostVariable::queueNC()
{
    wcsp->NC.push(&linkNCQueue, wcsp->nbNodes);
}

void CostVariable::queueAC()
{
    wcsp->AC.push(&linkACQueue, wcsp->nbNodes);
}

void CostVariable::changeNCBucket(int newBucket)
{
    if (NCBucket != newBucket) {
        if (ToulBar2::verbose >= 3) cout << "changeNCbucket " << getName() << ": " << NCBucket << " -> " << newBucket << endl;
        wcsp->changeNCBucket(NCBucket, newBucket, &linkNCBucket);
        NCBucket = newBucket;
    }
}

void CostVariable::setMaxUnaryCost(Value a)
{
    assert(canbe(a));
    maxCostValue = a;
    Cost cost = getCost(a);
    assert(cost >= 0);
    if (maxCost != cost) {
        maxCost = cost;
        int newbucket = min(cost2log2(cost), wcsp->NCBucketSize - 1);
        changeNCBucket(newbucket);
    }
}

void CostVariable::project(Value value, Cost cost)
{
    assert(cost >= 0);
    assert(enumerated);
    costs[toIndex(value)] += cost;
    Cost newcost = getCost(value);
    if (value == maxCostValue || newcost > maxCost) queueNC();
    if (newcost + wcsp->getLb() > wcsp->getUb()) remove(value);
}

void CostVariable::extend(Value value, Cost cost)
{
    assert(cost >= 0);
    assert(enumerated);
    assert(costs[toIndex(value)] >= cost);
    costs[toIndex(value)] -= cost;
    if (value == maxCostValue) queueNC();
}

void CostVariable::extendAll(Cost cost)
{
    assert(cost >= 0);
    deltaCost += cost;          // Warning! Possible overflow???
    queueNC();
}

void CostVariable::findSupport()
{
    assert(enumerated);
    if (cannotbe(support) || getCost(support) > 0) {
        Value newSupport = getInf();
        Cost minCost = getCost(newSupport);
        iterator iter = begin();
        for (++iter; minCost > 0 && iter != end(); ++iter) {
            if (getCost(*iter) < minCost) {
                minCost = getCost(*iter);
                newSupport = *iter;
            }
        }
        if (minCost > 0) {
            extendAll(minCost);
            if (ToulBar2::verbose >= 2) cout << "lower bound increased " << wcsp->getLb() << " -> " << wcsp->getLb()+minCost << endl;
            wcsp->getObjective()->increase(wcsp->getLb() + minCost);
        }
        assert(canbe(newSupport) && getCost(newSupport) == 0);
        support = newSupport;
    }
}

void CostVariable::propagateNC()
{
    if (ToulBar2::verbose >= 3) cout << "propagateNC for " << getName() << endl;
    Value maxcostvalue = getSup()+1;
    Cost maxcost = -1;
    // Warning! Avoid reinsertion into NC due to a value removal
    linkNCQueue.content.timeStamp = wcsp->nbNodes;
    // Warning! the first value must be visited because it may be removed
    for (iterator iter = begin(); iter != end(); ++iter) {
        if (getCost(*iter) + wcsp->getLb() > wcsp->getUb()) {
            remove(*iter);
        } else if (getCost(*iter) > maxcost) {
            maxcostvalue = *iter;
            maxcost = getCost(*iter);
        }
    }
    linkNCQueue.content.timeStamp = -1;
    setMaxUnaryCost(maxcostvalue);
}

bool CostVariable::verifyNC()
{
    bool supported = false;
    for (iterator iter = begin(); iter != end(); ++iter) {
        if (getCost(*iter) + wcsp->getLb() > wcsp->getUb()) return false;
        if (getCost(*iter) == 0) supported = true;
    }
    return supported;
}

void CostVariable::propagateAC()
{
    if (ToulBar2::verbose >= 3) cout << "propagateAC for " << getName() << endl;
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        (*iter).constr->propagate((*iter).scopeIndex);
    }
}

void CostVariable::increaseWCSP()
{
    if (getInf() > maxCostValue) queueNC();
    if (getInf() > getSupport()) findSupport();
    queueAC();
}

void CostVariable::decreaseWCSP()
{
    if (getSup() < maxCostValue) queueNC();
    if (getSup() < getSupport()) findSupport();
    queueAC();
}

void CostVariable::assignWCSP()
{
    changeNCBucket(-1);
    Cost cost = getCost(getValue());
    if (cost > 0) {
        if (ToulBar2::verbose >= 2) cout << "lower bound increased " << wcsp->getLb() << " -> " << wcsp->getLb()+cost << endl;
        wcsp->getObjective()->increase(wcsp->getLb() + cost);
    }
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        (*iter).constr->assign((*iter).scopeIndex);
    }
}
void CostVariable::removeWCSP(Value value)
{
    if (value == maxCostValue) queueNC();
    if (value == getSupport()) findSupport();
    queueAC();
}
