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
    linkIncDecQueue.content.var = this;
    linkIncDecQueue.content.timeStamp = -1;
    linkIncDecQueue.content.incdec = NOTHING_EVENT;
    linkACQueue.content.var = this;
    linkACQueue.content.timeStamp = -1;
    linkDACQueue.content.var = this;
    linkDACQueue.content.timeStamp = -1;
}

DLink<ConstraintLink> *CostVariable::postConstraint(Constraint *c, int index)
{
    ConstraintLink e;
    e.constr = c;
    e.scopeIndex = index;
    DLink<ConstraintLink> *elt = new DLink<ConstraintLink>;
    elt->content = e;
    constrs.push_back(elt,true);
    return elt;
}

int cmpConstraint(const void *p1, const void *p2)
{
    DLink<ConstraintLink> *c1 = *((DLink<ConstraintLink> **) p1);
    DLink<ConstraintLink> *c2 = *((DLink<ConstraintLink> **) p2);
    int v1 = c1->content.constr->getSmallestVarIndexInScope(c1->content.scopeIndex);
    int v2 = c2->content.constr->getSmallestVarIndexInScope(c2->content.scopeIndex);
    if (v1 < v2) return -1;
    else if (v1 > v2) return 1;
    else return 0;
}

void CostVariable::sortConstraints()
{
    int size = constrs.getSize();
    DLink<ConstraintLink> *sorted[size];
    int i=0;
    for (ConstraintList::iterator iter = constrs.begin(); iter != constrs.end(); ++iter) {
        sorted[i++] = iter.getElt();
    }
    qsort(sorted, size, sizeof(DLink<ConstraintLink> *), cmpConstraint);
    for (int i = 0; i < size; i++) {
        constrs.erase(sorted[i],true);
        constrs.push_back(sorted[i],true);
    }
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
            os << " > s:" << var.support;
 //        assert(var.canbe(var.support) && var.getCost(var.support)==0);
        } else {
            os << " " << var.getInfCost() << "," << var.getSupCost() << " >";
        }
    }
    if (ToulBar2::verbose >= 3) {
        for (ConstraintList::iterator iter=var.constrs.begin(); iter != var.constrs.end(); ++iter) {
            os << " (" << (*iter).constr << "," << (*iter).scopeIndex << ")";
        }
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

void CostVariable::queueInc()
{
    wcsp->IncDec.push(&linkIncDecQueue, INCREASE_EVENT, wcsp->nbNodes);
}

void CostVariable::queueDec()
{
    wcsp->IncDec.push(&linkIncDecQueue, DECREASE_EVENT, wcsp->nbNodes);
}

void CostVariable::queueAC()
{
    wcsp->AC.push(&linkACQueue, wcsp->nbNodes);
}

void CostVariable::queueDAC()
{
    wcsp->DAC.push(&linkDACQueue, wcsp->nbNodes);
}

void CostVariable::changeNCBucket(int newBucket)
{
    if (NCBucket != newBucket) {
        if (ToulBar2::verbose >= 3) cout << "changeNCbucket " << getName() << ": " << NCBucket << " -> " << newBucket << endl;
        wcsp->changeNCBucket(NCBucket, newBucket, &linkNCBucket);
        NCBucket = newBucket;
    }
}

void CostVariable::setMaxUnaryCost(Value a, Cost cost)
{
    assert(canbe(a));
    maxCostValue = a;
    assert(cost >= 0);
    if (maxCost != cost) {
        maxCost = cost;
        int newbucket = min(cost2log2(cost), wcsp->NCBucketSize - 1);
        changeNCBucket(newbucket);
    }
}

void CostVariable::project(Value value, Cost cost)
{
    assert(cost > 0);
    assert(enumerated);
    Cost oldcost = getCost(value);
    costs[toIndex(value)] += cost;
    Cost newcost = oldcost + cost;
    if (value == maxCostValue || newcost > maxCost) queueNC();
    if (oldcost == 0 && cost > 0) queueDAC();
    if (newcost + wcsp->getLb() > wcsp->getUb()) remove(true, value);     // Avoid any unary cost overflow
}

void CostVariable::projectInfCost(Cost cost)
{
    assert(cost > 0);
    if (enumerated) {
        Value value = getInf();
        project(value, cost);
        if (support == value) findSupport();
    } else {
        Cost oldcost = getInfCost();
        infCost += cost;
        Cost newcost = oldcost + cost;
        if (getInf() == maxCostValue || newcost > maxCost) queueNC();
        if (newcost + wcsp->getLb() > wcsp->getUb()) increase(true, getInf() + 1);     // Avoid any unary cost overflow
        if (getSup() == getInf() + 1 && getInfCost() > 0 && getSupCost() > 0) {
            Cost minCost = min(getInfCost(),getSupCost());
            extendAll(minCost);
            if (ToulBar2::verbose >= 2) cout << "lower bound increased " << wcsp->getLb() << " -> " << wcsp->getLb()+minCost << endl;
            wcsp->increaseLb(wcsp->getLb() + minCost);
        }
    }
}

void CostVariable::projectSupCost(Cost cost)
{
    assert(cost > 0);
    if (enumerated) {
        Value value = getSup();
        project(value, cost);
        if (support == value) findSupport();
    } else {
        Cost oldcost = getSupCost();
        supCost += cost;
        Cost newcost = oldcost + cost;
        if (getSup() == maxCostValue || newcost > maxCost) queueNC();
        if (newcost + wcsp->getLb() > wcsp->getUb()) decrease(true, getSup() - 1);     // Avoid any unary cost overflow
        if (getSup() == getInf() + 1 && getInfCost() > 0 && getSupCost() > 0) {
            Cost minCost = min(getInfCost(),getSupCost());
            extendAll(minCost);
            if (ToulBar2::verbose >= 2) cout << "lower bound increased " << wcsp->getLb() << " -> " << wcsp->getLb()+minCost << endl;
            wcsp->increaseLb(wcsp->getLb() + minCost);
        }
    }
}

void CostVariable::extend(Value value, Cost cost)
{
    assert(cost > 0);
    assert(enumerated);
    assert(costs[toIndex(value)] >= cost);
    costs[toIndex(value)] -= cost;
    if (value == maxCostValue) queueNC();
}

void CostVariable::extendAll(Cost cost)
{
    assert(cost > 0);
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
            Cost cost = getCost(*iter);
            if (cost < minCost) {
                minCost = cost;
                newSupport = *iter;
            }
        }
        if (minCost > 0) {
            extendAll(minCost);
            if (ToulBar2::verbose >= 2) cout << "lower bound increased " << wcsp->getLb() << " -> " << wcsp->getLb()+minCost << endl;
            wcsp->increaseLb(wcsp->getLb() + minCost);
        }
        assert(canbe(newSupport) && getCost(newSupport) == 0);
        support = newSupport;
    }
}

void CostVariable::propagateNC()
{
    if (ToulBar2::verbose >= 3) cout << "propagateNC for " << getName() << endl;
    if (enumerated) {
        Value maxcostvalue = getSup()+1;
        Cost maxcost = -1;
        // Warning! the first value must be visited because it may be removed
        for (iterator iter = begin(); iter != end(); ++iter) {
            Cost cost = getCost(*iter);
            if (cost + wcsp->getLb() > wcsp->getUb()) {
                remove(true, *iter);
            } else if (cost > maxcost) {
                maxcostvalue = *iter;
                maxcost = cost;
            }
        }
        setMaxUnaryCost(maxcostvalue, getCost(maxcostvalue));
    } else {
        if (getInfCost() + wcsp->getLb() > wcsp->getUb()) increase(true, getInf() + 1);
        if (getSupCost() + wcsp->getLb() > wcsp->getUb()) decrease(true, getSup() - 1);
        if (getInfCost() > getSupCost()) {
            setMaxUnaryCost(getInf(), getInfCost());
        } else {
            setMaxUnaryCost(getSup(), getSupCost());
        }
    }
}

bool CostVariable::verifyNC()
{
    bool supported = false;
    if (enumerated) {
        for (iterator iter = begin(); iter != end(); ++iter) {
            if (getCost(*iter) + wcsp->getLb() > wcsp->getUb()) {
                cout << *this << " not NC!" << endl;
                return false;
            }
            if (getCost(*iter) == 0) supported = true;
        }
        if (!supported) cout << *this << " not NC*!" << endl;
    } else {
        if (getInfCost() + wcsp->getLb() > wcsp->getUb()) {
            cout << *this << " has inf cost not NC!" << endl;
            return false;
        }
        if (getSupCost() + wcsp->getLb() > wcsp->getUb()) {
            cout << *this << " has sup cost not NC!" << endl;
            return false;
        }
        supported = (getDomainSize() > 2 || getInfCost() == 0 || getSupCost() == 0);
        if (!supported) cout << *this << " not NC*!" << endl;
    }
    return supported;
}

void CostVariable::propagateIncDec(int incdec)
{
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        if (incdec & INCREASE_EVENT) (*iter).constr->increase((*iter).scopeIndex);
        if (incdec & DECREASE_EVENT) (*iter).constr->decrease((*iter).scopeIndex);
    }
}

void CostVariable::propagateAC()
{
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        (*iter).constr->remove((*iter).scopeIndex);
    }
}

void CostVariable::propagateDAC()
{
    assert(enumerated);
    for (ConstraintList::iterator iter=constrs.rbegin(); iter != constrs.rend(); --iter) {
        (*iter).constr->projectFromZero((*iter).scopeIndex);
    }
}

// noZeroCostRemoved flag is used to avoid checking NC* and queueing in DAC if it was already done
void CostVariable::increaseFromOutside(bool noZeroCostRemoved)
{
    if (!enumerated) infCost = deltaCost;
    if (!noZeroCostRemoved) {
        if (getInf() > maxCostValue) queueNC();
        if (getEnumerated()) {
            if (getInf() > support) findSupport();
            queueDAC();
        }
    }
    queueInc();
}

void CostVariable::decreaseFromOutside(bool noZeroCostRemoved)
{
    if (!enumerated) supCost = deltaCost;
    if (!noZeroCostRemoved) {
        if (getSup() < maxCostValue) queueNC();
        if (getEnumerated()) {
            if (getSup() < support) findSupport();
            queueDAC();
        }
    }
    queueDec();
}

void CostVariable::removeFromOutside(bool noZeroCostRemoved, Value value)
{
    if (!noZeroCostRemoved) {
        if (value == maxCostValue) queueNC();
        if (getEnumerated()) {
            if (value == support) findSupport();
            queueDAC();
        }
    }
    queueAC();
}

void CostVariable::assignFromOutside(Value prevInf, Value prevSup)
{
    assert(prevInf != prevSup);
    changeNCBucket(-1);
    support = getValue();
    maxCostValue = getValue();
    maxCost = 0;
    Cost cost = ((enumerated)?getCost(maxCostValue):
            ((maxCostValue==prevInf)?getInfCost():((maxCostValue==prevSup)?getSupCost():0)));
    if (cost > 0) {
        deltaCost += cost;
        if (ToulBar2::verbose >= 2) cout << "lower bound increased " << wcsp->getLb() << " -> " << wcsp->getLb()+cost << endl;
        wcsp->increaseLb(wcsp->getLb() + cost);
    }
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        (*iter).constr->assign((*iter).scopeIndex);
    }
}
