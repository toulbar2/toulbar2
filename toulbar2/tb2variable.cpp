/*
 * ****** Variable with domain represented by an interval or an enumerated domain *******
 */
 
#include "tb2system.hpp"
#include "tb2variable.hpp"
#include "tb2wcsp.hpp"

/*
 * Constructors and misc.
 * 
 */

Variable::Variable(WCSP *w, string n, Value iinf, Value isup) : 
        WCSPLink(w,w->numberOfVariables()), name(n),
        inf(iinf, &w->getStore()->storeValue), sup(isup, &w->getStore()->storeValue), 
        constrs(&w->getStore()->storeConstraint), deltaCost(0, &w->getStore()->storeCost),
        maxCost(0, &w->getStore()->storeCost), maxCostValue(iinf, &w->getStore()->storeValue), 
        NCBucket(-1, &w->getStore()->storeValue)
        
{
    if (w->getStore()->getDepth() > 0) {
        cerr << "You cannot create a variable during the search!" << endl;
        exit(EXIT_FAILURE);
    }
    w->link(this);

    linkNCBucket.content = this;
    linkNCQueue.content.var = this;
    linkNCQueue.content.timeStamp = -1;
    linkIncDecQueue.content.var = this;
    linkIncDecQueue.content.timeStamp = -1;
    linkIncDecQueue.content.incdec = NOTHING_EVENT;
    linkEliminateQueue.content.var = this;
    linkEliminateQueue.content.timeStamp = -1;
}

DLink<ConstraintLink> *Variable::link(Constraint *c, int index)
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

void Variable::sortConstraints()
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

void Variable::deconnect(DLink<ConstraintLink> *link) {
    if (!link->removed) {
        getConstrs()->erase(link, true);
	if (getDegree() <= ToulBar2::elimLevel) queueEliminate();
    }
}

ostream& operator<<(ostream& os, Variable &var) {
    os << var.name;
    var.print(os);
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

void Variable::queueNC()
{
    wcsp->queueNC(&linkNCQueue);
}

void Variable::queueInc()
{
    wcsp->queueInc(&linkIncDecQueue);
}

void Variable::queueDec()
{
    wcsp->queueDec(&linkIncDecQueue);
}

void Variable::queueEliminate()
{
    wcsp->queueEliminate(&linkEliminateQueue);
}

void Variable::changeNCBucket(int newBucket)
{
    if (NCBucket != newBucket) {
        if (ToulBar2::verbose >= 3) cout << "changeNCbucket " << getName() << ": " << NCBucket << " -> " << newBucket << endl;
        wcsp->changeNCBucket(NCBucket, newBucket, &linkNCBucket);
        NCBucket = newBucket;
    }
}

void Variable::setMaxUnaryCost(Value a, Cost cost)
{
    assert(canbe(a));
    maxCostValue = a;
    assert(cost >= 0);
    if (maxCost != cost) {
        maxCost = cost;
        int newbucket = min(cost2log2(cost), wcsp->getNCBucketSize() - 1);
        changeNCBucket(newbucket);
    }
}

void Variable::extendAll(Cost cost)
{
    assert(cost > 0);
    deltaCost += cost;          // Warning! Possible overflow???
    queueNC();
}

void Variable::propagateIncDec(int incdec)
{
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        if (incdec & INCREASE_EVENT) (*iter).constr->increase((*iter).scopeIndex);
        if (incdec & DECREASE_EVENT) (*iter).constr->decrease((*iter).scopeIndex);
    }
}
