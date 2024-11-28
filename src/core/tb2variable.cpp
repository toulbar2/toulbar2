/*
 * ****** Variable with domain represented by an interval or an enumerated domain *******
 */

#include "tb2variable.hpp"
#include "tb2wcsp.hpp"
#include "tb2binconstr.hpp"
#include "tb2ternaryconstr.hpp"
#include "search/tb2clusters.hpp"

/*
 * Constructors and misc.
 *
 */

Variable::Variable(WCSP* w, string n, Value iinf, Value isup)
    : WCSPLink(w, w->numberOfVariables())
    , name(n)
    , dac(w->numberOfVariables())
    , timestamp(-1)
    , pos(-1)
    , inf(iinf)
    , sup(isup)
    , fulleac((isup > iinf) ? 0 : 1)
    , constrs(&Store::storeConstraint)
    // ,triangles(&Store::storeConstraint)
    , maxCost(MIN_COST)
    , maxCostValue(iinf)
    , NCBucket(-1)
{
    if (Store::getDepth() > 0) {
        cerr << "You cannot create a variable during the search!" << endl;
        throw InternalError();
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
    cluster = -1;
    isSep_ = false;
}

DLink<ConstraintLink>* Variable::link(Constraint* c, int index)
{
    ConstraintLink e;
    e.constr = c;
    e.scopeIndex = index;
    DLink<ConstraintLink>* elt = new DLink<ConstraintLink>;
    elt->content = e;
    //    if (c->isTriangle()) triangles.push_back(elt,true);
    //    else
    constrs.push_back(elt, true);
    return elt;
}

int Variable::getCurrentVarId()
{
    if (assigned())
        return -1;
    if (wcsp->getNbNodes() > timestamp)
        wcsp->updateCurrentVarsId();
    assert(pos >= 0);
    assert(wcsp->getNbNodes() == timestamp);
    return pos;
}

void Variable::setCurrentVarId(int idx)
{
    pos = idx;
    timestamp = wcsp->getNbNodes();
}

void Variable::sortConstraints()
{
    vector<DLink<ConstraintLink>*> sorted;
    for (ConstraintList::iterator iter = constrs.begin(); iter != constrs.end(); ++iter) {
        sorted.push_back(iter.getElt());
    }
    if (abs(ToulBar2::constrOrdering) == CONSTR_ORDER_RANDOM) {
        shuffle(sorted.begin(), sorted.end(), myrandom_generator);
    } else {
        stable_sort(sorted.begin(), sorted.end(), Constraint::cmpConstraintLink);
    }
    for (unsigned int i = 0; i < sorted.size(); i++) {
        constrs.erase(sorted[i], true);
        constrs.push_back(sorted[i], true);
    }
}

void Variable::deconnect(DLink<ConstraintLink>* link, bool reuse)
{
    if (!link->removed) {
        //        if (link->content.constr->isTriangle()) getTriangles()->erase(link, true);
        //        else
        getConstrs()->erase(link, true);

        if (getDegree() <= ToulBar2::elimDegree_ || (ToulBar2::elimDegree_preprocessing_ >= 0 && (getDegree() <= min(1, ToulBar2::elimDegree_preprocessing_) || getTrueDegree() <= ToulBar2::elimDegree_preprocessing_)))
            queueEliminate();
    }
    if (reuse) {
        //        assert(Store::getDepth() == 0);
        link->prev = NULL;
        link->next = NULL;
    }
}

void Variable::deconnect()
{
    for (ConstraintList::iterator iter = constrs.begin(); iter != constrs.end(); ++iter) {
        (*iter).constr->deconnect();
    }
    assign(getSupport());
}

int Variable::getTrueDegree()
{
    //	if (constrs.getSize() >= ToulBar2::weightedDegree) return getDegree(); ///\warning returns an approximate degree if the constraint list is too large!
    //    TSCOPE scope1,scope2,scope3;
    set<int> scope1;
    for (ConstraintList::iterator iter = constrs.begin(); iter != constrs.end(); ++iter) {
        if ((*iter).constr->isSep())
            continue;
        for (int k = 0; k < (*iter).constr->arity(); k++) {
            scope1.insert((*iter).constr->getVar(k)->wcspIndex);
        }
        //		(*iter).constr->getScope(scope2);
        //		for (TSCOPE::iterator iter2=scope2.begin(); iter2 != scope2.end(); ++iter2) {
        //		   scope1.insert( iter2->first );
        //		}
        //     	set_union( scope1.begin(), scope1.end(),
        //	  		   	   scope2.begin(), scope2.end(),
        //			  	   inserter(scope3, scope3.begin()) );
        //		scope1 = scope3;
        //		scope3.clear();
    }
    if (scope1.size() >= 1)
        return scope1.size() - 1;
    else
        return 0;
}

Double Variable::getMaxElimSize()
{
    if (getDegree() == 0)
        return getDomainSize();
    if (getDegree() == 1)
        return (*constrs.begin()).constr->size();
    //  if (constrs.getSize() >= ToulBar2::weightedDegree) return getDegree(); ///\warning returns an approximate degree if the constraint list is too large!
    //    TSCOPE scope1,scope2,scope3;
    map<int, unsigned int> scope1;
    for (ConstraintList::iterator iter = constrs.begin(); iter != constrs.end(); ++iter) {
        if ((*iter).constr->isSep())
            continue;
        for (int k = 0; k < (*iter).constr->arity(); k++) {
            scope1[(*iter).constr->getVar(k)->wcspIndex] = (*iter).constr->getVar(k)->getDomainSize();
        }
    }
    Double sz = 1.;
    for (map<int, unsigned int>::iterator iter = scope1.begin(); iter != scope1.end(); ++iter) {
        sz *= (*iter).second;
    }
    return sz;
}

Long Variable::getWeightedDegree()
{
    Long res = 0;
    for (ConstraintList::iterator iter = constrs.begin(); iter != constrs.end(); ++iter) {
        //    	if((*iter).constr->isSep()) continue;
        res += (*iter).constr->getConflictWeight((*iter).scopeIndex);
        if ((*iter).constr->isSep())
            res--; // do not count unused separators
    }
    return res;
}

void Variable::conflict()
{
    wcsp->conflict();
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
        if (ToulBar2::verbose >= 3)
            cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] changeNCbucket " << getName() << ": " << NCBucket << " -> " << newBucket << endl;
        wcsp->changeNCBucket(NCBucket, newBucket, &linkNCBucket);
        NCBucket = newBucket;
    }
}

void Variable::setMaxUnaryCost(Value a, Cost cost)
{
    assert(canbe(a));
    maxCostValue = a;
    assert(cost >= MIN_COST);
    if (maxCost != cost) {
        if (cost > maxCost)
            queueDEE();
        maxCost = cost;
        int newbucket = min(cost2log2gub(cost), wcsp->getNCBucketSize() - 1);
        changeNCBucket(newbucket);
    }
}

void Variable::projectLB(Cost cost)
{
    if (cost == MIN_COST)
        return;
    if (ToulBar2::verbose >= 2)
        cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] lower bound increased " << wcsp->getLb() << " -> " << wcsp->getLb() + cost << endl;
    if (cost < MIN_COST) {
        wcsp->decreaseLb(cost);
    } else {
        wcsp->increaseLb(cost); // done before cluster LB because of #CSP (assuming a contradiction will occur here)
    }
    if (wcsp->td) {
        if (ToulBar2::verbose >= 2)
            cout << " in cluster C" << getCluster() << " (from " << wcsp->td->getCluster(getCluster())->getLb() << " to " << wcsp->td->getCluster(getCluster())->getLb() + cost << ")" << endl;
        wcsp->td->getCluster(getCluster())->increaseLb(cost);
    }
}

void Variable::propagateIncDec(int incdec)
{
    for (ConstraintList::iterator iter = constrs.begin(); iter != constrs.end(); ++iter) {
        if (incdec & INCREASE_EVENT) {
            (*iter).constr->increase((*iter).scopeIndex);
        }
        if ((*iter).constr->connected() && (incdec & DECREASE_EVENT)) {
            (*iter).constr->decrease((*iter).scopeIndex);
        }
    }
}

// returns true if there is a global cost function currently linked to this variable
bool Variable::isGlobal()
{
    for (ConstraintList::iterator iter = constrs.begin(); iter != constrs.end(); ++iter) {
        if ((*iter).constr->isGlobal()) {
            return true;
        }
    }
    return false;
}

// Looks for the constraint that links this variable with x
BinaryConstraint* Variable::getConstr(Variable* x)
{
    BinaryConstraint* ctr2;
    TernaryConstraint* ctr3;
    for (ConstraintList::iterator iter = constrs.begin(); iter != constrs.end(); ++iter) {
        if ((*iter).constr->isSep() || (*iter).constr->isGlobal())
            continue;

        if ((*iter).constr->isBinary()) {
            ctr2 = (BinaryConstraint*)(*iter).constr;
            if (ctr2->getIndex(x) >= 0)
                return ctr2;
        } else if ((*iter).constr->isTernary()) {
            ctr3 = (TernaryConstraint*)(*iter).constr;
            int idx = ctr3->getIndex(x);
            if (idx >= 0) {
                int idt = (*iter).scopeIndex;
                if ((0 != idx) && (0 != idt))
                    return ctr3->yz;
                else if ((1 != idx) && (1 != idt))
                    return ctr3->xz;
                else
                    return ctr3->xy;
            }
        }
    }
    return NULL;
}

// Looks for the constraint that links this variable with x and which is not a duplicated constraint
BinaryConstraint* Variable::getConstrNotDuplicate(Variable* x)
{
    BinaryConstraint* ctr2;
    TernaryConstraint* ctr3;
    for (ConstraintList::iterator iter = constrs.begin(); iter != constrs.end(); ++iter) {
        if ((*iter).constr->isDuplicate() || (*iter).constr->isSep() || (*iter).constr->isGlobal())
            continue;

        if ((*iter).constr->isBinary()) {
            ctr2 = (BinaryConstraint*)(*iter).constr;
            if (ctr2->getIndex(x) >= 0)
                return ctr2;
        } else if ((*iter).constr->isTernary()) {
            ctr3 = (TernaryConstraint*)(*iter).constr;
            int idx = ctr3->getIndex(x);
            if (idx >= 0) {
                int idt = (*iter).scopeIndex;
                if ((0 != idx) && (0 != idt))
                    return ctr3->yz;
                else if ((1 != idx) && (1 != idt))
                    return ctr3->xz;
                else
                    return ctr3->xy;
            }
        }
    }
    return NULL;
}

BinaryConstraint* Variable::getConstr(Variable* x, int cid)
{
    BinaryConstraint* res;
    BinaryConstraint* ctr2;
    TernaryConstraint* ctr3;
    for (ConstraintList::iterator iter = constrs.begin(); iter != constrs.end(); ++iter) {
        if ((*iter).constr->isSep() || (*iter).constr->isGlobal())
            continue;

        if ((*iter).constr->isBinary()) {
            ctr2 = (BinaryConstraint*)(*iter).constr;
            if (ctr2->getIndex(x) >= 0) {
                res = ctr2;
                if (res->getCluster() == cid)
                    return res;
            }
        } else if ((*iter).constr->isTernary()) {
            ctr3 = (TernaryConstraint*)(*iter).constr;
            int idx = ctr3->getIndex(x);
            if (idx >= 0) {
                int idt = (*iter).scopeIndex;
                if ((0 != idx) && (0 != idt))
                    res = ctr3->yz;
                else if ((1 != idx) && (1 != idt))
                    res = ctr3->xz;
                else
                    res = ctr3->xy;

                if (res && res->getCluster() == cid)
                    return res;
            }
        }
    }
    return NULL;
}

// Looks for the ternary constraint that links this variable with x and y
TernaryConstraint* Variable::getConstr(Variable* x, Variable* y)
{
    TernaryConstraint* ctr;
    for (ConstraintList::iterator iter = constrs.begin(); iter != constrs.end(); ++iter) {
        if ((*iter).constr->isSep() || (*iter).constr->isGlobal())
            continue;

        if ((*iter).constr->isTernary()) {
            ctr = (TernaryConstraint*)(*iter).constr;
            if ((ctr->getIndex(x) >= 0) && (ctr->getIndex(y) >= 0))
                return ctr;
        }
    }
    return NULL;
}

TernaryConstraint* Variable::getConstr(Variable* x, Variable* y, int cid)
{
    TernaryConstraint* ctr = NULL;
    for (ConstraintList::iterator iter = constrs.begin(); iter != constrs.end(); ++iter) {
        if ((*iter).constr->isSep() || (*iter).constr->isGlobal())
            continue;

        if ((*iter).constr->isTernary()) {
            ctr = (TernaryConstraint*)(*iter).constr;
            if ((ctr->getIndex(x) >= 0) && (ctr->getIndex(y) >= 0)) {
                if (ctr->getCluster() == cid)
                    return ctr;
            }
        }
    }
    return NULL;
}

// returns a ternary constraint if the current variable is linked to one
TernaryConstraint* Variable::existTernary()
{
    TernaryConstraint* ctr;
    for (ConstraintList::iterator iter = constrs.begin(); iter != constrs.end(); ++iter) {
        if ((*iter).constr->isSep())
            continue;
        if ((*iter).constr->isTernary()) {
            ctr = (TernaryConstraint*)(*iter).constr;
            return ctr;
        }
    }
    return NULL;
}

double Variable::strongLinkedby(Variable*& strvar, TernaryConstraint*& tctr1max, TernaryConstraint*& tctr2max)
{
    double maxtight = -1;
    strvar = NULL;
    tctr1max = NULL;
    tctr2max = NULL;

    TernaryConstraint* tctr1 = NULL;

    for (ConstraintList::iterator iter = constrs.begin(); iter != constrs.end(); ++iter) {
        if ((*iter).constr->isSep() || (*iter).constr->isGlobal())
            continue;
        if ((*iter).constr->isBinary()) {
            BinaryConstraint* bctr = (BinaryConstraint*)(*iter).constr;
            double bintight = bctr->getTightness();
            if (bintight > maxtight) {
                maxtight = bintight;
                strvar = wcsp->getVar(bctr->getSmallestVarIndexInScope((*iter).scopeIndex));
                tctr1max = NULL;
                tctr2max = NULL;
            }
        } else if ((*iter).constr->isTernary()) {
            double terntight;
            tctr1 = (TernaryConstraint*)(*iter).constr;
            terntight = tctr1->getTightness() + tctr1->xy->getTightness() + tctr1->xz->getTightness() + tctr1->yz->getTightness();

            Variable *x1 = NULL, *x2 = NULL;
            switch ((*iter).scopeIndex) {
            case 0:
                x1 = tctr1->getVar(1);
                x2 = tctr1->getVar(2);
                break;
            case 1:
                x1 = tctr1->getVar(0);
                x2 = tctr1->getVar(2);
                break;
            case 2:
                x1 = tctr1->getVar(0);
                x2 = tctr1->getVar(1);
                break;
            default:;
            }

            if (terntight > maxtight) {
                maxtight = terntight;
                strvar = x1;
                tctr1max = tctr1;
                tctr1max = NULL;
            }

            for (ConstraintList::iterator iter2 = iter; iter2 != constrs.end(); ++iter2) {
                if ((*iter2).constr->isTernary()) {
                    TernaryConstraint* tctr2 = (TernaryConstraint*)(*iter2).constr;
                    Variable* commonvar = NULL;
                    if (tctr2->getIndex(x1) >= 0)
                        commonvar = x1;
                    else if (tctr2->getIndex(x2) >= 0)
                        commonvar = x2;

                    if (commonvar) {
                        terntight += tctr2->getTightness() + tctr2->xy->getTightness() + tctr2->xz->getTightness() + tctr2->yz->getTightness();

                        if (tctr1->xy->getIndex(commonvar) >= 0)
                            terntight -= tctr1->xy->getTightness();
                        else if (tctr1->xz->getIndex(commonvar) >= 0)
                            terntight -= tctr1->xz->getTightness();
                        else if (tctr1->yz->getIndex(commonvar) >= 0)
                            terntight -= tctr1->yz->getTightness();

                        if (terntight > maxtight) {
                            maxtight = terntight;
                            strvar = commonvar;
                            tctr1max = tctr1;
                            tctr2max = tctr2;
                        }
                    }
                }
            }
        }
    }

    return maxtight;
}

// take into account the current tree decomposition of adaptive BTD
bool Variable::isSep()
{
    if (ToulBar2::heuristicFreedom) {
        isSep_ = false;

        TSepLink::iterator it;

        it = clusters.begin();

        while (it != clusters.end()) {
            int c = (*it).first;
            if (wcsp->getTreeDec()->getCluster(c)->getIsCurrInTD()) {
                isSep_ = true;
                return true;
            } else {
                ++it;
            }
        }

        return false;
    } else {
        return isSep_;
    }
}

ostream& operator<<(ostream& os, Variable& var)
{
    os << var.name; // << " #" << var.dac;
    var.print(os);
    if (ToulBar2::verbose >= 3) {
        for (ConstraintList::iterator iter = var.constrs.begin(); iter != var.constrs.end(); ++iter) {
            os << " (" << (*iter).constr << "," << (*iter).scopeIndex << ")";
        }
    }
    return os;
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
