/*
 * **************** Abstract constraint **************
 */

#include "tb2constraint.hpp"
#include "tb2wcsp.hpp"
#include "search/tb2clusters.hpp"
#include "core/tb2knapsack.hpp"

/*
 * Constructor
 *
 */

Constraint::Constraint(WCSP* w)
    : WCSPLink(w, w->numberOfConstraints())
    , conflictWeight(1)
    , fromElim1(NULL)
    , fromElim2(NULL)
{
    w->link(this);
    tight = -1;
    isSep_ = false;
    isDuplicate_ = false;
    cluster = -1;
}

Constraint::Constraint(WCSP* w, int elimCtrIndex)
    : WCSPLink(w, elimCtrIndex)
    , conflictWeight(1)
    , fromElim1(NULL)
    , fromElim2(NULL)
{
    tight = -1;
    isSep_ = false;
    isDuplicate_ = false;
    cluster = -1;
}

bool Constraint::checkEACGreedySolution(int index, Value support)
{
    static Tuple t;
    int a = arity();
    t.resize(a);
    for (int i = 0; i < a; i++) {
        Variable* var = getVar(i);
        if (var->enumerated())
            t[i] = ((EnumeratedVariable*)var)->toIndex((i == index) ? support : var->getSupport());
        else
            t[i] = ((i == index) ? support : var->getSupport());
    }
    return (evalsubstr(t, this) == MIN_COST);
}

bool Constraint::reviseEACGreedySolution(int index, Value support)
{
    bool result = !isGlobal() && checkEACGreedySolution(index, support);
    if (!result) {
        if (index >= 0) {
            getVar(index)->unsetFullEAC();
        } else {
            int a = arity();
            for (int i = 0; i < a; i++) {
                getVar(i)->unsetFullEAC();
            }
        }
    }
    return result;
}

/// \return size of the cartesian product of all domains in the constraint scope.
/// \warning use deprecated MAX_DOMAIN_SIZE for performance.
Long Constraint::getDomainSizeProduct() const
{
    if (arity() == 0)
        return 1;
    Long cartesianProduct = 1;
    for (int i = 0; i < arity(); i++) {
// trap overflow numbers
#if __GNUC__ >= 5
        if (__builtin_smulll_overflow(cartesianProduct,
                getVar(i)->getDomainSize(),
                &cartesianProduct))
            return LONGLONG_MAX;
#else
        if (cartesianProduct > LONGLONG_MAX / MAX_DOMAIN_SIZE)
            return LONGLONG_MAX;
        cartesianProduct *= getVar(i)->getDomainSize();
#endif
    }
    return cartesianProduct;
}

/// \return size of the cartesian product of all initial domains in the constraint scope.
/// \warning use deprecated MAX_DOMAIN_SIZE for performance.
Long Constraint::getDomainInitSizeProduct() const
{
    if (arity() == 0)
        return 1;
    Long cartesianProduct = 1;
    for (int i = 0; i < arity(); i++) {
// trap overflow numbers
#if __GNUC__ >= 5
        if (__builtin_smulll_overflow(cartesianProduct,
                (getVar(i)->enumerated()) ? ((EnumeratedVariable*)getVar(i))->getDomainInitSize() : getVar(i)->getDomainSize(),
                &cartesianProduct))
            return LONGLONG_MAX;
#else
        if (cartesianProduct > LONGLONG_MAX / MAX_DOMAIN_SIZE)
            return LONGLONG_MAX;
        cartesianProduct *= (getVar(i)->enumerated()) ? ((EnumeratedVariable*)getVar(i))->getDomainInitSize() : getVar(i)->getDomainSize();
#endif
    }
    return cartesianProduct;
}

void Constraint::conflict()
{
    wcsp->conflict();
}

void Constraint::projectLB(Cost cost)
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

void Constraint::assigns()
{
    for (int i = 0; connected() && i < arity(); i++)
        if (getVar(i)->assigned())
            assign(i);
}

void Constraint::sumScopeIncluded(Constraint* ctr)
{
    Cost Top = wcsp->getUb();
    Cost c;
    Tuple t;

    if (getDefCost() < Top) { // enumeration case
        firstlex();
        while (nextlex(t, c)) {
            Cost cplus = ctr->evalsubstr(t, this);
            if (c + cplus < Top) {
                if (isNary() && getDefCost() > MIN_COST)
                    setTuple(t, c + cplus);
                else
                    addtoTuple(t, cplus);
            } else {
                setTuple(t, Top);
            }
        }
    } else {
        first();
        while (next(t, c)) {
            Cost cplus = ctr->evalsubstr(t, this);
            if (c + cplus < Top)
                addtoTuple(t, cplus);
            else
                setTuple(t, Top);
        }
    }
}

void Constraint::setCluster(int i)
{
    assert(cluster == -1 || Store::getDepth() == 0);
    if (ToulBar2::verbose >= 1 && cluster != -1 && i != cluster) {
        cout << *this << " change to cluster " << i << endl;
    }
    TreeDecomposition* td = wcsp->getTreeDec();
    if (td && cluster != -1) {
        td->getCluster(cluster)->removeCtr(this);
    }
    cluster = i;
    if (td) {
        td->getCluster(cluster)->addCtr(this);
    }
}

void Constraint::assignCluster()
{
    TreeDecomposition* td = wcsp->getTreeDec();
    assert(!td || !isDuplicate());
    if (!td)
        return;
    Cluster* lowest = td->getRoot();
    for (int i = 0; i < arity(); i++)
        if (getVar(i)->unassigned() || isSep()) { // keep separator Constraint::cluster unchanged
            Variable* x = getVar(i);
            Cluster* c = td->getCluster(x->getCluster());
            if (lowest->isDescendant(c))
                lowest = c;
        }
    setCluster(lowest->getId());
    if (isTernary() && connected() && !isSep()) { // side-effect: if it is a ternary cost function then its internal binary cost functions must belong to the same cluster
        TernaryConstraint* tctr = (TernaryConstraint*)this;
        tctr->xy->setCluster(tctr->getCluster());
        tctr->xz->setCluster(tctr->getCluster());
        tctr->yz->setCluster(tctr->getCluster());
    }
}

/// \warning always returns 0 for cost functions in intention
Cost Constraint::getMinCost()
{
    static Tuple tuple;
    if (!extension())
        return MIN_COST;

    // 	    Cost minc = MAX_COST;
    //         Tuple tuple;
    //         Cost cost;
    //         firstlex();
    //         while (nextlex(tuple,cost)) {
    //             if (cost < minc) minc = cost;
    //         }
    //         return minc;

    Cost minc = MAX_COST;
    Cost cost;
    Long nbtuples = 0;
    first();
    while (next(tuple, cost)) {
        nbtuples++;
        if (cost < minc)
            minc = cost;
    }
    if (getDefCost() < minc && nbtuples < getDomainSizeProduct())
        minc = getDefCost();
    return minc;
}

Cost Constraint::getCost()
{
    static Tuple t;
    int a = arity();
    t.resize(a);
    for (int i = 0; i < a; i++) {
        Variable* var = getVar(i);
        if (var->enumerated())
            t[i] = ((EnumeratedVariable*)var)->toIndex(var->getValue());
        else
            t[i] = var->getValue();
    }
    return evalsubstr(t, this);
}

/// \warning always returns false for cost functions in intention
bool Constraint::universal(Cost zero)
{
    static Tuple tuple;
    if (!extension())
        return false;

    //   Tuple tuple;
    //   Cost cost;
    //   firstlex();
    //   while (nextlex(tuple,cost)) {
    // 	if (cost > zero) return false;
    //   }
    //   return true;

    Cost cost;
    Long nbtuples = 0;
    first();
    while (next(tuple, cost)) {
        nbtuples++;
        if (cost > zero)
            return false;
    }
    if (getDefCost() > zero && nbtuples < getDomainSizeProduct())
        return false;
    return true;
}

/// \warning always returns MAX_COST for cost functions in intention
Cost Constraint::getMaxFiniteCost()
{
    static Tuple tuple;
    if (!extension())
        return MAX_COST;

    Cost maxcost = MIN_COST;
    Cost cost;
    Long nbtuples = 0;
    first();
    while (next(tuple, cost)) {
        nbtuples++;
        if (cost < wcsp->getUb() && cost > maxcost)
            maxcost = cost;
    }
    if (getDefCost() < wcsp->getUb() && getDefCost() > maxcost && nbtuples < getDomainSizeProduct())
        maxcost = getDefCost();
    return maxcost;
}

/// \warning always returns false for cost functions in intention
bool Constraint::ishard()
{
    static Tuple tuple;
    if (!extension())
        return false;

    Cost cost;
    firstlex();
    while (nextlex(tuple, cost)) {
        if (cost > MIN_COST && !CUT(cost, wcsp->getUb() - wcsp->getLb()))
            return false;
    }
    return true;
}

/// \warning always returns false for cost functions in intention
bool Constraint::isfinite()
{
    static Tuple tuple;
    if (!extension())
        return false;

    Cost cost;
    firstlex();
    while (nextlex(tuple, cost)) {
        if (CUT(wcsp->getLb() + cost, wcsp->getUb()))
            return false;
    }
    return true;
}

bool Constraint::verifySeparate(Constraint* ctr1, Constraint* ctr2)
{
    assert(scopeIncluded(ctr1));
    assert(scopeIncluded(ctr2));
    static Tuple tuple;
    Cost cost, c1, c2;
    firstlex();
    if (ToulBar2::verbose >= 3) {
        cout << "[ ";
        for (int i = 0; i < arity(); ++i)
            cout << getVar(i)->getName() << " ";
        cout << " ]\n";
    }
    while (nextlex(tuple, cost)) {
        c1 = ctr1->evalsubstr(tuple, this);
        c2 = ctr2->evalsubstr(tuple, this);
        if (ToulBar2::verbose >= 3) {
            for (int i = 0; i < arity(); ++i)
                cout << tuple[i] << " ";
            // cout << endl;
            cout << " : " << cost << " =? " << c1 << " + " << c2 << " : " << c1 + c2 << endl;
        }
        if (cost < wcsp->getUb() && c1 + c2 != cost)
            return false;
        if (cost >= wcsp->getUb() && c1 + c2 < wcsp->getUb())
            return false;
    }
    return true;
}

bool Constraint::findConditionalIndependences()
{
    bool sep = false;
    if (extension() && !universal() && arity() >= 3 && arity() <= ToulBar2::preprocessNary && ((isTernary() && ((getVar(0)->getTrueDegree() > (min(2, max(ToulBar2::elimDegree, ToulBar2::elimDegree_preprocessing))) && getVar(1)->getTrueDegree() > (min(2, max(ToulBar2::elimDegree, ToulBar2::elimDegree_preprocessing))) && getVar(2)->getTrueDegree() > (min(2, max(ToulBar2::elimDegree, ToulBar2::elimDegree_preprocessing)))))) || (isNary() && (getDomainSizeProduct() < MAX_NB_TUPLES) && (getDefCost() > MIN_COST || ((NaryConstraint*)this)->size() > 1)))) {
        TSCOPE scopeinv;
        getScope(scopeinv);
        EnumeratedVariable* vx = NULL;
        EnumeratedVariable* vz = NULL;
        for (TSCOPE::reverse_iterator it1 = scopeinv.rbegin(); it1 != scopeinv.rend() && !sep; ++it1) {
            TSCOPE::reverse_iterator it2 = it1;
            for (++it2; it2 != scopeinv.rend() && !sep; ++it2) {
                vx = (EnumeratedVariable*)wcsp->getVar((*it2).first);
                vz = (EnumeratedVariable*)wcsp->getVar((*it1).first);
                if (ToulBar2::verbose >= 1)
                    cout << /*"\n" <<*/ vx->getName() << " and " << vz->getName() << " are separable in ";
                sep = separability(vx, vz);
                if (sep && ToulBar2::verbose >= 1) {
                    cout << " YES";
#ifndef NDEBUG
                    if (!ishard())
                        cout << " with finite costs";
#endif
                    cout << endl;
                }
                if (!sep && ToulBar2::verbose >= 1)
                    cout << " NO" << endl;
            }
        }
        if (sep)
            separate(vx, vz);
        if (ToulBar2::verbose >= 3)
            cout << "=====================================================" << endl;
    }
    return sep;
}

Constraint* Constraint::copy()
{
    static Tuple t;
    int scope[arity()];
    for (int i = 0; i < arity(); i++)
        scope[i] = getVar(i)->wcspIndex;
    Cost defcost = getDefCost();
    int ctrIndex = wcsp->postNaryConstraintBegin(scope, arity(), defcost, size(), true); // be sure to not create a clause instead of NaryConstraint!
    Cost c;
    first();
    while (next(t, c)) {
        if (c != defcost)
            wcsp->postNaryConstraintTupleInternal(ctrIndex, t, c);
    }
    wcsp->getCtr(ctrIndex)->deconnect(true);
    return wcsp->getCtr(ctrIndex);
}

bool Constraint::cmpConstraintId(Constraint* c1, Constraint* c2)
{
    int v1 = c1->getSmallestVarIndexInScope();
    int v2 = c2->getSmallestVarIndexInScope();
    return (v1 < v2);
}

bool Constraint::cmpConstraintId(DLink<ConstraintLink>* c1, DLink<ConstraintLink>* c2)
{
    int v1 = c1->content.constr->getSmallestVarIndexInScope(c1->content.scopeIndex);
    int v2 = c2->content.constr->getSmallestVarIndexInScope(c2->content.scopeIndex);
    return (v1 < v2);
}

bool Constraint::cmpConstraintId(DLink<Constraint*>* c1, DLink<Constraint*>* c2)
{
    int v1 = c1->content->getSmallestVarIndexInScope();
    int v2 = c2->content->getSmallestVarIndexInScope();
    return (v1 < v2);
}

bool Constraint::cmpConstraintDAC(Constraint* c1, Constraint* c2)
{
    int v1 = c1->getDACVar(0)->getDACOrder();
    int v2 = c2->getDACVar(0)->getDACOrder();
    return (v1 > v2);
}

bool Constraint::cmpConstraintDAC(DLink<ConstraintLink>* c1, DLink<ConstraintLink>* c2)
{
    int v1 = c1->content.constr->getSmallestDACIndexInScope(c1->content.scopeIndex);
    int v2 = c2->content.constr->getSmallestDACIndexInScope(c2->content.scopeIndex);
    return (v1 > v2);
}

bool Constraint::cmpConstraintDAC(DLink<Constraint*>* c1, DLink<Constraint*>* c2)
{
    int v1 = c1->content->getDACVar(0)->getDACOrder();
    int v2 = c2->content->getDACVar(0)->getDACOrder();
    return (v1 > v2);
}

bool Constraint::cmpConstraintTightness(Constraint* c1, Constraint* c2)
{
    double v1 = c1->getTightness();
    double v2 = c2->getTightness();
    return (v1 > v2);
}

bool Constraint::cmpConstraintTightness(DLink<ConstraintLink>* c1, DLink<ConstraintLink>* c2)
{
    double v1 = c1->content.constr->getTightness();
    double v2 = c2->content.constr->getTightness();
    return (v1 > v2);
}

bool Constraint::cmpConstraintTightness(DLink<Constraint*>* c1, DLink<Constraint*>* c2)
{
    double v1 = c1->content->getTightness();
    double v2 = c2->content->getTightness();
    return (v1 > v2);
}

bool Constraint::cmpConstraintDACTightness(Constraint* c1, Constraint* c2)
{
    int v1 = c1->getDACVar(0)->getDACOrder();
    int v2 = c2->getDACVar(0)->getDACOrder();
    if (v1 != v2)
        return (v1 > v2);
    else {
        double v1 = c1->getTightness();
        double v2 = c2->getTightness();
        return (v1 > v2);
    }
}

bool Constraint::cmpConstraintDACTightness(DLink<ConstraintLink>* c1, DLink<ConstraintLink>* c2)
{
    int v1 = c1->content.constr->getSmallestDACIndexInScope(c1->content.scopeIndex);
    int v2 = c2->content.constr->getSmallestDACIndexInScope(c2->content.scopeIndex);
    if (v1 != v2)
        return (v1 > v2);
    else {
        double v1 = c1->content.constr->getTightness();
        double v2 = c2->content.constr->getTightness();
        return (v1 > v2);
    }
}

bool Constraint::cmpConstraintDACTightness(DLink<Constraint*>* c1, DLink<Constraint*>* c2)
{
    int v1 = c1->content->getDACVar(0)->getDACOrder();
    int v2 = c2->content->getDACVar(0)->getDACOrder();
    if (v1 != v2)
        return (v1 > v2);
    else {
        double v1 = c1->content->getTightness();
        double v2 = c2->content->getTightness();
        return (v1 > v2);
    }
}

bool Constraint::cmpConstraintTightnessDAC(Constraint* c1, Constraint* c2)
{
    double v1 = c1->getTightness();
    double v2 = c2->getTightness();
    if (v1 != v2)
        return (v1 > v2);
    else {
        int v1 = c1->getDACVar(0)->getDACOrder();
        int v2 = c2->getDACVar(0)->getDACOrder();
        return (v1 > v2);
    }
}

bool Constraint::cmpConstraintTightnessDAC(DLink<ConstraintLink>* c1, DLink<ConstraintLink>* c2)
{
    double v1 = c1->content.constr->getTightness();
    double v2 = c2->content.constr->getTightness();
    if (v1 != v2)
        return (v1 > v2);
    else {
        int v1 = c1->content.constr->getSmallestDACIndexInScope(c1->content.scopeIndex);
        int v2 = c2->content.constr->getSmallestDACIndexInScope(c2->content.scopeIndex);
        return (v1 > v2);
    }
}

bool Constraint::cmpConstraintTightnessDAC(DLink<Constraint*>* c1, DLink<Constraint*>* c2)
{
    double v1 = c1->content->getTightness();
    double v2 = c2->content->getTightness();
    if (v1 != v2)
        return (v1 > v2);
    else {
        int v1 = c1->content->getDACVar(0)->getDACOrder();
        int v2 = c2->content->getDACVar(0)->getDACOrder();
        return (v1 > v2);
    }
}

bool Constraint::cmpConstraintLAG(Constraint* c1, Constraint* c2)
{
    auto* k = dynamic_cast<KnapsackConstraint*>(c1);
    auto* k1 = dynamic_cast<KnapsackConstraint*>(c2);
    if (k && !k1)
        return false;
    if (!k && k1)
        return true;
    if (k && k1) {
        Double v1 = k->getLag();
        Double v2 = k1->getLag();
        if (v1 == v2)
            return k->computeTightness() < k1->computeTightness();
        return (v1 < v2);
    } else
        return cmpConstraintDAC(c1, c2);
}

bool Constraint::cmpConstraintLAG(DLink<ConstraintLink>* c1, DLink<ConstraintLink>* c2)
{
    auto* k = dynamic_cast<KnapsackConstraint*>(c1->content.constr);
    auto* k1 = dynamic_cast<KnapsackConstraint*>(c2->content.constr);
    if (k && !k1)
        return false;
    if (!k && k1)
        return true;
    if (k && k1) {
        Double v1 = k->getLag();
        Double v2 = k1->getLag();
        if (v1 == v2)
            return k->computeTightness() < k1->computeTightness();
        return (v1 < v2);
    } else
        return cmpConstraintDAC(c1, c2);
}

bool Constraint::cmpConstraintLAG(DLink<Constraint*>* c1, DLink<Constraint*>* c2)
{
    auto* k = dynamic_cast<KnapsackConstraint*>(c1->content);
    auto* k1 = dynamic_cast<KnapsackConstraint*>(c2->content);
    if (k && !k1)
        return false;
    if (!k && k1)
        return true;
    if (k && k1) {
        Double v1 = k->getLag();
        Double v2 = k1->getLag();
        if (v1 == v2)
            return k->computeTightness() < k1->computeTightness();
        return (v1 < v2);
    } else
        return cmpConstraintDAC(c1, c2);
}

bool Constraint::cmpConstraintArity(Constraint* c1, Constraint* c2)
{
    int v1 = c1->arity();
    int v2 = c2->arity();
    return (v1 < v2);
}

bool Constraint::cmpConstraintArity(DLink<ConstraintLink>* c1, DLink<ConstraintLink>* c2)
{
    int v1 = c1->content.constr->arity();
    int v2 = c2->content.constr->arity();
    return (v1 < v2);
}

bool Constraint::cmpConstraintArity(DLink<Constraint*>* c1, DLink<Constraint*>* c2)
{
    int v1 = c1->content->arity();
    int v2 = c2->content->arity();
    return (v1 < v2);
}

bool Constraint::cmpConstraintArityDAC(Constraint* c1, Constraint* c2)
{
    int v1 = c1->arity();
    int v2 = c2->arity();
    if (v1 != v2)
        return (v1 < v2);
    else {
        int v1 = c1->getDACVar(0)->getDACOrder();
        int v2 = c2->getDACVar(0)->getDACOrder();
        return (v1 > v2);
    }
}

bool Constraint::cmpConstraintArityDAC(DLink<ConstraintLink>* c1, DLink<ConstraintLink>* c2)
{
    int v1 = c1->content.constr->arity();
    int v2 = c2->content.constr->arity();
    if (v1 != v2)
        return (v1 < v2);
    else {
        int v1 = c1->content.constr->getSmallestDACIndexInScope(c1->content.scopeIndex);
        int v2 = c2->content.constr->getSmallestDACIndexInScope(c2->content.scopeIndex);
        return (v1 > v2);
    }
}

bool Constraint::cmpConstraintArityDAC(DLink<Constraint*>* c1, DLink<Constraint*>* c2)
{
    int v1 = c1->content->arity();
    int v2 = c2->content->arity();
    if (v1 != v2)
        return (v1 < v2);
    else {
        int v1 = c1->content->getDACVar(0)->getDACOrder();
        int v2 = c2->content->getDACVar(0)->getDACOrder();
        return (v1 > v2);
    }
}

// sort a list of constraints
int Constraint::cmpConstraint(Constraint* c1, Constraint* c2)
{
    assert(c1->arity() > 0 && c2->arity() > 0);
    bool result = false;
    switch (abs(ToulBar2::constrOrdering)) {
    case CONSTR_ORDER_ID:
        result = cmpConstraintId(c1, c2);
        break;
    case CONSTR_ORDER_DAC:
        result = cmpConstraintDAC(c1, c2);
        break;
    case CONSTR_ORDER_TIGHTNESS:
        result = cmpConstraintTightness(c1, c2);
        break;
    case CONSTR_ORDER_DAC_TIGHTNESS:
        result = cmpConstraintDACTightness(c1, c2);
        break;
    case CONSTR_ORDER_TIGHTNESS_DAC:
        result = cmpConstraintTightnessDAC(c1, c2);
        break;
    case CONSTR_ORDER_LAG:
        result = cmpConstraintTightnessDAC(c1, c2);
        break;
    case CONSTR_ORDER_ARITY:
        result = cmpConstraintArity(c1, c2);
        break;
    case CONSTR_ORDER_ARITY_DAC:
        result = cmpConstraintArityDAC(c1, c2);
        break;
    default:
        cerr << "Unknown constraint ordering value " << ToulBar2::constrOrdering << endl;
        throw BadConfiguration();
    }
    if (ToulBar2::constrOrdering >= 0) {
        return result;
    } else {
        return (!result);
    }
}

// sort a list of constraints related to a given variable
int Constraint::cmpConstraintLink(DLink<ConstraintLink>* c1, DLink<ConstraintLink>* c2)
{
    assert(c1->content.constr->arity() > 0 && c2->content.constr->arity() > 0);
    bool result = false;
    switch (abs(ToulBar2::constrOrdering)) {
    case CONSTR_ORDER_ID:
        result = cmpConstraintId(c1, c2);
        break;
    case CONSTR_ORDER_DAC:
        result = cmpConstraintDAC(c1, c2);
        break;
    case CONSTR_ORDER_TIGHTNESS:
        result = cmpConstraintTightness(c1, c2);
        break;
    case CONSTR_ORDER_DAC_TIGHTNESS:
        result = cmpConstraintDACTightness(c1, c2);
        break;
    case CONSTR_ORDER_TIGHTNESS_DAC:
        result = cmpConstraintTightnessDAC(c1, c2);
        break;
    case CONSTR_ORDER_LAG:
        result = cmpConstraintLAG(c1, c2);
        break;
    case CONSTR_ORDER_ARITY:
        result = cmpConstraintArity(c1, c2);
        break;
    case CONSTR_ORDER_ARITY_DAC:
        result = cmpConstraintArityDAC(c1, c2);
        break;
    default:
        cerr << "Unknown constraint ordering value " << ToulBar2::constrOrdering << endl;
        throw BadConfiguration();
    }
    if (ToulBar2::constrOrdering >= 0) {
        return result;
    } else {
        return (!result);
    }
}

// sort a bactrackable list of constraints
int Constraint::cmpConstraintLinkPointer(DLink<Constraint*>* c1, DLink<Constraint*>* c2)
{
    assert(c1->content->arity() > 0 && c2->content->arity() > 0);
    bool result = false;
    switch (abs(ToulBar2::constrOrdering)) {
    case CONSTR_ORDER_ID:
        result = cmpConstraintId(c1, c2);
        break;
    case CONSTR_ORDER_DAC:
        result = cmpConstraintDAC(c1, c2);
        break;
    case CONSTR_ORDER_TIGHTNESS:
        result = cmpConstraintTightness(c1, c2);
        break;
    case CONSTR_ORDER_DAC_TIGHTNESS:
        result = cmpConstraintDACTightness(c1, c2);
        break;
    case CONSTR_ORDER_TIGHTNESS_DAC:
        result = cmpConstraintTightnessDAC(c1, c2);
        break;
    case CONSTR_ORDER_LAG:
        result = cmpConstraintLAG(c1, c2);
        break;
    case CONSTR_ORDER_ARITY:
        result = cmpConstraintArity(c1, c2);
        break;
    case CONSTR_ORDER_ARITY_DAC:
        result = cmpConstraintArityDAC(c1, c2);
        break;
    default:
        cerr << "Unknown constraint ordering value " << ToulBar2::constrOrdering << endl;
        throw BadConfiguration();
    }
    if (ToulBar2::constrOrdering >= 0) {
        return result;
    } else {
        return (!result);
    }
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
