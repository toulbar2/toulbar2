/*
 * **************** Abstract constraints of predefined arities **************
 */

#include "tb2abstractconstr.hpp"
#include "tb2ternaryconstr.hpp"
#include "tb2binconstr.hpp"
#include "search/tb2clusters.hpp"

/*
 * Constructors and misc.
 *
 */

/// \return size of the cartesian product of all initial domains in the constraint scope.
/// \warning use deprecated MAX_DOMAIN_SIZE for performance.
Long AbstractNaryConstraint::getDomainInitSizeProduct()
{
    if (arity_ == 0)
        return 1;
    Long cartesianProduct = 1;
    for (int i = 0; i < arity_; i++) {
        // trap overflow numbers
        if (cartesianProduct > LONGLONG_MAX / MAX_DOMAIN_SIZE)
            return LONGLONG_MAX;
        cartesianProduct *= scope[i]->getDomainInitSize();
    }
    return cartesianProduct;
}

// sorts scope variables by increasing DAC order
int cmpDAC(const void* p1, const void* p2)
{
    EnumeratedVariable* var1 = *((EnumeratedVariable**)p1);
    EnumeratedVariable* var2 = *((EnumeratedVariable**)p2);
    int v1 = var1->getDACOrder();
    int v2 = var2->getDACOrder();
    if (v1 > v2)
        return 1;
    else if (v1 < v2)
        return -1;
    else
        return 0;
}

void AbstractNaryConstraint::firstlex()
{
    it_values.clear();
    EnumeratedVariable* var;
    for (int i = 0; i < arity_; i++) {
        var = (EnumeratedVariable*)getVar(i);
        it_values.push_back(var->begin());
    }
}

bool AbstractNaryConstraint::nextlex(Tuple& t, Cost& c)
{
    int i;
    int a = arity_;
    EnumeratedVariable* var = (EnumeratedVariable*)getVar(0);
    if (it_values[0] == var->end())
        return false;

    t.resize(a);
    for (i = 0; i < a; i++) {
        var = (EnumeratedVariable*)getVar(i);
        t[i] = var->toIndex(*it_values[i]);
    }
    c = eval(t);

    // and now increment
    bool finished = false;
    i = a - 1;
    while (!finished) {
        var = (EnumeratedVariable*)getVar(i);
        ++it_values[i];
        finished = it_values[i] != var->end();
        if (!finished) {
            if (i > 0) {
                it_values[i] = var->begin();
                i--;
            } else
                finished = true;
        }
    }
    return true;
}

template <class T>
Cost AbstractNaryConstraint::evalsubstrAny(const Tuple& s, T* ctr)
{
    int count = 0;

    for (int i = 0; i < arity_; i++) {
        int ind = ctr->getIndex(getVar(i));
        if (ind >= 0) {
            evalTuple[i] = s[ind];
            count++;
        }
    }
    assert(count <= arity_);

    Cost cost;
    if (count == arity_)
        cost = eval(evalTuple);
    else
        cost = MIN_COST;

    return cost;
}
Cost AbstractNaryConstraint::evalsubstr(const Tuple& s, Constraint* ctr) { return evalsubstrAny(s, ctr); }
Cost AbstractNaryConstraint::evalsubstr(const Tuple& s, NaryConstraint* ctr) { return evalsubstrAny(s, ctr); }

// projects n-ary cost function of arity less than 3 into a unary/binary/ternary cost function in extension before the search
void AbstractNaryConstraint::projectNaryBeforeSearch()
{
    Tuple t;
    assert(arity_ <= 3);
    deconnect(); // Warning! It assumes the default cost is not used if the cost function has zero arity
    if (arity_ == 3) {
        vector<Cost> costs;
        EnumeratedVariable* x = (EnumeratedVariable*)getVar(0);
        EnumeratedVariable* y = (EnumeratedVariable*)getVar(1);
        EnumeratedVariable* z = (EnumeratedVariable*)getVar(2);
        unsigned int sizeX = x->getDomainInitSize();
        unsigned int sizeY = y->getDomainInitSize();
        unsigned int sizeZ = z->getDomainInitSize();
        for (unsigned int a = 0; a < sizeX; a++) {
            for (unsigned int b = 0; b < sizeY; b++) {
                for (unsigned int c = 0; c < sizeZ; c++) {
                    costs.push_back(getDefCost());
                }
            }
        }
        Cost cost;
        first();
        while (next(t, cost)) {
            tValue a = t[0];
            tValue b = t[1];
            tValue c = t[2];
            costs[(size_t)(a * sizeY * sizeZ) + (size_t)(b * sizeZ) + (size_t)c] = cost;
        }
        wcsp->postTernaryConstraint(x->wcspIndex, y->wcspIndex, z->wcspIndex, costs);
    } else if (arity_ == 2) {
        vector<Cost> costs;
        EnumeratedVariable* x = (EnumeratedVariable*)getVar(0);
        EnumeratedVariable* y = (EnumeratedVariable*)getVar(1);
        unsigned int sizeX = x->getDomainInitSize();
        unsigned int sizeY = y->getDomainInitSize();
        for (unsigned int a = 0; a < sizeX; a++) {
            for (unsigned int b = 0; b < sizeY; b++) {
                costs.push_back(getDefCost());
            }
        }
        Cost cost;
        first();
        while (next(t, cost)) {
            tValue a = t[0];
            tValue b = t[1];
            costs[(a * sizeY) + b] = cost;
        }
        wcsp->postBinaryConstraint(x->wcspIndex, y->wcspIndex, costs);
    } else if (arity_ == 1) {
        vector<Cost> costs;
        EnumeratedVariable* x = (EnumeratedVariable*)getVar(0);
        unsigned int sizeX = x->getDomainInitSize();
        for (unsigned int a = 0; a < sizeX; a++) {
            costs.push_back(getDefCost());
        }
        Cost cost;
        first();
        while (next(t, cost)) {
            tValue a = t[0];
            costs[a] = cost;
        }
        wcsp->postUnaryConstraint(x->wcspIndex, costs);
    }
}

// USED ONLY DURING SEARCH
void AbstractNaryConstraint::projectNaryTernary(TernaryConstraint* xyz)
{
    TreeDecomposition* td = wcsp->getTreeDec();
    if (td && !ToulBar2::approximateCountingBTD)
        xyz->setCluster(cluster);
    EnumeratedVariable* x = (EnumeratedVariable*)xyz->getVar(0);
    EnumeratedVariable* y = (EnumeratedVariable*)xyz->getVar(1);
    EnumeratedVariable* z = (EnumeratedVariable*)xyz->getVar(2);
    TernaryConstraint* ctr = x->getConstr(y, z);
    if (ctr && td && !td->isSameCluster(ctr->getCluster(), getCluster())) {
        TernaryConstraint* ctr_ = x->getConstr(y, z, getCluster());
        if (ctr_)
            ctr = ctr_;
    }
    if (ToulBar2::verbose >= 2) {
        cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] project nary to ternary (" << x->wcspIndex << "," << y->wcspIndex << "," << z->wcspIndex << ") ";
        if (td)
            cout << "   cluster nary: " << getCluster() << endl;
        else
            cout << endl;
        if (ctr)
            cout << "ctr exists" << endl;
    }
    if (!ctr || (ctr && td && !td->isSameCluster(cluster, ctr->getCluster()))) {
        xyz->fillElimConstrBinaries();
        xyz->reconnect();
        if (ctr)
            xyz->setDuplicate();
    } else {
        ctr->addCosts(xyz);
        xyz = ctr;
    }
    xyz->propagate();
    assert(!td || (td->isSameCluster(xyz->getCluster(), xyz->xy->getCluster()) && td->isSameCluster(xyz->getCluster(), xyz->xz->getCluster()) && td->isSameCluster(xyz->getCluster(), xyz->yz->getCluster())));
}

void AbstractNaryConstraint::projectNaryBinary(BinaryConstraint* xy)
{
    TreeDecomposition* td = wcsp->getTreeDec();
    EnumeratedVariable* x = (EnumeratedVariable*)xy->getVar(0);
    EnumeratedVariable* y = (EnumeratedVariable*)xy->getVar(1);

    if (ToulBar2::verbose >= 2)
        cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] project nary to binary (" << x->wcspIndex << "," << y->wcspIndex << ")" << endl;

    BinaryConstraint* ctr = NULL;
    if (td)
        ctr = x->getConstr(y, getCluster());
    if (!ctr)
        ctr = x->getConstr(y);

    if ((ctr && !td) || (ctr && td && td->isSameCluster(getCluster(), ctr->getCluster()))) {
        if (ToulBar2::verbose >= 2)
            cout << " exists -> fusion" << endl;
        ctr->addCosts(xy);
        xy = ctr;
    } else {
        if (td && !ToulBar2::approximateCountingBTD) {
            if (ctr)
                xy->setDuplicate();
            xy->setCluster(getCluster());
        }
    }
    if (x->unassigned() && y->unassigned())
        xy->reconnect();
    xy->propagate();

    if (ToulBar2::verbose >= 2)
        cout << " and the result: " << *xy << endl;
}

// USED ONLY DURING SEARCH to project the nary constraint
void AbstractNaryConstraint::projectNary()
{
    wcsp->revise(this);
    int indexs[3];
    EnumeratedVariable* unassigned[3] = { NULL, NULL, NULL };
    bool flag = false;

    int i, nunassigned = 0;
    for (i = 0; i < arity_; i++) {
        EnumeratedVariable* var = (EnumeratedVariable*)getVar(i);
        if (var->unassigned()) {
            unassigned[nunassigned] = var;
            indexs[nunassigned] = i;
            evalTuple[i] = 0;
            nunassigned++;
        } else
            evalTuple[i] = var->toIndex(var->getValue());
    }

    EnumeratedVariable* x = unassigned[0];
    EnumeratedVariable* y = unassigned[1];
    EnumeratedVariable* z = unassigned[2];

    assert(nunassigned <= 3);
    if (nunassigned == 3) {
        TernaryConstraint* xyz = wcsp->newTernaryConstr(x, y, z, this);
        wcsp->elimTernOrderInc();
        for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
            for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
                for (EnumeratedVariable::iterator iterz = z->begin(); iterz != z->end(); ++iterz) {
                    Value xval = *iterx;
                    Value yval = *itery;
                    Value zval = *iterz;
                    evalTuple[indexs[0]] = x->toIndex(xval);
                    evalTuple[indexs[1]] = y->toIndex(yval);
                    evalTuple[indexs[2]] = z->toIndex(zval);
                    Cost curcost = eval(evalTuple);
                    if (curcost > MIN_COST)
                        flag = true;
                    xyz->setcost(x, y, z, xval, yval, zval, curcost);
                }
            }
        }
        if (flag)
            projectNaryTernary(xyz);
        // else cout << "ternary empty!" << endl;
    } else if (nunassigned == 2) {
        BinaryConstraint* xy = wcsp->newBinaryConstr(x, y, this);
        wcsp->elimBinOrderInc();
        for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
            for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
                Value xval = *iterx;
                Value yval = *itery;
                evalTuple[indexs[0]] = x->toIndex(xval);
                evalTuple[indexs[1]] = y->toIndex(yval);
                Cost curcost = eval(evalTuple);
                if (curcost > MIN_COST)
                    flag = true;
                xy->setcost(xval, yval, curcost);
                if (ToulBar2::verbose >= 5) {
                    cout << evalTuple;
                    cout << " " << curcost << endl;
                }
            }
        }
        if (flag)
            projectNaryBinary(xy);
        // else cout << "binary empty!" << endl;
    } else if (nunassigned == 1) {
        for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
            Value xval = *iterx;
            evalTuple[indexs[0]] = x->toIndex(xval);
            Cost c = eval(evalTuple);
            if (c > MIN_COST) {
                if (!CUT(c + wcsp->getLb(), wcsp->getUb())) {
                    TreeDecomposition* td = wcsp->getTreeDec();
                    if (td)
                        td->addDelta(cluster, x, xval, c);
                }
                x->project(xval, c);
            }
        }
        x->findSupport();
    } else {
        Cost c = eval(evalTuple);
        Constraint::projectLB(c); // warning! projectLB is not a virtual method
    }
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
