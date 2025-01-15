/*
 * ****** Binary constraint applied on variables with enumerated domains ******
 */

#include "tb2binconstr.hpp"
#include "tb2wcsp.hpp"
#include "search/tb2clusters.hpp"
#include "tb2vacutils.hpp"

/*
 * Constructors and misc.
 *
 */

BinaryConstraint::BinaryConstraint(WCSP* wcsp, EnumeratedVariable* xx, EnumeratedVariable* yy, vector<Cost>& tab)
    : AbstractBinaryConstraint<EnumeratedVariable, EnumeratedVariable>(wcsp, xx, yy)
    , sizeX(xx->getDomainInitSize())
    , sizeY(yy->getDomainInitSize())
{
    deltaCostsX = vector<StoreCost>(sizeX, StoreCost(MIN_COST));
    deltaCostsY = vector<StoreCost>(sizeY, StoreCost(MIN_COST));
    assert(tab.size() == sizeX * sizeY);
    supportX = vector<Value>(sizeX, y->getInf());
    supportY = vector<Value>(sizeY, x->getInf());
    trwsM = vector<Cost>(max(sizeX, sizeY), MIN_COST);

#ifdef NO_STORE_BINARY_COSTS
    costs = vector<Cost>((size_t)sizeX * (size_t)sizeY, MIN_COST);
#else
    costs = vector<StoreCost>((size_t)sizeX * (size_t)sizeY, StoreCost(MIN_COST));
#endif
    for (unsigned int a = 0; a < x->getDomainInitSize(); a++)
        for (unsigned int b = 0; b < y->getDomainInitSize(); b++)
            costs[a * sizeY + b] = tab[a * sizeY + b];
    if (ToulBar2::bilevel >= 2) {
        assert(ToulBar2::bilevel <= 3);
        deconnect(true);
        wcsp->delayedCtrBLP[ToulBar2::bilevel - 1].push_back(wcspIndex);
    } else {
        propagate();
    }
}

BinaryConstraint::BinaryConstraint(WCSP* wcsp)
    : AbstractBinaryConstraint<EnumeratedVariable, EnumeratedVariable>(wcsp)
    , sizeX(0)
    , sizeY(0)
{
    //	unsigned int maxdomainsize = wcsp->getMaxDomainSize();
    //    deltaCostsX = vector<StoreCost>(maxdomainsize,StoreCost(MIN_COST,storeCost));
    //    deltaCostsY = vector<StoreCost>(maxdomainsize,StoreCost(MIN_COST,storeCost));
    //    supportX = vector<Value>(maxdomainsize,0);
    //    supportY = vector<Value>(maxdomainsize,0);
    linkX = new DLink<ConstraintLink>;
    linkY = new DLink<ConstraintLink>;

    //    costs = vector<StoreCost>(maxdomainsize*maxdomainsize,StoreCost(MIN_COST,storeCost));
    //    for (unsigned int a = 0; a < maxdomainsize; a++)
    //         for (unsigned int b = 0; b < maxdomainsize; b++)
    //                costs[a * maxdomainsize + b] = MIN_COST;
}

void BinaryConstraint::print(ostream& os)
{
    os << this << " BinaryConstraint(" << x->getName() << "," << y->getName() << ")";
    if (ToulBar2::weightedDegree)
        os << "/" << getConflictWeight();
    if (wcsp->getTreeDec()) {
        os << "   cluster: " << getCluster();
        if (ToulBar2::heuristicFreedom && getCluster() != -1) {
            cout << "  freedom: " << wcsp->getTreeDec()->getCluster(getCluster())->getFreedom() << " isinTD: " << wcsp->getTreeDec()->getCluster(getCluster())->getIsCurrInTD();
        }
    }
    if (ToulBar2::verbose >= 8) {
        os << " supportX:[";
        for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
            os << " " << supportX[x->toIndex(*iterX)];
        }
        os << " ] supportY:[";
        for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
            os << " " << supportY[y->toIndex(*iterY)];
        }
        os << " ]";
        os << " deltaX:[";
        for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
            os << " " << deltaCostsX[x->toIndex(*iterX)];
        }
        os << " ] deltaY:[";
        for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
            os << " " << deltaCostsY[y->toIndex(*iterY)];
        }
        os << " ]";
    }
    os << endl;
    if (ToulBar2::verbose >= 5) {
        for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
            for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                os << " " << getCost(*iterX, *iterY);
            }
            os << endl;
        }
    }
}

void BinaryConstraint::dump(ostream& os, bool original)
{
    unsigned int tuples = 0;
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
            if (getCost(*iterX, *iterY) > MIN_COST) {
                tuples++;
            }
        }
    }
    os << "2 " << ((original) ? (x->wcspIndex) : x->getCurrentVarId()) << " " << ((original) ? (y->wcspIndex) : y->getCurrentVarId()) << " " << MIN_COST << " " << tuples << endl;
    int i = 0;
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX, i++) {
        int j = 0;
        for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY, j++) {
            if (getCost(*iterX, *iterY) > MIN_COST) {
                os << ((original) ? *iterX : i) << " " << ((original) ? *iterY : j) << " " << ((original) ? getCost(*iterX, *iterY) : min(wcsp->getUb(), getCost(*iterX, *iterY))) << endl;
            }
        }
    }
}

void BinaryConstraint::dump_CFN(ostream& os, bool original)
{
    bool printed = false;
    os << "\"F_" << ((original) ? (x->wcspIndex) : x->getCurrentVarId()) << "_" << ((original) ? (y->wcspIndex) : y->getCurrentVarId()) << "\":{\"scope\":[";
    os << "\"" << name2cfn(x->getName()) << "\",\"" << name2cfn(y->getName()) << "\"],";
    os << "\"defaultcost\":" << MIN_COST << ",\n\"costs\":[\n";
    int i = 0;
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX, i++) {
        int j = 0;
        for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY, j++) {
            if (getCost(*iterX, *iterY) != MIN_COST) {
                if (printed)
                    os << ",\n";
                os << ((original) ? x->toIndex(*iterX) : i) << "," << ((original) ? y->toIndex(*iterY) : j) << ","
                   << ((original) ? wcsp->DCost2Decimal(wcsp->Cost2RDCost(getCost(*iterX, *iterY))) : ((wcsp->getUb() > getCost(*iterX, *iterY)) ? wcsp->DCost2Decimal(wcsp->Cost2RDCost(getCost(*iterX, *iterY))) : "inf"));
                printed = true;
            }
        }
    }
    os << "\n]},\n";
}

/*
 * Propagation methods
 *
 */
bool BinaryConstraint::project(EnumeratedVariable* x, Value value, Cost cost, vector<StoreCost>& deltaCostsX)
{
    assert(ToulBar2::verbose < 4 || ((cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] project(C" << getVar(0)->getName() << "," << getVar(1)->getName() << ", (" << x->getName() << "," << value << "), " << cost << ")" << endl), true));

    // hard binary constraint costs are not changed
    if (!CUT(cost + wcsp->getLb(), wcsp->getUb())) {
        TreeDecomposition* td = wcsp->getTreeDec();
        if (td)
            td->addDelta(cluster, x, value, cost);
        deltaCostsX[x->toIndex(value)] += cost; // Warning! Possible overflow???
        assert(getCost(x, (EnumeratedVariable*)getVarDiffFrom(x), value, getVarDiffFrom(x)->getInf()) >= MIN_COST);
        assert(getCost(x, (EnumeratedVariable*)getVarDiffFrom(x), value, getVarDiffFrom(x)->getSup()) >= MIN_COST);
    }

    Cost oldcost = x->getCost(value);
    x->project(value, cost);
#ifdef DEECOMPLETE
    getVarDiffFrom(x)->queueDEE();
#endif
    return (x->getSupport() == value || SUPPORTTEST(oldcost, cost));
}

void BinaryConstraint::extend(EnumeratedVariable* x, Value value, Cost cost, vector<StoreCost>& deltaCostsX)
{
    assert(ToulBar2::verbose < 4 || ((cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] extend(C" << getVar(0)->getName() << "," << getVar(1)->getName() << ", (" << x->getName() << "," << value << "), " << cost << ")" << endl), true));

    TreeDecomposition* td = wcsp->getTreeDec();
    if (td)
        td->addDelta(cluster, x, value, -cost);

    deltaCostsX[x->toIndex(value)] -= cost; // Warning! Possible overflow???
    x->extend(value, cost);
}

void BinaryConstraint::permute(EnumeratedVariable* xin, Value a, Value b)
{
    EnumeratedVariable* yin = y;
    if (xin != x)
        yin = x;

    vector<Cost> aux;
    for (EnumeratedVariable::iterator ity = yin->begin(); ity != yin->end(); ++ity)
        aux.push_back(getCost(xin, yin, a, *ity));
    for (EnumeratedVariable::iterator ity = yin->begin(); ity != yin->end(); ++ity)
        setcost(xin, yin, a, *ity, getCost(xin, yin, b, *ity));

    vector<Cost>::iterator itc = aux.begin();
    for (EnumeratedVariable::iterator ity = yin->begin(); ity != yin->end(); ++ity) {
        setcost(xin, yin, b, *ity, *itc);
        ++itc;
    }
}

bool BinaryConstraint::isFunctional(EnumeratedVariable* xin, EnumeratedVariable* yin, map<Value, Value>& functional)
{
    assert(xin != yin);
    assert(getIndex(xin) >= 0);
    assert(getIndex(yin) >= 0);
    bool isfunctional = true;
    functional.clear();
    for (EnumeratedVariable::iterator itx = xin->begin(); isfunctional && itx != xin->end(); ++itx) {
        bool first = true;
        for (EnumeratedVariable::iterator ity = yin->begin(); isfunctional && ity != yin->end(); ++ity) {
            if (!CUT(getCost(xin, yin, *itx, *ity) + wcsp->getLb(), wcsp->getUb())) {
                if (first) {
                    functional[*itx] = *ity;
                    first = false;
                } else {
                    isfunctional = false;
                }
            }
        }
        assert(!first); // assumes it is SAC already
    }
    if (isfunctional)
        return true;
    functional.clear();
    return false;
}

pair<pair<Cost, Cost>, pair<Cost, Cost>> BinaryConstraint::getMaxCost(int varIndex, Value a, Value b)
{
    //    	cout << "getMaxCost(" << getVar(varIndex)->getName() << ") " << a << " <-> " << b << endl << *this << endl;
    Cost maxcosta = MIN_COST;
    Cost diffcosta = MIN_COST;
    Cost maxcostb = MIN_COST;
    Cost diffcostb = MIN_COST;
    if (varIndex == 0) {
        Cost ucosta = x->getCost(a);
        Cost ucostb = x->getCost(b);
        for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
            Cost costa = getCost(a, *iterY);
            Cost costb = getCost(b, *iterY);
            if (costa > maxcosta)
                maxcosta = costa;
            if (costb > maxcostb)
                maxcostb = costb;
            Cost ucosty = y->getCost(*iterY);
            if (!CUT(ucostb + costb + ucosty + wcsp->getLb(), wcsp->getUb())) {
                if (costa - costb > diffcosta)
                    diffcosta = costa - costb;
            }
            if (!CUT(ucosta + costa + ucosty + wcsp->getLb(), wcsp->getUb())) {
                if (costb - costa > diffcostb)
                    diffcostb = costb - costa;
            }
        }
    } else {
        Cost ucosta = y->getCost(a);
        Cost ucostb = y->getCost(b);
        for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
            Cost costa = getCost(*iterX, a);
            Cost costb = getCost(*iterX, b);
            if (costa > maxcosta)
                maxcosta = costa;
            if (costb > maxcostb)
                maxcostb = costb;
            Cost ucostx = x->getCost(*iterX);
            if (!CUT(ucostb + costb + ucostx + wcsp->getLb(), wcsp->getUb())) {
                if (costa - costb > diffcosta)
                    diffcosta = costa - costb;
            }
            if (!CUT(ucosta + costa + ucostx + wcsp->getLb(), wcsp->getUb())) {
                if (costb - costa > diffcostb)
                    diffcostb = costb - costa;
            }
        }
    }
    assert(maxcosta >= diffcosta);
    assert(maxcostb >= diffcostb);
    return make_pair(make_pair(maxcosta, diffcosta), make_pair(maxcostb, diffcostb));
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
