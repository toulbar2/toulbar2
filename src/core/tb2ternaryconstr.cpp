/*
 * ****** Binary constraint applied on variables with enumerated domains ******
 */

#include "tb2ternaryconstr.hpp"
#include "tb2wcsp.hpp"
#include "search/tb2clusters.hpp"

/*
 * Constructors and misc.
 *
 */

TernaryConstraint::TernaryConstraint(WCSP* wcsp,
    EnumeratedVariable* xx,
    EnumeratedVariable* yy,
    EnumeratedVariable* zz,
    BinaryConstraint* xy_,
    BinaryConstraint* xz_,
    BinaryConstraint* yz_,
    vector<Cost>& tab)
    : AbstractTernaryConstraint<EnumeratedVariable, EnumeratedVariable, EnumeratedVariable>(wcsp, xx, yy, zz)
    , sizeX(xx->getDomainInitSize())
    , sizeY(yy->getDomainInitSize())
    , sizeZ(zz->getDomainInitSize())
    , top(wcsp->getUb())
    , functionalX(true)
    , functionalY(true)
    , functionalZ(true)
{
    assert(tab.size() == (size_t)sizeX * (size_t)sizeY * (size_t)sizeZ);
    deltaCostsX = vector<StoreCost>(sizeX, StoreCost(MIN_COST));
    deltaCostsY = vector<StoreCost>(sizeY, StoreCost(MIN_COST));
    deltaCostsZ = vector<StoreCost>(sizeZ, StoreCost(MIN_COST));
    assert(getIndex(x) < getIndex(y) && getIndex(y) < getIndex(z));
    functionX = vector<Value>((size_t)sizeY * (size_t)sizeZ, WRONG_VAL);
    functionY = vector<Value>((size_t)sizeX * (size_t)sizeZ, WRONG_VAL);
    functionZ = vector<Value>((size_t)sizeX * (size_t)sizeY, WRONG_VAL);
    supportX = vector<pair<Value, Value>>(sizeX, make_pair(y->getInf(), z->getInf()));
    supportY = vector<pair<Value, Value>>(sizeY, make_pair(x->getInf(), z->getInf()));
    supportZ = vector<pair<Value, Value>>(sizeZ, make_pair(x->getInf(), y->getInf()));

    for (unsigned int a = 0; a < sizeX; a++) {
        for (unsigned int b = 0; b < sizeY; b++) {
            for (unsigned int c = 0; c < sizeZ; c++) {
                Cost cost = tab[(size_t)a * (size_t)sizeY * (size_t)sizeZ + (size_t)b * (size_t)sizeZ + (size_t)c];
                if (!CUT(cost, top) && x->canbe(x->toValue(a)) && y->canbe(y->toValue(b)) && z->canbe(z->toValue(c))) {
                    if (functionalX) {
                        if (functionX[b * sizeZ + c] == WRONG_VAL)
                            functionX[b * sizeZ + c] = x->toValue(a);
                        else
                            functionalX = false;
                    }
                    if (functionalY) {
                        if (functionY[a * sizeZ + c] == WRONG_VAL)
                            functionY[a * sizeZ + c] = y->toValue(b);
                        else
                            functionalY = false;
                    }
                    if (functionalZ) {
                        if (functionZ[a * sizeY + b] == WRONG_VAL)
                            functionZ[a * sizeY + b] = z->toValue(c);
                        else
                            functionalZ = false;
                    }
                }
            }
        }
    }

    setBinaries(xy_, xz_, yz_);

    if (functionalX) {
#ifdef NO_STORE_TERNARY_COSTS
        costsYZ = vector<Cost>((size_t)sizeY * (size_t)sizeZ, MIN_COST);
#else
        costsYZ = vector<StoreCost>((size_t)sizeY * (size_t)sizeZ, StoreCost(MIN_COST));
#endif
        for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
            for (unsigned int c = 0; c < z->getDomainInitSize(); c++) {
                if (functionX[b * sizeZ + c] != WRONG_VAL)
                    costsYZ[b * sizeZ + c] = tab[x->toIndex(functionX[b * sizeZ + c]) * sizeY * sizeZ + b * sizeZ + c];
            }
        }
        //    	costs.free_all();
    } else {
#ifdef NO_STORE_TERNARY_COSTS
        costs = vector<Cost>((size_t)sizeX * (size_t)sizeY * (size_t)sizeZ, MIN_COST);
#else
        costs = vector<StoreCost>((size_t)sizeX * (size_t)sizeY * (size_t)sizeZ, StoreCost(MIN_COST));
#endif
        for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
            for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
                for (unsigned int c = 0; c < z->getDomainInitSize(); c++) {
                    costs[(size_t)a * sizeY * sizeZ + (size_t)b * sizeZ + c] = tab[(size_t)a * sizeY * sizeZ + (size_t)b * sizeZ + c];
                }
            }
        }
    }

// Uncomment the following code if toulbar2 is used within numberjack
#ifdef NUMBERJACK
    vector<int>& vecX = wcsp->getListSuccessors()->at(x->wcspIndex);
    vector<int>& vecY = wcsp->getListSuccessors()->at(y->wcspIndex);
    vector<int>& vecZ = wcsp->getListSuccessors()->at(z->wcspIndex);
    // If variables x,y,z are not yet involved by a decomposable global cost function and there is a functional dependency between them
    // then suggests a good Berge DAC ordering:
    if ((std::find(vecX.begin(), vecX.end(), y->wcspIndex) == vecX.end()) && (std::find(vecX.begin(), vecX.end(), z->wcspIndex) == vecX.end()) && (std::find(vecY.begin(), vecY.end(), x->wcspIndex) == vecY.end()) && (std::find(vecY.begin(), vecY.end(), z->wcspIndex) == vecY.end()) && (std::find(vecZ.begin(), vecZ.end(), x->wcspIndex) == vecZ.end()) && (std::find(vecZ.begin(), vecZ.end(), y->wcspIndex) == vecZ.end())) {
        switch (functionalX + functionalY + functionalZ) {
        case 1:
            if (functionalX) {
                vecX.push_back(y->wcspIndex);
                vecX.push_back(z->wcspIndex);
            } else if (functionalY) {
                vecY.push_back(x->wcspIndex);
                vecY.push_back(z->wcspIndex);
            } else if (functionalZ) {
                vecZ.push_back(x->wcspIndex);
                vecZ.push_back(y->wcspIndex);
            }
            ToulBar2::Berge_Dec = true;
            break;
        case 2:
            if (functionalX && functionalY) {
                if (x->wcspIndex < y->wcspIndex) {
                    vecX.push_back(z->wcspIndex);
                    vecZ.push_back(y->wcspIndex);
                } else {
                    vecY.push_back(z->wcspIndex);
                    vecZ.push_back(x->wcspIndex);
                }
            } else if (functionalX && functionalZ) {
                if (x->wcspIndex < z->wcspIndex) {
                    vecX.push_back(y->wcspIndex);
                    vecY.push_back(z->wcspIndex);
                } else {
                    vecZ.push_back(y->wcspIndex);
                    vecY.push_back(x->wcspIndex);
                }
            } else if (functionalY && functionalZ) {
                if (y->wcspIndex < z->wcspIndex) {
                    vecY.push_back(x->wcspIndex);
                    vecX.push_back(z->wcspIndex);
                } else {
                    vecZ.push_back(x->wcspIndex);
                    vecX.push_back(y->wcspIndex);
                }
            }
            ToulBar2::Berge_Dec = true;
            break;
        case 3:
            if (x->wcspIndex < y->wcspIndex && x->wcspIndex < z->wcspIndex) {
                if (y->wcspIndex < z->wcspIndex) {
                    vecY.push_back(x->wcspIndex);
                    vecX.push_back(z->wcspIndex);
                } else {
                    vecZ.push_back(x->wcspIndex);
                    vecX.push_back(y->wcspIndex);
                }
            } else if (y->wcspIndex < x->wcspIndex && y->wcspIndex < z->wcspIndex) {
                if (x->wcspIndex < z->wcspIndex) {
                    vecX.push_back(y->wcspIndex);
                    vecY.push_back(z->wcspIndex);
                } else {
                    vecZ.push_back(y->wcspIndex);
                    vecY.push_back(x->wcspIndex);
                }
            } else if (z->wcspIndex < x->wcspIndex && z->wcspIndex < y->wcspIndex) {
                if (x->wcspIndex < y->wcspIndex) {
                    vecX.push_back(z->wcspIndex);
                    vecZ.push_back(y->wcspIndex);
                } else {
                    vecY.push_back(z->wcspIndex);
                    vecZ.push_back(x->wcspIndex);
                }
            }
            ToulBar2::Berge_Dec = true;
            break;
        default:
            break;
        }
    }
#endif

    if (ToulBar2::bilevel >= 2) {
        assert(ToulBar2::bilevel <= 3);
        deconnect(true);
        assert(xy->deconnected());
        assert(xz->deconnected());
        assert(yz->deconnected());
        assert(wcsp->delayedCtrBLP[ToulBar2::bilevel - 1].size() >= 3); // fresh new binary cost functions must have been automatically created and added
        assert(yz->wcspIndex == wcsp->delayedCtrBLP[ToulBar2::bilevel - 1].back());
        wcsp->delayedCtrBLP[ToulBar2::bilevel - 1].pop_back();
        assert(xz->wcspIndex == wcsp->delayedCtrBLP[ToulBar2::bilevel - 1].back());
        wcsp->delayedCtrBLP[ToulBar2::bilevel - 1].pop_back();
        assert(xy->wcspIndex == wcsp->delayedCtrBLP[ToulBar2::bilevel - 1].back());
        wcsp->delayedCtrBLP[ToulBar2::bilevel - 1].pop_back();
        wcsp->delayedCtrBLP[ToulBar2::bilevel - 1].push_back(wcspIndex);
    } else {
        propagate();
    }
}

TernaryConstraint::TernaryConstraint(WCSP* wcsp)
    : AbstractTernaryConstraint<EnumeratedVariable, EnumeratedVariable, EnumeratedVariable>(wcsp)
    , sizeX(0)
    , sizeY(0)
    , sizeZ(0)
    , functionalX(false)
    , functionalY(false)
    , functionalZ(false)
{
    //	unsigned int maxdom = wcsp->getMaxDomainSize();
    //    deltaCostsX = vector<StoreCost>(maxdom,StoreCost(MIN_COST,storeCost));
    //    deltaCostsY = vector<StoreCost>(maxdom,StoreCost(MIN_COST,storeCost));
    //    deltaCostsZ = vector<StoreCost>(maxdom,StoreCost(MIN_COST,storeCost));
    //    supportX = vector< pair<Value,Value> >(maxdom);
    //    supportY = vector< pair<Value,Value> >(maxdom);
    //    supportZ = vector< pair<Value,Value> >(maxdom);
    linkX = new DLink<ConstraintLink>;
    linkY = new DLink<ConstraintLink>;
    linkZ = new DLink<ConstraintLink>;

    //    costs = vector<StoreCost>(maxdom*maxdom*maxdom,StoreCost(MIN_COST,storeCost));
    //    for (unsigned int a = 0; a < maxdom; a++)
    //       for (unsigned int b = 0; b < maxdom; b++)
    //           for (unsigned int c = 0; c < maxdom; c++)
    //               costs[a * maxdom * maxdom + b * maxdom + c] = MIN_COST;
    xy = NULL;
    xz = NULL;
    yz = NULL;
}

double TernaryConstraint::computeTightness()
{
    int count = 0;
    double sum = 0;
    double tight = -1;
    Cost* costs = new Cost[x->getDomainSize() * y->getDomainSize() * z->getDomainSize()];
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
            for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                Cost c = getCost(*iterX, *iterY, *iterZ);
                sum += to_double(min(wcsp->getUb(), c));
                costs[count] = min(wcsp->getUb(), c);
                count++;
            }
        }
    }

    if (ToulBar2::weightedTightness == 2) {
        tight = to_double(stochastic_selection<Cost>(costs, 0, count - 1, count / 2));
    } else {
        tight = sum / (double)count;
    }
    delete[] costs;
    return tight;
}

void TernaryConstraint::print(ostream& os)
{
    os << this << " TernaryConstraint(" << x->getName() << ((functionalX) ? "!" : "") << "," << y->getName() << ((functionalY) ? "!" : "") << "," << z->getName() << ((functionalZ) ? "!" : "") << ")";
    if (ToulBar2::weightedDegree)
        os << "/" << getConflictWeight();
    if (wcsp->getTreeDec()) {
        cout << "   cluster: " << getCluster();
        if (ToulBar2::heuristicFreedom && getCluster() != -1) {
            cout << "  freedom: " << wcsp->getTreeDec()->getCluster(getCluster())->getFreedom() << " isinTD: " << wcsp->getTreeDec()->getCluster(getCluster())->getIsCurrInTD();
        }
    }
    cout << endl;

    if (ToulBar2::verbose >= 5) {
        for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
            for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                    os << " " << getCost(*iterX, *iterY, *iterZ);
                }
                os << " ; ";
            }
            os << endl;
        }
    }
}

void TernaryConstraint::dump(ostream& os, bool original)
{
    unsigned int tuples = 0;
    unsigned int tuplesNotTop = 0;
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
            for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                Cost cost = getCost(*iterX, *iterY, *iterZ);
                if (cost > MIN_COST) {
                    tuples++;
                }
                if (cost != top) {
                    tuplesNotTop++;
                }
            }
        }
    }
    Cost defaultCost = MIN_COST;
    if (tuplesNotTop < tuples) {
        defaultCost = top;
        tuples = tuplesNotTop;
    }
    os << "3 " << ((original) ? (x->wcspIndex) : x->getCurrentVarId()) << " " << ((original) ? (y->wcspIndex) : y->getCurrentVarId()) << " " << ((original) ? (z->wcspIndex) : z->getCurrentVarId()) << " " << defaultCost << " " << tuples << endl;
    int i = 0;
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX, i++) {
        int j = 0;
        for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY, j++) {
            int k = 0;
            for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ, k++) {
                Cost cost = getCost(*iterX, *iterY, *iterZ);
                if (cost != defaultCost) {
                    os << ((original) ? *iterX : i) << " " << ((original) ? *iterY : j) << " " << ((original) ? *iterZ : k) << " " << ((original) ? cost : min(wcsp->getUb(), cost)) << endl;
                }
            }
        }
    }
}

void TernaryConstraint::dump_CFN(ostream& os, bool original)
{
    unsigned int tuples = 0;
    unsigned int tuplesNotTop = 0;
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
            for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                Cost cost = getCost(*iterX, *iterY, *iterZ);
                if (cost > MIN_COST) {
                    tuples++;
                }
                if (cost != top) {
                    tuplesNotTop++;
                }
            }
        }
    }
    Cost defaultCost = MIN_COST;
    if (tuplesNotTop < tuples) {
        defaultCost = top;
        tuples = tuplesNotTop;
    }
    bool printed = false;
    os << "\"F_" << ((original) ? (x->wcspIndex) : x->getCurrentVarId()) << "_"
       << ((original) ? (y->wcspIndex) : y->getCurrentVarId()) << "_"
       << ((original) ? (z->wcspIndex) : z->getCurrentVarId()) << "\":{\"scope\":[\"";
    os << name2cfn(x->getName()) << "\",\""
       << name2cfn(y->getName()) << "\",\""
       << name2cfn(z->getName()) << "\"],";
    os << "\"defaultcost\":" << defaultCost << ",\n\"costs\":[\n";
    int i = 0;
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX, i++) {
        int j = 0;
        for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY, j++) {
            int k = 0;
            for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ, k++) {
                if (getCost(*iterX, *iterY, *iterZ) != defaultCost) {
                    if (printed)
                        os << ",\n";
                    os << ((original) ? x->toIndex(*iterX) : i) << "," << ((original) ? y->toIndex(*iterY) : j) << "," << ((original) ? z->toIndex(*iterZ) : k) << ","
                       << ((original) ? wcsp->DCost2Decimal(wcsp->Cost2RDCost(getCost(*iterX, *iterY, *iterZ))) : ((wcsp->getUb() > getCost(*iterX, *iterY, *iterZ)) ? wcsp->DCost2Decimal(wcsp->Cost2RDCost(getCost(*iterX, *iterY, *iterZ))) : "inf"));
                    printed = true;
                }
            }
        }
    }
    os << "]},\n";
}

/*
 * Propagation methods
 *
 */
bool TernaryConstraint::project(EnumeratedVariable* x, Value value, Cost cost, vector<StoreCost>& deltaCostsX)
{
    assert(ToulBar2::verbose < 4 || ((cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] project(C" << getVar(0)->getName() << "," << getVar(1)->getName() << "," << getVar(2)->getName() << ", (" << x->getName() << "," << value << "), " << cost << ")" << endl), true));

    // hard ternary constraint costs are not changed
    if (!CUT(cost + wcsp->getLb(), wcsp->getUb())) {
        TreeDecomposition* td = wcsp->getTreeDec();
        if (td)
            td->addDelta(cluster, x, value, cost);
        deltaCostsX[x->toIndex(value)] += cost; // Warning! Possible overflow???
    }

    Cost oldcost = x->getCost(value);
    x->project(value, cost);
#ifdef DEECOMPLETE
    int xindex = getIndex(x);
    getVar((xindex + 1) % 3)->queueDEE();
    getVar((xindex + 2) % 3)->queueDEE();
#endif
    return (x->getSupport() == value || SUPPORTTEST(oldcost, cost));
}

void TernaryConstraint::extend(EnumeratedVariable* x, Value value, Cost cost, vector<StoreCost>& deltaCostsX)
{
    assert(ToulBar2::verbose < 4 || ((cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] extend(C" << getVar(0)->getName() << "," << getVar(1)->getName() << "," << getVar(2)->getName() << ", (" << x->getName() << "," << value << "), " << cost << ")" << endl), true));
    TreeDecomposition* td = wcsp->getTreeDec();
    if (td)
        td->addDelta(cluster, x, value, -cost);
    deltaCostsX[x->toIndex(value)] -= cost; // Warning! Possible overflow???
    x->extend(value, cost);
}

pair<pair<Cost, Cost>, pair<Cost, Cost>> TernaryConstraint::getMaxCost(int varIndex, Value a, Value b)
{
    Cost maxcosta = MIN_COST;
    Cost diffcosta = MIN_COST;
    Cost maxcostb = MIN_COST;
    Cost diffcostb = MIN_COST;
    if (varIndex == 0) {
        Cost ucosta = x->getCost(a);
        Cost ucostb = x->getCost(b);
        for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
            Cost ucosty = y->getCost(*iterY);
            for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                Cost costa = getCost(a, *iterY, *iterZ);
                Cost costb = getCost(b, *iterY, *iterZ);
                if (costa > maxcosta)
                    maxcosta = costa;
                if (costb > maxcostb)
                    maxcostb = costb;
                Cost ucostz = z->getCost(*iterZ);
                if (!CUT(ucostb + getCostWithBinaries(b, *iterY, *iterZ) + ucosty + ucostz + wcsp->getLb(), wcsp->getUb())) {
                    if (costa - costb > diffcosta)
                        diffcosta = costa - costb;
                }
                if (!CUT(ucosta + getCostWithBinaries(a, *iterY, *iterZ) + ucosty + ucostz + wcsp->getLb(), wcsp->getUb())) {
                    if (costb - costa > diffcostb)
                        diffcostb = costb - costa;
                }
            }
        }
    } else if (varIndex == 1) {
        Cost ucosta = y->getCost(a);
        Cost ucostb = y->getCost(b);
        for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
            Cost ucostx = x->getCost(*iterX);
            for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                Cost costa = getCost(*iterX, a, *iterZ);
                Cost costb = getCost(*iterX, b, *iterZ);
                if (costa > maxcosta)
                    maxcosta = costa;
                if (costb > maxcostb)
                    maxcostb = costb;
                Cost ucostz = z->getCost(*iterZ);
                if (!CUT(ucostb + getCostWithBinaries(*iterX, b, *iterZ) + ucostx + ucostz + wcsp->getLb(), wcsp->getUb())) {
                    if (costa - costb > diffcosta)
                        diffcosta = costa - costb;
                }
                if (!CUT(ucosta + getCostWithBinaries(*iterX, a, *iterZ) + ucostx + ucostz + wcsp->getLb(), wcsp->getUb())) {
                    if (costb - costa > diffcostb)
                        diffcostb = costb - costa;
                }
            }
        }
    } else {
        assert(varIndex == 2);
        Cost ucosta = z->getCost(a);
        Cost ucostb = z->getCost(b);
        for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
            Cost ucostx = x->getCost(*iterX);
            for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                Cost costa = getCost(*iterX, *iterY, a);
                Cost costb = getCost(*iterX, *iterY, b);
                if (costa > maxcosta)
                    maxcosta = costa;
                if (costb > maxcostb)
                    maxcostb = costb;
                Cost ucosty = y->getCost(*iterY);
                if (!CUT(ucostb + getCostWithBinaries(*iterX, *iterY, b) + ucostx + ucosty + wcsp->getLb(), wcsp->getUb())) {
                    if (costa - costb > diffcosta)
                        diffcosta = costa - costb;
                }
                if (!CUT(ucosta + getCostWithBinaries(*iterX, *iterY, a) + ucostx + ucosty + wcsp->getLb(), wcsp->getUb())) {
                    if (costb - costa > diffcostb)
                        diffcostb = costb - costa;
                }
            }
        }
    }
    assert(maxcosta >= diffcosta);
    assert(maxcostb >= diffcostb);
    return make_pair(make_pair(maxcosta, diffcosta), make_pair(maxcostb, diffcostb));
}

bool TernaryConstraint::separability(EnumeratedVariable* vy, EnumeratedVariable* vz)
{
    Cost c1, c;
    Tuple tch(3, 0);
    bool neweq = true; // true if we have  not a difference value
    bool sep = true; // false if vy and vz are not separable
    Cost diff = 0;
    first(vy, vz);
    EnumeratedVariable::iterator itvyfirst = yvar->begin();
    if (ToulBar2::verbose >= 1)
        cout << " [ " << zvar->getName() << " " << xvar->getName() << " " << yvar->getName() << " ] ?"; // << endl;
    while (sep && itvyfirst != yvar->end()) {
        itvx = xvar->begin();
        itvy = yvar->begin();
        itvz = zvar->begin();
        while (sep && itvx != xvar->end()) {
            unsigned int ix = xvar->toIndex(*itvx);
            tch[0] = ix;
            while (sep && itvy != yvar->end()) {
                unsigned int iy = yvar->toIndex(*itvy);
                tch[1] = iy;
                while (sep && itvy != itvyfirst && itvz != zvar->end()) {
                    unsigned int iz = zvar->toIndex(*itvz);
                    tch[2] = iz;

                    c1 = getCost(xvar, yvar, zvar, *itvx, *(itvyfirst), *itvz);
                    c = getCost(xvar, yvar, zvar, *itvx, *itvy, *itvz);

                    if (!universe(c, c1, wcsp->getUb())) {
                        if (ToulBar2::verbose >= 3) {
                            if (!neweq)
                                cout << "= \n";
                            else
                                cout << endl;
                        }
                        if (ToulBar2::verbose >= 3)
                            cout << " C" << tch[2] << "." << tch[0] << "." << tch[1] << " -  C" << tch[2] << "." << tch[0] << "." << yvar->toIndex(*itvyfirst) << " = " << c << " - " << c1;

                        if (neweq) {
                            diff = squareminus(c, c1, wcsp->getUb());
                            neweq = false;
                        } else
                            sep = (diff == squareminus(c, c1, wcsp->getUb()));
                        if (ToulBar2::verbose >= 3)
                            cout << " = " << squareminus(c, c1, wcsp->getUb()) << endl;
                    } else {
                        if (ToulBar2::verbose >= 3)
                            cout << "universe\n";
                    }
                    ++itvz;
                }
                ++itvy;
                itvz = zvar->begin();
                neweq = true;
            }
            ++itvx;
            itvy = yvar->begin();
            neweq = true;
        }
        ++itvyfirst;
        itvx = xvar->begin();
        neweq = true;
        if (ToulBar2::verbose >= 3)
            cout << "---\n";
    }
    return sep;
}

void TernaryConstraint::separate(EnumeratedVariable* vy, EnumeratedVariable* vz)
{
    Cost cost, minCost = MAX_COST;
    // assert(separability(vy,vz));
    first(vy, vz);
    vector<Cost> costsZX((size_t)zvar->getDomainInitSize() * (size_t)xvar->getDomainInitSize(), MIN_COST);
    vector<Cost> costsXY((size_t)xvar->getDomainInitSize() * (size_t)yvar->getDomainInitSize(), MIN_COST);
    string xv(xvar->getName()), yv(yvar->getName()), zv(zvar->getName());
    if (ToulBar2::verbose == 1)
        cout << "\n";

    if (ToulBar2::verbose >= 3)
        cout << "[ " << zvar->getName() << " " << xvar->getName() << " ]" << endl;
    while (itvz != zvar->end()) {
        while (itvx != xvar->end()) {
            minCost = MAX_COST;
            while (itvy != yvar->end()) {
                cost = getCost(xvar, yvar, zvar, *itvx, *itvy, *itvz);
                if (cost < minCost)
                    minCost = cost;
                if (minCost >= wcsp->getUb())
                    minCost = wcsp->getUb();
                ++itvy;
            }
            if (ToulBar2::verbose >= 3)
                cout << *itvx << " " << *itvz << " : " << minCost << endl;
            costsZX[zvar->toIndex(*itvz) * xvar->getDomainInitSize() + xvar->toIndex(*itvx)] = minCost;
            ++itvx;
            itvy = yvar->begin();
        }
        ++itvz;
        itvx = xvar->begin();
    }
    BinaryConstraint* existZX = xvar->getConstr(zvar);
    assert(existZX);
    BinaryConstraint* zx = new BinaryConstraint(wcsp, zvar, xvar, costsZX);
    if (ToulBar2::verbose >= 3)
        cout << "-------------\n";
    if (ToulBar2::verbose >= 3)
        cout << "[ " << xvar->getName() << " " << yvar->getName() << " ]" << endl;

    first(vy, vz);
    Cost costzx;
    while (itvx != xvar->end()) {
        while (itvy != yvar->end()) {
            itvz = zvar->begin();
            do {
                cost = getCost(xvar, yvar, zvar, *itvx, *itvy, *itvz);
                costzx = zx->getCost(*itvz, *itvx);
                ++itvz;
            } while (itvz != zvar->end() && cost >= wcsp->getUb() && costzx >= wcsp->getUb());
            costsXY[xvar->toIndex(*itvx) * yvar->getDomainInitSize() + yvar->toIndex(*itvy)] = squareminus(cost, costzx, wcsp->getUb());
            if (ToulBar2::verbose >= 3)
                cout << *itvx << " " << *itvy << " : " << squareminus(cost, costzx, wcsp->getUb()) << endl;
            assert(squareminus(cost, costzx, wcsp->getUb()) >= MIN_COST);
            ++itvy;
        }
        ++itvx;
        itvy = yvar->begin();
    }
    BinaryConstraint* existXY = xvar->getConstr(yvar);
    assert(existXY);
    BinaryConstraint* xy = new BinaryConstraint(wcsp, xvar, yvar, costsXY);

    assert(verifySeparate(zx, xy));

    // fusion with the existing constraint (xz)
    if (!zx->universal()) {
        if (ToulBar2::verbose >= 1)
            cout << "[ " << zv << " " << xv << " ]" << endl;
        existZX->addCosts(zx);
        existZX->reconnect();
        existZX->propagate();
    }
    zx->deconnect(); //  unsafe to delete zx due to x and z lists of cost functions

    // fusion with the existing constraint (xy)
    if (!xy->universal()) {
        if (ToulBar2::verbose >= 1)
            cout << "[ " << xv << " " << yv << " ]" << endl;
        existXY->addCosts(xy);
        existXY->reconnect();
        existXY->propagate();
    }
    xy->deconnect(); //  unsafe to delete xy due to x and y lists of cost functions
    deconnect();
}

void TernaryConstraint::fillxy()
{
    TreeDecomposition* td = wcsp->getTreeDec();
    BinaryConstraint* xy_ = NULL;
    xy_ = x->getConstr(y);
    if (td && xy_ && !td->isSameCluster(getCluster(), xy_->getCluster())) {
        BinaryConstraint* xy__ = x->getConstr(y, getCluster());
        if (xy__)
            xy_ = xy__; // we have found another constraint of the same cluster
    }
    if (!xy_ || (xy_ && td && !td->isSameCluster(getCluster(), xy_->getCluster()))) {
        xy = wcsp->newBinaryConstr(x, y, this);
        if (td && !ToulBar2::approximateCountingBTD)
            xy->setCluster(getCluster());
        if (td && xy_ && !td->isSameCluster(getCluster(), xy_->getCluster()))
            xy->setDuplicate();
        wcsp->elimBinOrderInc();
    } else
        xy = xy_;
    if (xy->isDuplicate())
        setDuplicate();
}

void TernaryConstraint::fillxz()
{
    TreeDecomposition* td = wcsp->getTreeDec();
    BinaryConstraint* xz_ = NULL;
    xz_ = x->getConstr(z);
    if (td && xz_ && !td->isSameCluster(getCluster(), xz_->getCluster())) {
        BinaryConstraint* xz__ = x->getConstr(z, getCluster());
        if (xz__)
            xz_ = xz__; // we have found another constraint of the same cluster
    }
    if (!xz_ || (xz_ && td && !td->isSameCluster(getCluster(), xz_->getCluster()))) {
        xz = wcsp->newBinaryConstr(x, z, this);
        if (td && !ToulBar2::approximateCountingBTD)
            xz->setCluster(getCluster());
        if (td && xz_ && !td->isSameCluster(getCluster(), xz_->getCluster()))
            xz->setDuplicate();
        wcsp->elimBinOrderInc();
    } else
        xz = xz_;
    if (xz->isDuplicate())
        setDuplicate();
}

void TernaryConstraint::fillyz()
{
    TreeDecomposition* td = wcsp->getTreeDec();
    BinaryConstraint* yz_ = NULL;
    yz_ = y->getConstr(z);
    if (td && yz_ && !td->isSameCluster(getCluster(), yz_->getCluster())) {
        BinaryConstraint* yz__ = y->getConstr(z, getCluster());
        if (yz__)
            yz_ = yz__;
    }
    if (!yz_ || (yz_ && td && !td->isSameCluster(getCluster(), yz_->getCluster()))) {
        yz = wcsp->newBinaryConstr(y, z, this);
        if (td && !ToulBar2::approximateCountingBTD)
            yz->setCluster(getCluster());
        if (td && yz_ && !td->isSameCluster(getCluster(), yz_->getCluster()))
            yz->setDuplicate();
        wcsp->elimBinOrderInc();
    } else
        yz = yz_;
    if (yz->isDuplicate())
        setDuplicate();
}

void TernaryConstraint::fillElimConstrBinaries()
{
    fillxy();
    fillxz();
    fillyz();

    resetConflictWeight(); // if needed recompute tightness after having updated costs

    if (ToulBar2::verbose > 1)
        cout << "[" << Store::getDepth() << ",W" << wcsp->getIndex() << "] fillElimConstrBinaries (" << x->wcspIndex << "," << y->wcspIndex << "," << z->wcspIndex << ")  ";
}

void TernaryConstraint::setDuplicates()
{
    assert(wcsp->getTreeDec());
    if (xy->getCluster() != cluster) {
        BinaryConstraint* xy_ = x->getConstr(y, getCluster());
        if (xy_) {
            if (xy_->isDuplicate())
                setDuplicate();
            xy = xy_;
        } else {
            wcsp->initElimConstr();
            xy = wcsp->newBinaryConstr(x, y);
            xy->setCluster(getCluster());
            xy->setDuplicate();
            wcsp->elimBinOrderInc();
            setDuplicate();
        }
    }
    if (xz->getCluster() != cluster) {
        BinaryConstraint* xz_ = x->getConstr(z, getCluster());
        if (xz_) {
            xz = xz_;
            if (xz_->isDuplicate())
                setDuplicate();
        } else {
            wcsp->initElimConstr();
            xz = wcsp->newBinaryConstr(x, z);
            xz->setCluster(getCluster());
            xz->setDuplicate();
            wcsp->elimBinOrderInc();
            setDuplicate();
        }
    }
    if (yz->getCluster() != cluster) {
        BinaryConstraint* yz_ = y->getConstr(z, getCluster());
        if (yz_) {
            yz = yz_;
            if (yz_->isDuplicate())
                setDuplicate();
        } else {
            wcsp->initElimConstr();
            yz = wcsp->newBinaryConstr(y, z);
            yz->setCluster(getCluster());
            yz->setDuplicate();
            wcsp->elimBinOrderInc();
            setDuplicate();
        }
    }
    assert(xy->getCluster() == getCluster() && xz->getCluster() == getCluster() && yz->getCluster() == getCluster());
}

bool TernaryConstraint::verify(EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z)
{
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        Cost minCost = MAX_COST;
        for (EnumeratedVariable::iterator iterY = y->begin(); minCost > MIN_COST && iterY != y->end(); ++iterY) {
            for (EnumeratedVariable::iterator iterZ = z->begin(); minCost > MIN_COST && iterZ != z->end(); ++iterZ) {
                Cost cost = getCost(x, y, z, *iterX, *iterY, *iterZ);
#ifndef NO_STORE_TERNARY_COSTS
                if (ToulBar2::LcLevel >= LC_DAC && getIndex(x) == getDACScopeIndex())
                    cost += y->getCost(*iterY) + z->getCost(*iterZ);
#endif
                GLB(&minCost, cost);
            }
        }
        if (minCost > MIN_COST) {
            cout << "not FDAC: variable " << x->getName() << " value " << *iterX << " of " << *this;
            return false;
        }
    }
    return true;
}

bool TernaryConstraint::verify()
{
    TreeDecomposition* td = wcsp->getTreeDec();

    if (td) {
        if (!td->isSameCluster(cluster, xy->getCluster()) || !td->isSameCluster(cluster, xz->getCluster()) || !td->isSameCluster(cluster, yz->getCluster())) {
            if (ToulBar2::heuristicFreedom) {
                cout << " different cluster assignment for ternary: " << cluster << "(" << td->getCluster(cluster)->getFreedom() << ") xy: " << xy->getCluster() << "(" << td->getCluster(xy->getCluster())->getFreedom() << ") xz: " << xz->getCluster() << "(" << td->getCluster(xz->getCluster())->getFreedom() << ") yz: " << yz->getCluster() << "(" << td->getCluster(yz->getCluster())->getFreedom() << ")" << endl;
            } else {
                cout << " different cluster assignment for ternary: " << cluster << " xy: " << xy->getCluster() << " xz: " << xz->getCluster() << " yz: " << yz->getCluster() << endl;
            }
            cout << *this;
            cout << *xy;
            cout << *xz;
            cout << *yz;
            return false;
        }
    }

    if (ToulBar2::LcLevel == LC_DAC) {
        switch (getDACScopeIndex()) {
        case 0:
            return verifyX();
            break;
        case 1:
            return verifyY();
            break;
        case 2:
            return verifyZ();
            break;
        default:
            return false;
        }
    } else {
        return verifyX() && verifyY() && verifyZ();
    }
}

bool TernaryConstraint::checkTreeDecomposition()
{
    return (!wcsp->getTreeDec() || (wcsp->getTreeDec()->isSameCluster(cluster, xy->getCluster()) && wcsp->getTreeDec()->isSameCluster(cluster, xz->getCluster()) && wcsp->getTreeDec()->isSameCluster(cluster, yz->getCluster())));
}

// Triangle::Triangle(WCSP *wcsp,
//				  EnumeratedVariable *xx,
//				  EnumeratedVariable *yy,
//				  EnumeratedVariable *zz,
//				  BinaryConstraint* _xy,
//				  BinaryConstraint* _xz,
//				  BinaryConstraint* _yz,
//				  StoreStack<Cost, Cost> *storeCost)
//	: AbstractTernaryConstraint<EnumeratedVariable,EnumeratedVariable,EnumeratedVariable>(wcsp, xx, yy, zz),
//	  xy(_xy), xz(_xz), yz(_yz), xyz(xx->getConstr(yy,zz))
//						  {
//	// if xyz == NULL then create "empty" ternaryconstr like in NaryConstr.cpp
//	// 		BinaryConstraint* bctr;
////	TernaryConstraint* tctr = new TernaryConstraint(this, &storeData->storeCost);
////	elimTernConstrs.push_back(tctr);
////	for (int j = 0; j < 3; j++) {
////		if (!ToulBar2::vac) bctr = new BinaryConstraint(this, &storeData->storeCost);
////		else bctr = new VACBinaryConstraint(this, &storeData->storeCost);
////		elimBinConstrs.push_back(bctr);
////	}
//
//						  }

// activate {
//	xyz = wcsp->newTernaryConstr(x,y,z,this);
// }

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
