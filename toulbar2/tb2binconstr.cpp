/*
 * ****** Binary constraint applied on variables with enumerated domains ******
 */

#include "tb2binconstr.hpp"
#include "tb2wcsp.hpp"

// coding shorthand
#define GETCOST (this->*getBinaryCost)

/*
 * Constructors and misc.
 * 
 */

BinaryConstraint::BinaryConstraint(WCSP *wcsp, EnumeratedVariable *xx, EnumeratedVariable *yy, vector<Cost> &tab, StoreStack<Cost, Cost> *storeCost) : 
        AbstractBinaryConstraint<EnumeratedVariable,EnumeratedVariable>(wcsp, xx, yy), sizeX(xx->getDomainInitSize()), sizeY(yy->getDomainInitSize()), costs(tab)
{
    deltaCostsX = vector<StoreCost>(sizeX,StoreCost(0,storeCost));
    deltaCostsY = vector<StoreCost>(sizeY,StoreCost(0,storeCost));
    assert(tab.size() == sizeX * sizeY);
    supportX = vector<Value>(sizeX,y->getInf());
    supportY = vector<Value>(sizeY,x->getInf());
    xx->queueAC();
    xx->queueDAC();
    yy->queueAC();
    yy->queueDAC();
}

void BinaryConstraint::print(ostream& os)
{
    os << this << " BinaryConstraint(" << x->getName() << "," << y->getName() << ")" << endl;
    if (ToulBar2::verbose >= 5) {
        for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
            for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                os << " " << getCost(*iterX, *iterY);
            }
            os << endl;
        }
    }
}


/*
 * Propagation methods
 * 
 */

template <GetCostMember getBinaryCost>
void BinaryConstraint::findSupport(EnumeratedVariable *x, EnumeratedVariable *y,
        vector<Value> &supportX, vector<StoreCost> &deltaCostsX)
{
    if (ToulBar2::verbose >= 3) cout << "findSupport C" << x->getName() << "," << y->getName() << endl;
    bool supportBroken = false;
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        int xindex = x->toIndex(*iterX);
        Value support = supportX[xindex];
        if (y->cannotbe(support) || GETCOST(*iterX, support) > 0) {
            Value minCostValue = y->getInf();
            Cost minCost = GETCOST(*iterX, minCostValue);
            EnumeratedVariable::iterator iterY = y->begin();
            for (++iterY; minCost > 0 && iterY != y->end(); ++iterY) {
                Cost cost = GETCOST(*iterX, *iterY);
                if (cost < minCost) {
                    minCost = cost;
                    minCostValue = *iterY;
                }
            }
            if (minCost > 0) {
                // hard binary constraint costs are not changed
                if (minCost + wcsp->getLb() < wcsp->getUb()) deltaCostsX[xindex] += minCost;  // Warning! Possible overflow???
                if (x->getSupport() == *iterX) supportBroken = true;
                if (ToulBar2::verbose >= 2) cout << "binary projection of " << minCost << " from C" << x->getName() << "," << y->getName() << "(" << *iterX << "," << minCostValue << ")" << endl;
                x->project(*iterX, minCost);
            }
            supportX[xindex] = minCostValue;
        }
    }
    if (supportBroken) {
        x->findSupport();
    }
}

template <GetCostMember getBinaryCost>
void BinaryConstraint::findFullSupport(EnumeratedVariable *x, EnumeratedVariable *y, 
        vector<Value> &supportX, vector<StoreCost> &deltaCostsX, 
        vector<Value> &supportY, vector<StoreCost> &deltaCostsY)
{
    if (ToulBar2::verbose >= 3) cout << "findFullSupport C" << x->getName() << "," << y->getName() << endl;
    bool supportBroken = false;
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        int xindex = x->toIndex(*iterX);
        Value support = supportX[xindex];
        if (y->cannotbe(support) || GETCOST(*iterX, support) + y->getCost(support) > 0) {
            Value minCostValue = y->getInf();
            Cost minCost = GETCOST(*iterX, minCostValue) + y->getCost(minCostValue);
            EnumeratedVariable::iterator iterY = y->begin();
            for (++iterY; minCost > 0 && iterY != y->end(); ++iterY) {
                Cost cost = GETCOST(*iterX, *iterY) + y->getCost(*iterY);
                if (cost < minCost) {
                    minCost = cost;
                    minCostValue = *iterY;
                }
            }
            if (minCost > 0) {
                // extend unary to binary
                for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                    Cost cost = GETCOST(*iterX, *iterY);
                    if (minCost > cost) {
                        deltaCostsY[*iterY] -= (minCost - cost);  // Warning! Possible overflow???
                        y->extend(*iterY, minCost - cost);
                        int yindex = y->toIndex(*iterY);
                        supportY[yindex] = *iterX;
                     }
                }
                
                // hard binary constraint costs are not changed
                if (minCost + wcsp->getLb() < wcsp->getUb()) deltaCostsX[xindex] += minCost;  // Warning! Possible overflow???
                if (x->getSupport() == *iterX) supportBroken = true;
                if (ToulBar2::verbose >= 2) cout << "binary projection of " << minCost << " from C" << x->getName() << "," << y->getName() << "(" << *iterX << "," << minCostValue << ")" << endl;
                x->project(*iterX, minCost);
            }
            supportX[xindex] = minCostValue;
        }
    }
    if (supportBroken) {
        x->findSupport();
    }
}

template <GetCostMember getBinaryCost>
void BinaryConstraint::projection(EnumeratedVariable *x, Value valueY)
{
    bool supportBroken = false;
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        Cost cost = GETCOST(*iterX, valueY);
        if (cost > 0) {
            if (x->getSupport() == *iterX) supportBroken = true;
            if (ToulBar2::verbose >= 2) cout << "binary projection of " << cost << " on C(" << x->getName() << "," << *iterX << ")" << endl;
            x->project(*iterX, cost);
        }
    }
    if (supportBroken) {
        x->findSupport();
    }
}

template <GetCostMember getBinaryCost>
bool BinaryConstraint::verify(EnumeratedVariable *x, EnumeratedVariable *y)
{
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        Cost minCost = GETCOST(*iterX, y->getInf());
        if (x->wcspIndex < y->wcspIndex) minCost += y->getCost(y->getInf());
        EnumeratedVariable::iterator iterY = y->begin();
        for (++iterY; minCost > 0 && iterY != y->end(); ++iterY) {
            Cost cost = GETCOST(*iterX, *iterY);
            if (x->wcspIndex < y->wcspIndex) cost += y->getCost(*iterY);
            if (cost < minCost) {
                minCost = cost;
            }
        }
        if (minCost > 0) {
            cout << *this;
            return false;
        }
    }
    return true;
}
