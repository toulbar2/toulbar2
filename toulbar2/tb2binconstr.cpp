/*
 * ****** Binary constraint applied on variables with enumerated domains ******
 */

#include "tb2binconstr.hpp"


/*
 * Constructors and misc.
 * 
 */

BinaryConstraint::BinaryConstraint(CostVariable *xx, CostVariable *yy, vector<Cost> &tab, StoreStack<Cost,Cost> *store) : 
        x(xx), y(yy), sizeX(xx->getDomainInitSize()), sizeY(yy->getDomainInitSize()), costs(tab),
        linkX(NULL), linkY(NULL)
{
    assert(xx->getEnumerated() && yy->getEnumerated());
    deltaCostsX = vector<StoreCost>(sizeX,StoreCost(0,store));
    deltaCostsY = vector<StoreCost>(sizeY,StoreCost(0,store));
    linkX = xx->addConstraint(this,0);
    linkY = yy->addConstraint(this,1);
    assert(tab.size() == sizeX * sizeY);
    supportX = vector<Value>(sizeX,y->getInf());
    supportY = vector<Value>(sizeY,x->getInf());
}

void BinaryConstraint::print(ostream& os)
{
    os << this << " BinaryConstraint(" << x->getName() << "," << y->getName() << ")" << endl;
    if (ToulBar2::verbose >= 2) {
        for (Variable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
            for (Variable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
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

void findSupport(BinaryConstraint *c, CostVariable *x, CostVariable *y, 
        vector<Value> &supportX, vector<StoreCost> &deltaCostsX, GetCostFunc getCost)
{
    bool supportBroken = false;
    for (Variable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        int xindex = x->toIndex(*iterX);
        Value support = supportX[xindex];
        if (y->cannotbe(support) || getCost(c, *iterX, support) > 0) {
            Value minCostValue = y->getInf();
            Cost minCost = getCost(c, *iterX, minCostValue);
            Variable::iterator iterY = y->begin();
            for (++iterY; minCost > 0 && iterY != y->end(); ++iterY) {
                Cost cost = getCost(c, *iterX, *iterY);
                if (cost < minCost) {
                    minCost = cost;
                    minCostValue = *iterY;
                }
            }
            if (minCost > 0) {
                // hard binary constraint costs are not changed
                if (minCost <= c->wcsp->getUb()) deltaCostsX[xindex] += minCost;  // Warning! Possible overflow???
                if (x->getSupport() == *iterX) supportBroken = true;
                if (ToulBar2::verbose >= 2) cout << "binary projection of " << minCost << " from C(" << x->getName() << "," << y->getName() << "," << *iterX << "," << minCostValue << ")" << endl;
                x->project(*iterX, minCost);
            }
            supportX[xindex] = minCostValue;
        }
    }
    if (supportBroken) {
        x->findSupport();
    }
}

void projection(BinaryConstraint *c, CostVariable *x, Value valueY, GetCostFunc getCost)
{
    bool supportBroken = false;
    for (Variable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        Cost cost = getCost(c, *iterX, valueY);
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

bool verify(BinaryConstraint *c, CostVariable *x, CostVariable *y, GetCostFunc getCost)
{
    for (Variable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        Cost minCost = getCost(c, *iterX, y->getInf());
        Variable::iterator iterY = y->begin();
        for (++iterY; minCost > 0 && iterY != y->end(); ++iterY) {
            Cost cost = getCost(c, *iterX, *iterY);
            if (cost < minCost) {
                minCost = cost;
            }
        }
        if (minCost > 0) return false;
    }
    return true;
}

