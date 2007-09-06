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
        AbstractBinaryConstraint<EnumeratedVariable,EnumeratedVariable>(wcsp, xx, yy), sizeX(xx->getDomainInitSize()), sizeY(yy->getDomainInitSize())
{
    deltaCostsX = vector<StoreCost>(sizeX,StoreCost(0,storeCost));
    deltaCostsY = vector<StoreCost>(sizeY,StoreCost(0,storeCost));
    assert(tab.size() == sizeX * sizeY);
    supportX = vector<Value>(sizeX,y->getInf());
    supportY = vector<Value>(sizeY,x->getInf());

	maxCost = 0;
	costs = vector<StoreCost>(sizeX*sizeY,StoreCost(0,storeCost));
	
    for (unsigned int a = 0; a < x->getDomainInitSize(); a++) 
         for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
				Cost c = tab[a * sizeY + b];
                costs[a * sizeY + b] = c;
                if(maxCost < c) maxCost = c;
         }


    propagate();
}

BinaryConstraint::BinaryConstraint(WCSP *wcsp, StoreStack<Cost, Cost> *storeCost)
 	: AbstractBinaryConstraint<EnumeratedVariable,EnumeratedVariable>(wcsp), sizeX(wcsp->maxdomainsize), sizeY(wcsp->maxdomainsize)
{
	unsigned int maxdomainsize = wcsp->maxdomainsize;
    deltaCostsX = vector<StoreCost>(maxdomainsize,StoreCost(0,storeCost));
    deltaCostsY = vector<StoreCost>(maxdomainsize,StoreCost(0,storeCost));
    supportX = vector<Value>(maxdomainsize,0);
    supportY = vector<Value>(maxdomainsize,0);
    linkX = new DLink<ConstraintLink>;
    linkY = new DLink<ConstraintLink>;    

    costs = vector<StoreCost>(maxdomainsize*maxdomainsize,StoreCost(0,storeCost));         	
    for (unsigned int a = 0; a < maxdomainsize; a++) 
         for (unsigned int b = 0; b < maxdomainsize; b++) 
                costs[a * maxdomainsize + b] = 0;
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

void BinaryConstraint::dump(ostream& os)
{
    os << "2 " << x->wcspIndex << " " << y->wcspIndex << " 0 " << x->getDomainSize() * y->getDomainSize() << endl;
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
            os << *iterX << " " << *iterY << " " << getCost(*iterX, *iterY) << endl;
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
	assert(connected());
//    wcsp->revise(this);
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
                if (NOCUT(minCost + wcsp->getLb(),wcsp->getUb())) deltaCostsX[xindex] += minCost;  // Warning! Possible overflow???
                if (x->getSupport() == *iterX) supportBroken = true;
                if (ToulBar2::verbose >= 2) cout << "binary projection of " << minCost << " from C" << x->getName() << "," << y->getName() << "(" << *iterX << "," << minCostValue << ")" << endl;
                x->project(*iterX, minCost);
                if (deconnected()) return;
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
	assert(connected());
//    wcsp->revise(this);	
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
                        int yindex = y->toIndex(*iterY);
                        deltaCostsY[yindex] -= (minCost - cost);  // Warning! Possible overflow???
                        y->extend(*iterY, minCost - cost);
                        supportY[yindex] = *iterX;
                        if (ToulBar2::vacAlternative) {
                            x->queueVAC2();
                            y->queueVAC2();
                        }
                     }
                }
                
                // hard binary constraint costs are not changed
                if (NOCUT(minCost + wcsp->getLb(), wcsp->getUb())) deltaCostsX[xindex] += minCost;  // Warning! Possible overflow???
                if (x->getSupport() == *iterX) supportBroken = true;
                if (ToulBar2::verbose >= 2) cout << "binary projection of " << minCost << " from C" << x->getName() << "," << y->getName() << "(" << *iterX << "," << minCostValue << ")" << endl;
                x->project(*iterX, minCost);
                if (deconnected()) return;
            }
            supportX[xindex] = minCostValue;
        }
    }
    if (supportBroken) {
        x->findSupport();
    }
}

template <GetCostMember getBinaryCost>
void BinaryConstraint::projection(EnumeratedVariable *x, EnumeratedVariable *y, Value valueY)
{
    bool supportBroken = false;
//    wcsp->revise(this);
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        Cost cost = GETCOST(*iterX, valueY);
        if (cost > 0) {
            if (x->getSupport() == *iterX) supportBroken = true;
            if (ToulBar2::verbose >= 2) cout << "binary projection of " << cost << " on C(" << x->getName() << "," << *iterX << ")" << endl;
            addcost(x, y, *iterX, valueY, -cost);
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
        if (getDACScopeIndex() == getIndex(x)) minCost += y->getCost(y->getInf());
        EnumeratedVariable::iterator iterY = y->begin();
        for (++iterY; minCost > 0 && iterY != y->end(); ++iterY) {
            Cost cost = GETCOST(*iterX, *iterY);
            if (getDACScopeIndex() == getIndex(x)) cost += y->getCost(*iterY);
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
