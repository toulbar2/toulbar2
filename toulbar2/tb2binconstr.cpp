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
    assert(xx != yy);
    assert(xx->getEnumerated() && yy->getEnumerated());
    deltaCostsX = vector<StoreCost>(sizeX,StoreCost(0,store));
    deltaCostsY = vector<StoreCost>(sizeY,StoreCost(0,store));
    linkX = xx->postConstraint(this,0);
    linkY = yy->postConstraint(this,1);
    assert(tab.size() == sizeX * sizeY);
    supportX = vector<Value>(sizeX,y->getInf());
    supportY = vector<Value>(sizeY,x->getInf());
}

int BinaryConstraint::getSmallestVarIndexInScope(int forbiddenScopeIndex)
{
    return (forbiddenScopeIndex)?x->wcspIndex:y->wcspIndex;
}

void BinaryConstraint::print(ostream& os)
{
    os << this << " BinaryConstraint(" << x->getName() << "," << y->getName() << ")" << endl;
    if (ToulBar2::verbose >= 5) {
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

//void findSupport(BinaryConstraint *c, CostVariable *x, CostVariable *y, 
//        vector<Value> &supportX, vector<StoreCost> &deltaCostsX, GetCostFunc getCost)
//{
//    if (ToulBar2::verbose >= 3) cout << "findSupport C" << x->getName() << "," << y->getName() << endl;
//    bool supportBroken = false;
//    for (Variable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
//        int xindex = x->toIndex(*iterX);
//        Value support = supportX[xindex];
//        if (y->cannotbe(support) || getCost(c, *iterX, support) > 0) {
//            Value minCostValue = y->getInf();
//            Cost minCost = getCost(c, *iterX, minCostValue);
//            Variable::iterator iterY = y->begin();
//            for (++iterY; minCost > 0 && iterY != y->end(); ++iterY) {
//                Cost cost = getCost(c, *iterX, *iterY);
//                if (cost < minCost) {
//                    minCost = cost;
//                    minCostValue = *iterY;
//                }
//            }
//            if (minCost > 0) {
//                // hard binary constraint costs are not changed
//                if (minCost + c->wcsp->getLb() <= c->wcsp->getUb()) deltaCostsX[xindex] += minCost;  // Warning! Possible overflow???
//                if (x->getSupport() == *iterX) supportBroken = true;
//                if (ToulBar2::verbose >= 2) cout << "binary projection of " << minCost << " from C" << x->getName() << "," << y->getName() << "(" << *iterX << "," << minCostValue << ")" << endl;
//                x->project(*iterX, minCost);
//            }
//            supportX[xindex] = minCostValue;
//        }
//    }
//    if (supportBroken) {
//        x->findSupport();
//    }
//}

void BinaryConstraint::findSupportX()
{
    if (ToulBar2::verbose >= 3) cout << "findSupport C" << x->getName() << "," << y->getName() << endl;
    bool supportBroken = false;
    for (Variable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        int xindex = x->toIndex(*iterX);
        Value support = supportX[xindex];
        if (y->cannotbe(support) || getCost(*iterX, support) > 0) {
            Value minCostValue = y->getInf();
            Cost minCost = getCost(*iterX, minCostValue);
            Variable::iterator iterY = y->begin();
            for (++iterY; minCost > 0 && iterY != y->end(); ++iterY) {
                Cost cost = getCost(*iterX, *iterY);
                if (cost < minCost) {
                    minCost = cost;
                    minCostValue = *iterY;
                }
            }
            if (minCost > 0) {
                // hard binary constraint costs are not changed
                if (minCost + wcsp->getLb() <= wcsp->getUb()) deltaCostsX[xindex] += minCost;  // Warning! Possible overflow???
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

void BinaryConstraint::findSupportY()
{
    if (ToulBar2::verbose >= 3) cout << "findSupport C" << y->getName() << "," << x->getName() << endl;
    bool supportBroken = false;
    for (Variable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
        int yindex = y->toIndex(*iterY);
        Value support = supportY[yindex];
        if (x->cannotbe(support) || getCost(support,*iterY) > 0) {
            Value minCostValue = x->getInf();
            Cost minCost = getCost(minCostValue,*iterY);
            Variable::iterator iterX = x->begin();
            for (++iterX; minCost > 0 && iterX != x->end(); ++iterX) {
                Cost cost = getCost(*iterX, *iterY);
                if (cost < minCost) {
                    minCost = cost;
                    minCostValue = *iterX;
                }
            }
            if (minCost > 0) {
                // hard binary constraint costs are not changed
                if (minCost + wcsp->getLb() <= wcsp->getUb()) deltaCostsY[yindex] += minCost;  // Warning! Possible overflow???
                if (y->getSupport() == *iterY) supportBroken = true;
                if (ToulBar2::verbose >= 2) cout << "binary projection of " << minCost << " from C" << y->getName() << "," << x->getName() << "(" << *iterY << "," << minCostValue << ")" << endl;
                y->project(*iterY, minCost);
            }
            supportY[yindex] = minCostValue;
        }
    }
    if (supportBroken) {
        y->findSupport();
    }
}

//void findFullSupport(BinaryConstraint *c, CostVariable *x, CostVariable *y, 
//        vector<Value> &supportX, vector<StoreCost> &deltaCostsX, 
//        vector<Value> &supportY, vector<StoreCost> &deltaCostsY, GetCostFunc getCost)
//{
//    if (ToulBar2::verbose >= 3) cout << "findFullSupport C" << x->getName() << "," << y->getName() << endl;
//    bool supportBroken = false;
//    for (Variable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
//        int xindex = x->toIndex(*iterX);
//        Value support = supportX[xindex];
//        if (y->cannotbe(support) || getCost(c, *iterX, support) + y->getCost(support) > 0) {
//            Value minCostValue = y->getInf();
//            Cost minCost = getCost(c, *iterX, minCostValue) + y->getCost(minCostValue);
//            Variable::iterator iterY = y->begin();
//            for (++iterY; minCost > 0 && iterY != y->end(); ++iterY) {
//                Cost cost = getCost(c, *iterX, *iterY) + y->getCost(*iterY);
//                if (cost < minCost) {
//                    minCost = cost;
//                    minCostValue = *iterY;
//                }
//            }
//            if (minCost > 0) {
//                // extend unary to binary
//                for (Variable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
//                    Cost cost = getCost(c, *iterX, *iterY);
//                    if (minCost > cost) {
//                        deltaCostsY[*iterY] -= (minCost - cost);  // Warning! Possible overflow???
//                        y->extend(*iterY, minCost - cost);
//                        int yindex = y->toIndex(*iterY);
//                        supportY[yindex] = *iterX;
//                     }
//                }
//                
//                // hard binary constraint costs are not changed
//                if (minCost + c->wcsp->getLb() <= c->wcsp->getUb()) deltaCostsX[xindex] += minCost;  // Warning! Possible overflow???
//                if (x->getSupport() == *iterX) supportBroken = true;
//                if (ToulBar2::verbose >= 2) cout << "binary projection of " << minCost << " from C" << x->getName() << "," << y->getName() << "(" << *iterX << "," << minCostValue << ")" << endl;
//                x->project(*iterX, minCost);
//            }
//            supportX[xindex] = minCostValue;
//        }
//    }
//    if (supportBroken) {
//        x->findSupport();
//    }
//}

void BinaryConstraint::findFullSupportX()
{
    if (ToulBar2::verbose >= 3) cout << "findFullSupport C" << x->getName() << "," << y->getName() << endl;
    bool supportBroken = false;
    for (Variable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        int xindex = x->toIndex(*iterX);
        Value support = supportX[xindex];
        if (y->cannotbe(support) || getCost(*iterX, support) + y->getCost(support) > 0) {
            Value minCostValue = y->getInf();
            Cost minCost = getCost(*iterX, minCostValue) + y->getCost(minCostValue);
            Variable::iterator iterY = y->begin();
            for (++iterY; minCost > 0 && iterY != y->end(); ++iterY) {
                Cost cost = getCost(*iterX, *iterY) + y->getCost(*iterY);
                if (cost < minCost) {
                    minCost = cost;
                    minCostValue = *iterY;
                }
            }
            if (minCost > 0) {
                // extend unary to binary
                for (Variable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                    Cost cost = getCost(*iterX, *iterY);
                    if (minCost > cost) {
                        deltaCostsY[*iterY] -= (minCost - cost);  // Warning! Possible overflow???
                        y->extend(*iterY, minCost - cost);
                        int yindex = y->toIndex(*iterY);
                        supportY[yindex] = *iterX;
                     }
                }
                
                // hard binary constraint costs are not changed
                if (minCost + wcsp->getLb() <= wcsp->getUb()) deltaCostsX[xindex] += minCost;  // Warning! Possible overflow???
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

void BinaryConstraint::findFullSupportY()
{
    if (ToulBar2::verbose >= 3) cout << "findFullSupport C" << y->getName() << "," << x->getName() << endl;
    bool supportBroken = false;
    for (Variable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
        int yindex = y->toIndex(*iterY);
        Value support = supportY[yindex];
        if (x->cannotbe(support) || getCost(support, *iterY) + x->getCost(support) > 0) {
            Value minCostValue = x->getInf();
            Cost minCost = getCost(minCostValue, *iterY) + x->getCost(minCostValue);
            Variable::iterator iterX = x->begin();
            for (++iterX; minCost > 0 && iterX != x->end(); ++iterX) {
                Cost cost = getCost(*iterX, *iterY) + x->getCost(*iterX);
                if (cost < minCost) {
                    minCost = cost;
                    minCostValue = *iterX;
                }
            }
            if (minCost > 0) {
                // extend unary to binary
                for (Variable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
                    Cost cost = getCost(*iterX, *iterY);
                    if (minCost > cost) {
                        deltaCostsX[*iterX] -= (minCost - cost);  // Warning! Possible overflow???
                        x->extend(*iterX, minCost - cost);
                        int xindex = x->toIndex(*iterX);
                        supportX[xindex] = *iterY;
                     }
                }
                
                // hard binary constraint costs are not changed
                if (minCost + wcsp->getLb() <= wcsp->getUb()) deltaCostsY[yindex] += minCost;  // Warning! Possible overflow???
                if (y->getSupport() == *iterY) supportBroken = true;
                if (ToulBar2::verbose >= 2) cout << "binary projection of " << minCost << " from C" << y->getName() << "," << x->getName() << "(" << *iterY << "," << minCostValue << ")" << endl;
                y->project(*iterY, minCost);
            }
            supportY[yindex] = minCostValue;
        }
    }
    if (supportBroken) {
        y->findSupport();
    }
}

//void projection(BinaryConstraint *c, CostVariable *x, Value valueY, GetCostFunc getCost)
//{
//    bool supportBroken = false;
//    for (Variable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
//        Cost cost = getCost(c, *iterX, valueY);
//        if (cost > 0) {
//            if (x->getSupport() == *iterX) supportBroken = true;
//            if (ToulBar2::verbose >= 2) cout << "binary projection of " << cost << " on C(" << x->getName() << "," << *iterX << ")" << endl;
//            x->project(*iterX, cost);
//        }
//    }
//    if (supportBroken) {
//        x->findSupport();
//    }
//}

void BinaryConstraint::projectX()
{
    bool supportBroken = false;
    for (Variable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        Cost cost = getCost(*iterX, y->getValue());
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

void BinaryConstraint::projectY()
{
    bool supportBroken = false;
    for (Variable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
        Cost cost = getCost(x->getValue(), *iterY);
        if (cost > 0) {
            if (y->getSupport() == *iterY) supportBroken = true;
            if (ToulBar2::verbose >= 2) cout << "binary projection of " << cost << " on C(" << y->getName() << "," << *iterY << ")" << endl;
            y->project(*iterY, cost);
        }
    }
    if (supportBroken) {
        y->findSupport();
    }
}

//bool verify(BinaryConstraint *c, CostVariable *x, CostVariable *y, GetCostFunc getCost)
//{
//    for (Variable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
//        Cost minCost = getCost(c, *iterX, y->getInf());
//        if (x->wcspIndex < y->wcspIndex) minCost += y->getCost(y->getInf());
//        Variable::iterator iterY = y->begin();
//        for (++iterY; minCost > 0 && iterY != y->end(); ++iterY) {
//            Cost cost = getCost(c, *iterX, *iterY);
//            if (x->wcspIndex < y->wcspIndex) cost += y->getCost(*iterY);
//            if (cost < minCost) {
//                minCost = cost;
//            }
//        }
//        if (minCost > 0) {
//            cout << *c;
//            return false;
//        }
//    }
//    return true;
//}

bool BinaryConstraint::verifyX()
{
    for (Variable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        Cost minCost = getCost(*iterX, y->getInf());
        if (x->wcspIndex < y->wcspIndex) minCost += y->getCost(y->getInf());
        Variable::iterator iterY = y->begin();
        for (++iterY; minCost > 0 && iterY != y->end(); ++iterY) {
            Cost cost = getCost(*iterX, *iterY);
            if (x->wcspIndex < y->wcspIndex) cost += y->getCost(*iterY);
            if (cost < minCost) {
                minCost = cost;
            }
        }
        if (minCost > 0) {
            cout << this;
            return false;
        }
    }
    return true;
}

bool BinaryConstraint::verifyY()
{
    for (Variable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
        Cost minCost = getCost(x->getInf(), *iterY);
        if (y->wcspIndex < x->wcspIndex) minCost += x->getCost(x->getInf());
        Variable::iterator iterX = x->begin();
        for (++iterX; minCost > 0 && iterX != x->end(); ++iterX) {
            Cost cost = getCost(*iterX, *iterY);
            if (y->wcspIndex < x->wcspIndex) cost += x->getCost(*iterX);
            if (cost < minCost) {
                minCost = cost;
            }
        }
        if (minCost > 0) {
            cout << this;
            return false;
        }
    }
    return true;
}
