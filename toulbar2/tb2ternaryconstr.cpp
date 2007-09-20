/*
 * ****** Binary constraint applied on variables with enumerated domains ******
 */

#include "tb2ternaryconstr.hpp"
#include "tb2wcsp.hpp"


/*
 * Constructors and misc.
 * 
 */

TernaryConstraint::TernaryConstraint(WCSP *wcsp, 
									 EnumeratedVariable *xx, 
									 EnumeratedVariable *yy,
									 EnumeratedVariable *zz,
									 BinaryConstraint *xy_,
									 BinaryConstraint *xz_,
									 BinaryConstraint *yz_,
									 vector<Cost> &tab, 
									 StoreStack<Cost, Cost> *storeCost) 
									 
					: AbstractTernaryConstraint<EnumeratedVariable,EnumeratedVariable,EnumeratedVariable>(wcsp, xx, yy, zz), 
					  sizeX(xx->getDomainInitSize()), sizeY(yy->getDomainInitSize()), sizeZ(zz->getDomainInitSize())
{
    deltaCostsX = vector<StoreCost>(sizeX,StoreCost(MIN_COST,storeCost));
    deltaCostsY = vector<StoreCost>(sizeY,StoreCost(MIN_COST,storeCost));
    deltaCostsZ = vector<StoreCost>(sizeZ,StoreCost(MIN_COST,storeCost));
    assert(tab.size() == sizeX * sizeY * sizeZ);
    supportX = vector< pair<Value,Value> >(sizeX,make_pair(y->getInf(),z->getInf()));
    supportY = vector< pair<Value,Value> >(sizeY,make_pair(x->getInf(),z->getInf()));
    supportZ = vector< pair<Value,Value> >(sizeZ,make_pair(x->getInf(),y->getInf()));

    costs = vector<StoreCost>(sizeX*sizeY*sizeZ,StoreCost(MIN_COST,storeCost));
    
    for (unsigned int a = 0; a < x->getDomainInitSize(); a++) 
         for (unsigned int b = 0; b < y->getDomainInitSize(); b++) 
            for (unsigned int c = 0; c < z->getDomainInitSize(); c++) 
                costs[a * sizeY * sizeZ + b * sizeZ + c] = tab[a * sizeY * sizeZ + b * sizeZ + c];
 
    xy = xy_;
    xz = xz_;
    yz = yz_;
                
    propagate();
}


TernaryConstraint::TernaryConstraint(WCSP *wcsp, StoreStack<Cost, Cost> *storeCost) 
					: AbstractTernaryConstraint<EnumeratedVariable,EnumeratedVariable,EnumeratedVariable>(wcsp), 
					  sizeX(wcsp->maxdomainsize), sizeY(wcsp->maxdomainsize), sizeZ(wcsp->maxdomainsize)
{
	unsigned int maxdom = wcsp->maxdomainsize;
    deltaCostsX = vector<StoreCost>(maxdom,StoreCost(MIN_COST,storeCost));
    deltaCostsY = vector<StoreCost>(maxdom,StoreCost(MIN_COST,storeCost));
    deltaCostsZ = vector<StoreCost>(maxdom,StoreCost(MIN_COST,storeCost));
    supportX = vector< pair<Value,Value> >(maxdom);
    supportY = vector< pair<Value,Value> >(maxdom);
    supportZ = vector< pair<Value,Value> >(maxdom);

    costs = vector<StoreCost>(maxdom*maxdom*maxdom,StoreCost(MIN_COST,storeCost));
    
    for (unsigned int a = 0; a < maxdom; a++) 
       for (unsigned int b = 0; b < maxdom; b++) 
           for (unsigned int c = 0; c < maxdom; c++) 
               costs[a * maxdom * maxdom + b * maxdom + c] = MIN_COST;
 
    xy = NULL;
    xz = NULL;
    yz = NULL ;
}


double TernaryConstraint::computeTightness()
{
   int count = 0;
   double sum = 0;
   for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
      for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
	      for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
			sum += to_double(getCost(*iterX, *iterY, *iterZ));
			count++;
       }
     }
   }
   
   tight = sum / (double) count;
   return tight;
}


void TernaryConstraint::print(ostream& os)
{
    os << this << " TernaryConstraint(" << x->getName() << "," << y->getName() << "," << z->getName() << ")" << endl;
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

void TernaryConstraint::dump(ostream& os)
{
    os << "3 " << x->wcspIndex << " " << y->wcspIndex << " " << z->wcspIndex << " " << MIN_COST << " " << x->getDomainSize() * y->getDomainSize() * z->getDomainSize() << endl;
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
           for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
             os << *iterX << " " << *iterY << " " << *iterZ << " " << getCost(*iterX, *iterY, *iterZ) << endl;
           }
        }
    }
}

/*
 * Propagation methods
 * 
 */

void TernaryConstraint::findSupport(EnumeratedVariable *x, EnumeratedVariable *y, EnumeratedVariable *z,
            vector< pair<Value,Value> > &supportX, vector<StoreCost> &deltaCostsX, 
            vector< pair<Value,Value> > &supportY, vector< pair<Value,Value> > &supportZ)
{
    assert(getIndex(y) < getIndex(z));  // check that support.first/.second is consistent with y/z parameters
    assert(connected());
//    wcsp->revise(this);
    if (ToulBar2::verbose >= 3) cout << "findSupport C" << x->getName() << "," << y->getName() << "," << z->getName() << endl;
    bool supportBroken = false;
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        int xindex = x->toIndex(*iterX);
        pair<Value,Value> support = supportX[xindex];
        if (y->cannotbe(support.first) || z->cannotbe(support.second) || 
            getCost(x,y,z,*iterX,support.first,support.second) > MIN_COST) {

            Cost minCost = MAX_COST;
            for (EnumeratedVariable::iterator iterY = y->begin(); minCost > MIN_COST && iterY != y->end(); ++iterY) {
                for (EnumeratedVariable::iterator iterZ = z->begin(); minCost > MIN_COST && iterZ != z->end(); ++iterZ) {
                    Cost cost = getCost(x,y,z,*iterX,*iterY,*iterZ);
                    if (GLB(&minCost, cost)) {
                        support = make_pair(*iterY,*iterZ);
                    }
                }
            }
            if (minCost > MIN_COST) {
                // hard ternary constraint costs are not changed
                if (!CUT(minCost + wcsp->getLb(), wcsp->getUb())) deltaCostsX[xindex] += minCost;  // Warning! Possible overflow???
                if (x->getSupport() == *iterX || SUPPORTTEST(x->getCost(*iterX), minCost)) supportBroken = true;
                if (ToulBar2::verbose >= 2) cout << "ternary projection of " << minCost << " from C" << x->getName() << "," << y->getName()<< "," << z->getName() << "(" << *iterX << "," << support.first << "," << support.second << ")" << endl;
                x->project(*iterX, minCost);
                if (deconnected()) return;
            }
            supportX[xindex] = support;
            assert(getIndex(y) < getIndex(z));
            int yindex = y->toIndex(support.first);
            int zindex = z->toIndex(support.second);
            // warning! do not break DAC support for the variable in the constraint scope having the smallest wcspIndex
            if (getIndex(y)!=getDACScopeIndex()) {
                if (getIndex(x) < getIndex(z)) supportY[yindex] = make_pair(*iterX, support.second);
                else supportY[yindex] = make_pair(support.second, *iterX);
            }
            if (getIndex(z)!=getDACScopeIndex()) {
                if (getIndex(x) < getIndex(y)) supportZ[zindex] = make_pair(*iterX, support.first);
                else supportZ[zindex] = make_pair(support.first, *iterX);
            }
        }
    }
    if (supportBroken) {
        x->findSupport();
    }
}

// take into account associated binary constraints and perform unary extension to binary instead of ternary
void TernaryConstraint::findFullSupport(EnumeratedVariable *x, EnumeratedVariable *y, EnumeratedVariable *z,
            vector< pair<Value,Value> > &supportX, vector<StoreCost> &deltaCostsX, 
            vector< pair<Value,Value> > &supportY, vector<StoreCost> &deltaCostsY,
            vector< pair<Value,Value> > &supportZ, vector<StoreCost> &deltaCostsZ,
            BinaryConstraint* xy, BinaryConstraint* xz, BinaryConstraint* yz)
{
    assert(connected());
//    wcsp->revise(this);   
    if (ToulBar2::verbose >= 3) cout << "findFullSupportEAC C" << x->getName() << "," << y->getName()<< "," << z->getName() << endl;
    bool supportBroken = false;
    bool supportReversed = (getIndex(y) > getIndex(z));
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        int xindex = x->toIndex(*iterX);
        pair<Value,Value> support = (supportReversed)?make_pair(supportX[xindex].second,supportX[xindex].first):make_pair(supportX[xindex].first,supportX[xindex].second);
        if (y->cannotbe(support.first) || z->cannotbe(support.second) || 
            getCostWithBinaries(x,y,z,*iterX,support.first,support.second) + y->getCost(support.first) + z->getCost(support.second) > MIN_COST) {

            Cost minCost = MAX_COST;
            for (EnumeratedVariable::iterator iterY = y->begin(); minCost > MIN_COST && iterY != y->end(); ++iterY) {
                for (EnumeratedVariable::iterator iterZ = z->begin(); minCost > MIN_COST && iterZ != z->end(); ++iterZ) {
                    Cost cost = getCostWithBinaries(x,y,z,*iterX,*iterY,*iterZ) + y->getCost(*iterY) + z->getCost(*iterZ);
                    if (GLB(&minCost, cost)) {
                        support = make_pair(*iterY,*iterZ);
                    }
                }
            }
            assert(minCost < MAX_COST);
            
            if (minCost > MIN_COST) {
                // extend unary to binary
                for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                    if (y->getCost(*iterY) > MIN_COST) {
                        Cost costfromy = MIN_COST;
                        for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                            Cost cost = getCost(x,y,z,*iterX,*iterY,*iterZ);
                            if (LUBTEST(cost, minCost)) {
                                Cost zcost = z->getCost(*iterZ);
                                Cost xycost = (xy->connected())?xy->getCost(x,y,*iterX,*iterY):MIN_COST;
                                Cost xzcost = (xz->connected())?xz->getCost(x,z,*iterX,*iterZ):MIN_COST;
                                Cost yzcost = (yz->connected())?yz->getCost(y,z,*iterY,*iterZ):MIN_COST;
                                Cost remain = minCost - (cost + xycost + xzcost + yzcost + zcost);
                                LUB(&costfromy, remain);
                            }
                        }
                        assert(costfromy <= y->getCost(*iterY));
                        if (costfromy > MIN_COST) {
                            xy->reconnect();    // must be done before using the constraint
                            xy->extend(xy->getIndex(y), *iterY, costfromy);
                        }
                    }
                }
                for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                    if (z->getCost(*iterZ) > MIN_COST) {
                        Cost costfromz = MIN_COST;
                        for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                            Cost cost = getCost(x,y,z,*iterX,*iterY,*iterZ);
                            if (LUBTEST(cost, minCost)) {
                                Cost xycost = (xy->connected())?xy->getCost(x,y,*iterX,*iterY):MIN_COST;
                                Cost xzcost = (xz->connected())?xz->getCost(x,z,*iterX,*iterZ):MIN_COST;
                                Cost yzcost = (yz->connected())?yz->getCost(y,z,*iterY,*iterZ):MIN_COST;
                                Cost remain = minCost - (cost + xycost + xzcost + yzcost);
                                LUB(&costfromz, remain);
                            }
                        }
                        assert(costfromz <= z->getCost(*iterZ));
                        if (costfromz > MIN_COST) {
                            xz->reconnect();    // must be done before using the constraint
                            xz->extend(xz->getIndex(z), *iterZ, costfromz);
                        }
                    }
                }
                // extend binary to ternary
                if (xy->connected()) {
                    for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                        Cost xycost = xy->getCost(x,y,*iterX,*iterY);
                        if (xycost > MIN_COST) {
                            Cost costfromxy = MIN_COST;
                            for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                                Cost cost = getCost(x,y,z,*iterX,*iterY,*iterZ);
                                if (LUBTEST(cost, minCost)) {
                                    Cost xzcost = (xz->connected())?xz->getCost(x,z,*iterX,*iterZ):MIN_COST;
                                    Cost yzcost = (yz->connected())?yz->getCost(y,z,*iterY,*iterZ):MIN_COST;
                                    Cost remain = minCost - (cost + xzcost + yzcost);
                                    LUB(&costfromxy, remain);
                                }
                            }
                            assert(costfromxy <= xycost);
                            if (costfromxy > MIN_COST) {
                                extend(xy,x,y,z,*iterX,*iterY,costfromxy);
                            }
                        }
                    }
                }
                if (xz->connected()) {
                    for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                        Cost xzcost = xz->getCost(x,z,*iterX,*iterZ);
                        if (xzcost > MIN_COST) {
                            Cost costfromxz = MIN_COST;
                            for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                                Cost cost = getCost(x,y,z,*iterX,*iterY,*iterZ);
                                if (LUBTEST(cost, minCost)) {
                                    Cost yzcost = (yz->connected())?yz->getCost(y,z,*iterY,*iterZ):MIN_COST;
                                    Cost remain = minCost - (cost + yzcost);
                                    LUB(&costfromxz, remain);
                                }
                            }
                            assert(costfromxz <= xzcost);
                            if (costfromxz > MIN_COST) {
                                extend(xz,x,z,y,*iterX,*iterZ,costfromxz);
                            }
                        }
                    }
                }
                if (yz->connected()) {
                    for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                        for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                            Cost yzcost = yz->getCost(y,z,*iterY,*iterZ);
                            if (yzcost > MIN_COST) {
                                Cost cost = getCost(x,y,z,*iterX,*iterY,*iterZ);
                                if (LUBTEST(cost, minCost)) {
                                    assert(yzcost >= minCost - cost);
                                    extend(yz,y,z,x,*iterY,*iterZ,minCost - cost);
                                }
                            }
                        }
                    }
                }
                
                // hard ternary constraint costs are not changed
                if (!CUT(minCost + wcsp->getLb(), wcsp->getUb())) deltaCostsX[xindex] += minCost;  // Warning! Possible overflow???
                if (x->getSupport() == *iterX || SUPPORTTEST(x->getCost(*iterX), minCost)) supportBroken = true;
                if (ToulBar2::verbose >= 2) cout << "ternary projection of " << minCost << " from C" << x->getName() << "," << y->getName()<< "," << z->getName() << "(" << *iterX << "," << support.first << "," << support.second << ")" << endl;
                x->project(*iterX, minCost);
                if (deconnected()) return;
            }
            
            supportX[xindex] = (supportReversed)?make_pair(support.second, support.first):make_pair(support.first, support.second);
            
            int yindex = y->toIndex(support.first);
            int zindex = z->toIndex(support.second);
            if (getIndex(x) < getIndex(z)) supportY[yindex] = make_pair(*iterX, support.second);
            else supportY[yindex] = make_pair(support.second, *iterX);
            if (getIndex(x) < getIndex(y)) supportZ[zindex] = make_pair(*iterX, support.first);
            else supportZ[zindex] = make_pair(support.first, *iterX);
        }
    }
    if (supportBroken) {
        x->findSupport();
    }
}

void TernaryConstraint::projectTernaryBinary( BinaryConstraint* yzin )
{
    bool flag = false;
    EnumeratedVariable* xin = NULL;
    if(yzin->getIndex(x) < 0)      xin = x;
    else if(yzin->getIndex(y) < 0) xin = y;
    else if(yzin->getIndex(z) < 0) xin = z;
    
    EnumeratedVariable* xx = xin;
    EnumeratedVariable* yy = (EnumeratedVariable*) yzin->getVar(0);
    EnumeratedVariable* zz = (EnumeratedVariable*) yzin->getVar(1);

    for (EnumeratedVariable::iterator itery = yy->begin(); itery != yy->end(); ++itery) {
    for (EnumeratedVariable::iterator iterz = zz->begin(); iterz != zz->end(); ++iterz) {
        Cost mincost = wcsp->getUb();
        for (EnumeratedVariable::iterator iterx = xx->begin(); iterx != xx->end(); ++iterx) {
            Cost curcost = getCost(xx, yy, zz, *iterx, *itery, *iterz);                            
            GLB(&mincost, curcost);
        }
        if (mincost > MIN_COST) {
            flag = true;
            yzin->addcost(*itery,*iterz,mincost);   // project mincost to binary
            if (ToulBar2::verbose >= 2) cout << "ternary projection of " << mincost << " from C" << xx->getName() << "," << yy->getName() << "," << zz->getName() << " to C" << yy->getName() << "," << zz->getName() << "(" << *itery << "," << *iterz << ")" << endl;

			//BUG: connected() incompatible with preproject ternary ???            
			//if (connected() && mincost + wcsp->getLb() < wcsp->getUb()) {  
            // substract mincost from ternary constraint
                for (EnumeratedVariable::iterator iterx = xx->begin(); iterx != xx->end(); ++iterx) {
                    addcost(xx, yy, zz, *iterx, *itery, *iterz, -mincost);                            
                }
			//}               
        }
    }}
    
    if (flag) {
        yzin->reconnect();
        yzin->propagate();
    }
}

bool TernaryConstraint::isEAC(EnumeratedVariable *x, Value a, EnumeratedVariable *y, EnumeratedVariable *z,
            vector< pair<Value,Value> > &supportX)
{
    int xindex = x->toIndex(a);
    assert(getIndex(y) < getIndex(z));
    pair<Value,Value> support = supportX[xindex];
    if (y->cannotbe(support.first) || z->cannotbe(support.second) || 
        getCostWithBinaries(x,y,z,a,support.first,support.second) + y->getCost(support.first) + z->getCost(support.second) > MIN_COST) {
        for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
            for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                if (getCostWithBinaries(x,y,z,a,*iterY,*iterZ) + y->getCost(*iterY) + z->getCost(*iterZ) == MIN_COST) {
                    supportX[xindex] = make_pair(*iterY,*iterZ);
                    return true;
                }
            }
        }
//        cout << x->getName() << " = " << a << " not EAC due to constraint " << *this << endl; 
        return false;
    }
    return true;
}

bool TernaryConstraint::verify(EnumeratedVariable *x, EnumeratedVariable *y, EnumeratedVariable *z)
{
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        Cost minCost = MAX_COST;
        for (EnumeratedVariable::iterator iterY = y->begin(); minCost > MIN_COST && iterY != y->end(); ++iterY) {
            for (EnumeratedVariable::iterator iterZ = z->begin(); minCost > MIN_COST && iterZ != z->end(); ++iterZ) {
                Cost cost = getCost(x,y,z,*iterX,*iterY,*iterZ);
                if (getIndex(x)==getDACScopeIndex()) cost += y->getCost(*iterY) + z->getCost(*iterZ);
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
