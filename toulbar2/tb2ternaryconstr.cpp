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
    deltaCostsX = vector<StoreCost>(sizeX,StoreCost(0,storeCost));
    deltaCostsY = vector<StoreCost>(sizeY,StoreCost(0,storeCost));
    deltaCostsZ = vector<StoreCost>(sizeZ,StoreCost(0,storeCost));
    assert(tab.size() == sizeX * sizeY * sizeZ);
    supportX = vector< pair<Value,Value> >(sizeX,make_pair(y->getInf(),z->getInf()));
    supportY = vector< pair<Value,Value> >(sizeY,make_pair(x->getInf(),z->getInf()));
    supportZ = vector< pair<Value,Value> >(sizeZ,make_pair(x->getInf(),y->getInf()));

    costs = vector<StoreCost>(sizeX*sizeY*sizeZ,StoreCost(0,storeCost));
    
    for (unsigned int a = 0; a < x->getDomainInitSize(); a++) 
         for (unsigned int b = 0; b < y->getDomainInitSize(); b++) 
            for (unsigned int c = 0; c < z->getDomainInitSize(); c++) 
                costs[a * sizeY * sizeZ + b * sizeZ + c] = tab[a * sizeY * sizeZ + b * sizeZ + c];
 
    xy = xy_;
    xz = xz_;
    yz = yz_;
                
    propagate();
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
    os << "3 " << x->wcspIndex << " " << y->wcspIndex << " " << z->wcspIndex << " 0 " << x->getDomainSize() * y->getDomainSize() * z->getDomainSize() << endl;
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
    assert(connected());
//    wcsp->revise(this);
    if (ToulBar2::verbose >= 3) cout << "findSupport C" << x->getName() << "," << y->getName() << "," << z->getName() << endl;
    bool supportBroken = false;
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        int xindex = x->toIndex(*iterX);
        pair<Value,Value> support = supportX[xindex];
        if (y->cannotbe(support.first) || z->cannotbe(support.second) || 
            getCost(x,y,z,*iterX,support.first,support.second) > 0) {

            Cost minCost = MAX_COST;
            for (EnumeratedVariable::iterator iterY = y->begin(); minCost > 0 && iterY != y->end(); ++iterY) {
                for (EnumeratedVariable::iterator iterZ = z->begin(); minCost > 0 && iterZ != z->end(); ++iterZ) {
                    Cost cost = getCost(x,y,z,*iterX,*iterY,*iterZ);
                    if (cost < minCost) {
                        minCost = cost;
                        support = make_pair(*iterY,*iterZ);
                    }
                }
            }
            if (minCost > 0) {
                // hard ternary constraint costs are not changed
                if (minCost + wcsp->getLb() < wcsp->getUb()) deltaCostsX[xindex] += minCost;  // Warning! Possible overflow???
                if (x->getSupport() == *iterX) supportBroken = true;
                if (ToulBar2::verbose >= 2) cout << "ternary projection of " << minCost << " from C" << x->getName() << "," << y->getName()<< "," << z->getName() << "(" << *iterX << "," << support.first << "," << support.second << ")" << endl;
                x->project(*iterX, minCost);
                if (deconnected()) return;
            }
            supportX[xindex] = support;
            assert(getIndex(y) < getIndex(z));
            int yindex = y->toIndex(support.first);
            int zindex = z->toIndex(support.second);
            // warning! do not break DAC for the variable in the constraint scope having the smallest wcspIndex
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

void TernaryConstraint::findFullSupport(EnumeratedVariable *x, EnumeratedVariable *y, EnumeratedVariable *z,
            vector< pair<Value,Value> > &supportX, vector<StoreCost> &deltaCostsX, 
            vector< pair<Value,Value> > &supportY, vector<StoreCost> &deltaCostsY,
            vector< pair<Value,Value> > &supportZ, vector<StoreCost> &deltaCostsZ)
{
    assert(connected());
//    wcsp->revise(this);   
    if (ToulBar2::verbose >= 3) cout << "findFullSupport C" << x->getName() << "," << y->getName()<< "," << z->getName() << endl;
    bool supportBroken = false;
    bool yACRevise = false;
    bool zACRevise = false;
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        int xindex = x->toIndex(*iterX);
        pair<Value,Value> support = supportX[xindex];
        if (y->cannotbe(support.first) || z->cannotbe(support.second) || 
            getCost(x,y,z,*iterX,support.first,support.second) + y->getCost(support.first) + z->getCost(support.second) > 0) {

            Cost minCost = MAX_COST;
            for (EnumeratedVariable::iterator iterY = y->begin(); minCost > 0 && iterY != y->end(); ++iterY) {
                for (EnumeratedVariable::iterator iterZ = z->begin(); minCost > 0 && iterZ != z->end(); ++iterZ) {
                    Cost cost = getCost(x,y,z,*iterX,*iterY,*iterZ) + y->getCost(*iterY) + z->getCost(*iterZ);
                    if (cost < minCost) {
                        minCost = cost;
                        support = make_pair(*iterY,*iterZ);
                    }
                }
            }
            if (minCost > 0) {
                // extend unary to ternary
                for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                    for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                        Cost cost = getCost(x,y,z,*iterX,*iterY,*iterZ);
                        if (minCost > cost) {
                            int yindex = y->toIndex(*iterY);
                            int zindex = z->toIndex(*iterZ);
                            Cost ycost = y->getCost(*iterY);
                            Cost zcost = z->getCost(*iterZ);
                            Cost remain = minCost - cost;
                            // Try to extend unary cost from highest variable first
                            if (z->wcspIndex > y->wcspIndex) {
                                if (zcost > 0) {
                                    Cost costfromz = min(remain,zcost);
                                    deltaCostsZ[zindex] -= costfromz;  // Warning! Possible overflow???
                                    z->extend(*iterZ, costfromz);
                                    remain = remain - costfromz;
                                    yACRevise = true;   // AC may be broken by unary extension to ternary
                                }
                                if (remain > 0) {
                                    assert(remain <= ycost);
                                    deltaCostsY[yindex] -= remain;  // Warning! Possible overflow???
                                    y->extend(*iterY, remain);
                                    zACRevise = true;   // AC may be broken by unary extension to ternary
                                }
                            } else {
                                if (ycost > 0) {
                                    Cost costfromy = min(remain,ycost);
                                    deltaCostsY[yindex] -= costfromy;  // Warning! Possible overflow???
                                    y->extend(*iterY, costfromy);
                                    remain = remain - costfromy;
                                    zACRevise = true;   // AC may be broken by unary extension to ternary
                                }
                                if (remain > 0) {
                                    assert(remain <= zcost);
                                    deltaCostsZ[zindex] -= remain;  // Warning! Possible overflow???
                                    z->extend(*iterZ, remain);
                                    yACRevise = true;   // AC may be broken by unary extension to ternary
                                }
                            }
                            if (getIndex(x) < getIndex(z)) supportY[yindex] = make_pair(*iterX, *iterZ);
                            else supportY[yindex] = make_pair(*iterZ, *iterX);
                            if (getIndex(x) < getIndex(y)) supportZ[zindex] = make_pair(*iterX, *iterY);
                            else supportZ[zindex] = make_pair(*iterY, *iterX);
                         }
                    }
                }
                // hard ternary constraint costs are not changed
                if (minCost + wcsp->getLb() < wcsp->getUb()) deltaCostsX[xindex] += minCost;  // Warning! Possible overflow???
                if (x->getSupport() == *iterX) supportBroken = true;
                if (ToulBar2::verbose >= 2) cout << "ternary projection of " << minCost << " from C" << x->getName() << "," << y->getName()<< "," << z->getName() << "(" << *iterX << "," << support.first << "," << support.second << ")" << endl;
                x->project(*iterX, minCost);
                if (deconnected()) return;
            }
            supportX[xindex] = support;
            assert(getIndex(y) < getIndex(z));
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
    if (yACRevise && connected()) {
        if (getIndex(x) < getIndex(z)) findSupport(y,x,z,supportY,deltaCostsY,supportX,supportZ);
        else findSupport(y,z,x,supportY,deltaCostsY,supportZ,supportX);
    }
    if (zACRevise && connected()) {
        if (getIndex(x) < getIndex(y)) findSupport(z,x,y,supportZ,deltaCostsZ,supportX,supportY);
        else findSupport(z,y,x,supportZ,deltaCostsZ,supportY,supportX);
    }
}

// take into account associated binary constraints and perform unary extension to binary instead of ternary
void TernaryConstraint::findFullSupportEAC(EnumeratedVariable *x, EnumeratedVariable *y, EnumeratedVariable *z,
            vector< pair<Value,Value> > &supportX, vector<StoreCost> &deltaCostsX, 
            vector< pair<Value,Value> > &supportY, vector<StoreCost> &deltaCostsY,
            vector< pair<Value,Value> > &supportZ, vector<StoreCost> &deltaCostsZ,
            BinaryConstraint* xy, BinaryConstraint* xz, BinaryConstraint* yz)
{
    assert(connected());
//    wcsp->revise(this);   
    if (ToulBar2::verbose >= 3) cout << "findFullSupportEAC C" << x->getName() << "," << y->getName()<< "," << z->getName() << endl;
    bool supportBroken = false;
    bool yACRevise = false;
    bool zACRevise = false;
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        int xindex = x->toIndex(*iterX);
        pair<Value,Value> support = supportX[xindex];
        if (y->cannotbe(support.first) || z->cannotbe(support.second) || 
            getCostWithBinaries(x,y,z,*iterX,support.first,support.second) + y->getCost(support.first) + z->getCost(support.second) > 0) {

            Cost minCost = MAX_COST;
            for (EnumeratedVariable::iterator iterY = y->begin(); minCost > 0 && iterY != y->end(); ++iterY) {
                for (EnumeratedVariable::iterator iterZ = z->begin(); minCost > 0 && iterZ != z->end(); ++iterZ) {
                    Cost cost = getCostWithBinaries(x,y,z,*iterX,*iterY,*iterZ) + y->getCost(*iterY) + z->getCost(*iterZ);
                    if (cost < minCost) {
                        minCost = cost;
                        support = make_pair(*iterY,*iterZ);
                    }
                }
            }
            if (minCost > 0) {
                // extend unary to binary and extend binary to ternary
                for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                    for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                        Cost cost = getCost(x,y,z,*iterX,*iterY,*iterZ);
                        if (minCost > cost) {
                            int yindex = y->toIndex(*iterY);
                            int zindex = z->toIndex(*iterZ);
                            Cost ycost = y->getCost(*iterY);
                            Cost zcost = z->getCost(*iterZ);
                            Cost xycost = (xy->connected())?xy->getCost(x,y,*iterX,*iterY):0;
                            Cost xzcost = (xz->connected())?xz->getCost(x,z,*iterX,*iterZ):0;
                            Cost yzcost = (yz->connected())?yz->getCost(y,z,*iterY,*iterZ):0;
                            Cost remain = minCost - cost - xycost - xzcost - yzcost;
                            if (remain > 0) {
                                // extend unary to binary
                                // Try to extend unary cost from highest variable first
                                if (z->wcspIndex > y->wcspIndex) {
                                    if (zcost > 0) {
                                        Cost costfromz = min(remain,zcost);
                                        xz->reconnect();    // must be done before using the constraint
                                        xz->extend(xz->getIndex(z), *iterZ, costfromz);
                                        remain = remain - costfromz;
                                    }
                                    if (remain > 0) {
                                        assert(remain <= ycost);
                                        xy->reconnect();    // must be done before using the constraint
                                        xy->extend(xy->getIndex(y), *iterY, remain);
                                    }
                                } else {
                                    if (ycost > 0) {
                                        Cost costfromy = min(remain,ycost);
                                        xy->reconnect();    // must be done before using the constraint
                                        xy->extend(xy->getIndex(y), *iterY, costfromy);
                                        remain = remain - costfromy;
                                    }
                                    if (remain > 0) {
                                        assert(remain <= zcost);
                                        xz->reconnect();    // must be done before using the constraint
                                        xz->extend(xz->getIndex(z), *iterZ, remain);
                                    }
                                }
                            }
                            remain = minCost - cost;
                            xycost = (xy->connected())?xy->getCost(x,y,*iterX,*iterY):0;
                            xzcost = (xz->connected())?xz->getCost(x,z,*iterX,*iterZ):0;
                            yzcost = (yz->connected())?yz->getCost(y,z,*iterY,*iterZ):0;
                            assert(xycost + xzcost + yzcost >= remain);
                            // extend binary to ternary
                            if (xycost > 0) {
                                Cost costfromxy = min(remain,xycost);
                                extend(xy,x,y,z,*iterX,*iterY,costfromxy);
                                zACRevise = true;   // AC may be broken by binary extension to ternary
                                remain = remain - costfromxy;
                            }
                            if (remain > 0 && xzcost > 0) {
                                Cost costfromxz = min(remain,xzcost);
                                extend(xz,x,z,y,*iterX,*iterZ,costfromxz);
                                yACRevise = true;   // AC may be broken by binary extension to ternary
                                remain = remain - costfromxz;
                            }
                            if (remain > 0 && yzcost > 0) {
                                Cost costfromyz = min(remain,yzcost);
                                extend(yz,y,z,x,*iterY,*iterZ,costfromyz);
                                remain = remain - costfromyz;
                            }
                            assert(remain==0);
                            
                            if (getIndex(x) < getIndex(z)) supportY[yindex] = make_pair(*iterX, *iterZ);
                            else supportY[yindex] = make_pair(*iterZ, *iterX);
                            if (getIndex(x) < getIndex(y)) supportZ[zindex] = make_pair(*iterX, *iterY);
                            else supportZ[zindex] = make_pair(*iterY, *iterX);
                         }
                    }
                }
                // hard ternary constraint costs are not changed
                if (minCost + wcsp->getLb() < wcsp->getUb()) deltaCostsX[xindex] += minCost;  // Warning! Possible overflow???
                if (x->getSupport() == *iterX) supportBroken = true;
                if (ToulBar2::verbose >= 2) cout << "ternary projection of " << minCost << " from C" << x->getName() << "," << y->getName()<< "," << z->getName() << "(" << *iterX << "," << support.first << "," << support.second << ")" << endl;
                x->project(*iterX, minCost);
                if (deconnected()) return;
            }
            
            supportX[xindex] = support;
            
            assert(getIndex(y) < getIndex(z));
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
    if (yACRevise && connected()) {
        if (getIndex(x) < getIndex(z)) findSupport(y,x,z,supportY,deltaCostsY,supportX,supportZ);
        else findSupport(y,z,x,supportY,deltaCostsY,supportZ,supportX);
    }
    if (zACRevise && connected()) {
        if (getIndex(x) < getIndex(y)) findSupport(z,x,y,supportZ,deltaCostsZ,supportX,supportY);
        else findSupport(z,y,x,supportZ,deltaCostsZ,supportY,supportX);
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
        Cost mincost = MAX_COST;
        for (EnumeratedVariable::iterator iterx = xx->begin(); iterx != xx->end(); ++iterx) {
            Cost curcost = getCost(xx, yy, zz, *iterx, *itery, *iterz);                            
            if (curcost < mincost) mincost = curcost;
        }
        assert(mincost != MAX_COST);
        if (mincost > 0) {
            flag = true;
            yzin->addcost(*itery,*iterz,mincost);   // project mincost to binary
            if (ToulBar2::verbose >= 2) cout << "ternary projection of " << mincost << " from C" << xx->getName() << "," << yy->getName() << "," << zz->getName() << " to C" << yy->getName() << "," << zz->getName() << "(" << *itery << "," << *iterz << ")" << endl;
            
            if (connected() && mincost + wcsp->getLb() < wcsp->getUb()) {  
                // substract mincost from ternary constraint
                for (EnumeratedVariable::iterator iterx = xx->begin(); iterx != xx->end(); ++iterx) {
                    addcost(xx, yy, zz, *iterx, *itery, *iterz, -mincost);                            
                }
            }               
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
        getCostWithBinaries(x,y,z,a,support.first,support.second) + y->getCost(support.first) + z->getCost(support.second) > 0) {
        for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
            for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                if (getCostWithBinaries(x,y,z,a,*iterY,*iterZ) + y->getCost(*iterY) + z->getCost(*iterZ) == 0) {
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
        for (EnumeratedVariable::iterator iterY = y->begin(); minCost > 0 && iterY != y->end(); ++iterY) {
            for (EnumeratedVariable::iterator iterZ = z->begin(); minCost > 0 && iterZ != z->end(); ++iterZ) {
                Cost cost = getCost(x,y,z,*iterX,*iterY,*iterZ);
                if (getIndex(x)==getDACScopeIndex()) cost += y->getCost(*iterY) + z->getCost(*iterZ);
                if (cost < minCost) {
                    minCost = cost;
                }
            }
        }
        if (minCost > 0) {
            cout << "not FDAC: variable " << x->getName() << " value " << *iterX << " of " << *this;
            return false;
        }
    }
    return true;
}
