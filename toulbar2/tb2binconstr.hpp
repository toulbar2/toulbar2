/** \file tb2binconstr.hpp
 *  \brief Binary constraint applied on variables with enumerated domains.
 *
 */

#ifndef TB2BINCONSTR_HPP_
#define TB2BINCONSTR_HPP_

#include "tb2abstractconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2wcsp.hpp"

class BinaryConstraint;
typedef Cost (BinaryConstraint::*GetCostMember)(Value vx, Value vy);

class BinaryConstraint : public AbstractBinaryConstraint<EnumeratedVariable,EnumeratedVariable>
{
    unsigned int sizeX;
    unsigned int sizeY;
    vector<StoreCost> deltaCostsX;
    vector<StoreCost> deltaCostsY;
    vector<StoreCost> costs;
    
    vector<Value> supportX;
    vector<Value> supportY;

    Cost getCostReverse(Value vy, Value vx) {return getCost(vx,vy);}
    
    template <GetCostMember getBinaryCost> void findSupport(EnumeratedVariable *x, EnumeratedVariable *y, 
            vector<Value> &supportX, vector<StoreCost> &deltaCostsX);
    template <GetCostMember getBinaryCost> void findFullSupport(EnumeratedVariable *x, EnumeratedVariable *y, 
            vector<Value> &supportX, vector<StoreCost> &deltaCostsX, 
            vector<Value> &supportY, vector<StoreCost> &deltaCostsY);
    template <GetCostMember getBinaryCost> void projection(EnumeratedVariable *x, EnumeratedVariable *y, Value valueY);
    template <GetCostMember getBinaryCost> bool verify(EnumeratedVariable *x, EnumeratedVariable *y);
    // return true if unary support of x is broken
    bool project(EnumeratedVariable *x, Value value, Cost cost, vector<StoreCost> &deltaCostsX)
    {
        // hard binary constraint costs are not changed
        if (cost + wcsp->getLb() < wcsp->getUb()) deltaCostsX[x->toIndex(value)] += cost;  // Warning! Possible overflow???
        x->project(value, cost);
        return (x->getSupport() == value);
    }
    
    void extend(EnumeratedVariable *x, Value value, Cost cost, vector<StoreCost> &deltaCostsX)
    {
        assert( ToulBar2::verbose < 4 || ((cout << "extend " << x->getName() << " (" << value << ") to C" << getVar(0)->getName() << "," << getVar(1)->getName() << " += " << cost << endl), true) );
        deltaCostsX[x->toIndex(value)] -= cost;  // Warning! Possible overflow???
        x->extend(value, cost);
    }
    
    void findSupportX() {findSupport<&BinaryConstraint::getCost>(x,y,supportX,deltaCostsX);}
    void findSupportY() {findSupport<&BinaryConstraint::getCostReverse>(y,x,supportY,deltaCostsY);}
    void findFullSupportX() {findFullSupport<&BinaryConstraint::getCost>(x,y,supportX,deltaCostsX,supportY,deltaCostsY);}
    void findFullSupportY() {findFullSupport<&BinaryConstraint::getCostReverse>(y,x,supportY,deltaCostsY,supportX,deltaCostsX);}
    void projectX() {projection<&BinaryConstraint::getCost>(x,y,y->getValue());}
    void projectY() {projection<&BinaryConstraint::getCostReverse>(y,x,x->getValue());}
    bool verifyX() {return verify<&BinaryConstraint::getCost>(x,y);}
    bool verifyY() {return verify<&BinaryConstraint::getCostReverse>(y,x);}
    bool projectX(Value value, Cost cost) {return project(x,value,cost,deltaCostsX);}
    bool projectY(Value value, Cost cost) {return project(y,value,cost,deltaCostsY);}
    void extendX(Value value, Cost cost) {extend(x,value,cost,deltaCostsX);}
    void extendY(Value value, Cost cost) {extend(y,value,cost,deltaCostsY);}
    
public:
    BinaryConstraint(WCSP *wcsp, EnumeratedVariable *xx, EnumeratedVariable *yy, vector<Cost> &tab, StoreStack<Cost, Cost> *storeCost);

    BinaryConstraint(WCSP *wcsp, StoreStack<Cost, Cost> *storeCost);

    ~BinaryConstraint() {}
    
    Cost getCost(Value vx, Value vy) {
        int ix = x->toIndex(vx);
        int iy = y->toIndex(vy);
        Cost res = costs[ix * sizeY + iy];
// BUG: incompatible with ternary preproject heurisitc
//        if (res + wcsp->getLb() < wcsp->getUb()) res -= deltaCostsX[ix] + deltaCostsY[iy];
        res -= deltaCostsX[ix] + deltaCostsY[iy];
        assert(res >= 0);
        return res;
    }
    
    Cost getCost(EnumeratedVariable *xx, EnumeratedVariable *yy, Value vx, Value vy) {
        int vindex[2];
        vindex[ getIndex(xx) ] = xx->toIndex(vx);
        vindex[ getIndex(yy) ] = yy->toIndex(vy);
        Cost res = costs[vindex[0] * sizeY + vindex[1]] - deltaCostsX[vindex[0]] - deltaCostsY[vindex[1]];
        assert(res >= 0);
        return res;
    }

    Cost getCostNoDelta(Value vx, Value vy) {
        int ix = x->toIndex(vx);
        int iy = y->toIndex(vy);
        Cost res = costs[ix * sizeY + iy];
        assert(res >= 0);
        return res;
    }


   void addcost( int vx, int vy, Cost mincost ) {
   	        int ix = x->toIndex(vx);
            int iy = y->toIndex(vy);
   	       	costs[ix * sizeY + iy] += mincost;
   }

   void addcost( EnumeratedVariable* xin, EnumeratedVariable* yin, int vx, int vy, Cost mincost ) {
      if (xin==x) costs[x->toIndex(vx) * sizeY + y->toIndex(vy)] += mincost;
      else costs[x->toIndex(vy) * sizeY + y->toIndex(vx)] += mincost;
   }

   void setCost( Cost c ) {   	
	    for (unsigned int a = 0; a < sizeX; a++) 
	         for (unsigned int b = 0; b < sizeY; b++) 
	                costs[a * sizeY + b] = c;
    }

   void setcost( int vx, int vy, Cost mincost ) {
   	        int ix = x->toIndex(vx);
            int iy = y->toIndex(vy);
   	       	costs[ix * sizeY + iy] = mincost;
    }

   void addCosts( EnumeratedVariable* xin, EnumeratedVariable* yin, vector<Cost>& costsin ) {
		assert(costsin.size() == costs.size());
		int ix, iy;
		for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
		for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
        	ix = xin->toIndex(*iterx);  	iy = yin->toIndex(*itery);
        	if(xin == x) costs[ix * sizeY + iy] += costsin[ix * sizeY + iy];
			else	     costs[ix * sizeY + iy] += costsin[iy * sizeX + ix];
	    }}
    }
    
   void addCosts( BinaryConstraint* xy ) {
		assert( ((x == xy->x) && (y == xy->y)) || ((x == xy->y) && (y == xy->x)) );
		int ix, iy;
		for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
		for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
        	ix = x->toIndex(*iterx); iy = y->toIndex(*itery);
        	if(xy->x == x) costs[ix * sizeY + iy] += xy->getCost(*iterx,*itery);
			else	       costs[ix * sizeY + iy] += xy->getCost(*itery,*iterx);
	    }}
    }
    
    void fillElimConstr( EnumeratedVariable* xin, EnumeratedVariable* yin)
	{
		x = xin;
		y = yin;
		sizeX = x->getDomainInitSize();
		sizeY = y->getDomainInitSize();
		linkX->removed = true;
		linkY->removed = true;
		linkX->content.constr = this;
		linkY->content.constr = this; 
		linkX->content.scopeIndex = 0;
		linkY->content.scopeIndex = 1;
        if (xin->wcspIndex < yin->wcspIndex) dacvar = 0; else dacvar = 1;
	}

    bool project(int varIndex, Value value, Cost cost) {
        if (varIndex == 0) return projectX(value, cost);
        else return projectY(value, cost);
    }
    void extend(int varIndex, Value value, Cost cost) {
        if (varIndex == 0) return extendX(value, cost);
        else return extendY(value, cost);
    }
    
    void propagate() {
        if (x->assigned()) {
            assign(0);
            return;
        }
        if (y->assigned()) {
            assign(1);
            return;
        }
        // delay true propagation in order to not interfer with ternary findFullSupportEAC
        if (getDACScopeIndex()==0) {
//            findSupportY();             // must do AC before DAC
//            if(connected()) findFullSupportX();
            x->queueAC(); 
            y->queueDAC();
            if (wcsp->isternary) {
                x->queueEAC1();
                y->queueEAC1();
            } else {
                if(connected()) {
                    int yindex = y->toIndex(y->getSupport());
                    if (y->cannotbe(y->getSupport()) || x->cannotbe(supportY[yindex]) ||
                        x->getCost(supportY[yindex]) > 0 || getCost(supportY[yindex], y->getSupport()) > 0) {
                        y->queueEAC1();
                    }
                }
            }
        } else {
//            findSupportX();             // must do AC before DAC
//            if(connected()) findFullSupportY();
            y->queueAC(); 
            x->queueDAC();
            if (wcsp->isternary) {
                x->queueEAC1();
                y->queueEAC1();
            } else {
                if(connected()) {
                    int xindex = x->toIndex(x->getSupport());
                    if (x->cannotbe(x->getSupport()) || y->cannotbe(supportX[xindex]) ||
                        y->getCost(supportX[xindex]) > 0 || getCost(x->getSupport(), supportX[xindex]) > 0) {
                        x->queueEAC1();
                    }
                }
            }
        }
    }
    void remove(int varIndex) {
        if (getDACScopeIndex()==0) {
            if (varIndex == 0) findSupportY();
        } else {
            if (varIndex == 1) findSupportX();
        }
    }
    void projectFromZero(int varIndex) {
        if (getDACScopeIndex()==0) {
            if (varIndex == 1) findFullSupportX();
        } else {
            if (varIndex == 0) findFullSupportY();
        }
    } 
    //Trick! instead of doing remove(index) now, let AC queue do the job. 
    //So several incdec events on the same constraint can be merged into one AC event
    void increase(int index) {if (index==0) x->queueAC(); else y->queueAC();}
    void decrease(int index) {if (index==0) x->queueAC(); else y->queueAC();}  // Trick! instead of remove(index);
    void assign(int varIndex) {
        deconnect();                    // Warning! deconnection has to be done before the projection
        if (varIndex == 0) projectY(); else projectX();
    }

  void fillEAC2(int varIndex) {
    if (getDACScopeIndex()==0) {
      if (varIndex==0) {
	   assert(y->canbe(y->getSupport()));
	   int yindex = y->toIndex(y->getSupport());
	   if (x->cannotbe(supportY[yindex]) || x->getCost(supportY[yindex]) > 0 || getCost(supportY[yindex],y->getSupport()) > 0) {
	       y->queueEAC2();
	   }
      }
    } else {
      if (varIndex==1) {
	   assert(x->canbe(x->getSupport()));
	   int xindex = x->toIndex(x->getSupport());
	   if (y->cannotbe(supportX[xindex]) || y->getCost(supportX[xindex]) > 0 || getCost(x->getSupport(),supportX[xindex]) > 0) {
	       x->queueEAC2();
	   }
      }
    }
  }
  
  bool isEAC(int varIndex, Value a) {
    if (varIndex==0) {
        int xindex = x->toIndex(a);
        if (y->cannotbe(supportX[xindex]) || y->getCost(supportX[xindex]) > 0 || getCost(a, supportX[xindex]) > 0) {
            for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                if (y->getCost(*iterY) == 0 && getCost(a,*iterY) == 0) {
                    supportX[xindex] = *iterY;
                    return true;
                }
            }
            return false;
        }
    } else {
        int yindex = y->toIndex(a);
        if (x->cannotbe(supportY[yindex]) || x->getCost(supportY[yindex]) > 0 || getCost(supportY[yindex], a) > 0) {
            for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
                if (x->getCost(*iterX) == 0 && getCost(*iterX, a) == 0) {
                    supportY[yindex] = *iterX;
                    return true;
                }
            }
            return false;
        }
    }
    return true;
  }
  
    void findFullSupportEAC(int varIndex) {
        if (varIndex == 0) findFullSupportX();
        else findFullSupportY();
    } 
    
    bool verify() {return verifyX() && verifyY();}
    
	double computeTightness() {
	   int count = 0;
	   double sum = 0;
	   for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
	      for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                sum += to_double(getCost(*iterX, *iterY));
				count++;
	       }
	    }
	    tight = sum / (double) count;
	    return tight;
	}

    void print(ostream& os);
};

#endif /*TB2BINCONSTR_HPP_*/
