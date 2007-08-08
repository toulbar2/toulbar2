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
protected:
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
        if (NOCUT(cost + wcsp->getLb(),wcsp->getUb())) deltaCostsX[x->toIndex(value)] += cost;  // Warning! Possible overflow???
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
		// BUG: incompatible with ternary projection ???
		//if (res >= wcsp->getUb() || res - deltaCostsX[ix] - deltaCostsY[iy] + wcsp->getLb() >= wcsp->getUb()) return wcsp->getUb();
		res -= deltaCostsX[ix] + deltaCostsY[iy];
        assert(res >= 0);
        return res;
    }
    
    Cost getCost(EnumeratedVariable *xx, EnumeratedVariable *yy, Value vx, Value vy) {
        int vindex[2];
        vindex[ getIndex(xx) ] = xx->toIndex(vx);
        vindex[ getIndex(yy) ] = yy->toIndex(vy);
        Cost res = costs[vindex[0] * sizeY + vindex[1]];
		//if (res >= wcsp->getUb() || res - deltaCostsX[vindex[0]] - deltaCostsY[vindex[1]] + wcsp->getLb() >= wcsp->getUb()) return wcsp->getUb();
		res -= deltaCostsX[vindex[0]] + deltaCostsY[vindex[1]];
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

   void setcost( EnumeratedVariable* xin, EnumeratedVariable* yin, int vx, int vy, Cost mincost ) {
      if (xin==x) costs[x->toIndex(vx) * sizeY + y->toIndex(vy)] = mincost;
      else costs[x->toIndex(vy) * sizeY + y->toIndex(vx)] = mincost;
   }

   void setcost( int vx, int vy, Cost mincost ) {
   	   costs[x->toIndex(vx) * sizeY + y->toIndex(vy)] = mincost;
   }


   void addCosts( EnumeratedVariable* xin, EnumeratedVariable* yin, vector<Cost>& costsin ) {
		assert(costsin.size() == costs.size());
		int ix, iy;
		for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
		for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
        	ix = x->toIndex(*iterx);  	iy = y->toIndex(*itery);
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
 
    void clearCosts() {
        for (unsigned int i=0; i<sizeX; i++) deltaCostsX[i] = 0;
        for (unsigned int j=0; j<sizeY; j++) deltaCostsY[j] = 0;
        for (unsigned int i=0; i<sizeX; i++) {
          for (unsigned int j=0; j<sizeY; j++) {
            costs[i * sizeY + j] = 0;
          }
        }
    }


	Cost evalsubstr( string& s, Constraint* ctr )
	{
		Value vals[2];
		int count = 0;
		
		for(int i=0;i<arity();i++) {
			EnumeratedVariable* var = (EnumeratedVariable*) getVar(i);
			int ind = ctr->getIndex( var );
			if(ind >= 0) { vals[i] = var->toValue(s[ind] - CHAR_FIRST); count++; }		
		}
		if(count == 2) return getCost(vals[0], vals[1]);
		else return 0;
	}    
    
    Cost getDefCost() { return wcsp->getUb(); }
	void setDefCost( Cost df ) {}
    
   
    EnumeratedVariable* xvar;
    EnumeratedVariable* yvar;
    EnumeratedVariable::iterator itvx;
    EnumeratedVariable::iterator itvy;

    void first() {
    	itvx = x->begin(); 
    	itvy = y->begin(); 
    	xvar = x;
    	yvar = y;
    }
   
    bool next( string& t, Cost& c) 
    { 
    	char tch[3];
    	if(itvx != xvar->end()) {
    		int ix = xvar->toIndex(*itvx);
	    	tch[0] = ix + CHAR_FIRST;
	    	if(itvy != yvar->end()) {
	    		int iy = yvar->toIndex(*itvy);
		    	tch[1] = iy + CHAR_FIRST;
		    	tch[2] = '\0';
		    	t = tch;
		    	c = getCost(xvar,yvar,*itvx, *itvy); 
				++itvy;
		    	return true;
	    	} else {
	    		++itvx;
	    		itvy = yvar->begin();
	    		return next(t,c);
	    	}
    	}
    	return false; 
    }
    
    void firstlex() { first(); }
    bool nextlex( string& t, Cost& c) { return next(t,c); } 
    

	void setTuple( string& t, Cost c, EnumeratedVariable** scope_in )  {
		Value v0 = scope_in[0]->toValue(t[0]-CHAR_FIRST);
		Value v1 = scope_in[1]->toValue(t[1]-CHAR_FIRST);
		setcost( scope_in[0], scope_in[1], v0, v1, c );				
	}
	
	void addtoTuple( string& t, Cost c, EnumeratedVariable** scope_in )  {
		Value v0 = scope_in[0]->toValue(t[0]-CHAR_FIRST);
		Value v1 = scope_in[1]->toValue(t[1]-CHAR_FIRST);
		addcost( scope_in[0], scope_in[1], v0, v1, c );				
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
    if (varIndex==getDACScopeIndex()) return true;
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
        if (varIndex==getDACScopeIndex()) return;
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

	EnumeratedVariable* commonVar( BinaryConstraint* bctr ) {
		if(getIndex(bctr->getVar(0)) >= 0) return (EnumeratedVariable*) bctr->getVar(0);
		else if (getIndex(bctr->getVar(1)) >= 0) return (EnumeratedVariable*) bctr->getVar(1);
		else return NULL;
	}

    void print(ostream& os);
    void dump(ostream& os);
};

#endif /*TB2BINCONSTR_HPP_*/
