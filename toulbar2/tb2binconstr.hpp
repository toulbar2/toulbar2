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
    template <GetCostMember getBinaryCost> void projection(EnumeratedVariable *x, Value valueY);
    template <GetCostMember getBinaryCost> bool verify(EnumeratedVariable *x, EnumeratedVariable *y);

    void findSupportX() {findSupport<&BinaryConstraint::getCost>(x,y,supportX,deltaCostsX);}
    void findSupportY() {findSupport<&BinaryConstraint::getCostReverse>(y,x,supportY,deltaCostsY);}
    void findFullSupportX() {findFullSupport<&BinaryConstraint::getCost>(x,y,supportX,deltaCostsX,supportY,deltaCostsY);}
    void findFullSupportY() {findFullSupport<&BinaryConstraint::getCostReverse>(y,x,supportY,deltaCostsY,supportX,deltaCostsX);}
    void projectX() {projection<&BinaryConstraint::getCost>(x,y->getValue());}
    void projectY() {projection<&BinaryConstraint::getCostReverse>(y,x->getValue());}
    bool verifyX() {return verify<&BinaryConstraint::getCost>(x,y);}
    bool verifyY() {return verify<&BinaryConstraint::getCostReverse>(y,x);}
    
public:
    BinaryConstraint(WCSP *wcsp, EnumeratedVariable *xx, EnumeratedVariable *yy, vector<Cost> &tab, StoreStack<Cost, Cost> *storeCost);

    BinaryConstraint(WCSP *wcsp, StoreStack<Cost, Cost> *storeCost);

    ~BinaryConstraint() {}
    
    Cost getCost(Value vx, Value vy) {
        int ix = x->toIndex(vx);
        int iy = y->toIndex(vy);
        Cost res = costs[ix * sizeY + iy];
        if (res < wcsp->getUb()) res -= deltaCostsX[ix] + deltaCostsY[iy];
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
	}
    
    void propagate() {
        if(connected()) {
            if (x->wcspIndex < y->wcspIndex) {
                findSupportY();             // must do AC before DAC
                if(connected()) findFullSupportX();
            } else {
                findSupportX();             // must do AC before DAC
                if(connected()) findFullSupportY();
            }
        }
    }
    
    void remove(int varIndex) {
        if (x->wcspIndex < y->wcspIndex) {
            if (varIndex == 0) findSupportY();
        } else {
            if (varIndex == 1) findSupportX();
        }
    }
    void projectFromZero(int varIndex) {
        if (x->wcspIndex < y->wcspIndex) {
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
        
    bool verify() {return verifyX() && verifyY();}
    

	double computeTightness() {
	   int count = 0;
	   double sum = 0;
	   for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
	      for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
				sum += getCost(*iterX, *iterY);
				count++;
	       }
	    }
	    tight = sum / (double) count;
	    return tight;
	}

    
    void print(ostream& os);
};

#endif /*TB2BINCONSTR_HPP_*/
