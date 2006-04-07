/** \file tb2binconstr.hpp
 *  \brief Binary constraint applied on variables with enumerated domains.
 *
 */

#ifndef TB2BINCONSTR_HPP_
#define TB2BINCONSTR_HPP_

#include "tb2abstractconstr.hpp"
#include "tb2enumvar.hpp"

class BinaryConstraint;
typedef Cost (BinaryConstraint::*GetCostMember)(Value vx, Value vy);
#define GETCOST (this->*getCost)

class BinaryConstraint : public AbstractBinaryConstraint<EnumeratedVariable,EnumeratedVariable>
{
    unsigned int sizeX;
    unsigned int sizeY;
    vector<StoreCost> deltaCostsX;
    vector<StoreCost> deltaCostsY;
    vector<Cost> costs;
    vector<Value> supportX;
    vector<Value> supportY;

    Cost getCost(Value vx, Value vy) {
        int ix = x->toIndex(vx);
        int iy = y->toIndex(vy);
        Cost res = costs[ix * sizeY + iy] - deltaCostsX[ix] - deltaCostsY[iy];
        assert(res >= 0);
        return res;
    }
    Cost getCostReverse(Value vy, Value vx) {return getCost(vx,vy);}
    
    void findSupport(EnumeratedVariable *x, EnumeratedVariable *y, 
            vector<Value> &supportX, vector<StoreCost> &deltaCostsX, GetCostMember getCost);
    void findFullSupport(EnumeratedVariable *x, EnumeratedVariable *y, 
            vector<Value> &supportX, vector<StoreCost> &deltaCostsX, 
            vector<Value> &supportY, vector<StoreCost> &deltaCostsY, GetCostMember getCost);
    void projection(EnumeratedVariable *x, Value valueY, GetCostMember getCost);
    bool verify(EnumeratedVariable *x, EnumeratedVariable *y, GetCostMember getCost);

    void findSupportX() {findSupport(x,y,supportX,deltaCostsX,&BinaryConstraint::getCost);}
    void findSupportY() {findSupport(y,x,supportY,deltaCostsY,&BinaryConstraint::getCostReverse);}
    void findFullSupportX() {findFullSupport(x,y,supportX,deltaCostsX,supportY,deltaCostsY,&BinaryConstraint::getCost);}
    void findFullSupportY() {findFullSupport(y,x,supportY,deltaCostsY,supportX,deltaCostsX,&BinaryConstraint::getCostReverse);}
    void projectX() {projection(x,y->getValue(),&BinaryConstraint::getCost);}
    void projectY() {projection(y,x->getValue(),&BinaryConstraint::getCostReverse);}
    bool verifyX() {return verify(x,y,&BinaryConstraint::getCost);}
    bool verifyY() {return verify(y,x,&BinaryConstraint::getCostReverse);}
    
public:
    BinaryConstraint(WCSP *wcsp, EnumeratedVariable *xx, EnumeratedVariable *yy, vector<Cost> &tab, StoreStack<Cost, Cost> *storeCost);

    ~BinaryConstraint() {}
    
    void propagate() {
        if (x->wcspIndex < y->wcspIndex) {
            findSupportY();             // must do AC before DAC
            findFullSupportX();
        } else {
            findSupportX();             // must do AC before DAC
            findFullSupportY();
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
    
    void print(ostream& os);
};

#endif /*TB2BINCONSTR_HPP_*/
