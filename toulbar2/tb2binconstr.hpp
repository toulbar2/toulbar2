/** \file tb2binconstr.hpp
 *  \brief Binary constraint applied on variables with enumerated domains.
 *
 */

#ifndef TB2BINCONSTR_HPP_
#define TB2BINCONSTR_HPP_

#include "tb2abstractconstr.hpp"

class BinaryConstraint;
typedef Cost (*GetCostFunc)(BinaryConstraint *c, Value vx, Value vy);

class BinaryConstraint : public AbstractBinaryConstraint
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

//    friend Cost getCost(BinaryConstraint *c, Value vx, Value vy) {return c->getCost(vx,vy);}
//    friend Cost getCostReverse(BinaryConstraint *c, Value vx, Value vy) {return c->getCost(vy,vx);}
//    friend void findSupport(BinaryConstraint *c, CostVariable *x, CostVariable *y, 
//            vector<Value> &supportX, vector<StoreCost> &deltaCostsX, GetCostFunc getCost);
//    friend void findFullSupport(BinaryConstraint *c, CostVariable *x, CostVariable *y, 
//            vector<Value> &supportX, vector<StoreCost> &deltaCostsX, 
//            vector<Value> &supportY, vector<StoreCost> &deltaCostsY, GetCostFunc getCost);
//    friend void projection(BinaryConstraint *c, CostVariable *x, Value valueY, GetCostFunc getCost);
//    friend bool verify(BinaryConstraint *c, CostVariable *x, CostVariable *y, GetCostFunc getCost);

//    void findSupportX() {::findSupport(this,x,y,supportX,deltaCostsX,::getCost);}
//    void findSupportY() {::findSupport(this,y,x,supportY,deltaCostsY,::getCostReverse);}
//    void findFullSupportX() {::findFullSupport(this,x,y,supportX,deltaCostsX,supportY,deltaCostsY,::getCost);}
//    void findFullSupportY() {::findFullSupport(this,y,x,supportY,deltaCostsY,supportX,deltaCostsX,::getCostReverse);}
//    void projectX() {::projection(this,x,y->getValue(),::getCost);}
//    void projectY() {::projection(this,y,x->getValue(),::getCostReverse);}
//    bool verifyX() {return ::verify(this,x,y,::getCost);}
//    bool verifyY() {return ::verify(this,y,x,::getCostReverse);}

    void findSupportX();
    void findSupportY();
    void findFullSupportX();
    void findFullSupportY();
    void projectX();
    void projectY();
    bool verifyX();
    bool verifyY();
    
public:
    BinaryConstraint(CostVariable *xx, CostVariable *yy, vector<Cost> &tab, StoreStack<Cost,Cost> *store);

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
