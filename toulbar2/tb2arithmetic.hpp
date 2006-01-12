/** \file tb2arithmetic.hpp
 *  \brief Binary arithmetic soft constraints.
 *  Warning! Value and Cost must be the same.
 *
 */

#ifndef TB2ARITHMETIC_HPP_
#define TB2ARITHMETIC_HPP_

#include "tb2abstractconstr.hpp"

/** semantic: soft(x >= y + cst) : max( y + cst - x , 0 )
 */
class Supxyc : public AbstractBinaryConstraint
{
    Value cst;
    StoreCost deltaCost;
    StoreValue deltaValueXinf;
    StoreValue deltaValueYsup;
    StoreCost deltaCostXinf;
    StoreCost deltaCostYsup;
    
public:
    Supxyc(CostVariable *xx, CostVariable *yy, Value c, StoreStack<Cost,Cost> *storeCost, StoreStack<Value,Value> *storeValue);

    ~Supxyc() {}
    
    void propagate();
    
    void remove(int varIndex) {}

//    void increase(int index) {
//        if (index==0) 
//        else 
//    }
//    
//    void decrease(int index) {
//        if (index==0) 
//        else 
//    }
//    
    void assign(int varIndex) {
        if (getX()->assigned() && getY()->assigned()) deconnect();  // Warning! deconnection has to be done before the projection
        propagate();
    }
        
    bool verify();
    
    void print(ostream& os);
};

#endif /*TB2ARITHMETIC_HPP_*/
