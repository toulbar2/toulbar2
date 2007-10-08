/** \file tb2arithmetic.hpp
 *  \brief Binary arithmetic soft constraints.
 *  Warning! Value and Cost must be the same.
 *
 */

#ifndef TB2ARITHMETIC_HPP_
#define TB2ARITHMETIC_HPP_

#include "tb2abstractconstr.hpp"
#include "tb2intervar.hpp"

#include <set>

/** semantic: soft(penalty, x in SetOfValues) : (x in Set)?0:penalty
 */
class Unary : public AbstractUnaryConstraint<IntervalVariable>
{
    set<Value> permitted;
    Cost penalty;
    StoreValue deltaValueXinf;
    StoreValue deltaValueXsup;
    
public:
    Unary(WCSP *wcsp, IntervalVariable *xx, Value *d, int dsize, Cost penalty, StoreStack<Value, Value> *storeValue);

    ~Unary() {}
    
    void propagate();
    
    void remove(int varIndex) {}

    void assign(int varIndex) {
        assert(connected());
        deconnect();  // Warning! deconnection has to be done before the projection
        if (permitted.find(x->getValue()) == permitted.end()) {
		  wcsp->increaseLb(wcsp->getLb() + penalty);
		}
    }
        
    bool verify();
    
    double  computeTightness() { return (double) permitted.size() / x->getDomainSize(); }
    
    void print(ostream& os);
    void dump(ostream& os);
};

/** semantic: soft(x >= y + cst) : max( (y + cst - x <= deltamax)?(y + cst - x):top , 0 )
 */
class Supxyc : public AbstractBinaryConstraint<Variable, Variable>
{
    Value cst;
    Value deltamax;
    StoreCost deltaCost;
    StoreValue deltaValueXinf;
    StoreValue deltaValueYsup;
    StoreCost deltaCostXinf;
    StoreCost deltaCostYsup;
    
public:
    Supxyc(WCSP *wcsp, Variable *xx, Variable *yy, Value c, Value delta,
        StoreStack<Cost, Cost> *storeCost, StoreStack<Value, Value> *storeValue);

    ~Supxyc() {}
    
    void propagate();
    
    void remove(int varIndex) {}

    void assign(int varIndex) {
        if (x->assigned() && y->assigned()) deconnect();
        propagate();
    }
        
    bool verify();
    
    double  computeTightness() { return 0; }
    
    void print(ostream& os);
    void dump(ostream& os);
};

/** semantic: soft(penalty, x >= y + csty or y >= x + cstx) : (x >= y + csty || y >= x + cstx)?0:penalty
 */
class Disjunction : public AbstractBinaryConstraint<Variable, Variable>
{
    Value cstx;
    Value csty;
    Cost penalty;
    StoreValue deltaValueXinf;
    StoreValue deltaValueYinf;
    StoreValue deltaValueXsup;
    StoreValue deltaValueYsup;
    
public:
    Disjunction(WCSP *wcsp, Variable *xx, Variable *yy, Value cxx, Value cyy, Cost penalty, 
        StoreStack<Value, Value> *storeValue);

    ~Disjunction() {}
    
    void propagate();
    
    void remove(int varIndex) {}

    bool verify();
    
    double  computeTightness() { return 0; }
    
    void print(ostream& os);
    void dump(ostream& os);
};

#endif /*TB2ARITHMETIC_HPP_*/
