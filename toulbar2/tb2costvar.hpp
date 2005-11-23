/** \file tb2costvar.hpp
 *  \brief Variable extended with unary costs.
 * 
 */
 
#ifndef TB2COSTVARIABLE_HPP_
#define TB2COSTVARIABLE_HPP_

#include "tb2variable.hpp"
#include "tb2queue.hpp"

class CostVariable : public WCSPLink
{
    bool enumerated;            // should be a constant
    Variable *var;
    ConstraintList constrs;

    StoreCost infCost;
    StoreCost supCost;
    StoreCost deltaCost;
    vector<StoreCost> costs;
    StoreValue support;     // Warning! unary support needs to be backtrackable in order to maintain findSupport
    
    // incremental NC data
    StoreCost maxCost;
    StoreValue maxCostValue;
    StoreInt NCBucket;
    DLink< CostVariable * > linkNCBucket;
    
    DLink<CostVariableWithTimeStamp> linkNCQueue;
    DLink<CostVariableWithTimeStamp> linkACQueue;

    // Do not call theses methods directly except in the corresponding functions in the WCSP class
    void increaseWCSP();
    void decreaseWCSP();
    void assignWCSP();
    void removeWCSP(Value val);

    void changeNCBucket(int newBucket);
        
    // make it private because we don't want copy nor assignment
    CostVariable(const CostVariable &x);
    CostVariable& operator=(const CostVariable &x);

public:    
    CostVariable(Variable *v, StoreStack<Cost,Cost> *storeCost, 
            StoreStack<ConstraintList, DLink<ConstraintLink> *> *storeConstraint,
            StoreStack<Value,Value> *storeValue, StoreStack<int,int> *storeInt);

    bool getEnumerated() const {return enumerated;}
    
    // Warning! Only valid if the variable is represented by an enumerated domain
    unsigned int getDomainInitSize() const {return var->getDomainInitSize();}
    unsigned int toIndex(Value v) const {return var->toIndex(v);}
    Value toValue(int index) const {return var->toValue(index);}

    Variable *getVar() const {return var;}
    ConstraintList *getConstrs() {return &constrs;}
    string getName() const {return var->getName();}        
    Value getInf() const {return var->getInf();}
    Value getSup() const {return var->getSup();}
    Value getValue() const {return var->getValue();}
    unsigned int getDomainSize() const {return var->getDomainSize();}

    bool assigned() const {return var->assigned();}
    bool unassigned() const {return var->unassigned();}
    bool canbe(Value v) const {return var->canbe(v);}
    bool cannotbe(Value v) const {return var->cannotbe(v);}

    void increase(Value newInf) {var->increase(newInf);}
    void decrease(Value newSup) {var->decrease(newSup);}
    void assign(Value newValue) {var->assign(newValue);}
    void remove(Value value) {var->remove(value);}

    // Warning! Only if the variable is represented by an interval
    Cost getInfCost() const {assert(!enumerated); return infCost - deltaCost;}
    Cost getSupCost() const {assert(!enumerated); return supCost - deltaCost;}
    void projectInfCost(Cost cost);
    void projectSupCost(Cost cost);

    // Warning! Only if the variable is represented by an enumerated domain
    void project(Value value, Cost cost);
    void extend(Value value, Cost cost);
#ifdef FASTWCSP
    void projectFast(Value value, Cost cost) {
        assert(cost >= 0);
        assert(enumerated);
        costs[toIndex(value)] += cost;
    }
    void extendFast(Value value, Cost cost) {
        assert(cost >= 0);
        assert(enumerated);
        assert(costs[toIndex(value)] >= cost);
        costs[toIndex(value)] -= cost;
    }
#endif
    Value getSupport() const {assert(enumerated);return support;}
    void setSupport(Value val) {assert(enumerated);support = val;}
    
    Cost getCost(Value value) const {
        if (enumerated) return costs[toIndex(value)] - deltaCost;
        else if (value == getInf()) return getInfCost();
        else if (value == getSup()) return getSupCost();
        else return 0;
    }
    
    void extendAll(Cost cost);
    void propagateAC();
    void propagateNC();
    void findSupport();
    bool verifyNC();

#ifdef FASTWCSP
    void extendAllFast(Cost cost) {
        assert(cost >= 0);
        deltaCost += cost;
    }
    void propagateNCFast();
    bool removeFast(Value val) {
        if (var->removeFast(val)) {
            assignWCSP();
            return true;
        } else return false;
    }
#endif
    
    Cost getMaxCostValue() const {return maxCostValue;}
    Cost getMaxCost() const {return maxCost;}
    void setMaxUnaryCost(Value a);
    void queueNC();
    void queueAC();

    typedef Variable::iterator iterator;
    iterator begin() {return var->begin();}
    iterator end() {return var->end();}
    iterator rbegin() {return var->rbegin();}
    iterator rend() {return var->rend();}
    iterator lower_bound(Value v) {return var->lower_bound(v);}
    iterator upper_bound(Value v) {return var->upper_bound(v);}

    DLink<ConstraintLink> *addConstraint(Constraint *c, int index);

    int getDegree() {return constrs.getSize();}

    friend ostream& operator<<(ostream& os, CostVariable &var);
    
    friend class WCSP;
};

#endif /*TB2COSTVARIABLE_HPP_*/
