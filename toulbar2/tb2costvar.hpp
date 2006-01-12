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
    StoreValue support;     // Warning! the unary support has to be backtrackable 
    
    // incremental NC data
    StoreCost maxCost;
    StoreValue maxCostValue;
    StoreInt NCBucket;
    DLink< CostVariable * > linkNCBucket;
    
    DLink<CostVariableWithTimeStamp> linkNCQueue;
    DLink<CostVariableWithTimeStamp> linkIncDecQueue;
    DLink<CostVariableWithTimeStamp> linkACQueue;
    DLink<CostVariableWithTimeStamp> linkDACQueue;

    // Do not call theses methods directly except in the corresponding functions in the WCSP class
    void increaseFromOutside(bool noZeroCostRemoved);
    void decreaseFromOutside(bool noZeroCostRemoved);
    void assignFromOutside(Value prevInf, Value prevSup);
    void removeFromOutside(bool noZeroCostRemoved, Value val);

    void changeNCBucket(int newBucket);
    Cost getMaxCostValue() const {return maxCostValue;}
    Cost getMaxCost() const {return maxCost;}
    void setMaxUnaryCost(Value a, Cost cost);
    void queueNC();
    void queueInc();
    void queueDec();
    void queueDAC();
        
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

    void increase(bool noZeroCostRemoved, Value newInf) {if (noZeroCostRemoved) var->increase(wcsp, newInf); else var->increase(NULL, newInf);}
    void decrease(bool noZeroCostRemoved, Value newSup) {if (noZeroCostRemoved) var->decrease(wcsp, newSup); else var->decrease(NULL, newSup);}
    void remove(bool noZeroCostRemoved, Value value) {if (noZeroCostRemoved) var->remove(wcsp, value); else var->remove(NULL, value);}
    void increase(Value newInf) {increase(false, newInf);}
    void decrease(Value newSup) {decrease(false, newSup);}
    void remove(Value value) {remove(false, value);}
    void assign(Value newValue) {var->assign(newValue);}

    Cost getInfCost() const {if (enumerated) return costs[toIndex(getInf())] - deltaCost; else return infCost - deltaCost;}
    Cost getSupCost() const {if (enumerated) return costs[toIndex(getSup())] - deltaCost; else return supCost - deltaCost;}
    void projectInfCost(Cost cost);
    void projectSupCost(Cost cost);

    // Warning! Only if the variable is represented by an enumerated domain
    void project(Value value, Cost cost);
    void extend(Value value, Cost cost);
    Value getSupport() const {assert(enumerated);return support;}
    void setSupport(Value val) {assert(enumerated);support = val;}    
    inline Cost getCost(const Value value) const {
        assert(enumerated);
        return costs[toIndex(value)] - deltaCost;
    }

    // this method can be applied to interval or enumerated domain
    Cost getUnaryCost(Value value) const {
        if (enumerated) return costs[toIndex(value)] - deltaCost;
        else if (value == getInf()) return getInfCost();
        else if (value == getSup()) return getSupCost();
        else return 0;
    }
    
    void extendAll(Cost cost);
    void queueAC();                     // public method used also by tb2binconstr.hpp
    void propagateIncDec(int incdec);
    void propagateAC();
    void propagateNC();

    // Warning! Only if the variable is represented by an enumerated domain
    void propagateDAC();
    void findSupport();
    
    bool verifyNC();

    typedef Variable::iterator iterator;
    iterator begin() {return var->begin();}
    iterator end() {return var->end();}
    iterator rbegin() {return var->rbegin();}
    iterator rend() {return var->rend();}
    iterator lower_bound(Value v) {return var->lower_bound(v);}
    iterator upper_bound(Value v) {return var->upper_bound(v);}

    DLink<ConstraintLink> *postConstraint(Constraint *c, int index);
    void sortConstraints();

    int getDegree() {return constrs.getSize();}

    friend ostream& operator<<(ostream& os, CostVariable &var);
    
    friend class WCSP;
};

#endif /*TB2COSTVARIABLE_HPP_*/
