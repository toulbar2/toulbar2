/** \file tb2variable.hpp
 *  \brief Abstract Variable class extended with unary costs.
 * 
 */
 
#ifndef TB2VARIABLE_HPP_
#define TB2VARIABLE_HPP_

#include "tb2btlist.hpp"
#include "tb2queue.hpp"
#include "tb2domain.hpp"

/*
 * Main class
 * 
 */

class Variable : public WCSPLink
{
protected:
    string name;
    
    StoreValue inf;
    StoreValue sup;

    ConstraintList constrs;

    StoreCost deltaCost;
    
    // incremental NC data
    StoreCost maxCost;
    StoreValue maxCostValue;
    StoreInt NCBucket;
    DLink< Variable * > linkNCBucket;
    
    DLink<VariableWithTimeStamp> linkNCQueue;
    DLink<VariableWithTimeStamp> linkIncDecQueue;

    void setMaxUnaryCost(Value a, Cost cost);
    void changeNCBucket(int newBucket);
        
    // make it private because we don't want copy nor assignment
    Variable(const Variable &x);
    Variable& operator=(const Variable &x);
    
public:    
    Variable(WCSP *w, string n, Value iinf, Value isup);
        
    virtual ~Variable() {}
    
    virtual bool enumerated() const =0;
    
    string getName() const {return name;}
    Value getInf() const {return inf;}
    Value getSup() const {return sup;}
    Value getValue() const {assert(assigned()); return inf;}
    virtual unsigned int getDomainSize() const =0;

    bool assigned() const {return inf == sup;}
    bool unassigned() const {return inf != sup;}
    virtual bool canbe(Value v) const {return v >= inf && v <= sup;}
    virtual bool cannotbe(Value v) const {return v < inf || v > sup;}
    
    virtual void increase(Value newInf) =0;
    virtual void decrease(Value newSup) =0;
    virtual void remove(Value remValue) =0;
    virtual void assign(Value newValue) =0;
    
    ConstraintList *getConstrs() {return &constrs;}
    int getDegree() {return constrs.getSize();}
    DLink<ConstraintLink> *link(Constraint *c, int index);
    void sortConstraints();

    virtual Cost getInfCost() const =0;
    virtual Cost getSupCost() const =0;
    virtual void projectInfCost(Cost cost) =0;
    virtual void projectSupCost(Cost cost) =0;
    virtual Cost getCost(const Value value) const =0;

    virtual Value getSupport() const {return inf;}      // If there is no defined support then return inf
    
    Cost getMaxCost() const {return maxCost;}
    Cost getMaxCostValue() const {return maxCostValue;}
    
    void extendAll(Cost cost);
    virtual void propagateNC() =0;    
    virtual bool verifyNC() =0;

    void queueNC();
    void queueInc();
    void queueDec();

    void propagateIncDec(int incdec);

/*
    class iterator;
    friend class iterator;
    class iterator {
    public:
        virtual Value operator*() const =0;
        
        virtual iterator &operator++() =0;     // Prefix form
        virtual iterator &operator--() =0;     // Prefix form

        // To see if you're at the end:
        virtual bool operator==(const iterator &iter) const =0;
        virtual bool operator!=(const iterator &iter) const =0;
    };
    virtual iterator begin() =0;
    virtual iterator end() =0;
    virtual iterator rbegin() =0;
    virtual iterator rend() =0;

    //Finds the first available element whose value is greater or equal to v
    virtual iterator lower_bound(Value v) =0;

    //Finds the first available element whose value is lower or equal to v
    virtual iterator upper_bound(Value v) =0;
*/

    virtual void print(ostream& os) =0;
    
    friend ostream& operator<<(ostream& os, Variable &var);
};

#endif /*TB2VARIABLE_HPP_*/
