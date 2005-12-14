/** \file tb2wcsp.hpp
 *  \brief Global soft constraint representing a weighted CSP
 * 
 */
 
#ifndef TB2WCSP_HPP_
#define TB2WCSP_HPP_

#include "tb2constraint.hpp"
#include "tb2costvar.hpp"

class WCSP {
    Variable *objective;
    Store *storeData;
    vector<CostVariable *> vars;
    vector<Constraint *> constrs;
    int NCBucketSize;
    vector< CostVariableList > NCBuckets;     // vector of backtrackable lists
    Queue NC;                                 // non backtrackable list
    Queue IncDec;                                 // non backtrackable list
    Queue AC;                                 // non backtrackable list
    Queue DAC;                                // non backtrackable list
    bool objectiveChanged;
    long long nbNodes;                        // used as a time-stamp by Queue methods
    
    int link(Variable *x);
    int link(Constraint *c);

    void changeNCBucket(int oldBucket, int newBucket, DLink<CostVariable *> *elt) {
        if (oldBucket >= 0) NCBuckets[oldBucket].erase(elt, true);
        if (newBucket >= 0) NCBuckets[newBucket].push_back(elt, true);
    }
    void printNCBuckets();

    void propagateNC();
    void propagateIncDec();
    void propagateAC();
    void propagateDAC();

    // make it private because we don't want copy nor assignment
    WCSP(const WCSP &wcsp);
    WCSP& operator=(const WCSP &wcsp);
    
public:
    WCSP(Variable *objective, Store *s);
    
    ~WCSP();

    Cost getLb() const {return objective->getInf();}
    Cost getUb() const {return objective->getSup();}

    Cost getUnaryCost(int varIndex, Value v) const {return (varIndex>=0)?vars[varIndex]->getUnaryCost(v):0;}

    Value getSupport(int varIndex) const {return (varIndex>=0)?vars[varIndex]->getSupport():getLb();}

    int getDegree(int varIndex) const {return (varIndex>=0)?vars[varIndex]->getDegree():0;}

    Value getDomainSizeSum();

    // avoid cost overflow
    Value add(const Value a, const Value b) const {return (a + b > getUb())?(getUb()+1):(a+b);}
    // avoid weakning hard costs
    Value sub(const Value a, const Value b) const {assert(b <= a); return (a > getUb())?a:(a-b);}

    void increaseLb(Value newlb) {objective->increase(this, newlb);}
    void decreaseUb(Value newub) {objective->decrease(this, newub);}
    
    void increase(bool noZeroCostRemoved, int varIndex) {if (varIndex >= 0) vars[varIndex]->increaseFromOutside(noZeroCostRemoved); else objectiveChanged=true;}
    void decrease(bool noZeroCostRemoved, int varIndex) {if (varIndex >= 0) vars[varIndex]->decreaseFromOutside(noZeroCostRemoved); else objectiveChanged=true;}
    void assign(int varIndex, Value prevInf, Value prevSup) {if (varIndex >= 0) vars[varIndex]->assignFromOutside(prevInf,prevSup); else objectiveChanged=true;}
    void remove(bool noZeroCostRemoved, int varIndex, Value val) {if (varIndex >= 0) vars[varIndex]->removeFromOutside(noZeroCostRemoved, val);}

    void whenContradiction();       // after a contradiction, reset propagation queues and increase nbNodes
    void propagate();               // propagate until a fix point and increase nbNodes
    bool verify();
    
    int postBinaryConstraint(Variable *x, Variable *y, vector<Cost> &tab);

    void sortConstraints();
    
    void read_wcsp(const char *fileName, Solver *solver);
    
    friend ostream& operator<<(ostream& os, WCSP &wcsp);
    
    friend class CostVariable;
};

#endif /*TB2WCSP_HPP_*/
