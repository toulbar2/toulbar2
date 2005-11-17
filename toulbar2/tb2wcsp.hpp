/** \file tb2wcsp.hpp
 *  \brief General ToulBar2 header.
 * 
 */
 
#ifndef WCSP_HPP_
#define WCSP_HPP_

#include "tb2constraint.hpp"
#include "tb2costvar.hpp"

class WCSP {
    Variable *objective;
    Store *storeData;
    vector<CostVariable *> vars;
    vector<Constraint *> constrs;
    int NCBucketSize;
    vector< CostVariableList > NCBuckets;     // vector of backtrackable lists
    Queue NC;         // non backtrackable list
    Queue AC;         // non backtrackable list
    Cost upperBound;
    bool objectiveChanged;
    long long nbNodes;
    long long nbBacktracks;
    vector<Variable *> readVars;
    
    int link(Variable *x);
    int link(Constraint *c);

    void changeNCBucket(int oldBucket, int newBucket, DLink<CostVariable *> *elt) {
        if (oldBucket >= 0) NCBuckets[oldBucket].erase(elt, true);
        if (newBucket >= 0) NCBuckets[newBucket].push_back(elt, true);
    }

    // Search methods
    CostVariable *getVarMinDomainDivMaxDegree();
    void binaryChoicePoint(CostVariable *x, Value value);
    void recursiveSolve();

    // make it private because we don't want copy nor assignment
    WCSP(const WCSP &wcsp);
    WCSP& operator=(const WCSP &wcsp);
    
public:
    WCSP(Variable *objective, Store *s);

    ~WCSP() {
        for (unsigned int i=0; i<vars.size(); i++) delete vars[i];
        for (unsigned int i=0; i<constrs.size(); i++) delete constrs[i];
        for (unsigned int i=0; i<readVars.size(); i++) delete readVars[i];
    }

    // Warning! Can return NULL if index corresponds to the objective
    CostVariable *getVar(int index) {return (index >= 0)?vars[index]:NULL;}
    Constraint *getConstraint(int index) {return constrs[index];}
    Variable *getObjective() {return objective;}
    Cost getLb() const {return objective->getInf();}
    Cost getUb() const {return objective->getSup();}
    Cost getUnaryCost(int varIndex, Value v) const {return vars[varIndex]->getCost(v);}

    // avoid cost overflow
    Value add(const Value a, const Value b) const {return (a + b > getUb())?(getUb()+1):(a+b);}
    // avoid weakning hard costs
    Value sub(const Value a, const Value b) const {assert(b <= a); return (a > getUb())?a:(a-b);}

    void store();
    void restore();
    void propagateNC();
    void propagateAC();
    void propagate();
    bool verify();
    
    void increase(int varIndex) {if (varIndex >= 0) vars[varIndex]->increaseWCSP(); else objectiveChanged=true;}
    void decrease(int varIndex) {if (varIndex >= 0) vars[varIndex]->decreaseWCSP(); else objectiveChanged=true;}
    void assign(int varIndex) {if (varIndex >= 0) vars[varIndex]->assignWCSP(); else objectiveChanged=true;}
    void remove(int varIndex, Value val) {if (varIndex >= 0) vars[varIndex]->removeWCSP(val);}
    
    int addBinaryConstraint(Variable *x, Variable *y, vector<Cost> &tab);
    void read_wcsp(const char *fileName);

    void printNCBuckets();
    
    bool solve();
    
    friend ostream& operator<<(ostream& os, WCSP &wcsp);
    
    friend class CostVariable;
};

#endif /*TB2WCSP_HPP_*/
