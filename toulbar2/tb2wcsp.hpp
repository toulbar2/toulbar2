/** \file tb2wcsp.hpp
 *  \brief Global soft constraint representing a weighted CSP
 * 
 */
 
#ifndef TB2WCSP_HPP_
#define TB2WCSP_HPP_

#include "toulbar2.hpp"
#include "tb2constraint.hpp"
#include "tb2variable.hpp"
#include "tb2enumvar.hpp"
#include "tb2intervar.hpp"

class WCSP : public WeightedCSP {
    static int wcspCounter; // count the number of instantiations of WCSP
    int instance; // instantiation occurence number
    Store *storeData;
    StoreCost lb;
    Cost ub;
    vector<Variable *> vars;
    vector<Constraint *> constrs;
    int NCBucketSize;
    vector< VariableList > NCBuckets;         // vector of backtrackable lists
    Queue NC;                                 // non backtrackable list
    Queue IncDec;                             // non backtrackable list
    Queue AC;                                 // non backtrackable list
    Queue DAC;                                // non backtrackable list
    bool objectiveChanged;
    long long nbNodes;                        // used as a time-stamp by Queue methods

    Store *getStore() {return storeData;}
    
    void link(Variable *x) {vars.push_back(x);}
    void link(Constraint *c) {constrs.push_back(c);}

    void changeNCBucket(int oldBucket, int newBucket, DLink<Variable *> *elt) {
        if (oldBucket >= 0) NCBuckets[oldBucket].erase(elt, true);
        if (newBucket >= 0) NCBuckets[newBucket].push_back(elt, true);
    }
    void printNCBuckets();

    void propagateNC();
    void propagateIncDec();
    void propagateAC();
    void propagateDAC();

    void sortConstraints();

    // make it private because we don't want copy nor assignment
    WCSP(const WCSP &wcsp);
    WCSP& operator=(const WCSP &wcsp);
    
public:
    WCSP(Store *s, Cost upperBound);
    
    virtual ~WCSP();
    
    int getIndex() const {return instance;} // instantiation occurence number of current WCSP object
    
    Cost getLb() const {return lb;}
    Cost getUb() const {return ub;}

    // avoid cost overflow
    Cost add(const Cost a, const Cost b) const {return (a + b >= getUb())?getUb():(a+b);}
    // avoid weakning hard costs
    Cost sub(const Cost a, const Cost b) const {assert(b <= a); return (a >= getUb())?a:(a-b);}

    void updateUb(Cost newUb) {
        ub = min(ub,newUb);
    }
    void enforceUb() {
        if (lb >= ub) throw Contradiction();
        objectiveChanged=true;
    }
    void increaseLb(Cost newLb) {
        if (newLb > lb) {
            if (newLb >= ub) throw Contradiction();
            lb = newLb;
            objectiveChanged=true;
        }
    }
    void decreaseUb(Cost newUb) {
        if (newUb < ub) {
            if (newUb <= lb) throw Contradiction();
            ub = newUb;
            objectiveChanged=true;
        }
    }

    bool enumerated(int varIndex) const {return vars[varIndex]->enumerated();}
    
    string getName(int varIndex) const {return vars[varIndex]->getName();}
    Value getInf(int varIndex) const {return vars[varIndex]->getInf();}
    Value getSup(int varIndex) const {return vars[varIndex]->getSup();}
    Value getValue(int varIndex) const {return vars[varIndex]->getValue();}
    unsigned int getDomainSize(int varIndex) const {return vars[varIndex]->getDomainSize();}

    bool assigned(int varIndex) const {return vars[varIndex]->assigned();}
    bool unassigned(int varIndex) const {return vars[varIndex]->unassigned();}
    bool canbe(int varIndex, Value v) const {return vars[varIndex]->canbe(v);}
    bool cannotbe(int varIndex, Value v) const {return vars[varIndex]->cannotbe(v);}
       
    void increase(int varIndex, Value newInf) {vars[varIndex]->increase(newInf);}
    void decrease(int varIndex, Value newSup) {vars[varIndex]->decrease(newSup);}
    void assign(int varIndex, Value newValue) {vars[varIndex]->assign(newValue);}
    void remove(int varIndex, Value remValue) {vars[varIndex]->remove(remValue);}
        
    Cost getUnaryCost(int varIndex, Value v) const {return vars[varIndex]->getCost(v);}
    Value getSupport(int varIndex) const {return vars[varIndex]->getSupport();}
    
    int getDegree(int varIndex) const {return vars[varIndex]->getDegree();}

    void whenContradiction();       // after a contradiction, reset propagation queues and increase nbNodes
    void propagate();               // propagate until a fix point and increase nbNodes
    bool verify();

    unsigned int numberOfVariables() const {return vars.size();};
    unsigned int numberOfConstraints() const {return constrs.size();};
    Value getDomainSizeSum();       // total current number of values

    int makeEnumeratedVariable(string n, Value iinf, Value isup);
    int makeEnumeratedVariable(string n, Value *d, int dsize);
    int makeIntervalVariable(string n, Value iinf, Value isup);
    
    void postBinaryConstraint(int xIndex, int yIndex, vector<Cost> &costs);
    void postSupxyc(int xIndex, int yIndex, Value cste);
    
    void read_wcsp(const char *fileName);
    
    friend ostream& operator<<(ostream& os, WCSP &wcsp);
    
    friend class Variable;
    friend class EnumeratedVariable;
    friend class IntervalVariable;
    friend class Constraint;
    template<class T1, class T2> friend class AbstractBinaryConstraint;
    friend class Supxyc;
    friend class BinaryConstraint;
};

#endif /*TB2WCSP_HPP_*/
