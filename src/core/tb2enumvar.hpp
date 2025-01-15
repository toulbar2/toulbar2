/** \file tb2enumvar.hpp
 *  \brief Variable with domain represented by an enumerated domain.
 *
 */

#ifndef TB2ENUMVAR_HPP_
#define TB2ENUMVAR_HPP_

#include "tb2variable.hpp"
#include "tb2domain.hpp"

class EnumeratedVariable : public Variable {
protected:
    Domain domain;
    vector<StoreCost> costs;
    StoreCost deltaCost;
    StoreValue support; // Warning! the unary support has to be backtrackable
    Double trwsGamma; // The gamma factor used in TRW-S
    vector<string> valueNames;

    DLink<VariableWithTimeStamp> linkACQueue;
    DLink<VariableWithTimeStamp> linkDACQueue;
    DLink<VariableWithTimeStamp> linkEAC1Queue;
    DLink<VariableWithTimeStamp> linkEAC2Queue;
    DLink<VariableWithTimeStamp> linkDEEQueue;
    DLink<VariableWithTimeStamp> linkFEACQueue;

    bool watchForIncrease; ///< \warning should be true if there exists a cost function on this variable watching for increase events
    bool watchForDecrease; ///< \warning should be true if there exists a cost function on this variable watching for decrease events
    ConstraintLink DEE; ///< \brief residue for dead-end elimination
    vector<ConstraintLink> DEE2; ///< \brief residue for generalized dead-end elimination

    void init();

    void increaseFast(Value newInf); // Do not check for a support nor insert in NC and DAC queue
    void decreaseFast(Value newSup); // Do not check for a support nor insert in NC and DAC queue
    void removeFast(Value val); // Do not check for a support nor insert in NC and DAC queue

public:
    EnumeratedVariable(WCSP* wcsp, string n, Value iinf, Value isup);
    EnumeratedVariable(WCSP* wcsp, string n, vector<Value>& dom);

    bool enumerated() const FINAL { return true; }

    bool isValueNames() const { return valueNames.size() == getDomainInitSize(); }
    void addValueName(const string& vname) { valueNames.push_back(vname); }
    const string& getValueName(int index) const
    {
        static const string None = std::string("");
        if (isValueNames())
            return valueNames[index];
        else
            return None;
    }
    string getValueNameOrGenerate(int index) const
    {
        if (isValueNames())
            return valueNames[index];
        else
            return to_string("v") + to_string(toValue((unsigned int)index));
    }
    unsigned int toIndex(const string& vname)
    {
        vector<string>::iterator iter = find_if(valueNames.begin(), valueNames.end(), [&vname](const string& val) { return (val == vname); });
        return (unsigned int)std::distance(valueNames.begin(), iter);
    }

    unsigned int getDomainInitSize() const { return domain.getInitSize(); }
#if defined(WCSPFORMATONLY) && !defined(NUMBERJACK)
    unsigned int toIndex(Value v) const
    {
        return (unsigned int)v;
    }
    Value toValue(unsigned int idx) const { return idx; }
#else
    unsigned int toIndex(Value v) const
    {
        return domain.toIndex(v);
    }
    Value toValue(unsigned int idx) const { return domain.toValue(idx); }
#endif
    unsigned int toCurrentIndex(Value v)
    {
        return domain.toCurrentIndex(v);
    } // return value position in current domain
    unsigned int getDomainSize() const FINAL
    {
        if (assigned())
            return 1;
        else
            return domain.getSize(); ///< \warning can return a negative size in the case of a wrong list utilization
    }
    void getDomain(set<Value>& array);
    void getDomain(Value* array);
    void getDomainAndCost(ValueCost* array);

    bool canbe(Value v) const FINAL { return v >= inf && v <= sup && domain.canbe(v); }
    bool canbeAfterElim(Value v) const { return domain.canbe(v); }
    bool cannotbe(Value v) const FINAL { return v < inf || v > sup || domain.cannotbe(v); }

    void increase(Value newInf, bool isDecision = false) FINAL;
    void decrease(Value newSup, bool isDecision = false) FINAL;
    void remove(Value value, bool isDecision = false) FINAL;
    void assign(Value newValue, bool isDecision = false) FINAL;
    void assignWhenEliminated(Value newValue);
    void restoreInitialDomainWhenEliminated();
    void assignLS(Value newValue, ConstraintSet& delayedCtrs, bool force = false) FINAL;

    void project(Value value, Cost cost, bool delayed = false); ///< \param delayed if true, it does not check for forbidden cost/value and let node consistency do the job later
    void extend(Value value, Cost cost);
    void extendAll(Cost cost);
    Value getSupport() const FINAL { return support; }
    void setSupport(Value val)
    {
        if (support != val) {
            if (ToulBar2::verbose >= 8)
                cout << "change support for " << getName() << " from " << support << " to " << val << endl;
            support = val;
            if (ToulBar2::FullEAC)
                queueFEAC();
        }
    }
    inline Cost getCost(const Value value) const FINAL
    {
        return costs[toIndex(value)] - deltaCost;
    }
    Cost getBinaryCost(ConstraintLink c, Value myvalue, Value itsvalue);
    Cost getBinaryCost(BinaryConstraint* c, Value myvalue, Value itsvalue);

    vector<Cost> getCosts()
    {
        vector<Cost> costs_;
        for (Cost cost : costs) {
            costs_.push_back(cost);
        }
        return costs_;
    }
    void setCosts(const vector<Cost>& update)
    {
        assert(costs.size() == update.size());
        for (unsigned int a = 0; a < costs.size(); a++) {
            costs[a] = update[a];
        }
    }
    Cost getDeltaCost()
    {
        return deltaCost;
    }
    void setDeltaCost(Cost cost)
    {
        deltaCost = cost;
    }

    void setTRWSGamma(Double g) { trwsGamma = g; }
    Double getTRWSGamma() const { return trwsGamma; }

    Cost getInfCost() const FINAL { return costs[toIndex(getInf())] - deltaCost; }
    Cost getSupCost() const FINAL { return costs[toIndex(getSup())] - deltaCost; }
    void projectInfCost(Cost cost) FINAL;
    void projectSupCost(Cost cost) FINAL;

    void propagateNC() FINAL;
    bool verifyNC() FINAL;
    void queueAC(); // public method used also by tb2binconstr.hpp
    void queueDAC();
    void propagateAC();
    void propagateDAC();
    void findSupport();
    Cost normalizeTRWS();
    bool verify();

    void queueEAC1();
    void queueEAC2();
    void queueFEAC();
    void fillEAC2(bool self);
    bool isEAC(Value a);
    bool isEAC() FINAL;
    void propagateEAC();
    void setCostProvidingPartition();
    bool checkEACGreedySolution();
    bool reviseEACGreedySolution();
    void shrink() FINAL;

    void eliminate() FINAL;
    bool elimVar(BinaryConstraint* xy, BinaryConstraint* xy_duplicate = NULL);
    bool elimVar(ConstraintLink xylink, ConstraintLink xzlink);
    bool elimVar(TernaryConstraint* xyz);

    void queueDEE() FINAL;
    void propagateDEE(Value a, Value b, bool dee = true);
    bool verifyDEE(Value a, Value b);
    bool verifyDEE() FINAL;

    // merge current cost functions to x's list by replacing current variable y by x thanks to functional constraint xy (i.e., y := functional[x])
    void mergeTo(BinaryConstraint* xy, map<Value, Value>& functional);
    bool canbeMerged(EnumeratedVariable* x);

    class iterator;
    friend class iterator;
    class iterator {
        EnumeratedVariable* var;
        Domain::iterator diter;

    public:
        iterator() { var = NULL; }
        iterator(EnumeratedVariable* v, Domain::iterator iter)
            : var(v)
            , diter(iter)
        {
        }

        Value operator*() const { return *diter; }

        iterator& operator++()
        { // Prefix form //TODO: add a const_iterator to speed-up iterations (should be inlined?)
            if (var->unassigned())
                ++diter;
            else {
                if (*diter < var->getValue())
                    diter = var->domain.lower_bound(var->getValue());
                else
                    diter = var->domain.end();
            }
            return *this;
        }

        iterator& operator--()
        { // Prefix form
            if (var->unassigned())
                --diter;
            else {
                if (*diter > var->getValue())
                    diter = var->domain.lower_bound(var->getValue());
                else
                    diter = var->domain.end();
            }
            return *this;
        }

        // To see if you're at the end:
        bool operator==(const iterator& iter) const { return diter == iter.diter; }
        bool operator!=(const iterator& iter) const { return diter != iter.diter; }
    };
    iterator begin()
    {
        if (assigned())
            return iterator(this, domain.lower_bound(getValue()));
        else
            return iterator(this, domain.begin());
    }
    iterator end() { return iterator(this, domain.end()); }
    iterator rbegin()
    {
        if (assigned())
            return iterator(this, domain.upper_bound(getValue()));
        else
            return iterator(this, domain.rbegin());
    }
    iterator rend() { return end(); }

    // Finds the first available element whose value is greater or equal to v
    iterator lower_bound(Value v)
    {
        if (assigned()) {
            if (v <= getValue())
                return iterator(this, domain.lower_bound(getValue()));
            else
                return end();
        } else if (v > sup) {
            return end();
        } else
            return iterator(this, domain.lower_bound(max(getInf(), v)));
    }

    // Finds the first available element whose value is lower or equal to v
    iterator upper_bound(Value v)
    {
        if (assigned()) {
            if (v >= getValue())
                return iterator(this, domain.upper_bound(getValue()));
            else
                return end();
        } else if (v < inf) {
            return end();
        } else
            return iterator(this, domain.upper_bound(min(getSup(), v)));
    }

    void permuteDomain(int numberOfPermutations);
    void permuteDomain(Value a, Value b);
    ValueCost* sortDomain(vector<Cost>& costs);

    virtual void print(ostream& os);
};

#endif /*TB2ENUMVAR_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
