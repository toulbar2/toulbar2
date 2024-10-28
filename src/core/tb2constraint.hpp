/** \file tb2constraint.hpp
 *  \brief Abstract constraint.
 *
 */

#ifndef TB2CONSTRAINT_HPP_
#define TB2CONSTRAINT_HPP_

#include "tb2types.hpp"
#include "utils/tb2btlist.hpp"

class Constraint : public WCSPLink {
    Long conflictWeight;
    string name;
    Constraint* fromElim1; // remember the original constraint(s) from which this constraint is derived
    Constraint* fromElim2; // it can be from variable elimination during search or n-ary constraint projection

    // make it private because we don't want copy nor assignment
    Constraint(const Constraint& c);
    Constraint& operator=(const Constraint& c);

public:
    Constraint(WCSP* wcsp);
    Constraint(WCSP* wcsp, int elimCtrIndex);

    virtual ~Constraint() {}

    virtual bool extension() const { return false; } // return true if the cost function is defined in extension (with explicit tuples)
    virtual bool isBinary() const { return false; } // return true if the cost function class is a BinaryConstraint
    virtual bool isTernary() const { return false; } // return true if the cost function class is a TernaryConstraint
    virtual bool isNary() const { return false; } // return true if the cost function class is a NaryConstraint
    virtual bool isKnapsack() const { return false; } // return true if the cost function class is a KnapsackConstraint
    virtual bool isGlobal() const { return false; } // return true if it is a global cost function (flow-based monolithic propagation)
    //    virtual bool isTriangle() const {return false;} // return true if it is a triangle of three binary cost functions (maxRPC/PIC)

    virtual bool connected() const
    {
        cout << "dummy connected on (" << this << ")!" << endl;
        return true;
    }
    virtual bool deconnected() const
    {
        cout << "dummy deconnected on (" << this << ")!" << endl;
        return false;
    }
    // remove a constraint from the set of active constraints
    // (reuse=true if the constraint is empty and canbe reuse immediately)
    virtual void deconnect(bool reuse = false) { cout << "dummy deconnect on (" << this << ")!" << endl; }
    virtual void reconnect() { cout << "dummy reconnect on (" << this << ")!" << endl; }

    virtual int arity() const = 0;
    virtual Variable* getVar(int scopeIndex) const = 0;
    virtual int getIndex(Variable* var) const = 0;
    string getName() { return this->name; };
    void setName(const string& name) { this->name = name; };

    void conflict();
    virtual Long getConflictWeight() const { return conflictWeight; }
    virtual Long getConflictWeight(int varIndex) const { return conflictWeight; }
    virtual void incConflictWeight(Constraint* from)
    {
        if (from == this || deconnected())
            conflictWeight++;
        if (fromElim1)
            fromElim1->incConflictWeight(from);
        if (fromElim2)
            fromElim2->incConflictWeight(from);
    }
    void incConflictWeight(Long incval) { conflictWeight += incval; }
    virtual void resetConflictWeight() { conflictWeight = 1 + ((ToulBar2::weightedTightness) ? getTightness() : 0); } // side-effect: recompute tightness if needed
    void elimFrom(Constraint* from1, Constraint* from2 = NULL)
    {
        fromElim1 = from1;
        fromElim2 = from2;
    }

    double tight;
    virtual void resetTightness() { tight = -1; }
    double getTightness()
    {
        if (tight < 0)
            tight = computeTightness();
        assert(tight >= 0);
        return tight;
    }
    virtual double computeTightness() = 0;

    // return the smallest wcsp index in the constraint scope except for one variable having a forbidden scope index
    virtual int getSmallestVarIndexInScope(int forbiddenScopeIndex) = 0;
    virtual int getSmallestVarIndexInScope() = 0;
    virtual int getDACScopeIndex() const
    {
        cout << "dummy getDACScopeIndex on (" << this << ")!" << endl;
        return 0;
    }
    virtual void setDACScopeIndex() {}
    // return the smallest DAC ordering index in the constraint scope except for one variable having a forbidden scope index
    virtual int getSmallestDACIndexInScope(int forbiddenScopeIndex) = 0;
    virtual Variable* getDACVar(int scopeDACIndex) const = 0; // return scope variable with associated index in the sorted scope by DAC ordering

    virtual void propagate() = 0;
    virtual void increase(int index) { propagate(); }
    virtual void decrease(int index) { propagate(); }
    virtual void remove(int index) { propagate(); }
    virtual void projectFromZero(int index) {}
    virtual void assign(int index) { propagate(); }
    void assigns();

    virtual bool checkEACGreedySolution(int index = -1, Value a = 0); // if index is non-negative use value a instead of its support for the corresponding variable
    virtual bool reviseEACGreedySolution(int index = -1, Value a = 0); // if index is negative revise all variables else only one
    virtual void fillEAC2(int index) {}
    virtual bool isEAC(int index, Value a)
    {
        if (ToulBar2::FullEAC && !isSep()) {
            bool res = reviseEACGreedySolution(index, a);
            if (ToulBar2::verbose >= 4 && !res)
                cout << "Current greedy solution with variable index " << index << " assigned to " << a << " violates this constraint: " << *this << endl;
        }
        return true;
    }
    virtual void findFullSupportEAC(int index) {}

    virtual void linkCostProvidingPartition(int index, Variable* support) {}
    virtual void showCostProvidingPartition(int index) {}

    void projectLB(Cost cost);

    virtual bool verify() { return true; };

    virtual void print(ostream& os) { os << this << " Unknown constraint!"; }

    virtual void dump(ostream& os, bool original = true) { os << this << " Unknown constraint!"; }
    virtual void dump_CFN(ostream& os, bool original = true) { os << this << " Unknown constraint!"; }

    virtual Long getDomainSizeProduct() const; // warning! return LONGLONG_MAX if overflow occurs
    virtual Long getDomainInitSizeProduct() const; // warning! return LONGLONG_MAX if overflow occurs
    virtual Long size() const { return getDomainInitSizeProduct(); } ///< \brief number of tuples stored by the cost function
    virtual Long space() const
    {
        Long sz = size();
        Long eltsz = sizeof(Cost) + arity() * sizeof(tValue);
        if (sz < LONGLONG_MAX / eltsz)
            return (sz * eltsz);
        else
            return LONGLONG_MAX;
    } ///< \brief estimate of the cost function memory space size

    virtual void firstlex() {} ///< \brief enumerate all **valid** tuples of the cost function in lexicographic order (initialization call)
    virtual bool nextlex(Tuple& t, Cost& c)
    {
        cout << "dummy nextlex on (" << this << ")!" << endl;
        return false;
    } ///< \brief enumerate all **valid** tuples of the cost function in lexicographic order

    virtual void first() { firstlex(); } ///< \brief enumerate **valid** tuples of the cost function in undefined order, possibly skipping some valid tuples with a default cost (initialization call)
    virtual bool next(Tuple& t, Cost& c) { return nextlex(t, c); } ///< \brief enumerate **valid** tuples of the cost function in undefined order, possibly skipping some valid tuples with a default cost

    virtual void first(EnumeratedVariable* alpha, EnumeratedVariable* beta) {}
    virtual bool separability(EnumeratedVariable* alpha, EnumeratedVariable* beta) { return false; }
    virtual void separate(EnumeratedVariable* a, EnumeratedVariable* c) {}
    bool findConditionalIndependences();
    Cost squareminus(Cost c1, Cost c2, Cost top)
    {
        Cost c;
        if (c1 >= top && c2 >= top)
            c = top; // c = 0;
        else if (c1 >= top)
            c = MEDIUM_COST * top;
        else if (c2 >= top)
            c = -MEDIUM_COST * top;
        else
            c = c1 - c2;
        return c;
    }
    bool universe(Cost c1, Cost c2, Cost top)
    {
        if (c1 >= top && c2 >= top)
            return true;
        else
            return false;
    }
    bool verifySeparate(Constraint* ctr1, Constraint* ctr2);

    virtual void setTuple(const Tuple& t, Cost c) {}
    virtual void addtoTuple(const Tuple& t, Cost c) {}
    virtual void clearCosts() { cout << "dummy clearcosts on (" << this << ")!" << endl; }
    virtual void clearFiniteCosts() { cout << "dummy clearfinitecosts on (" << this << ")!" << endl; }

    virtual void getScope(TSCOPE& scope_inv) {}
    virtual Cost evalsubstr(const Tuple& s, Constraint* ctr)
    {
        cerr << "dummy evalsubstr call on:" << *this << endl;
        return MIN_COST;
    }
    virtual Cost evalsubstr(const Tuple& s, NaryConstraint* ctr)
    {
        cerr << "dummy evalsubstr call on:" << *this << endl;
        return MIN_COST;
    }
    virtual Cost getDefCost() { return MIN_COST; }
    virtual Cost getCost(); /// \warning all variables must be properly assigned and the returned cost might be different than the original cost due to propagation

    virtual bool universal(Cost zero = MIN_COST);
    virtual bool ishard(); ///< \brief returns true if all valid tuples have costs either zero or greater than the current upper bound
    virtual bool isfinite(); ///< \brief returns true if all valid tuples have costs strictly lower than the current upper bound minus the current lower bound

    virtual Cost getMinCost();
    virtual pair<pair<Cost, Cost>, pair<Cost, Cost>> getMaxCost(int index, Value a, Value b) { return std::make_pair(std::make_pair(MAX_COST, MAX_COST), std::make_pair(MAX_COST, MAX_COST)); }
    virtual Cost getMaxFiniteCost();
    virtual void setInfiniteCost(Cost ub) {}

    Constraint* copy(); ///< \brief returns a copy of itself as a new deconnected NaryConstraint (DO NOT USE DURING SEARCH!)

    void sumScopeIncluded(Constraint* ctr); ///< \brief adds a cost function into ourself ; scope of ctr must be included in our scope

    bool scopeIncluded(Constraint* ctr)
    {
        bool isincluded = true;
        int a_in = ctr->arity();
        if (a_in >= arity())
            return false;
        for (int i = 0; isincluded && i < a_in; i++)
            isincluded = isincluded && (getIndex(ctr->getVar(i)) >= 0);
        return isincluded;
    }

    void scopeCommon(TSCOPE& scope_out, Constraint* ctr)
    {
        TSCOPE scope1, scope2;
        getScope(scope1);
        ctr->getScope(scope2);

        TSCOPE::iterator it1 = scope1.begin();
        TSCOPE::iterator it2 = scope2.begin();
        while (it1 != scope1.end()) {
            it1->second = 0;
            ++it1;
        }
        while (it2 != scope2.end()) {
            it2->second = 0;
            ++it2;
        }
        set_intersection(scope1.begin(), scope1.end(),
            scope2.begin(), scope2.end(),
            inserter(scope_out, scope_out.begin()));
    }

    void scopeUnion(TSCOPE& scope_out, Constraint* ctr)
    {
        TSCOPE scope1, scope2;
        getScope(scope1);
        ctr->getScope(scope2);

        assert(arity() == (int)scope1.size());
        assert(ctr->arity() == (int)scope2.size());

        set_union(scope1.begin(), scope1.end(),
            scope2.begin(), scope2.end(),
            inserter(scope_out, scope_out.begin()));
    }

    void scopeDifference(TSCOPE& scope_out, Constraint* ctr)
    {
        TSCOPE scope1, scope2;
        getScope(scope1);
        ctr->getScope(scope2);
        set_difference(scope1.begin(), scope1.end(),
            scope2.begin(), scope2.end(),
            inserter(scope_out, scope_out.begin()));
    }

    int order(Constraint* ctr)
    {
        if (arity() < ctr->arity())
            return 1;
        else if (arity() > ctr->arity())
            return -1;
        TSCOPE scope1, scope2;
        getScope(scope1);
        ctr->getScope(scope2);
        TSCOPE::iterator it1 = scope1.begin();
        TSCOPE::iterator it2 = scope2.begin();
        while (it1 != scope1.end()) {
            if (it1->first < it2->first)
                return 1;
            else if (it1->first > it2->first)
                return -1;
            ++it1;
            ++it2;
        }
        return 0;
    }

    //   added for tree decomposition
    int cluster; ///< \brief warning! this value should not change during search (except for initialization in binary/ternary elimination pool)
    int getCluster() { return cluster; }
    void setCluster(int i);
    void assignCluster();

    bool isSep_; ///< \brief true if the constraint is used to record nogoods associated to a Separator
    void setSep() { isSep_ = true; }
    bool isSep() { return isSep_; }

    bool isDuplicate_; ///< \brief true if it exists another binary constraint with the same scope or if it is a ternary constraint with a duplicated binary constraint inside
    void setDuplicate()
    {
        isDuplicate_ = true;
        if (ToulBar2::verbose >= 1) {
            cout << *this << " set duplicate" << endl;
        }
    }
    bool isDuplicate() { return isDuplicate_; }

    virtual ConstraintSet subConstraint()
    {
        ConstraintSet s;
        return s;
    };

    friend ostream& operator<<(ostream& os, Constraint& c)
    {
        c.print(os);
        return os;
    }

    static bool cmpConstraintId(Constraint* c1, Constraint* c2);
    static bool cmpConstraintId(DLink<ConstraintLink>* c1, DLink<ConstraintLink>* c2);
    static bool cmpConstraintId(DLink<Constraint*>* c1, DLink<Constraint*>* c2);
    static bool cmpConstraintDAC(Constraint* c1, Constraint* c2);
    static bool cmpConstraintDAC(DLink<ConstraintLink>* c1, DLink<ConstraintLink>* c2);
    static bool cmpConstraintDAC(DLink<Constraint*>* c1, DLink<Constraint*>* c2);
    static bool cmpConstraintTightness(Constraint* c1, Constraint* c2);
    static bool cmpConstraintTightness(DLink<ConstraintLink>* c1, DLink<ConstraintLink>* c2);
    static bool cmpConstraintTightness(DLink<Constraint*>* c1, DLink<Constraint*>* c2);
    static bool cmpConstraintDACTightness(Constraint* c1, Constraint* c2);
    static bool cmpConstraintDACTightness(DLink<ConstraintLink>* c1, DLink<ConstraintLink>* c2);
    static bool cmpConstraintDACTightness(DLink<Constraint*>* c1, DLink<Constraint*>* c2);
    static bool cmpConstraintTightnessDAC(Constraint* c1, Constraint* c2);
    static bool cmpConstraintTightnessDAC(DLink<ConstraintLink>* c1, DLink<ConstraintLink>* c2);
    static bool cmpConstraintTightnessDAC(DLink<Constraint*>* c1, DLink<Constraint*>* c2);
    static bool cmpConstraintLAG(Constraint* c1, Constraint* c2);
    static bool cmpConstraintLAG(DLink<ConstraintLink>* c1, DLink<ConstraintLink>* c2);
    static bool cmpConstraintLAG(DLink<Constraint*>* c1, DLink<Constraint*>* c2);
    static bool cmpConstraintArity(Constraint* c1, Constraint* c2);
    static bool cmpConstraintArity(DLink<ConstraintLink>* c1, DLink<ConstraintLink>* c2);
    static bool cmpConstraintArity(DLink<Constraint*>* c1, DLink<Constraint*>* c2);
    static bool cmpConstraintArityDAC(Constraint* c1, Constraint* c2);
    static bool cmpConstraintArityDAC(DLink<ConstraintLink>* c1, DLink<ConstraintLink>* c2);
    static bool cmpConstraintArityDAC(DLink<Constraint*>* c1, DLink<Constraint*>* c2);

    // sort a list of constraints
    static int cmpConstraint(Constraint* c1, Constraint* c2);

    // sort a list of constraints related to a given variable
    static int cmpConstraintLink(DLink<ConstraintLink>* c1, DLink<ConstraintLink>* c2);

    // sort a backtrackable list of constraints
    static int cmpConstraintLinkPointer(DLink<Constraint*>* c1, DLink<Constraint*>* c2);
};

#endif /*TB2CONSTRAINT_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
