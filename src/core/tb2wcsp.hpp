/** \file tb2wcsp.hpp
 *  \brief Weighted constraint satisfaction problem modeling and local reasoning
 *
 */

#ifndef TB2WCSP_HPP_
#define TB2WCSP_HPP_

#include "toulbar2lib.hpp"
#include "tb2variable.hpp"
#include "tb2constraint.hpp"
#include "tb2enumvar.hpp"
#include "tb2intervar.hpp"
#include "search/tb2solver.hpp"

class NaryConstraint;
class GlobalConstraint;
class FlowBasedGlobalConstraint;
class AllDiffConstraint;
class GlobalCardinalityConstraint;
class SameConstraint;
class RegularFlowConstraint;

/** Concrete class WCSP containing a weighted constraint satisfaction problem
 *	- problem lower and upper bound
 *	- list of variables
 *	- list of cost functions (created before and during search by variable elimination) composed of constrs, elimBinConstrs, and elimTernConstrs
 *	- propagation queues (variable-based except for separators)
 *	- data for bounded variable elimination during search (limited to degree at most two)
 *
 * \note In most WCSP member functions, variables are referenced by their lexicographic index number (position in WCSP::vars vector, as returned by \e eg WCSP::makeEnumeratedVariable)
 *
 */

class WCSP FINAL : public WeightedCSP {
    static int wcspCounter; ///< count the number of instances of WCSP class
    int instance; ///< instance number
    string name; ///< problem name
    void* solver; ///< special hook to access solver information
    StoreCost lb; ///< current problem lower bound
    Cost ub; ///< current problem upper bound
    StoreCost negCost; ///< shifting value to be added to problem lowerbound when computing the partition function
    vector<Variable*> vars; ///< list of all variables
    vector<Variable*> divVariables; ///< list of variables submitted to diversity requirements
    vector<map<int, int>> divVarsId; // vector[j][idx] = index of the dual variable that encodes the diversity constraint on sol j at position idx
    vector<map<int, int>> divHVarsId; // vector[j][idx] = index of the hidden variable that encodes the diversity constraint on sol j between idx/idx+1
    vector<Value> bestValues; ///< hint for some value ordering heuristics (ONLY used by RDS)
    vector<Value> solution; ///< remember last solution found
    vector<pair<Double, vector<Value>>> solutions; ///< remember all solutions found
    Cost solutionCost; ///< and its cost
    vector<Constraint*> constrs; ///< list of original cost functions
    int NCBucketSize; ///< number of buckets for NC bucket sort
    vector<VariableList> NCBuckets; ///< NC buckets: vector of backtrackable lists of variables
    Queue NC; ///< NC queue (non backtrackable list)
    Queue IncDec; ///< BAC queue (non backtrackable list)
    Queue AC; ///< AC queue (non backtrackable list)
    Queue DAC; ///< DAC queue (non backtrackable list)
    Queue EAC1; ///< EAC intermediate queue (non backtrackable list)
    Queue EAC2; ///< EAC queue (non backtrackable list)
    Queue Eliminate; ///< Variable Elimination queue (non backtrackable list)
    Queue FEAC; ///< FullEAC queue (non backtrackable list)
    SeparatorList PendingSeparator; ///< List of pending separators for BTD-like methods (backtrackable list)
    Queue DEE; ///< Dead-End Elimination queue (non backtrackable list)
    bool objectiveChanged; ///< flag if lb or ub has changed (NC propagation needs to be done)
    Long nbNodes; ///< current number of calls to propagate method (roughly equal to number of search nodes), used as a time-stamp by Queue methods
    Long nbDEE; ///< number of value removals due to DEE
    Constraint* lastConflictConstr; ///< hook for last conflict variable heuristic
    int maxdomainsize; ///< maximum initial domain size found in all variables
    vector<GlobalConstraint*> globalconstrs; ///< a list of all original global constraints (also inserted in constrs)
    vector<int> delayedNaryCtr; ///< a list of all original nary constraints in extension (also inserted in constrs)
    bool isDelayedNaryCtr; ///< postpone naryctr propagation after all variables have been created
    vector<vector<int>> listofsuccessors; ///< list of topologic order of var used when q variables are  added for decomposing global constraint (berge acyclic)
    StoreInt isPartOfOptimalSolution; ///< true if the current assignment belongs to an optimal solution recorded into bestValues

    // make it private because we don't want copy nor assignment
    WCSP(const WCSP& wcsp);
    WCSP& operator=(const WCSP& wcsp);
    friend class CFNStreamReader;
    friend class VACExtension;

public:
    /// \brief variable elimination information used in backward phase to get a solution during search
    /// \warning restricted to at most two neighbor variables
    typedef struct {
        EnumeratedVariable* x; ///< eliminated variable
        EnumeratedVariable* y; ///< first neighbor variable if any
        EnumeratedVariable* z; ///< second neighbor variable if any
        BinaryConstraint* xy; ///< corresponding binary cost function if any
        BinaryConstraint* xz; ///< corresponding binary cost function if any
        TernaryConstraint* xyz; ///< corresponding ternary cost function if any
        Constraint* ctr; ///< corresponding sum/join cost function if using generic variable elimination
    } elimInfo;

    /// \brief used for preprocesssTernary to make a list of potential ternary constraints with assciated information
    struct TripleVarCostSize {
        EnumeratedVariable* x;
        EnumeratedVariable* y;
        EnumeratedVariable* z;
        float meancost;
        long unsigned int size;

        bool operator<(const TripleVarCostSize& a) const
        {
            return (meancost) > (a.meancost);
        }
    };

    StoreInt elimOrder; ///< current number of eliminated variables
    vector<elimInfo> elimInfos; ///< variable elimination information used in backward phase to get a solution
    StoreInt elimBinOrder; ///< current number of extra binary cost functions consumed in the corresponding pool
    StoreInt elimTernOrder; ///< current number of extra ternary cost functions consumed in the corresponding pool
    vector<Constraint*> elimBinConstrs; ///< pool of (fresh) binary cost functions
    vector<Constraint*> elimTernConstrs; ///< pool of (fresh) ternary cost functions
    int maxDegree; ///< maximum degree of eliminated variables found in preprocessing
    Long elimSpace; ///< estimate of total space required for generic variable elimination

    VACExtension* vac; ///< link to VAC management system

#ifdef XMLFLAG
    map<int, int> varsDom; ///< structures for solution translation: we don't have to parse the XML file again
    vector<vector<int>> Doms; ///< structures for solution translation: we don't have to parse the XML file again
#endif

    WCSP(Cost upperBound, void* solver = NULL);

    virtual ~WCSP();

    // -----------------------------------------------------------
    // General API for weighted CSP global constraint

    int getIndex() const { return instance; } ///< \brief instantiation occurrence number of current WCSP object
    string getName() const { return (name.size()>0)?name:"problem"; }
    void setName(const string& problem) { name = problem; }
    void* getSolver() const { return solver; }

    Cost getLb() const { return lb; } ///< \brief gets problem internal lower bound
    Cost getUb() const { return ub; } ///< \brief gets problem internal upper bound

    Double getDDualBound() const { return Cost2ADCost(lb); } ///< \brief gets problem dual bound as a Double representing a decimal cost (upper resp. lower bound for minimization resp. maximization)
    Double getDPrimalBound() const { return Cost2ADCost(ub); } ///< \brief gets problem primal bound as a Double representing a decimal cost (lower resp. upper bound for minimization resp. maximization)

    Double getDUb() const { return (ToulBar2::costMultiplier < 0 ? Cost2ADCost(lb) : Cost2ADCost(ub)); } ///< \brief gets problem upper bound as a Double representing a decimal cost
    Double getDLb() const { return (ToulBar2::costMultiplier < 0 ? Cost2ADCost(ub) : Cost2ADCost(lb)); } ///< \brief gets problem lower bound as a Double representing a decimal cost

    void setLb(Cost newLb) { lb = newLb; } ///< \internal sets problem lower bound
    void setUb(Cost newUb) { ub = newUb; } ///< \internal sets problem upper bound

    /// \brief sets problem upper bound when a new solution is found
    /// \warning side-effect: adjusts maximum number of buckets (see \ref ncbucket) if called before adding variables to a problem
    void updateUb(Cost newUb)
    {
        if (newUb < ub) {
            ub = newUb;
            if (vars.size() == 0)
                NCBucketSize = cost2log2gub(ub) + 1;
        }
    }

    /// \brief enforces problem upper bound when exploring an alternative search node
    void enforceUb()
    {
        if (CUT((Cost)lb, ub))
            THROWCONTRADICTION;
        objectiveChanged = true;
    }

    /// \brief sets problem upper bound and asks for propagation
    /// \deprecated
    void decreaseUb(Cost newUb)
    {
        if (newUb < ub) {
            if (CUT((Cost)lb, newUb))
                THROWCONTRADICTION;
            ub = newUb;
            objectiveChanged = true;
        }
    }

    /// \brief increases problem lower bound thanks to \e eg soft local consistencies
    /// \param addLb increment value to be \b added to the problem lower bound
    void increaseLb(Cost addLb)
    {
        assert(addLb >= MIN_COST);
        if (addLb > MIN_COST) {
            //		   incWeightedDegree(addLb);
            Cost newLb = lb + addLb;
            if (CUT(newLb, ub))
                THROWCONTRADICTION;
            lb = newLb;
            objectiveChanged = true;
            if (ToulBar2::setminobj)
                (*ToulBar2::setminobj)(getIndex(), -1, newLb, getSolver());
        }
    }

    /// \brief computes the worst-case assignment finite cost (sum of maximum finite cost over all cost functions)
    /// \return the worst-case assignment finite cost plus one
    /// \warning current problem should be completely loaded before calling this function
    Cost finiteUb() const;

    /// \brief updates infinite costs in all cost functions accordingly to the problem global lower and upper bounds
    /// \warning to be used in preprocessing only
    void setInfiniteCost();

    void decreaseLb(Cost cost)
    {
        negCost += cost;
    } ///< \internal manages negative costs in probabilistic inference
    Cost getNegativeLb() const { return negCost; } ///< \internal manages negative costs in probabilistic inference

    bool enumerated(int varIndex) const { return vars[varIndex]->enumerated(); } ///< \brief true if the variable has an enumerated domain

    string getName(int varIndex) const { return vars[varIndex]->getName(); } ///< \note by default, variables names are integers, starting at zero
    unsigned int getVarIndex(const string& s) const { int i = std::distance(vars.begin(), find_if(vars.begin(), vars.end(), [&s](const Variable *var){return (var->getName()==s);})); assert (i >= 0); return static_cast<unsigned>(i); }
    Value getInf(int varIndex) const { return vars[varIndex]->getInf(); } ///< \brief minimum current domain value
    Value getSup(int varIndex) const { return vars[varIndex]->getSup(); } ///< \brief maximum current domain value
    Value getValue(int varIndex) const { return vars[varIndex]->getValue(); } ///< \brief current assigned value \warning undefined if not assigned yet
    unsigned int getDomainSize(int varIndex) const { return vars[varIndex]->getDomainSize(); } ///< \brief current domain size
    vector<Value> getEnumDomain(int varIndex)
    {
        vector<Value> array(getDomainSize(varIndex));
        assert(enumerated(varIndex));
        getEnumDomain(varIndex, array.data());
        return array;
    }
    bool getEnumDomain(int varIndex, Value* array);
    vector<pair<Value, Cost>> getEnumDomainAndCost(int varIndex)
    {
        vector<pair<Value, Cost>> array(getDomainSize(varIndex));
        assert(enumerated(varIndex));
        getEnumDomainAndCost(varIndex, (ValueCost*)array.data());
        return array;
    }
    bool getEnumDomainAndCost(int varIndex, ValueCost* array);
    unsigned int getDomainInitSize(int varIndex) const
    {
        assert(vars[varIndex]->enumerated());
        return ((EnumeratedVariable*)vars[varIndex])->getDomainInitSize();
    } ///< \brief gets initial domain size (warning! assumes EnumeratedVariable)
    Value toValue(int varIndex, unsigned int idx)
    {
        assert(vars[varIndex]->enumerated());
        return ((EnumeratedVariable*)vars[varIndex])->toValue(idx);
    } ///< \brief gets value from index (warning! assumes EnumeratedVariable)
    unsigned int toIndex(int varIndex, Value value)
    {
        assert(vars[varIndex]->enumerated());
        return ((EnumeratedVariable*)vars[varIndex])->toIndex(value);
    } ///< \brief gets index from value (warning! assumes EnumeratedVariable)
    unsigned int toIndex(int varIndex, const string& valueName)
    {
        assert(vars[varIndex]->enumerated());
        return ((EnumeratedVariable*)vars[varIndex])->toIndex(valueName);
    }///< \brief gets index from value name (warning! assumes EnumeratedVariable)
    int getDACOrder(int varIndex) const { return vars[varIndex]->getDACOrder(); } ///< \brief index of the variable in the DAC variable ordering
    void updateCurrentVarsId(); ///< \brief determines the position of each variable in the current list of unassigned variables (see \ref WCSP::dump)

    bool assigned(int varIndex) const { return vars[varIndex]->assigned(); }
    bool unassigned(int varIndex) const { return vars[varIndex]->unassigned(); }
    bool canbe(int varIndex, Value v) const { return vars[varIndex]->canbe(v); }
    bool cannotbe(int varIndex, Value v) const { return vars[varIndex]->cannotbe(v); }
    Value nextValue(int varIndex, Value v) const
    {
        if (enumerated(varIndex)) {
            EnumeratedVariable::iterator iter = ((EnumeratedVariable*)vars[varIndex])->lower_bound(v + 1);
            if (iter != ((EnumeratedVariable*)vars[varIndex])->end())
                return *iter;
            else
                return v;
        } else {
            IntervalVariable::iterator iter = ((IntervalVariable*)vars[varIndex])->lower_bound(v + 1);
            if (iter != ((IntervalVariable*)vars[varIndex])->end())
                return *iter;
            else
                return v;
        }
    }

    void increase(int varIndex, Value newInf) { vars[varIndex]->increase(newInf, true); } ///< \brief changes domain lower bound
    void decrease(int varIndex, Value newSup) { vars[varIndex]->decrease(newSup, true); } ///< \brief changes domain upper bound
    void assign(int varIndex, Value newValue) { vars[varIndex]->assign(newValue, true); } ///< \brief assigns a variable and immediately propagates this assignment
    void remove(int varIndex, Value remValue) { vars[varIndex]->remove(remValue, true); } ///< \brief removes a domain value

    /// \brief assigns a set of variables at once and propagates
    /// \param varIndexes vector of variable indexes as returned by makeXXXVariable
    /// \param newValues vector of values to be assigned to the corresponding variables
    /// \param force boolean if true then apply assignLS even if the variable is already assigned
    /// \note this function is equivalent but faster than a sequence of \ref WCSP::assign. it is particularly useful for Local Search methods such as Large Neighborhood Search.
    void assignLS(vector<int>& varIndexes, vector<Value>& newValues, bool force = false)
    {
        assert(varIndexes.size() == newValues.size());
        unsigned int size = varIndexes.size();
        assignLS((size > 0) ? &varIndexes[0] : NULL, (size > 0) ? &newValues[0] : NULL, size, true, force);
    }

    void assignLS(int* varIndexes, Value* newValues, unsigned int size, bool dopropagate, bool force = false)
    {
        ConstraintSet delayedctrs;
        for (unsigned int i = 0; i < size; i++)
            vars[varIndexes[i]]->assignLS(newValues[i], delayedctrs, force);
        for (ConstraintSet::iterator it = delayedctrs.begin(); it != delayedctrs.end(); ++it)
            if (!(*it)->isGlobal()) {
                if ((*it)->isSep())
                    (*it)->assigns();
                else
                    (*it)->propagate();
            }
        if (dopropagate)
            propagate();
    }

    Cost getUnaryCost(int varIndex, Value v) const { return vars[varIndex]->getCost(v); } ///< \brief unary cost associated to a domain value
    Cost getMaxUnaryCost(int varIndex) const { return vars[varIndex]->getMaxCost(); } ///< \brief maximum unary cost in the domain
    Value getMaxUnaryCostValue(int varIndex) const { return vars[varIndex]->getMaxCostValue(); } ///< \brief a value having the maximum unary cost in the domain
    Value getSupport(int varIndex) const { return vars[varIndex]->getSupport(); } ///< \brief unary (NC/EAC) support value
    Value getBestValue(int varIndex) const { return bestValues[varIndex]; } ///< \brief hint for some value ordering heuristics (ONLY used by RDS)
    void setBestValue(int varIndex, Value v) { bestValues[varIndex] = v; } ///< \brief hint for some value ordering heuristics (ONLY used by RDS)
    bool getIsPartOfOptimalSolution() { return (isPartOfOptimalSolution != 0); } ///< \brief special flag used for debugging purposes only
    void setIsPartOfOptimalSolution(bool v) { isPartOfOptimalSolution = (v ? 1 : 0); } ///< \brief special flag used for debugging purposes only

    int getDegree(int varIndex) const { return vars[varIndex]->getDegree(); } ///< \brief approximate degree of a variable (\e ie number of active cost functions, see \ref varelim)
    int getTrueDegree(int varIndex) const { return vars[varIndex]->getTrueDegree(); } ///< \brief degree of a variable
    Long getWeightedDegree(int varIndex) const { return vars[varIndex]->getWeightedDegree(); } ///< \brief weighted degree heuristic
    void resetWeightedDegree(int varIndex) { vars[varIndex]->resetWeightedDegree(); } ///< \brief initialize weighted degree heuristic
    void revise(Constraint* c) { lastConflictConstr = c; } ///< \internal last conflict heuristic
    /// \internal last conflict heuristic
    void conflict()
    {
        if (lastConflictConstr) {
            if (ToulBar2::verbose >= 2)
                cout << "Last conflict on " << *lastConflictConstr << endl;
            lastConflictConstr->incConflictWeight(lastConflictConstr);
            lastConflictConstr = NULL;
        }
    }
    /// \internal \deprecated
    void incWeightedDegree(Long incval)
    {
        if (lastConflictConstr) {
            lastConflictConstr->incConflictWeight(incval);
        }
    }

    //  set<int> lastConflictSet;
    //  set<int> getLastConflicts(){
    //      return lastConflictSet;
    //  };
    //  void registerConflicts(){
    //      lastConflictSet=getConflictVars();
    //  };
    //  set<int> getConflictVars(){
    //      set<int> conflictvar;
    //      for(vector<Constraint*>::iterator it = constrs.begin();it!=constrs.end();++it)
    //      {
    //          if((*it)->getCost() > MIN_COST) // Warning!!! It should test initial costs (or without propagation)
    //          {
    //              TSCOPE s;
    //              (*it)->getScope(s);
    //              for(TSCOPE::iterator it2=s.begin(); it2!= s.end();++it2)
    //                  conflictvar.insert((*it2).first);
    //          }
    //      }
    //      return conflictvar;
    //  }

    void whenContradiction(); ///< \brief after a contradiction, resets propagation queues and increases \ref WCSP::nbNodes
    void propagate(); ///< \brief propagates until a fix point is reached (or throws a contradiction) and then increases \ref WCSP::nbNodes
    bool verify(); ///< \brief checks the propagation fix point is reached \warning might change EAC supports

    unsigned int numberOfVariables() const { return vars.size(); } ///< \brief current number of created variables
    /// \brief returns current number of unassigned variables
    unsigned int numberOfUnassignedVariables() const
    {
        int res = 0;
        for (unsigned int i = 0; i < vars.size(); i++)
            if (unassigned(i))
                res++;
        return res;
    }
    unsigned int numberOfConstraints() const { return constrs.size(); } ///< \brief initial number of cost functions
    unsigned int numberOfConnectedConstraints() const; ///< \brief current number of cost functions
    unsigned int numberOfConnectedBinaryConstraints() const; ///< \brief current number of binary cost functions
    unsigned int medianDomainSize() const; ///< \brief median current domain size of variables
    unsigned int medianDegree() const; ///< \brief median current degree of variables
    unsigned int medianArity() const; ///< \brief median arity of current cost functions
    int getMaxDomainSize() const { return maxdomainsize; } ///< \brief maximum initial domain size found in all variables
    int getMaxCurrentDomainSize() const; ///< \brief maximum current domain size found in all variables
    unsigned int getDomainSizeSum() const; ///< \brief total sum of current domain sizes
    /// \brief Cartesian product of current domain sizes
    /// \param cartesianProduct result obtained by the GNU Multiple Precision Arithmetic Library GMP
    void cartProd(BigInteger& cartesianProduct)
    {
        for (vector<Variable*>::iterator it = vars.begin(); it != vars.end(); it++) {
            Variable* x = *it;
            mpz_mul_si(cartesianProduct.integer, cartesianProduct.integer, x->getDomainSize());
        }
    }
#ifdef BOOST
    int diameter();
    int connectedComponents();
    int biConnectedComponents();
    void minimumDegreeOrderingBGL(vector<int>& order);
    void spanningTreeOrderingBGL(vector<int>& order);
    void reverseCuthillMcKeeOrderingBGL(vector<int>& order);
    void maximumCardinalitySearch(vector<int>& order);
    void minimumFillInOrdering(vector<int>& order);
    void minimumDegreeOrdering(vector<int>& order);
#endif

    int makeEnumeratedVariable(string n, Value iinf, Value isup);
    int makeEnumeratedVariable(string n, Value* d, int dsize);
    void addValueName(int xIndex, const string& name);
    int makeIntervalVariable(string n, Value iinf, Value isup);

    void postNullaryConstraint(Double cost);
    void postNullaryConstraint(Cost cost);
    void postUnary(int xIndex, vector<Cost>& costs);
    int postUnary(int xIndex, Value* d, int dsize, Cost penalty);
    void postUnaryConstraint(int xIndex, vector<Double>& costs, bool incremental = false);
    void postUnaryConstraint(int xIndex, vector<Cost>& costs) { postUnary(xIndex, costs); }
    void postIncrementalUnaryConstraint(int xIndex, vector<Cost>& costs) { postUnary(xIndex, costs); }
    int postUnaryConstraint(int xIndex, Value* d, int dsize, Cost penalty) { return postUnary(xIndex, d, dsize, penalty); }
    int postSupxyc(int xIndex, int yIndex, Value cst, Value deltamax = MAX_VAL - MIN_VAL);
    int postDisjunction(int xIndex, int yIndex, Value cstx, Value csty, Cost penalty);
    int postSpecialDisjunction(int xIndex, int yIndex, Value cstx, Value csty, Value xinfty, Value yinfty, Cost costx, Cost costy);
    int postBinaryConstraint(int xIndex, int yIndex, vector<Double>& costs, bool incremental = false);
    int postBinaryConstraint(int xIndex, int yIndex, vector<Cost>& costs);
    int postIncrementalBinaryConstraint(int xIndex, int yIndex, vector<Cost>& costs);
    int postTernaryConstraint(int xIndex, int yIndex, int zIndex, vector<Double>& costs, bool incremental = false);
    int postTernaryConstraint(int xIndex, int yIndex, int zIndex, vector<Cost>& costs);
    int postIncrementalTernaryConstraint(int xIndex, int yIndex, int zIndex, vector<Cost>& costs);
    int postNaryConstraintBegin(vector<int>& scope, Cost defval, Long nbtuples = 0, bool forcenary = false) { return postNaryConstraintBegin(scope.data(), scope.size(), defval, nbtuples, forcenary); }
    int postNaryConstraintBegin(int* scopeIndex, int arity, Cost defval, Long nbtuples = 0, bool forcenary = false); /// \warning must call postNaryConstraintEnd after giving cost tuples ; \warning it may create a WeightedClause instead of NaryConstraint
    void postNaryConstraintTuple(int ctrindex, vector<Value>& tuple, Cost cost) { postNaryConstraintTuple(ctrindex, tuple.data(), tuple.size(), cost); }
    void postNaryConstraintTuple(int ctrindex, Value* tuple, int arity, Cost cost);
    void postNaryConstraintTuple(int ctrindex, const Tuple& tuple, Cost cost);
    void postNaryConstraintEnd(int ctrindex);

    // -----------------------------------------------------------
    // Methods for diverse solutions
    // -----------------------------------------------------------

    void addDivConstraint(const vector<Value> solution, int sol_id, Cost cost); // to look for the (j+1)-th solution, with j = sol_id
    void addHDivConstraint(const vector<Value> solution, int sol_id, Cost cost);
    void addTDivConstraint(const vector<Value> solution, int sol_id, Cost cost);
    void addMDDConstraint(Mdd mdd, int relaxed);
    void addHMDDConstraint(Mdd mdd, int relaxed);
    void addTMDDConstraint(Mdd mdd, int relaxed);

    const vector<Variable*>& getDivVariables()
    {
        return divVariables;
    }

    int postCliqueConstraint(vector<int>& scope, const string& arguments)
    {
        istringstream file(arguments);
        return postCliqueConstraint(scope.data(), scope.size(), file);
    }
    int postCliqueConstraint(int* scopeIndex, int arity, istream& file);

    int postKnapsackConstraint(vector<int>& scope, const string& arguments)
    {
        istringstream file(arguments);
        return postKnapsackConstraint(scope.data(), scope.size(), file);
    }
    int postKnapsackConstraint(int* scopeIndex, int arity, istream& file);
    int postGlobalConstraint(int* scopeIndex, int arity, const string& gcname, istream& file, int* constrcounter = NULL, bool mult = true); ///< \deprecated should use WCSP::postGlobalCostFunction instead \warning does not work for arity below 4 (use binary or ternary cost functions instead)

    GlobalConstraint* postGlobalCostFunction(int* scopeIndex, int arity, const string& name, int* constrcounter = NULL);

    int postWAmong(vector<int>& scope, const string& semantics, const string& propagator, Cost baseCost, const vector<Value>& values, int lb, int ub) { return postWAmong(scope.data(), scope.size(), semantics, propagator, baseCost, values, lb, ub); } ///< \brief post a soft among cost function
    int postWAmong(int* scopeIndex, int arity, const string& semantics, const string& propagator, Cost baseCost, const vector<Value>& values, int lb, int ub); ///< \deprecated
    void postWAmong(int* scopeIndex, int arity, string semantics, Cost baseCost, Value* values, int nbValues, int lb, int ub); ///< \deprecated post a weighted among cost function decomposed as a cost function network
    void postWVarAmong(vector<int>& scope, const string& semantics, Cost baseCost, vector<Value>& values, int varIndex) { postWVarAmong(scope.data(), scope.size(), semantics, baseCost, values.data(), values.size(), varIndex); } ///< \brief post a weighted among cost function with the number of values encoded as a variable with index \a varIndex (\e network-based propagator only)
    void postWVarAmong(int* scopeIndex, int arity, const string& semantics, Cost baseCost, Value* values, int nbValues, int varIndex); ///< \deprecated
    int postWRegular(vector<int>& scope, const string& semantics, const string& propagator, Cost baseCost,
        int nbStates,
        const vector<WeightedObjInt>& initial_States,
        const vector<WeightedObjInt>& accepting_States,
        const vector<DFATransition>& Wtransitions) { return postWRegular(scope.data(), scope.size(), semantics, propagator, baseCost, nbStates, initial_States, accepting_States, Wtransitions); } ///< \brief post a soft or weighted regular cost function
    int postWRegular(int* scopeIndex, int arity, const string& semantics, const string& propagator, Cost baseCost,
        int nbStates,
        const vector<WeightedObjInt>& initial_States,
        const vector<WeightedObjInt>& accepting_States,
        const vector<DFATransition>& Wtransitions); ///< \deprecated
    void postWRegular(int* scopeIndex, int arity, int nbStates, vector<pair<int, Cost>> initial_States, vector<pair<int, Cost>> accepting_States, int** Wtransitions, vector<Cost> transitionsCosts); ///< \deprecated post a weighted regular cost function decomposed as a cost function network
    int postWAllDiff(int* scopeIndex, int arity, const string& semantics, const string& propagator, Cost baseCost); ///< \brief post a soft alldifferent cost function
    void postWAllDiff(int* scopeIndex, int arity, string semantics, Cost baseCost); ///< \deprecated post a soft alldifferent cost function decomposed as a cost function network
    int postWGcc(int* scopeIndex, int arity, const string& semantics, const string& propagator, Cost baseCost,
        const vector<BoundedObjValue>& values); ///< \brief post a soft global cardinality cost function
    void postWGcc(int* scopeIndex, int arity, string semantics, Cost baseCost, Value* values, int nbValues, int* lb, int* ub); ///< \deprecated post a soft global cardinality cost function decomposed as a cost function network
    int postWSame(int* scopeIndexG1, int arityG1, int* scopeIndexG2, int arityG2, const string& semantics, const string& propagator, Cost baseCost); ///< \brief post a soft same cost function (a group of variables being a permutation of another group with the same size)
    void postWSame(int* scopeIndex, int arity, string semantics, Cost baseCost); ///< \deprecated post a soft same cost function
    void postWSameGcc(int* scopeIndex, int arity, string semantics, Cost baseCost, Value* values, int nbValues, int* lb, int* ub); ///< \brief post a combination of a same and gcc cost function decomposed as a cost function network
    int postWGrammarCNF(int* scopeIndex, int arity, const string& semantics, const string& propagator, Cost baseCost,
        int nbSymbols,
        int startSymbol,
        const vector<CFGProductionRule> WRuleToTerminal); ///< \brief post a soft/weighted grammar cost function with the dynamic programming propagator and grammar in Chomsky normal form
    int postMST(int* scopeIndex, int arity, const string& semantics, const string& propagator, Cost baseCost); ///< \brief post a Spanning Tree hard constraint
    int postMaxWeight(int* scopeIndex, int arity, const string& semantics, const string& propagator, Cost baseCost,
        const vector<WeightedVarValPair> weightFunction); ///< \brief post a weighted max cost function (maximum cost of a set of unary cost functions associated to a set of variables)
    void postWSum(int* scopeIndex, int arity, string semantics, Cost baseCost, string comparator, int rightRes); ///< \brief post a soft linear constraint with unit coefficients
    void postWVarSum(int* scopeIndex, int arity, string semantics, Cost baseCost, string comparator, int varIndex); ///< \brief post a soft linear constraint with unit coefficients and variable right-hand side
    void postWOverlap(int* scopeIndex, int arity, string semantics, Cost baseCost, string comparator, int rightRes); /// \brief post a soft overlap cost function (a group of variables being point-wise equivalent -- and not equal to zero -- to another group with the same size)

    bool isGlobal() { return (globalconstrs.size() > 0); } ///< \brief true if there are soft global cost functions defined in the problem

    Cost read_wcsp(const char* fileName); ///< \brief load problem in any of all formats managed by tb2. Return the global UB known to the solver at start (file and command line).
    void read_uai2008(const char* fileName); ///< \brief load problem in UAI 2008 format (see http://graphmod.ics.uci.edu/uai08/FileFormat and http://www.cs.huji.ac.il/project/UAI10/fileFormat.php) \warning UAI10 evidence file format not recognized by toulbar2 as it does not allow multiple evidence (you should remove the first value in the file)
    void read_random(int n, int m, vector<int>& p, int seed, bool forceSubModular = false, string globalname = ""); ///< \brief create a random WCSP with \e n variables, domain size \e m, array \e p where the first element is a percentage of tuples with a nonzero cost and next elements are the number of random cost functions for each different arity (starting with arity two), random seed, a flag to have a percentage (last element in the array \e p) of the binary cost functions being permutated submodular, and a string to use a specific global cost function instead of random cost functions in extension
    void read_wcnf(const char* fileName); ///< \brief load problem in (w)cnf format (see http://www.maxsat.udl.cat/08/index.php?disp=requirements)
    void read_qpbo(const char* fileName); ///< \brief load quadratic pseudo-Boolean optimization problem in unconstrained quadratic programming text format (first text line with n, number of variables and m, number of triplets, followed by the m triplets (x,y,cost) describing the sparse symmetric nXn cost matrix with variable indexes such that x <= y and any positive or negative real numbers for costs)
    void read_opb(const char* fileName); ///< \brief load pseudo-Boolean optimization problem
    void read_legacy(const char* fileName); ///< \brief common ending section for all readers

    void read_XML(const char* fileName); ///< \brief load problem in XML format (see http://www.cril.univ-artois.fr/~lecoutre/benchmarks.html)
    void solution_XML(bool opt = false); ///< \brief output solution in Max-CSP 2008 output format
    void solution_UAI(Cost res); ///< \brief output solution in UAI 2008 output format

    const vector<Value> getSolution() { return solution; }
    Double getSolutionValue() const { return Cost2ADCost(solutionCost); }
    Cost getSolutionCost() const { return solutionCost; }
    const vector<Value> getSolution(Cost* cost_ptr)
    {
        if (cost_ptr != NULL)
            *cost_ptr = solutionCost;
        return solution;
    }
    vector<pair<Double, vector<Value>>> getSolutions() const { return solutions; }
    void initSolutionCost() { solutionCost = MAX_COST; }
    void setSolution(Cost cost, TAssign* sol = NULL)
    {
        solutionCost = cost;
        for (unsigned int i = 0; i < numberOfVariables(); i++) {
            Value v = ((sol != NULL) ? (*sol)[i] : getValue(i));
            if (!ToulBar2::verifyOpt && ToulBar2::solutionBasedPhaseSaving)
                setBestValue(i, v);
            solution[i] = ((ToulBar2::sortDomains && ToulBar2::sortedDomains.find(i) != ToulBar2::sortedDomains.end()) ? ToulBar2::sortedDomains[i][toIndex(i, v)].value : v);
        }
        solutions.push_back(make_pair(Cost2ADCost(solutionCost), solution));
    }
    void printSolution()
    {
        for (unsigned int i = 0; i < numberOfVariables(); i++) {
            if (enumerated(i) && ((EnumeratedVariable*)getVar(i))->isValueNames()) {
                EnumeratedVariable* myvar = (EnumeratedVariable*)getVar(i);
                Value myvalue = solution[i];
                string valuelabel = myvar->getValueName(myvar->toIndex(myvalue));
                string varlabel = myvar->getName();

                switch (ToulBar2::showSolutions) {
                case 1:
                    cout << myvalue;
                    break;
                case 2:
                    cout << valuelabel;
                    break;
                case 3:
                    cout << varlabel << "=" << valuelabel;
                    break;
                default:
                    break;
                }
            } else {
                cout << solution[i];
            }
            cout << (i < numberOfVariables() - 1 ? " " : "");
        }
    }
    void printSolution(ostream& os)
    {
        for (unsigned int i = 0; i < numberOfVariables(); i++) {
            if (enumerated(i) && ((EnumeratedVariable*)getVar(i))->isValueNames()) {
                EnumeratedVariable* myvar = (EnumeratedVariable*)getVar(i);
                Value myvalue = solution[i];
                string valuelabel = myvar->getValueName(myvar->toIndex(myvalue));
                string varlabel = myvar->getName();

                switch (ToulBar2::writeSolution) {
                case 1:
                    os << myvalue;
                    break;
                case 2:
                    os << valuelabel;
                    break;
                case 3:
                    os << varlabel << "=" << valuelabel;
                    break;
                default:
                    break;
                }
            } else {
                os << solution[i];
            }
            os << (i < numberOfVariables() - 1 ? " " : "");
        }
    }
    void printSolution(FILE* f)
    {
        for (unsigned int i = 0; i < numberOfVariables(); i++) {
            if (enumerated(i) && ((EnumeratedVariable*)getVar(i))->isValueNames()) {
                EnumeratedVariable* myvar = (EnumeratedVariable*)getVar(i);
                Value myvalue = solution[i];
                string valuelabel = myvar->getValueName(myvar->toIndex(myvalue));
                string varlabel = myvar->getName();

                switch (ToulBar2::writeSolution) {
                case 1:
                    fprintf(f, "%d", myvalue);
                    break;
                case 2:
                    fprintf(f, "%s", valuelabel.c_str());
                    break;
                case 3:
                    fprintf(f, "%s=%s", varlabel.c_str(), valuelabel.c_str());
                    break;
                default:
                    break;
                }
            } else {
                fprintf(f, "%d", solution[i]);
            }
            if (i < numberOfVariables() - 1)
                fprintf(f, " ");
        }
    }
    void printSolutionMaxSAT(ostream& os)
    {
        os << "v";
        for (unsigned int i = 0; i < numberOfVariables(); i++) {
            os << " " << ((solution[i]) ? ((int)i + 1) : -((int)i + 1));
        };
        os << endl;
    }

    void print(ostream& os); ///< \brief print current domains and active cost functions (see \ref verbosity)
    void dump(ostream& os, bool original = true); ///< \brief output the current WCSP into a file in wcsp format \param os output file \param original if true then keeps all variables with their original domain size else uses unassigned variables and current domains recoding variable indexes
    void dump_CFN(ostream& os, bool original = true); ///< \brief output the current WCSP into a file in CFN format \param os output file \param original if true then keeps all variables with their original domain size else uses unassigned variables and current domains recoding variable indexes
    friend ostream& operator<<(ostream& os, WCSP& wcsp); ///< \relates WCSP::print

    // -----------------------------------------------------------
    // Specific API for Variable and Constraint classes

    Variable* getVar(int varIndex) const { return vars[varIndex]; }
    vector<vector<int>>* getListSuccessors() { return &listofsuccessors; }
    Constraint* getCtr(int ctrIndex) const
    {
        if (ctrIndex >= 0) {
            return constrs[ctrIndex];
        } else {
            if (-ctrIndex - 1 >= MAX_ELIM_BIN) {
                return elimTernConstrs[-ctrIndex - 1 - MAX_ELIM_BIN];
            } else {
                return elimBinConstrs[-ctrIndex - 1];
            }
        }
    }

    void link(Variable* x)
    {
        vars.push_back(x);
        bestValues.push_back(x->getSup() + 1);
        solution.push_back(x->getSup() + 1);
    }
    void link(Constraint* c) { constrs.push_back(c); }

    VariableList* getNCBucket(int ibucket) { return &NCBuckets[ibucket]; }
    int getNCBucketSize() const { return NCBucketSize; }
    void changeNCBucket(int oldBucket, int newBucket, DLink<Variable*>* elt)
    {
        assert(newBucket < NCBucketSize);
        if (oldBucket >= 0)
            NCBuckets[oldBucket].erase(elt, true);
        if (newBucket >= 0)
            NCBuckets[newBucket].push_back(elt, true);
    }
    void printNCBuckets();

    Long getNbNodes() const { return nbNodes; }
    Long getNbDEE() const { return nbDEE; }
    void incNbDEE(Long v = 1LL) { nbDEE += v; }

    void queueNC(DLink<VariableWithTimeStamp>* link) { NC.push(link, nbNodes); }
    void queueInc(DLink<VariableWithTimeStamp>* link) { IncDec.push(link, INCREASE_EVENT, nbNodes); }
    void queueDec(DLink<VariableWithTimeStamp>* link) { IncDec.push(link, DECREASE_EVENT, nbNodes); }
    void queueAC(DLink<VariableWithTimeStamp>* link) { AC.push(link, nbNodes); }
    void queueDAC(DLink<VariableWithTimeStamp>* link) { DAC.push(link, nbNodes); }
    void queueEAC1(DLink<VariableWithTimeStamp>* link) { EAC1.push(link, nbNodes); }
    void queueEAC2(DLink<VariableWithTimeStamp>* link) { EAC2.push(link, nbNodes); }
    void queueEliminate(DLink<VariableWithTimeStamp>* link) { Eliminate.push(link, nbNodes); }
    void queueSeparator(DLink<Separator*>* link) { PendingSeparator.push_back(link, true); }
    void unqueueSeparator(DLink<Separator*>* link) { PendingSeparator.erase(link, true); }
    void queueDEE(DLink<VariableWithTimeStamp>* link) { DEE.push(link, nbNodes); }
    void queueFEAC(DLink<VariableWithTimeStamp>* link) { FEAC.push(link, nbNodes); }

    void propagateNC(); ///< \brief removes forbidden values
    void propagateIncDec(); ///< \brief ensures unary bound arc consistency supports (remove forbidden domain bounds)
    void propagateAC(); ///< \brief ensures unary and binary and ternary arc consistency supports
    void propagateDAC(); ///< \brief ensures unary and binary and ternary directed arc consistency supports
    void propagateTRWS(); ///< \brief iterates TRW-S until convergence
    void fillEAC2();
    Queue* getQueueEAC1() { return &EAC1; }
    void propagateEAC(); ///< \brief ensures unary existential arc consistency supports
    void propagateSeparator(); ///< \brief exploits graph-based learning
    void propagateDEE(); ///< \brief removes dominated values (dead-end elimination and possibly soft neighborhood substitutability)
    void propagateFEAC(); ///< \brief seek if new EAC support is also FullEAC support (i.e., compatible with all its EAC value neighbors)

    /// \brief sorts the list of constraints associated to each variable based on smallest problem variable indexes
    /// \warning side-effect: updates DAC order according to an existing variable elimination order
    void sortConstraints();

    /// \brief applies preprocessing techniques before the search (depending on toublar2 options)
    /// \warning to be done \b only before the search
    void preprocessing();

    // -----------------------------------------------------------
    // Methods for Variable Elimination

    void initElimConstr();
    void initElimConstrs();

    int getElimOrder() { return (int)elimOrder; }
    int getElimBinOrder() { return (int)elimBinOrder; }
    int getElimTernOrder() { return (int)elimTernOrder; }
    void elimOrderInc() { elimOrder = elimOrder + 1; }
    void elimBinOrderInc() { elimBinOrder = elimBinOrder + 1; }
    void elimTernOrderInc() { elimTernOrder = elimTernOrder + 1; }
    Constraint* getElimBinCtr(int elimBinIndex) const { return elimBinConstrs[elimBinIndex]; }
    Constraint* getElimTernCtr(int elimTernIndex) const { return elimTernConstrs[elimTernIndex]; }

    BinaryConstraint* newBinaryConstr(EnumeratedVariable* x, EnumeratedVariable* y, Constraint* from1 = NULL, Constraint* from2 = NULL);
    BinaryConstraint* newBinaryConstr(EnumeratedVariable* x, EnumeratedVariable* y, vector<Cost>& costs);
    TernaryConstraint* newTernaryConstr(EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z, Constraint* from1 = NULL);
    TernaryConstraint* newTernaryConstr(EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z, vector<Cost>& costs);

    void eliminate();
    void restoreSolution(Cluster* c = NULL);

    Constraint* sum(Constraint* ctr1, Constraint* ctr2);
    void project(Constraint*& ctr_inout, EnumeratedVariable* var, Constraint* ctr_copy = NULL);
    void variableElimination(EnumeratedVariable* var);

    void processTernary(); ///< \brief projects&subtracts ternary cost functions (see \ref preprocessing)
    void ternaryCompletion();
    bool kconsistency(int xIndex, int yIndex, int zIndex, BinaryConstraint* xy, BinaryConstraint* yz, BinaryConstraint* xz);

    // -----------------------------------------------------------
    // Data and methods for Virtual Arc Consistency

    void histogram(Cost c); /// \brief initializes histogram of costs used by Virtual Arc Consistency to speed up its convergence (Bool\f$_\theta\f$ of P)
    void iniSingleton();
    void updateSingleton();
    void removeSingleton();
    void printVACStat();

    // -----------------------------------------------------------
    // Data and methods for Cluster Tree Decomposition

    TreeDecomposition* td;
    TreeDecomposition* getTreeDec() { return td; }
    static bool isAlreadyTreeDec(char* filename); ///< \brief finds if the given file is a variable ordering or a tree decomposition
    void buildTreeDecomposition();
    void elimOrderFile2Vector(char* elimVarOrderFilename, vector<int>& elimVarOrder); ///< \brief returns a reverse topological order from a variable elimination order
    void treeDecFile2Vector(char* treeDecFilename, vector<int>& elimVarOrder); ///< \brief returns a reverse topological order from a tree decomposition
    void setDACOrder(vector<int>& elimVarOrder); ///< \brief change DAC order and propagate from scratch

    // dac order reordering when Berge acyclic gobal constraint are present in the wcsp
    //
    void visit(int i, vector<int>& revdac, vector<bool>& marked, const vector<vector<int>>& listofsuccessors);

    // -----------------------------------------------------------
    // Functions dealing with all representations of Costs
    // warning: ToulBar2::NormFactor has to be initialized

    Cost decimalToCost(const string& decimalToken, const unsigned int lineNumber) const;
    Cost DoubletoCost(const Double& c) const { return Round(c * powl(10.0, ToulBar2::decimalPoint)) + negCost; }
    Double Cost2ADCost(const Cost& c) const { return Cost2RDCost(c - negCost); } // Absolute costs
    Double Cost2RDCost(const Cost& c) const { return ((Double)(c) / Exp10(ToulBar2::decimalPoint) / ToulBar2::costMultiplier); } //Relative costs
    Cost Prob2Cost(TProb p) const;
    TProb Cost2Prob(Cost c) const;
    TLogProb Cost2LogProb(Cost c) const;
    Cost LogProb2Cost(TLogProb p) const;
    Cost LogSumExp(Cost c1, Cost c2) const;
    TLogProb LogSumExp(TLogProb logc1, Cost c2) const;
    TLogProb LogSumExp(TLogProb logc1, TLogProb logc2) const;
};

#endif /*TB2WCSP_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
