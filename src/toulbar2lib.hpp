/** \file toulbar2lib.hpp
 *  \brief Main protocol class of a global soft constraint representing a weighted CSP and a generic WCSP complete tree-search-based solver
 *
<pre>
    Copyright (c) 2006-2022, toulbar2 team

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

    toulbar2 is currently maintained by Simon de Givry, INRAE - MIAT, Toulouse, France (simon.de-givry@inrae.fr)
</pre>
 */

/* Warning in case of 'mainpage' modifications :
 * The text of 'mainpage' Doxygen keyword is not automatically generated
 * into Sphinx documentation, where it has been manually copied.
 */

/*! @mainpage

    <table>
        <tr><th>Cost Function Network Solver     <td>toulbar2
        <tr><th>Copyright               <td>toulbar2 team
        <tr><th>Source                  <td>https://github.com/toulbar2/toulbar2
    </table>

    See the @link md_README.html README @endlink for more details.

    toulbar2 can be used as a stand-alone solver reading various problem file formats (wcsp, uai, wcnf, qpbo) or as a C++ library.\n
    This document describes the wcsp native file format and the toulbar2 C++ library API.
    \note Use cmake flags LIBTB2=ON and TOULBAR2_ONLY=OFF to get the toulbar2 C++ library libtb2.so and toulbar2test executable example.
    \see ./src/toulbar2test.cpp

    \defgroup wcspformat
    \defgroup modeling
    \defgroup solving
    \defgroup verbosity
    \defgroup preprocessing
    \defgroup heuristics
    \defgroup softac
    \defgroup VAC
    \defgroup ncbucket
    \defgroup varelim
    \defgroup propagation
    \defgroup backtrack

*/

#ifndef TOULBAR2LIB_HPP_
#define TOULBAR2LIB_HPP_

#include "core/tb2types.hpp"

/** Abstract class WeightedCSP representing a weighted constraint satisfaction problem
 *	- problem lower and upper bounds
 *	- list of variables with their finite domains (either represented by an enumerated list of values, or by a single interval)
 *	- list of cost functions (created before and during search by variable elimination of variables with small degree)
 *	- local consistency propagation (variable-based propagation) including cluster tree decomposition caching (separator-based cache)
 *
 * \note Variables are referenced by their lexicographic index number (as returned by \e eg WeightedCSP::makeEnumeratedVariable)
 * \note Cost functions are referenced by their lexicographic index number (as returned by \e eg WeightedCSP::postBinaryConstraint)
 *
 */

class WeightedCSP {
public:
    static WeightedCSP* makeWeightedCSP(Cost upperBound, void* solver = NULL); ///< \brief Weighted CSP factory

    virtual ~WeightedCSP() {}

    virtual int getIndex() const = 0; ///< \brief instantiation occurrence number of current WCSP object
    virtual string getName() const = 0; ///< \brief get WCSP problem name (defaults to filename with no extension)
    virtual void setName(const string& problem) = 0; ///< \brief set WCSP problem name
    virtual void* getSolver() const = 0; ///< \brief special hook to access solver information

    virtual Cost getLb() const = 0; ///< \brief gets internal dual lower bound
    virtual Cost getUb() const = 0; ///< \brief gets internal primal upper bound

    virtual Double getDPrimalBound() const = 0; ///< \brief gets problem primal bound as a Double representing a decimal cost (upper resp. lower bound for minimization resp. maximization)
    virtual Double getDDualBound() const = 0; ///< \brief gets problem dual bound as a Double representing a decimal cost (lower resp. upper bound for minimization resp. maximization)

    virtual Double getDLb() const = 0; ///< \brief gets problem lower bound as a Double representing a decimal cost
    virtual Double getDUb() const = 0; ///< \brief gets problem upper bound as a Double representing a decimal cost

    /// \brief sets problem upper bound as a Double representing a decimal cost
    virtual void updateDUb(Double newDUb) = 0;

    /// \brief sets initial problem upper bound and each time a new solution is found
    virtual void updateUb(Cost newUb) = 0;
    /// \brief enforces problem upper bound when exploring an alternative search node
    virtual void enforceUb() = 0;
    /// \brief increases problem lower bound thanks to \e eg soft local consistencies
    /// \param addLb increment value to be \b added to the problem lower bound
    virtual void increaseLb(Cost addLb) = 0;
    /// \brief shift problem optimum toward negative costs
    /// \param shift positive shifting value to be subtracted to the problem optimum when printing the solutions
    virtual void decreaseLb(Cost shift) = 0;
    virtual Cost getNegativeLb() const = 0; ///< \brief gets constant term used to subtract to the problem optimum when printing the solutions
    /// \brief computes the worst-case assignment finite cost (sum of maximum finite cost over all cost functions plus one)
    /// \return the worst-case assignment finite cost
    /// \warning current problem should be completely loaded and propagated before calling this function
    virtual Cost finiteUb() const = 0;
    /// \brief updates infinite costs in all cost functions accordingly to the problem global lower and upper bounds
    /// \warning to be used in preprocessing only
    virtual void setInfiniteCost() = 0;
    /// \brief returns true if any complete assignment using current domains is a valid tuple with finite cost (i.e., cost strictly less than the problem upper bound minus the lower bound)
    virtual bool isfinite() const = 0;

    virtual bool enumerated(int varIndex) const = 0; ///< \brief true if the variable has an enumerated domain

    virtual string getName(int varIndex) const = 0; ///< \note by default, variables names are integers, starting at zero
    virtual unsigned int getVarIndex(const string& s) const = 0; ///< return variable index from its name, or numberOfVariables() if not found
    virtual Value getInf(int varIndex) const = 0; ///< \brief minimum current domain value
    virtual Value getSup(int varIndex) const = 0; ///< \brief maximum current domain value
    virtual Value getValue(int varIndex) const = 0; ///< \brief current assigned value \warning undefined if not assigned yet
    virtual unsigned int getDomainSize(int varIndex) const = 0; ///< \brief current domain size
    virtual vector<Value> getEnumDomain(int varIndex) = 0; ///< \brief gets current domain values in an array
    virtual bool getEnumDomain(int varIndex, Value* array) = 0; ///< \deprecated
    virtual vector<pair<Value, Cost>> getEnumDomainAndCost(int varIndex) = 0; ///< \brief gets current domain values and unary costs in an array
    virtual bool getEnumDomainAndCost(int varIndex, ValueCost* array) = 0; ///< \deprecated
    virtual unsigned int getDomainInitSize(int varIndex) const = 0; ///< \brief gets initial domain size (warning! assumes EnumeratedVariable)
    virtual Value toValue(int varIndex, unsigned int idx) = 0; ///< \brief gets value from index (warning! assumes EnumeratedVariable)
    virtual unsigned int toIndex(int varIndex, Value value) = 0; ///< \brief gets index from value (warning! assumes EnumeratedVariable)
    virtual unsigned int toIndex(int varIndex, const string& valueName) = 0; ///< \brief gets index from value name (warning! assumes EnumeratedVariable with value names)
    virtual int getDACOrder(int varIndex) const = 0; ///< \brief index of the variable in the DAC variable ordering

    virtual bool assigned(int varIndex) const = 0;
    virtual bool unassigned(int varIndex) const = 0;
    virtual bool canbe(int varIndex, Value v) const = 0;
    virtual bool cannotbe(int varIndex, Value v) const = 0;
    virtual Value nextValue(int varIndex, Value v) const = 0; ///< \brief first value after v in the current domain or v if there is no value

    virtual void increase(int varIndex, Value newInf) = 0; ///< \brief changes domain lower bound
    virtual void decrease(int varIndex, Value newSup) = 0; ///< \brief changes domain upper bound
    virtual void assign(int varIndex, Value newValue) = 0; ///< \brief assigns a variable and immediately propagates this assignment
    virtual void remove(int varIndex, Value remValue) = 0; ///< \brief removes a domain value (valid if done for an enumerated variable or on its domain bounds)

    /// \brief assigns a set of variables at once and propagates (used by Local Search methods such as Large Neighborhood Search)
    /// \param varIndexes vector of variable indexes as returned by makeXXXVariable
    /// \param newValues vector of values to be assigned to the corresponding variables
    /// \param force boolean if true then apply assignLS even if the variable is already assigned
    /// Note this function is equivalent but faster than a sequence of assign.
    virtual void assignLS(vector<int>& varIndexes, vector<Value>& newValues, bool force = false) = 0;
    virtual void assignLS(int* varIndexes, Value* newValues, unsigned int size, bool dopropagate, bool force = false) = 0;

    /// \brief deconnects a set of variables from the rest of the problem and assigns them to their support value (used by Incremental Search)
    /// \param varIndexes vector of variable indexes as returned by makeXXXVariable
    virtual void deconnect(vector<int>& varIndexes) = 0;

    virtual Cost getUnaryCost(int varIndex, Value v) const = 0; ///< \brief unary cost associated to a domain value
    virtual Cost getMaxUnaryCost(int varIndex) const = 0; ///< \brief maximum unary cost in the domain
    virtual Value getMaxUnaryCostValue(int varIndex) const = 0; ///< \brief a value having the maximum unary cost in the domain
    virtual Value getSupport(int varIndex) const = 0; ///< \brief NC/EAC unary support value
    virtual Value getBestValue(int varIndex) const = 0; ///< \brief hint for some value ordering heuristics (only used by RDS)
    virtual void setBestValue(int varIndex, Value v) = 0; ///< \brief hint for some value ordering heuristics (only used by RDS)
    virtual bool getIsPartOfOptimalSolution() = 0; ///< \brief special flag used for debugging purposes only
    virtual void setIsPartOfOptimalSolution(bool v) = 0; ///< \brief special flag used for debugging purposes only

    virtual int getDegree(int varIndex) const = 0; ///< \brief approximate degree of a variable (\e ie number of active cost functions, see \ref varelim)
    virtual int getTrueDegree(int varIndex) const = 0; ///< \brief degree of a variable
    virtual Long getWeightedDegree(int varIndex) const = 0; ///< \brief weighted degree heuristic
    virtual void resetWeightedDegree() = 0; ///< \brief initialize weighted degree heuristic
    virtual void resetTightness() = 0; ///< \brief initialize constraint tightness used by some heuristics (including weighted degree)
    virtual void resetTightnessAndWeightedDegree() = 0; ///< \brief initialize tightness and weighted degree heuristics

    virtual void preprocessing() = 0; ///< \brief applies various preprocessing techniques to simplify the current problem
    /// \brief sorts the list of cost functions associated to each variable based on smallest problem variable indexes
    /// \warning side-effect: updates DAC order according to an existing variable elimination order
    /// \note must be called after creating all the cost functions and before solving the problem
    virtual void sortConstraints() = 0;

    virtual void whenContradiction() = 0; ///< \brief after a contradiction, resets propagation queues
    virtual void deactivatePropagate() = 0; ///< \brief forbids propagate calls
    virtual bool isactivatePropagate() = 0; ///< \brief are propagate calls authorized?
    virtual void reactivatePropagate() = 0; ///< \brief re-authorizes propagate calls
    virtual void propagate(bool fromscratch = false) = 0; ///< \brief (if authorized) propagates until a fix point is reached (or throws a contradiction). If fromscratch is true then propagates every cost function at least once.
    virtual bool verify() = 0; ///< \brief checks the propagation fix point is reached
#ifdef BOOST
    virtual void addAMOConstraints() = 0;
#endif
    virtual unsigned int numberOfVariables() const = 0; ///< \brief number of created variables
    virtual unsigned int numberOfUnassignedVariables() const = 0; ///< \brief current number of unassigned variables
    virtual unsigned int numberOfConstraints() const = 0; ///< \brief initial number of cost functions (before variable elimination)
    virtual unsigned int numberOfConnectedConstraints() const = 0; ///< \brief current number of cost functions
    virtual unsigned int numberOfConnectedBinaryConstraints() const = 0; ///< \brief current number of binary cost functions
    virtual unsigned int numberOfConnectedKnapsackConstraints() const = 0; ///< \brief current number of knapsack cost functions
    virtual unsigned int medianDomainSize() const = 0; ///< \brief median current domain size of variables
    virtual unsigned int medianDegree() const = 0; ///< \brief median current degree of variables
    virtual unsigned int medianArity() const = 0; ///< \brief median arity of current cost functions
    virtual unsigned int getMaxDomainSize() const = 0; ///< \brief maximum initial domain size found in all variables
    virtual unsigned int getMaxCurrentDomainSize() const = 0; ///< \brief maximum current domain size found in all variables
    virtual unsigned int getDomainSizeSum() const = 0; ///< \brief total sum of current domain sizes
    /// \brief Cartesian product of current domain sizes
    /// \param cartesianProduct result obtained by the GNU Multiple Precision Arithmetic Library GMP
    virtual void cartProd(BigInteger& cartesianProduct) = 0;
    virtual Long getNbDEE() const = 0; ///< \brief number of value removals due to dead-end elimination

    /// \defgroup modeling Variable and cost function modeling
    /// Modeling a Weighted CSP consists in creating variables and cost functions.\n
    /// Domains of variables can be of two different types:
    /// - enumerated domain allowing direct access to each value (array) and iteration on current domain in times proportional to the current number of values (double-linked list)
    /// - interval domain represented by a lower value and an upper value only (useful for large domains)
    /// .
    /// Warning : Current implementation of toulbar2 has limited modeling and solving facilities for interval domains.
    /// There is no cost functions accepting both interval and enumerated variables for the moment, which means all the variables should have the same type.

    /// \addtogroup modeling
    /// Cost functions can be defined in extension (table or maps) or having a specific semantic.\n
    /// Cost functions in extension depend on their arity:
    /// - unary cost function (directly associated to an enumerated variable)
    /// - binary and ternary cost functions (table of costs)
    /// - n-ary cost functions (n >= 4) defined by a list of tuples with associated costs and a default cost for missing tuples (allows for a compact representation)
    /// .
    ///
    /// Cost functions having a specific semantic (see \ref  wcspformat) are:
    /// - simple arithmetic and scheduling (temporal disjunction) cost functions on interval variables
    /// - global cost functions (\e eg soft alldifferent, soft global cardinality constraint, soft same, soft regular, etc) with three different propagator keywords:
    ///   - \e flow propagator based on flow algorithms with "s" prefix in the keyword (\e salldiff, \e sgcc, \e ssame, \e sregular)
    ///   - \e DAG propagator based on dynamic programming algorithms with "s" prefix and "dp" postfix (\e samongdp, salldiffdp, sgccdp, sregulardp, sgrammardp, smstdp, smaxdp)
    ///   - \e network propagator based on cost function network decomposition with "w" prefix (\e wsum, \e wvarsum, \e walldiff, \e wgcc, \e wsame, \e wsamegcc, \e wregular, \e wamong, \e wvaramong, \e woverlap)
    ///   .
    /// .
    ///
    /// Note : The default semantics (using \e var keyword) of monolithic (flow and DAG-based propagators) global cost functions is to count the number of variables to change in order to restore consistency and to multiply it by the basecost. Other particular semantics may be used in conjunction with the flow-based propagator
    ///
    /// Note : The semantics of the network-based propagator approach is either a hard constraint ("hard" keyword) or a soft constraint by multiplying the number of changes by the basecost ("lin" or "var" keyword) or by multiplying the square value of the number of changes by the basecost ("quad" keyword)
    ///
    /// Note : A decomposable version exists for each monolithic global cost function, except grammar and MST. The decomposable ones may propagate less than their monolithic counterpart and they introduce extra variables but they can be much faster in practice
    ///
    /// Warning : Each global cost function may have less than three propagators implemented
    ///
    /// Warning : Current implementation of toulbar2 has limited solving facilities for monolithic global cost functions (no BTD-like methods nor variable elimination)
    ///
    /// Warning : Current implementation of toulbar2 disallows global cost functions with less than or equal to three variables in their scope (use cost functions in extension instead)
    ///
    /// Warning : Before modeling the problem using make and post, call ::tb2init method to initialize toulbar2 global variables
    ///
    /// Warning : After modeling the problem using make and post, call WeightedCSP::sortConstraints method to initialize correctly the model before solving it

    virtual int makeEnumeratedVariable(string n, Value iinf, Value isup) = 0; ///< \brief create an enumerated variable with its domain bounds
    virtual int makeEnumeratedVariable(string n, vector<Value>& dom) = 0; ///< \brief create an enumerated variable with its domain values
    virtual void addValueName(int xIndex, const string& valuename) = 0; ///< \brief add next value name \warning should be called on EnumeratedVariable object as many times as its number of initial domain values
    virtual const string& getValueName(int xIndex, Value value) = 0; ///< \brief return the name associated to a value as defined by addValueName or an empty string if no name found
    virtual int makeIntervalVariable(string n, Value iinf, Value isup) = 0; ///< \brief create an interval variable with its domain bounds

    virtual void postNullaryConstraint(Double cost) = 0; ///< \brief add a zero-arity cost function with floating-point cost
    virtual void postUnaryConstraint(int xIndex, vector<Double>& costs, bool incremental = false) = 0; ///< \brief add a unary cost function with floating-point costs (if incremental is true then it disappears upon backtrack)
    virtual int postBinaryConstraint(int xIndex, int yIndex, vector<Double>& costs, bool incremental = false) = 0; ///< \brief add a binary cost function with floating-point costs (if incremental is true then it disappears upon backtrack)
    virtual int postTernaryConstraint(int xIndex, int yIndex, int zIndex, vector<Double>& costs, bool incremental = false) = 0; ///< \brief add a ternary cost function with floating-point costs (if incremental is true then it disappears upon backtrack)

    virtual void postNullaryConstraint(Cost cost) = 0;
    virtual void postUnary(int xIndex, vector<Cost>& costs) = 0; ///< \deprecated Please use the postUnaryConstraint method instead
    virtual void postUnaryConstraint(int xIndex, vector<Cost>& costs) = 0;
    virtual void postIncrementalUnaryConstraint(int xIndex, vector<Cost>& costs) = 0;
    virtual int postBinaryConstraint(int xIndex, int yIndex, vector<Cost>& costs) = 0;
    virtual int postIncrementalBinaryConstraint(int xIndex, int yIndex, vector<Cost>& costs) = 0;
    virtual int postTernaryConstraint(int xIndex, int yIndex, int zIndex, vector<Cost>& costs) = 0;
    virtual int postIncrementalTernaryConstraint(int xIndex, int yIndex, int zIndex, vector<Cost>& costs) = 0;
    virtual int postNaryConstraintBegin(vector<int> scope, Cost defval, Long nbtuples = 0, bool forcenary = !NARY2CLAUSE) = 0; /// \warning must call WeightedCSP::postNaryConstraintEnd after giving cost tuples
    virtual int postNaryConstraintBegin(int* scope, int arity, Cost defval, Long nbtuples = 0, bool forcenary = !NARY2CLAUSE) = 0; /// \deprecated
    virtual void postNaryConstraintTuple(int ctrindex, vector<Value>& tuple, Cost cost) = 0;
    virtual void postNaryConstraintTuple(int ctrindex, Value* tuple, int arity, Cost cost) = 0; /// \deprecated
    virtual void postNaryConstraintEnd(int ctrindex) = 0; /// \warning must call WeightedCSP::sortConstraints after all cost functions have been posted (see WeightedCSP::sortConstraints)
    virtual int postUnary(int xIndex, Value* d, int dsize, Cost penalty) = 0; ///< \deprecated Please use the postUnaryConstraint method instead
    virtual int postUnaryConstraint(int xIndex, Value* d, int dsize, Cost penalty) = 0;
    virtual int postSupxyc(int xIndex, int yIndex, Value cst, Value deltamax = MAX_VAL - MIN_VAL) = 0;
    virtual int postDisjunction(int xIndex, int yIndex, Value cstx, Value csty, Cost penalty) = 0;
    virtual int postSpecialDisjunction(int xIndex, int yIndex, Value cstx, Value csty, Value xinfty, Value yinfty, Cost costx, Cost costy) = 0;
    virtual int postCliqueConstraint(vector<int> scope, const string& arguments) = 0;
    virtual int postCliqueConstraint(int* scopeIndex, int arity, istream& file) = 0; /// \deprecated
    virtual int postKnapsackConstraint(vector<int> scope, const string& arguments, bool isclique = false, int kp = 0, bool conflict = false, Tuple wcnf = {}) = 0;
    virtual int postKnapsackConstraint(int* scopeIndex, int arity, istream& file, bool isclique = false, int kp = 0, bool conflict = false, Tuple wcnf = {}) = 0; /// \deprecated
    virtual int postWeightedCSPConstraint(vector<int> scope, WeightedCSP* problem, WeightedCSP* negproblem, Cost lb = MIN_COST, Cost ub = MAX_COST, bool duplicateHard = false, bool strongDuality = false) = 0; ///< \brief create a hard constraint such that the input cost function network (problem) must have its optimum cost in [lb,ub[ interval. \warning The input scope must contain all variables in problem in the same order. \warning if duplicateHard is true it assumes any forbidden tuple in the original input problem is also forbidden by another constraint in the main model (you must duplicate any hard constraints in your input model into the main model). \warning if strongDuality is true then it assumes the propagation is complete when all channeling variables in the scope are assigned and the semantic of the constraint enforces that the optimum on the remaining variables is between lb and ub.
    virtual int postGlobalConstraint(int* scopeIndex, int arity, const string& gcname, istream& file, int* constrcounter = NULL, bool mult = true) = 0; ///< \deprecated Please use the postWxxx methods instead
    virtual void postGlobalFunction(vector<int> scope, const string& gcname, const string& arguments) = 0; ///< \brief generic function to post any global cost function

    /// \brief post a soft among cost function
    /// \param scopeIndex an array of variable indexes as returned by WeightedCSP::makeEnumeratedVariable
    /// \param arity the size of the array
    /// \param semantics the semantics of the global cost function: "var" or -- "hard" or "lin" or "quad" (network-based propagator only)--
    /// \param propagator the propagation method (only "DAG" or "network")
    /// \param baseCost the scaling factor of the violation
    /// \param values a vector of values to be restricted
    /// \param lb a fixed lower bound for the number variables to be assigned to the values in \a values
    /// \param ub a fixed upper bound for the number variables to be assigned to the values in \a values
    virtual int postWAmong(vector<int> scope, const string& semantics, const string& propagator, Cost baseCost, const vector<Value>& values, int lb, int ub) = 0; ///< post a soft weighted among cost function
    virtual int postWAmong(int* scopeIndex, int arity, const string& semantics, const string& propagator, Cost baseCost, const vector<Value>& values, int lb, int ub) = 0; ///< \deprecated
    virtual void postWAmong(int* scopeIndex, int arity, string semantics, Cost baseCost, Value* values, int nbValues, int lb, int ub) = 0; ///< \deprecated post a weighted among cost function decomposed as a cost function network
    virtual void postWVarAmong(vector<int> scope, const string& semantics, Cost baseCost, vector<Value>& values, int varIndex) = 0; ///< \brief post a weighted among cost function with the number of values encoded as a variable with index \a varIndex (\e network-based propagator only)
    virtual void postWVarAmong(int* scopeIndex, int arity, const string& semantics, Cost baseCost, Value* values, int nbValues, int varIndex) = 0; ///< \deprecated

    /// \brief post a soft or weighted regular cost function
    /// \param scopeIndex an array of variable indexes as returned by WeightedCSP::makeEnumeratedVariable
    /// \param arity the size of the array
    /// \param semantics the semantics of the soft global cost function: "var" or "edit" (flow-based propagator) or -- "var" (DAG-based propagator)-- (unused parameter for network-based propagator)
    /// \param propagator the propagation method ("flow", "DAG", "network")
    /// \param baseCost the scaling factor of the violation ("flow", "DAG")
    /// \param nbStates the number of the states in the corresponding DFA. The states are indexed as 0, 1, ..., nbStates-1
    /// \param initial_States a vector of WeightedObjInt specifying the starting states with weight
    /// \param accepting_States a vector of WeightedObjInt specifying the final states
    /// \param Wtransitions a vector of (weighted) transitions
    /// \warning Weights are ignored in the current implementation of DAG and flow-based propagators
    virtual int postWRegular(vector<int> scope, const string& semantics, const string& propagator, Cost baseCost, int nbStates, const vector<WeightedObjInt>& initial_States, const vector<WeightedObjInt>& accepting_States, const vector<DFATransition>& Wtransitions) = 0; ///< post a soft weighted regular cost function
    virtual int postWRegular(int* scopeIndex, int arity, const string& semantics, const string& propagator, Cost baseCost, int nbStates, const vector<WeightedObjInt>& initial_States, const vector<WeightedObjInt>& accepting_States, const vector<DFATransition>& Wtransitions) = 0; ///< \deprecated
    virtual void postWRegular(int* scopeIndex, int arity, int nbStates, vector<pair<int, Cost>> initial_States, vector<pair<int, Cost>> accepting_States, int** Wtransitions, vector<Cost> transitionsCosts) = 0; ///< \deprecated post a weighted regular cost function decomposed as a cost function network

    /// \brief post a soft alldifferent cost function
    /// \param scopeIndex an array of variable indexes as returned by WeightedCSP::makeEnumeratedVariable
    /// \param arity the size of the array
    /// \param semantics the semantics of the global cost function: for flow-based propagator: "var" or "dec" or "decbi" (decomposed into a binary cost function complete network), for DAG-based propagator: "var", for network-based propagator: "hard" or "lin" or "quad" (decomposed based on wamong)
    /// \param propagator the propagation method ("flow", "DAG", "network")
    /// \param baseCost the scaling factor of the violation
    virtual int postWAllDiff(vector<int> scope, const string& semantics, const string& propagator, Cost baseCost) = 0; ///< post a soft alldifferent cost function
    virtual int postWAllDiff(int* scopeIndex, int arity, const string& semantics, const string& propagator, Cost baseCost) = 0; ///< \deprecated
    virtual void postWAllDiff(int* scopeIndex, int arity, string semantics, Cost baseCost) = 0; ///< \deprecated post a soft alldifferent cost function decomposed as a cost function network

    /// \brief post a soft global cardinality cost function
    /// \param scopeIndex an array of variable indexes as returned by WeightedCSP::makeEnumeratedVariable
    /// \param arity the size of the array
    /// \param semantics the semantics of the global cost function: "var" (DAG-based propagator only) or -- "var" or "dec" or "wdec" (flow-based propagator only) or -- "hard" or "lin" or "quad" (network-based propagator only)--
    /// \param propagator the propagation method ("flow", "DAG", "network")
    /// \param baseCost the scaling factor of the violation
    /// \param values a vector of BoundedObjValue, specifying the lower and upper bounds of each value, restricting the number of variables can be assigned to them
    virtual int postWGcc(int* scopeIndex, int arity, const string& semantics, const string& propagator, Cost baseCost,
        const vector<BoundedObjValue>& values)
        = 0;
    virtual void postWGcc(int* scopeIndex, int arity, string semantics, Cost baseCost, Value* values, int nbValues, int* lb, int* ub) = 0; ///< \deprecated post a soft global cardinality cost function decomposed as a cost function network

    /// \brief post a soft same cost function (a group of variables being a permutation of another group with the same size)
    /// \param scopeIndexG1 an array of the first group of variable indexes as returned by WeightedCSP::makeEnumeratedVariable
    /// \param arityG1 the size of \a scopeIndexG1
    /// \param scopeIndexG2 an array of the second group of variable indexes as returned by WeightedCSP::makeEnumeratedVariable
    /// \param arityG2 the size of \a scopeIndexG2
    /// \param semantics the semantics of the global cost function: "var" or -- "hard" or "lin" or "quad" (network-based propagator only)--
    /// \param propagator the propagation method ("flow" or "network")
    /// \param baseCost the scaling factor of the violation.
    virtual int postWSame(int* scopeIndexG1, int arityG1, int* scopeIndexG2, int arityG2, const string& semantics, const string& propagator, Cost baseCost) = 0;
    virtual void postWSame(int* scopeIndex, int arity, string semantics, Cost baseCost) = 0; ///< \deprecated post a soft same cost function
    virtual void postWSameGcc(int* scopeIndex, int arity, string semantics, Cost baseCost, Value* values, int nbValues, int* lb, int* ub) = 0; ///< \brief post a combination of a same and gcc cost function decomposed as a cost function network

    /// \brief post a soft/weighted grammar cost function with the dynamic programming propagator and grammar in Chomsky normal form
    /// \param scopeIndex an array of the first group of variable indexes as returned by WeightedCSP::makeEnumeratedVariable
    /// \param arity the size of \a scopeIndex
    /// \param semantics the semantics of the global cost function: "var" or "weight"
    /// \param propagator the propagation method ("DAG" only)
    /// \param baseCost the scaling factor of the violation
    /// \param nbSymbols the number of symbols in the corresponding grammar. Symbols are indexed as 0, 1, ..., nbSymbols-1
    /// \param startSymbol the index of the starting symbol
    /// \param WRuleToTerminal a vector of \a ::CFGProductionRule. Note that:
    ///  - if \a order in \a CFGProductionRule is set to 0, it is classified as A -> v, where A is the index of the terminal symbol and v is the value.
    ///  - if \a order in \a CFGProductionRule is set to 1, it is classified as A -> BC, where A,B,C  the index of the nonterminal symbols.
    ///  - if \a order in \a CFGProductionRule is set to 2, it is classified as weighted A -> v, where A is the index of the terminal symbol and v is the value.
    ///  - if \a order in \a CFGProductionRule is set to 3, it is classified as weighted A -> BC, where A,B,C  the index of the nonterminal symbols.
    ///  - if \a order in \a CFGProductionRule is set to values greater than 3, it is ignored.
    virtual int postWGrammarCNF(int* scopeIndex, int arity, const string& semantics, const string& propagator, Cost baseCost,
        int nbSymbols,
        int startSymbol,
        const vector<CFGProductionRule> WRuleToTerminal)
        = 0;

    /// \brief post a Spanning Tree hard constraint
    /// \param scopeIndex an array of variable indexes as returned by WeightedCSP::makeEnumeratedVariable
    /// \param arity the size of \a scopeIndex
    /// \param semantics the semantics of the global cost function: "hard"
    /// \param propagator the propagation method ("DAG" only)
    /// \param baseCost unused in the current implementation (MAX_COST)
    virtual int postMST(int* scopeIndex, int arity, const string& semantics, const string& propagator, Cost baseCost) = 0;

    /// \brief post a weighted max cost function (maximum cost of a set of unary cost functions associated to a set of variables)
    /// \param scopeIndex an array of variable indexes as returned by WeightedCSP::makeEnumeratedVariable
    /// \param arity the size of \a scopeIndex
    /// \param semantics the semantics of the global cost function: "val"
    /// \param propagator the propagation method ("DAG" only)
    /// \param baseCost if a variable-value pair does not exist in \a weightFunction, its weight will be mapped to baseCost.
    /// \param weightFunction a vector of WeightedVarValPair containing a mapping from variable-value pairs to their weights.
    virtual int postMaxWeight(int* scopeIndex, int arity, const string& semantics, const string& propagator, Cost baseCost,
        const vector<WeightedVarValPair> weightFunction)
        = 0;

    /// \brief post a soft linear constraint with unit coefficients
    /// \param scopeIndex an array of variable indexes as returned by WeightedCSP::makeEnumeratedVariable
    /// \param arity the size of \a scopeIndex
    /// \param semantics the semantics of the global cost function: "hard" or "lin" or "quad" (network-based propagator only)
    /// \param propagator the propagation method ("network" only)
    /// \param baseCost the scaling factor of the violation
    /// \param comparator the comparison operator of the linear constraint ("==", "!=", "<", "<=", ">,", ">=")
    /// \param rightRes right-hand side value of the linear constraint
    virtual void postWSum(int* scopeIndex, int arity, string semantics, Cost baseCost, string comparator, int rightRes) = 0;
    virtual void postWVarSum(int* scopeIndex, int arity, string semantics, Cost baseCost, string comparator, int varIndex) = 0; ///< \brief post a soft linear constraint with unit coefficients and variable right-hand side

    /// \brief post a soft overlap cost function (a group of variables being point-wise equivalent -- and not equal to zero -- to another group with the same size)
    /// \param scopeIndex an array of variable indexes as returned by WeightedCSP::makeEnumeratedVariable
    /// \param arity the size of \a scopeIndex (should be an even value)
    /// \param semantics the semantics of the global cost function: "hard" or "lin" or "quad" (network-based propagator only)
    /// \param propagator the propagation method ("network" only)
    /// \param baseCost the scaling factor of the violation.
    /// \param comparator the point-wise comparison operator applied to the number of equivalent variables ("==", "!=", "<", "<=", ">,", ">=")
    /// \param rightRes right-hand side value of the comparison
    virtual void postWOverlap(int* scopeIndex, int arity, string semantics, Cost baseCost, string comparator, int rightRes) = 0;

    /// \brief post a diversity Hamming distance constraint between a list of variables and a given fixed assignment
    /// \param scope a vector of variable indexes as returned by WeightedCSP::makeEnumeratedVariable
    /// \param distance the Hamming distance minimum bound
    /// \param values a vector of values (same size as scope)
    /// \param method the network decomposition method (0:Dual, 1:Hidden, 2:Ternary)
    /// \note depending on the decomposition method, it adds dual and/or hidden variables
    virtual void postWDivConstraint(vector<int> scope, unsigned int distance, vector<Value>& values, int method = 0) = 0;

    virtual vector<vector<int>>* getListSuccessors() = 0; ///< \brief generating additional variables vector created when berge decomposition are included in the WCSP
    virtual vector<int> getBergeDecElimOrder() = 0; ///< \brief return an elimination order compatible with Berge acyclic decomposition of global decomposable cost functions (if possible keep reverse of previous DAC order)
    virtual void setDACOrder(vector<int>& elimVarOrder) = 0; ///< \brief change DAC order and propagate from scratch

    virtual bool isKnapsack() = 0; ///< \brief true if there are knapsack constraints defined in the problem
    virtual bool isGlobal() = 0; ///< \brief true if there are soft global constraints defined in the problem
#ifdef ILOGCPLEX
    virtual bool isPLPS() = 0; ///< \brief true if there are Polytime Linear Projection-Safe global cost functions (slinear)
#endif

    virtual Cost read_wcsp(const char* fileName) = 0; ///< \brief load problem in all format supported by toulbar2. Returns the UB known to the solver before solving (file and command line).
    virtual void read_legacy(const char* fileName) = 0; ///< \brief load problem in wcsp legacy format
    virtual void read_uai2008(const char* fileName) = 0; ///< \brief load problem in UAI 2008 format (see http://graphmod.ics.uci.edu/uai08/FileFormat and http://www.cs.huji.ac.il/project/UAI10/fileFormat.php) \warning UAI10 evidence file format not recognized by toulbar2 as it does not allow multiple evidence (you should remove the first value in the file)
    virtual void read_random(int n, int m, vector<int>& p, int seed, bool forceSubModular = false, string globalname = "") = 0; ///< \brief create a random WCSP with \e n variables, domain size \e m, array \e p where the first element is a percentage of tuples with a nonzero cost and next elements are the number of random cost functions for each different arity (starting with arity two), random seed, a flag to have a percentage (last element in the array \e p) of the binary cost functions being permutated submodular, and a string to use a specific global cost function instead of random cost functions in extension
    virtual void read_wcnf(const char* fileName) = 0; ///< \brief load problem in (w)cnf format (see http://www.maxsat.udl.cat/08/index.php?disp=requirements)
    virtual void read_qpbo(const char* fileName) = 0; ///< \brief load quadratic pseudo-Boolean optimization problem in unconstrained quadratic programming text format (first text line with n, number of variables and m, number of triplets, followed by the m triplets (x,y,cost) describing the sparse symmetric nXn cost matrix with variable indexes such that x <= y and any positive or negative real numbers for costs)
    virtual void read_opb(const char* fileName) = 0; ///< \brief load pseudo-Boolean optimization problem
    virtual void read_lp(const char* fileName) = 0; ///< \brief load integer linear programming problem

    virtual const vector<Value> getSolution() = 0; ///< \brief after solving the problem, return the optimal solution (warning! do not use it if doing solution counting or if there is no solution, see WeightedCSPSolver::solve output for that)
    virtual Double getSolutionValue() const = 0; ///< \brief returns current best solution cost or MAX_COST if no solution found
    virtual Cost getSolutionCost() const = 0; ///< \brief returns current best solution cost or MAX_COST if no solution found
    virtual const vector<Value> getSolution(Cost* cost_ptr) = 0; ///< \deprecated \brief returns current best solution and its cost
    virtual vector<pair<Double, vector<Value>>> getSolutions() const = 0; ///< \brief returns all solutions found
    virtual void initSolutionCost() = 0; ///< \brief invalidate best solution by changing its cost to MAX_COST
    virtual void setSolution(Cost cost, TAssign* sol = NULL) = 0; ///< \brief set best solution from current assigned values or from a given assignment (for BTD-like methods)
    virtual void printSolution() = 0; ///< \brief prints current best solution on standard output (using variable and value names if cfn format and ToulBar2::showSolution>1)
    virtual void printSolution(ostream& os) = 0; ///< \brief prints current best solution (using variable and value names if cfn format and ToulBar2::writeSolution>1)
    virtual void printSolution(FILE* f) = 0; ///< \brief prints current best solution (using variable and value names if cfn format and ToulBar2::writeSolution>1)

    virtual void print(ostream& os) = 0; ///< \brief print current domains and active cost functions (see \ref verbosity)
    virtual void dump(ostream& os, bool original = true) = 0; ///< \brief output the current WCSP into a file in wcsp format \param os output file \param original if true then keeps all variables with their original domain size else uses unassigned variables and current domains recoding variable indexes
    virtual void dump_CFN(ostream& os, bool original = true) = 0; ///< \brief output the current WCSP into a file in wcsp format \param os output file \param original if true then keeps all variables with their original domain size else uses unassigned variables and current domains recoding variable indexes

    // -----------------------------------------------------------
    // Functions dealing with all representations of Costs
    // warning: ToulBar2::NormFactor has to be initialized

    virtual Cost decimalToCost(const string& decimalToken, const unsigned int lineNumber) const = 0;
    virtual Cost DoubletoCost(const Double& c) const = 0;
    virtual Double Cost2ADCost(const Cost& c) const = 0; ///< \brief converts an integer cost from a lower or upper bound of the whole problem into a real value
    virtual Double Cost2RDCost(const Cost& c) const = 0; ///< \brief converts an integer cost from a local cost function into a real value
    virtual Cost Prob2Cost(TProb p) const = 0;
    virtual TProb Cost2Prob(Cost c) const = 0;
    virtual TLogProb Cost2LogProb(Cost c) const = 0;
    virtual Cost LogProb2Cost(TLogProb p) const = 0;
    virtual Cost LogSumExp(Cost c1, Cost c2) const = 0;
    virtual TLogProb LogSumExp(TLogProb logc1, Cost c2) const = 0;
    virtual TLogProb LogSumExp(TLogProb logc1, TLogProb logc2) const = 0;

    // -----------------------------------------------------------
    // Internal WCSP functions DO NOT USE THEM

    virtual void setLb(Cost newLb) = 0; ///< \internal sets problem lower bound
    virtual void setUb(Cost newUb) = 0; ///< \internal sets problem upper bound
    virtual void restoreSolution(Cluster* c = NULL) = 0; ///< \internal restores correct values to eliminated variables when all the variables have been assigned

    virtual void buildTreeDecomposition() = 0;
    virtual TreeDecomposition* getTreeDec() = 0;

    virtual vector<Variable*>& getDivVariables() = 0; ///< \brief returns all variables on which a diversity request exists
    virtual void initDivVariables() = 0; ///< \brief initializes diversity variables with all decision variables in the problem

    virtual void iniSingleton() = 0;
    virtual void updateSingleton() = 0;
    virtual void removeSingleton() = 0;
    virtual void printVACStat() = 0;
};

ostream& operator<<(ostream& os, WeightedCSP& wcsp); ///< \see WeightedCSP::print

/** Abstract class WeightedCSPSolver representing a WCSP solver
 *	- link to a WeightedCSP
 *	- generic complete solving method configurable through global variables (see ::ToulBar2 class and command line options)
 *	- optimal solution available after problem solving
 *	- elementary decision operations on domains of variables
 *	- statistics information (number of nodes and backtracks)
 *	- problem file format reader (multiple formats, see \ref wcspformat)
 *	- solution checker (output the cost of a given solution)
 *
 */

class WeightedCSPSolver {
public:
#ifdef OPENMPI
    static const int MASTER = 0; // Master MPI rank number
    static const int WORKTAG = 1; // MPI tag value for still working
    static const int DIETAG = 2; // MPI tag value for stop working
    static const int IDLETAG = 3; // MPI tag value for no more working
#endif
    static WeightedCSPSolver* makeWeightedCSPSolver(Cost initUpperBound, WeightedCSP* wcsp = NULL); ///< \brief WeightedCSP Solver factory

    virtual ~WeightedCSPSolver() {}

    virtual WeightedCSP* getWCSP() = 0; ///< \brief access to its associated Weighted CSP

    virtual Long getNbNodes() const = 0; ///< \brief number of search nodes (see WeightedCSPSolver::increase, WeightedCSPSolver::decrease, WeightedCSPSolver::assign, WeightedCSPSolver::remove)
    virtual Long getNbBacktracks() const = 0; ///< \brief number of backtracks

    virtual void increase(int varIndex, Value value, bool reverse = false) = 0; ///< \brief changes domain lower bound and propagates
    virtual void decrease(int varIndex, Value value, bool reverse = false) = 0; ///< \brief changes domain upper bound and propagates
    virtual void assign(int varIndex, Value value, bool reverse = false) = 0; ///< \brief assigns a variable and propagates
    virtual void remove(int varIndex, Value value, bool reverse = false) = 0; ///< \brief removes a domain value and propagates (valid if done for an enumerated variable or on its domain bounds)

    /** \defgroup solving Solving cost function networks
     * After creating a Weighted CSP, it can be solved using a local search method like INCOP or PILS (see WeightedCSPSolver::narycsp or WeightedCSPSolver::pils) and/or an exact search method (see WeightedCSPSolver::solve).
     *
     * Various options of the solving methods are controlled by ::Toulbar2 static class members (see files ./src/core/tb2types.hpp and ./src/tb2main.cpp).\n
     * A brief code example reading a wcsp problem given as a single command-line parameter and solving it:
     * \code
    #include "toulbar2lib.hpp"
    #include <string.h>
    #include <stdio.h>
    #include <stdlib.h>
    #include <unistd.h>
    int main(int argc, char **argv) {

        tb2init(); // must be call before setting specific ToulBar2 options and creating a model

        // Create a solver object
        initCosts(); // last check for compatibility issues between ToulBar2 options and Cost data-type
        WeightedCSPSolver *solver = WeightedCSPSolver::makeWeightedCSPSolver(MAX_COST);

        // Read a problem file in wcsp format
        solver->read_wcsp(argv[1]);

        ToulBar2::verbose = -1;  // change to 0 or higher values to see more trace information

        // Uncomment if solved using INCOP local search followed by a partial Limited Discrepancy Search with a maximum discrepancy of one
        //  ToulBar2::incop_cmd = "0 1 3 idwa 100000 cv v 0 200 1 0 0";
        //  ToulBar2::lds = -1;  // remove it or change to a positive value then the search continues by a complete B&B search method
        // Uncomment the following lines if solved using Decomposition Guided Variable Neighborhood Search with min-fill cluster decomposition and absorption
        // ToulBar2::lds = 4;
        // ToulBar2::restart = 10000;
        // ToulBar2::searchMethod = DGVNS;
        // ToulBar2::vnsNeighborVarHeur = CLUSTERRAND;
        // ToulBar2::boostingBTD = 0.7;
        // ToulBar2::varOrder = reinterpret_cast<char*>(-3);

        if (solver->solve()) {
            // show (sub-)optimal solution
            vector<Value> sol;
            Cost ub = solver->getSolution(sol);
            cout << "Best solution found cost: " << ub << endl;
            cout << "Best solution found:";
            for (unsigned int i=0; i<sol.size(); i++) cout << ((i>0)?",":"") << " x" << i << " = " << sol[i];
            cout << endl;
        } else {
            cout << "No solution found!" << endl;
        }
        delete solver;
    }
    \endcode
     *
     * See : another code example in ./src/toulbar2test.cpp
     *
     * Warning : variable domains must start at zero, otherwise recompile libtb2.so without flag WCSPFORMATONLY
    **/

    virtual Cost read_wcsp(const char* fileName) = 0; ///< \brief reads a Cost function network from a file (format as indicated by ToulBar2:: global variables)
    virtual void read_random(int n, int m, vector<int>& p, int seed, bool forceSubModular = false, string globalname = "") = 0; ///< \brief create a random WCSP, see WeightedCSP::read_random

    /// \brief simplifies and solves to optimality the problem
    /// \return false if there is no solution found
    /// \warning after solving, the current problem has been modified by various preprocessing techniques
    /// \warning DO NOT READ VALUES OF ASSIGNED VARIABLES USING WeightedCSP::getValue (temporally wrong assignments due to variable elimination in preprocessing) BUT USE WeightedCSPSolver::getSolution INSTEAD
    virtual bool solve(bool first = true) = 0;

    // internal methods called by solve, for advanced programmers only!!!
    virtual void beginSolve(Cost ub) = 0;
    virtual Cost preprocessing(Cost ub) = 0;
    virtual void recursiveSolve(Cost lb = MIN_COST) = 0;
    virtual void recursiveSolveLDS(int discrepancy) = 0;
    virtual pair<Cost, Cost> hybridSolve() = 0;
    virtual void endSolve(bool isSolution, Cost cost, bool isComplete) = 0;
    // end of internal solve methods

    /// \brief solves the current problem using INCOP local search solver by Bertrand Neveu
    /// \return best solution cost found
    /// \param cmd command line argument for narycsp INCOP local search solver (cmd format: lowerbound randomseed nbiterations method nbmoves neighborhoodchoice neighborhoodchoice2 minnbneighbors maxnbneighbors  neighborhoodchoice3 autotuning tracemode)
    /// \param solution best solution assignment found (MUST BE INITIALIZED WITH A DEFAULT COMPLETE ASSIGNMENT)
    /// \warning cannot solve problems with global cost functions
    /// \note side-effects: updates current problem upper bound and propagates, best solution saved (using WCSP::setBestValue)
    virtual Cost narycsp(string cmd, vector<Value>& solution) = 0;

    // PILS local search
    /// \brief solves the current problem using PILS local search @ Francois Beuvin, David Simoncini, Sebastien Verel
    /// \return best solution cost found
    /// \param cmd command line argument for PILS local search solver (cmd format: nbruns perturb_mode perturb_strength flatMaxIter nbEvalHC nbEvalMax strengthMin strengthMax incrFactor decrFactor)
    /// \warning cannot solve problems with non-binary cost functions
    virtual Cost pils(string cmd, vector<Value>& solution) = 0;

    // LR-BCD local search
    /// \brief solves the current problem using LR-BCD local search @ Valentin Durante, George Katsirelos, Thomas Schiex
    /// \return best solution cost found
    /// \param cmd command line argument for LR-BCD local search solver (cmd format: maxiter rank nbroundings)
    /// \warning cannot solve problems with non-binary cost functions
    virtual Cost lrBCD(string cmd, vector<Value>& solution) = 0;

    /// \brief quadratic unconstrained pseudo-Boolean optimization
    /// Maximize \f$h' \times W \times h\f$ where \f$W\f$ is expressed by all its
    /// non-zero half squared matrix costs (can be positive or negative, with \f$\forall i, posx[i] \leq posy[i]\f$)
    /// \note costs for \f$posx \neq posy\f$ are multiplied by 2 by this method
    /// \note by convention: \f$h = 1 \equiv x = 0\f$ and \f$h = -1 \equiv x = 1\f$
    /// \warning does not allow infinite costs (no forbidden assignments, unconstrained optimization)
    /// \return true if at least one solution has been found (array \e sol being filled with the best solution)
    /// \see ::solvesymmax2sat_ for Fortran call
    virtual bool solve_symmax2sat(int n, int m, int* posx, int* posy, double* cost, int* sol) = 0;

    virtual void dump_wcsp(const char* fileName, bool original = true, ProblemFormat format = WCSP_FORMAT) = 0; ///< \brief output current problem in a file \see WeightedCSP::dump
    virtual void read_solution(const char* fileName, bool updateValueHeuristic = true) = 0; ///< \brief read a solution from a file
    virtual void parse_solution(const char* certificate, bool updateValueHeuristic = true) = 0; ///< \brief read a solution from a string (see ToulBar2 option \e -x)

    virtual const vector<Value> getSolution() = 0; ///< \brief after solving the problem, return the optimal solution (warning! do not use it if doing solution counting or if there is no solution, see WeightedCSPSolver::solve output for that)
    virtual Double getSolutionValue() const = 0; ///< \brief after solving the problem, return the optimal solution value (can be an arbitrary real cost in minimization or preference in maximization, see CFN format) (warning! do not use it if doing solution counting or if there is no solution, see WeightedCSPSolver::solve output for that)
    virtual Cost getSolutionCost() const = 0; ///< \brief after solving the problem, return the optimal solution nonnegative integer cost (warning! do not use it if doing solution counting or if there is no solution, see WeightedCSPSolver::solve output for that)
    virtual Cost getSolution(vector<Value>& solution) const = 0; ///< \deprecated \brief after solving the problem, add the optimal solution in the input/output vector and returns its optimum cost (warning! do not use it if doing solution counting or if there is no solution, see WeightedCSPSolver::solve output for that)
    virtual vector<pair<Double, vector<Value>>> getSolutions() const = 0; ///< \brief after solving the problem, return all solutions found with their corresponding value
    virtual Double getDDualBound() const = 0; ///< \brief after (partially) solving the problem (possibly interrupted before the search is complete), return a global problem dual bound as a Double representing a decimal cost (lower resp. upper bound for minimization resp. maximization)

    // -----------------------------------------------------------
    // Internal Solver functions DO NOT USE THEM

    virtual set<int> getUnassignedVars() const = 0; ///< \internal returns the set of unassigned variable indexes \warning not valid before the search (see WeightedCSPSolver::solve)
    virtual int numberOfUnassignedVariables() const = 0; ///< \internal returns the number of unassigned variables \warning not valid before the search (returns -1)
};

/// \brief initialization of ToulBar2 global variables (needed by numberjack/toulbar2)
extern void tb2init();
/// \brief reinitialization of ToulBar2 global variables before next call to solving methods
extern void tb2reinit();
/// \brief checks compatibility between selected options of ToulBar2 (needed by numberjack/toulbar2)
extern void tb2checkOptions();
#endif /*TOULBAR2LIB_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
