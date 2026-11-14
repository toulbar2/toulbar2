/** \file tb2alldifferent.hpp
 *  \brief Propagates an All-Different or a Permutation Constraint.
 *
 *  It maintains different levels of consistency:
 *  - tests if the constraint is always satisfied or violated when enough variables are assigned (NC level)
 *  TODO - enforces domain/bound arc consistency considering the hard constraint only (AC level)
 *  TODO - approximates full zero-inverse soft consistency (FDAC or EDAC or VAC/AC)
 *  TODO - enforces virtual arc consistency (VAC-lin)
 *
 *  TOBEVERIFIED full zero-inverse soft consistency may be not enforced because it requires propagating at every unary cost move and not only when it increases/projects from zero???
 *
 */

#ifndef TB2ALLDIFFERENT_HPP_
#define TB2ALLDIFFERENT_HPP_

#pragma once
#include <utility>
#include <variant>
#include <functional>
#include <numeric>
#include "tb2abstractconstr.hpp"
#include "tb2ternaryconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2wcsp.hpp"
#include "search/tb2clusters.hpp"
#include "utils/tb2lapjv.hpp"
#include <unordered_set>
#include <unordered_map>
#include <boost/heap/pairing_heap.hpp>
#include <random>
#include <algorithm>

using namespace std;
using PQ = boost::heap::pairing_heap<pair<Cost, int>, boost::heap::compare<greater<>>>;

class AllDifferentConstraint : public AbstractNaryConstraint {
    Cost Original_ub; // Initial upper bound when creating the constraint
    StoreCost lb; // Projected cost to problem lower bound (if zero, all deltaCosts must be zero)
    StoreCost assigneddeltas; // Accumulated deltas from assigned values (used in cost propagation)
    vector<Long> conflictWeights; // Used by weighted degree heuristics to prioritize variables
    vector<StoreCost> deltaCosts; // Extended unary costs to all values
    vector<int> ValuesCapacity; // Maximum number of variables allowed to take each value (1, or arity_ for excepted values)
    vector<StoreValue> storeLastAssignment; // Stores the last optimal assignment of variables
    StoreInt storeAssignment; // True if an optimal assignment is currently stored
    int* rowSol = nullptr; // Row solution array: rowSol[var] = value index assigned to variable var by the LAP solver
    int* colSol = nullptr; // Column solution array: colSol[val] = variable index assigned to value val by the LAP solver
    Cost* ReduceCostRow = nullptr; // Reduced costs per row (per variable) after solving the LAP
    Cost* ReduceCostCol = nullptr; // Reduced costs per column (per value) after solving the LAP
    vector<Cost> costMatrix; // Flattened cost matrix passed to the Jonker LAP solver (layout depends on the case, see the index lambdas)
    int NbValues; // Total number of distinct values across all variable domains (size of UnionVarDomain)
    vector<int> AssignedVar; // Indices (in scope) of currently assigned variables
    vector<int> AssignedVal; // Value indices (in UnionVarDomain) assigned to each variable in AssignedVar
    vector<int> NoAssignedVar; // Indices (in scope) of currently unassigned variables
    int NbAssigned; // Number of currently assigned variables
    int NbNoAssigned; // Number of currently unassigned variables
    int NbNoAssignedVal; // Number of values not yet taken by any assigned variable (columns in the reduced cost matrix)
    vector<bool> isAssignedValue; // isAssignedValue[v] = true if value v (index in UnionVarDomain) is already taken by an assigned variable
    vector<bool> varAlreadyProcessed; // varAlreadyProcessed[i] = true if variable i has already been processed during domain filtering
    vector<Value> exceptedValues; // Values exempt from the AllDifferent constraint (may be shared by multiple variables)
    vector<int> exceptedValIndex; // Indices in UnionVarDomain of each excepted value
    vector<uint8_t> isExceptedVal; // isExceptedVal[v] = 1 if value v (index in UnionVarDomain) is an excepted value
    bool excepted; // True if at least one excepted value exists, false otherwise
    bool isSquare; // True if the cost matrix is square (number of variables == number of values)
    vector<string> UnionVarDomain; // Sorted union of all variable domains (unique value names)
    vector<int> VarDomainSize; // Initial domain size for each variable (indexed by variable position in scope)
    unordered_map<string, int> mapDomainValToIndex; // Maps each value name to its index in UnionVarDomain
    bool SameDomain; // True if all variables share the same domain
    vector<Cost> ReduceCostMatrix; // Reduced cost matrix after LAP solving: ReduceCostMatrix[var * NbNoAssignedVal + val] = reduced cost of assigning value val to variable var
    vector<uint8_t> visited; // visited[v] = 1 if variable/value node v has been finalized by Dijkstra
    vector<uint8_t> inHeap; // inHeap[v] = 1 if variable/value node v is currently in the priority queue
    vector<Cost> distanceToVar; // Shortest distances from the Dijkstra source variable to each other variable node in the residual graph
    vector<Cost> distanceToVal; // Shortest distances from the Dijkstra source variable to each value node in the residual graph (reverse graph)
    int Q; // Number of Dijkstra runs to perform during reduced-cost filtering (derived from FiltLevel and NbNoAssigned)
    double FiltLevel; // Fraction of variables used as Dijkstra sources during filtering (0 = disabled, 1 = all variables)
    vector<int> trackList; // Variables not yet finalized by Dijkstra (dense array, used by the bimodal traversal strategy)
    vector<int> trackPos; // trackPos[var] = position of var in trackList, or -1 if finalized
    unordered_set<int> trackingList; // Set of variable/value nodes not yet finalized by Dijkstra, used by the bimodal traversal strategy
    vector<PQ::handle_type> dijkstraHandles; // Handles of the priority queue, reused between two Dijkstra runs
    Cost MaxReducedCost; // Maximum reduced cost observed over all arcs in the current residual graph (used as filtering threshold)
    Cost ReducedCost; // Temporary reduced cost of a single arc computed during cost projection
    Cost VarMaxReducedCost; // Maximum reduced cost among all arcs incident to a given variable (used to guide source selection)
    int findConflict; // Flag returned by the LAP solver: non-zero if a conflict (infeasible sub-problem) was detected, value encodes the number of conflicting variables
    int lower, upper; // Lower and upper thresholds on NbNoAssigned controlling when reduced-cost filtering is activated (lower) and when Dijkstra-based filtering is applied (upper)
    vector<int> varValIdx; // varValIdx[var * NbValues + val] = index of value val (in UnionVarDomain) inside the initial domain of variable var, or -1 if var cannot take it
    vector<int> capacity; // Capacity of each column/row of the LAP in the slack formulation (the dummy "slack" variable absorbs the unused values)
    StoreInt isSlack = 0; // 1 while the slack (dummy variable) formulation is in use, backtrackable
    bool slack = false; // True if the slack formulation was selected when reading the constraint
    bool relativeInterior = false; // If true, the dual solution of the (square or slack) LAP is moved to the relative interior of the dual optimizers before projecting costs

    // Projects cost c to the problem lower bound (only if it is strictly positive)
    void projectLB(Cost c)
    {
        if (c > MIN_COST) {
            Constraint::projectLB(c);
        }
    }

    // returns true if the constraint can be projected to small local cost function in extension
    // bool canbeProjectedInExtension();

    // Depending of the value and the cost, extend or project the cost on the index value of the variable var
    void ExtOrProJ(int var, Value value, Cost C)
    {
        // int value_idx = scope[var]->toIndex(value);
        TreeDecomposition* td = wcsp->getTreeDec();
        if (C > MIN_COST) {
            if (!isSquare) {
                if (td && scope[var]->canbe(value))
                    td->addDelta(cluster, scope[var], value, -C);
            }
            assert(scope[var]->getCost(value) >= C);
            scope[var]->extend(value, C);
        } else if (C < MIN_COST) {
            if (!isSquare) {
                if (td && scope[var]->canbe(value))
                    td->addDelta(cluster, scope[var], value, -C);
            }
            scope[var]->project(value, -C, true);
        }
    }

    // Converts a value index of UnionVarDomain into the actual Value of variable varIndex (WRONG_VAL if the variable does not have this value in its initial domain)
    Value toVarValue(int varIndex, int valIndex) const
    {
        int idx = varValIdx[varIndex * NbValues + valIndex];
        return (idx < 0) ? WRONG_VAL : scope[varIndex]->toValue((unsigned int)idx);
    }

    // Returns the index in UnionVarDomain of a value name, or -1 if the name is unknown or if the value has been removed from the union
    int unionIndex(const string& name) const
    {
        auto it = mapDomainValToIndex.find(name);
        if (it == mapDomainValToIndex.end() || it->second >= NbValues)
            return -1;
        return it->second;
    }

    // Rebuilds varValIdx and VarDomainSize from mapDomainValToIndex and NbValues
    void rebuildVarMapping()
    {
        varValIdx.assign((size_t)arity_ * NbValues, -1);
        for (int varIndex = 0; varIndex < arity_; varIndex++) {
            auto* variable = scope[varIndex];
            int count = 0;
            for (int valIndex = 0; valIndex < (int)variable->getDomainInitSize(); valIndex++) {
                auto it = mapDomainValToIndex.find(variable->getValueName(valIndex));
                if (it == mapDomainValToIndex.end() || it->second >= NbValues)
                    continue; // value removed from the union
                varValIdx[varIndex * NbValues + it->second] = valIndex;
                count++;
            }
            VarDomainSize[varIndex] = count;
            if (count == 0)
                THROWCONTRADICTION;
        }
    }


public:
    // Constructor for the AllDifferentConstraint class
    AllDifferentConstraint(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in)
        : AbstractNaryConstraint(wcsp, scope_in, arity_in)
        , Original_ub(wcsp->getUb())
        , lb(0)
        , assigneddeltas(0)
        , storeAssignment(false)
        , NbValues(0)
        , NbAssigned(0)
        , NbNoAssigned(arity_in)
        , excepted(false)
        , isSquare(false)
        , SameDomain(true)
        , FiltLevel(0.0)

    {
        if (arity_in > 0) {
            isSquare = true;
            unordered_set<string> seen;
            vector<string> TempUnionVarDomain;
            NbValues = 0;
            int DomainSize;
            int varIndex = 0;
            // add dummy value names needed by GCC
            for (int var_ind = 0; var_ind < arity_in; var_ind++) {
                auto* variable = scope[var_ind];
                if (variable->isValueNames())
                    continue;
                for (int val_ind = 0; val_ind < (int)variable->getDomainInitSize(); val_ind++) {
                    variable->addValueName(string("v") + to_string(variable->toValue(val_ind)));
                }
            }
            conflictWeights.push_back(0);

            // Get the domain size of the first variable
            auto* variable = scope[varIndex];
            DomainSize = (int)variable->getDomainInitSize();
            VarDomainSize.push_back(DomainSize);

            // Collect values from the first variable's domain
            for (int valIndex = 0; valIndex < DomainSize; valIndex++) {
                if (seen.insert(variable->getValueName(valIndex)).second) {
                    TempUnionVarDomain.push_back(variable->getValueName(valIndex));
                    NbValues++;
                }
            }

            // Repeat for the rest of the variables
            for (int varIndex = 1; varIndex < arity_in; varIndex++) {
                conflictWeights.push_back(0);
                auto* variable = scope[varIndex];
                DomainSize = (int)variable->getDomainInitSize();
                VarDomainSize.push_back(DomainSize);
                if (VarDomainSize[varIndex - 1] != VarDomainSize[varIndex])
                    SameDomain = false;

                // Add only new unique values to the union
                for (int valIndex = 0; valIndex < DomainSize; valIndex++) {
                    if (seen.insert(variable->getValueName(valIndex)).second) {
                        TempUnionVarDomain.push_back(variable->getValueName(valIndex));
                        NbValues++;
                        SameDomain = false; // Different domains detected
                    }
                }
            }

            // Value names really present in the initial domain of each variable
            vector<unordered_set<string>> varNames(arity_in);
            for (int vi = 0; vi < arity_in; vi++) {
                for (int k = 0; k < (int)scope[vi]->getDomainInitSize(); k++) {
                    varNames[vi].insert(scope[vi]->getValueName(k));
                }
            }

            // Keep a value only if at least one variable can still take it with a finite unary cost.
            // Values that nobody can take (removed, or unary cost >= ub) are erased from the union.
            int tempNbValues = 0;
            Cost Max_cost = wcsp->getUb();
            for (int valIndex = 0; valIndex < NbValues; valIndex++) {
                const string& name = TempUnionVarDomain[valIndex];
                for (int vi = 0; vi < arity_in; vi++) {
                    if (!varNames[vi].count(name))
                        continue;
                    auto* var = scope[vi];
                    Value value = var->toValue(var->toIndex(name));
                    if (var->canbe(value) && (var->getCost(value) < Max_cost)) {
                        UnionVarDomain.push_back(name);
                        mapDomainValToIndex[name] = tempNbValues;
                        ++tempNbValues;
                        break;
                    }
                }
            }
            NbValues = tempNbValues;

            // Map each (variable, unified value) to the index of the value in the variable's own domain
            // (also recomputes VarDomainSize as the number of kept values of each variable)
            rebuildVarMapping();

            // If more values than variables, cost matrix is not square
            if (NbValues > arity_in) {
                isSquare = false;
                deltaCosts.assign(NbValues, MIN_COST);
            }

            // Percentage of variables used as Dijkstra sources in reduced-cost filtering (0 = disabled)
            int rate = ToulBar2::ReducedCostsFiltering;

            if (rate > 0) {
                FiltLevel = rate / 100.0;
            }

            // Test value symmetries
            //            for (unsigned int a = 0; a < NbValues; ++a) {
            //                for (unsigned int b = a+1; b < NbValues; ++b) {
            //                    if (valueSymmetry(UnionVarDomain[a], UnionVarDomain[b])) {
            //                        if (ToulBar2::verbose >= 1) {
            //                            cout << "detect value symmetry between " <<  UnionVarDomain[a] << " and " << UnionVarDomain[b] << endl;
            //                        }
            //                    }
            //                }
            //            }

            // Initialize

            // lower = static_cast<int>(arity_in * 0.5);
            // upper = static_cast<int>(arity_in * 0.8);
            //if (FiltLevel >= 1) {
                upper = arity_in;
            //}
            storeLastAssignment = vector<StoreValue>(arity_in, StoreValue(WRONG_VAL));
            NoAssignedVar = vector<int>(arity_in, -1);
            isExceptedVal = vector<uint8_t>(NbValues, 0);
            ValuesCapacity = vector<int>(NbValues, 1);
            AssignedVar = vector<int>(arity_in, -1);
            AssignedVal = vector<int>(arity_in, -1);
            costMatrix = vector<Cost>((arity_ + 1) * NbValues, MAX_COST);
            ReduceCostMatrix = vector<Cost>((arity_ + 1) * NbValues, MAX_COST);
            capacity = vector<int>(arity_ + 1, 1);

            // Allocate memory
            rowSol = new int[arity_ + 1];
            colSol = new int[NbValues];
            ReduceCostRow = new Cost[arity_ + 1];
            ReduceCostCol = new Cost[NbValues];
        } else {
            deconnect();
        }
    }

    // Destructor: frees the arrays used by the LAP solver
    virtual ~AllDifferentConstraint()
    {
        delete[] rowSol;
        delete[] colSol;
        delete[] ReduceCostRow;
        delete[] ReduceCostCol;
    }

    // Reads the parameters of the constraint from a WCSP file: excepted values, removed values and (non-square case) extended unary costs
    void read(istream& files) // TODO: add a parameter for controlling the level of propagation if necessary
    {
        string line;
        getline(files, line);
        vector<Value> TempexceptedValues;
        stringstream file(line);
        Cost Max_cost = wcsp->getUb();
        int nbExcepted = 0;

        // Excepted values: number followed by the list
        file >> nbExcepted;
        int tempNbExcepted = nbExcepted;
        for (int v = 0; v < nbExcepted; v++) {
            Value except;
            file >> except;
            TempexceptedValues.push_back(except);
        }

        // Values removed from the domains: number followed by the list
        int Nb_removeValues = 0;
        file >> Nb_removeValues;
        vector<Value> ListremoveVal;
        ListremoveVal.reserve(Nb_removeValues);
        for (int val = 0; val < Nb_removeValues; val++) {
            Value rmval;
            file >> rmval;
            ListremoveVal.push_back(rmval);
        }

        // Register the excepted values (their index is shifted by the number of removed values below them)
        for (int v = 0; v < (int)TempexceptedValues.size(); v++) {
            bool findval = false;
            Value except = TempexceptedValues[v];
            if (Nb_removeValues > 0) {
                int distance = 0;
                for (int r = 0; r < Nb_removeValues; r++) {
                    if (ListremoveVal[r] < except)
                        distance++;
                }
                except = except - distance;
            }

            // The excepted value is kept only if a variable can still take it
            for (int varIndex = 0; varIndex < arity_; ++varIndex) {
                auto* variable = scope[varIndex];
                if (variable->canbe(except) && (variable->getCost(except) < Max_cost)) {
                    int idx = unionIndex(variable->getValueName(variable->toIndex(except)));
                    if (idx < 0)
                        continue;
                    exceptedValIndex.push_back(idx);
                    exceptedValues.push_back(except);
                    isExceptedVal[idx] = 1;
                    ValuesCapacity[idx] = arity_;
                    findval = true;
                    break;
                }
            }
            if (!findval)
                tempNbExcepted--;
        }
        nbExcepted = tempNbExcepted;
        excepted = nbExcepted > 0;
        if (excepted) {
            isSquare = false;
            deltaCosts.assign(NbValues, MIN_COST);
            /*for (int varIndex = 0; varIndex < arity_; varIndex++) {
                deltaCosts.emplace_back(NbValues, MIN_COST);
            }*/
        } else {
            // contradiction
            if (NbValues < arity_)
                THROWCONTRADICTION;
        }

        // Choose the slack (dummy variable) formulation of the LAP if requested
        int slackVar = ToulBar2::AlldiffSlackVar;
        int ri = ToulBar2::RelativeInteriorSolution;
        if (excepted)
            FiltLevel = 0;


        if ((slackVar > 0 && !isSquare && !excepted) || (!isSquare && FiltLevel > 0) || (!isSquare && ri > 0) ) {
            isSlack = 1;
            slack = true;
        }

        // Relative-interior dual solution: square case, or slack formulation
        if (ri > 0 && (isSquare || slack)) {
            relativeInterior = true;
        }

        // Extended unary costs (non-square case): number followed by (value, cost) pairs
        if (!isSquare) {
            int nbDelta = 0;
            file >> nbDelta;
            if (nbDelta) {
                for (int v = 0; v < nbDelta; v++) {
                    Value value;
                    Cost delta;
                    file >> value;
                    file >> delta;
                    if (Nb_removeValues > 0) {
                        int distance = 0;
                        for (int r = 0; r < Nb_removeValues; r++) {
                            if (ListremoveVal[r] < value)
                                distance++;
                        }
                        value = value - distance;
                    }
                    for (int varIndex = 0; varIndex < arity_; ++varIndex) {
                        auto* variable = scope[varIndex];
                        if (variable->canbe(value)) {
                            int idx = unionIndex(variable->getValueName(variable->toIndex(value)));
                            if (idx < 0)
                                continue;
                            deltaCosts[idx] = delta;
                            break;
                        }
                    }
                }
            }
        }
    }

    bool extension() const FINAL { return false; } // TODO: allows functional variable elimination but no other preprocessing
    bool isAllDiff() const FINAL { return true; } // this constraint is an alldifferent
    bool isAllDiffSquare() const FINAL { return isSquare; } // true if #variables == #values (permutation)

    // Returns the values exempt from the alldifferent constraint
    vector<Value> getExceptedValues() const FINAL { return exceptedValues; }

    // Cost of a forbidden tuple
    Cost getDefCost() FINAL { return MAX_COST; }

    // Global conflict weight (weighted degree heuristics)
    Long getConflictWeight() const override { return Constraint::getConflictWeight(); }
    // Conflict weight of the variable at position varIndex (its own weight + the constraint weight)
    Long getConflictWeight(int varIndex) const override
    {
        assert(varIndex >= 0);
        assert(varIndex < arity_);
        return conflictWeights[varIndex] + Constraint::getConflictWeight();
    }
    // Increases the conflict weights after a contradiction caused by constraint 'from'
    void incConflictWeight(Constraint* from) override
    {
        assert(from != NULL);
        if (from == this) {
            // The contradiction comes from this constraint
            if (getNonAssigned() == arity_ || deconnected()) {

                if (arity_ < (int)wcsp->numberOfVariables()) { // SdG: no need to increase weightedDegree of all problem variables
                    Constraint::incConflictWeight(1);
                }
            } else {
                // Increase the weight of the assigned variables
                for (int i = 0; i < arity_; i++) {
                    // for (int i = 0; i < numConflictVars; i++) {
                    auto* variable = scope[i];

                    if (!variable->unassigned()) {
                        // if (connected(i)) {
                        // conflictWeights[lastConflictVars[i]]++;
                        conflictWeights[i]++;
                    }
                }
            }
        } else if (deconnected()) {
            // The contradiction comes from another constraint: increase the weight of the variables in common
            for (int i = 0; i < from->arity(); i++) {
                int index = getIndex(from->getVar(i));
                if (index >= 0) { // the last conflict constraint may be derived from two binary constraints (boosting search), each one derived from an n-ary constraint with a scope which does not include parameter constraint from
                    assert(index < arity_);
                    conflictWeights[index]++;
                }
            }
        }
    }
    // Resets all conflict weights to zero
    void resetConflictWeight() override
    {
        conflictWeights.assign(conflictWeights.size(), 0);
        Constraint::resetConflictWeight();
    }

    /// \brief returns true if constraint always satisfied and has (less than) zero cost only
    // bool universal(Cost zero = MIN_COST) override;

    // Returns true if the binary constraint ctr is implied by this alldifferent
    // (i.e. ctr has zero cost on every pair of different values, or when one of them is excepted)
    bool implies(Constraint* ctr) FINAL
    {
        if (scopeIncluded(ctr) && ctr->isBinary()) {
            BinaryConstraint* bctr = (BinaryConstraint*)ctr;
            EnumeratedVariable* x = (EnumeratedVariable*)(bctr->getVar(0));
            EnumeratedVariable* y = (EnumeratedVariable*)(bctr->getVar(1));
            if (x->isValueNames() && y->isValueNames()) {
                // Values are compared by name
                for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
                    for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
                        if ((x->getValueName(x->toIndex(*iterx)) != y->getValueName(y->toIndex(*itery)) || (excepted && (find(exceptedValues.begin(), exceptedValues.end(), *iterx) != exceptedValues.end() || find(exceptedValues.begin(), exceptedValues.end(), *itery) != exceptedValues.end()))) && bctr->getCost(*iterx, *itery) > MIN_COST) {
                            return false;
                        }
                    }
                }
            } else {
                // Values are compared directly
                for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
                    for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
                        if ((*iterx != *itery || (excepted && find(exceptedValues.begin(), exceptedValues.end(), *iterx) != exceptedValues.end())) && bctr->getCost(*iterx, *itery) > MIN_COST) {
                            return false;
                        }
                    }
                }
            }
            return true;
        }
        return false;
    }

    // Projects the alldifferent (forbidden equal values) onto the binary constraint ctr included in its scope
    void projects(Constraint* ctr) FINAL
    {
        if (scopeIncluded(ctr) && ctr->isBinary()) {
            BinaryConstraint* bctr = (BinaryConstraint*)ctr;
            EnumeratedVariable* x = (EnumeratedVariable*)(bctr->getVar(0));
            EnumeratedVariable* y = (EnumeratedVariable*)(bctr->getVar(1));
            Cost mult_ub = (wcsp->getUb() < (MAX_COST / MEDIUM_COST)) ? (max(LARGE_COST, wcsp->getUb() * MEDIUM_COST)) : wcsp->getUb();
            if (x->isValueNames() && y->isValueNames()) {
                // Values are matched by name
                for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
                    string s = x->getValueName(x->toIndex(*iterx));
                    unsigned int yindex = y->toIndex(s);
                    Value yval = y->toValue(yindex);
                    if ((!excepted || (find(exceptedValues.begin(), exceptedValues.end(), *iterx) == exceptedValues.end() && find(exceptedValues.begin(), exceptedValues.end(), yval) == exceptedValues.end())) && y->canbe(yval) && !CUT(bctr->getCost(*iterx, yval), wcsp->getUb())) {
                        bctr->addcost(*iterx, yval, mult_ub - bctr->getCost(*iterx, yval));
                    }
                }
            } else {
                // Values are matched directly
                for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
                    if ((!excepted || find(exceptedValues.begin(), exceptedValues.end(), *iterx) == exceptedValues.end()) && y->canbe(*iterx) && !CUT(bctr->getCost(*iterx, *iterx), wcsp->getUb())) {
                        bctr->addcost(*iterx, *iterx, mult_ub - bctr->getCost(*iterx, *iterx));
                    }
                }
            }
        }
    }

    // Returns true if the two values (given by name) are symmetric for the problem
    // (same unary costs and same binary costs for all variables of the scope). Only for the square case.
    bool valueSymmetry(string stra, string strb)
    {
        if (!isSquare) {
            return false;
        }
        for (int j = 0; j < arity_; j++) {
            EnumeratedVariable* xj = (EnumeratedVariable*)getVar(j);
            unsigned int ida = xj->toIndex(stra);
            Value xja = xj->toValue(ida);
            unsigned int idb = xj->toIndex(strb);
            Value xjb = xj->toValue(idb);
            // Both values must be available and have the same unary cost
            if (xj->cannotbe(xja) || xj->cannotbe(xjb) || xj->getCost(xja) != xj->getCost(xjb)) {
                return false;
            }
            // Check the binary constraints of xj
            if (xj->unassigned() && xj->getDegree() > 1) {
                ConstraintList* constrsj = xj->getConstrs();
                for (ConstraintList::iterator it = constrsj->begin(); it != constrsj->end(); ++it) {
                    Constraint* ctr = (*it).constr;
                    if (ctr->isBinary() && !ctr->isSep() && scopeIncluded(ctr)) {
                        BinaryConstraint* cjk = (BinaryConstraint*)ctr;
                        EnumeratedVariable* xk = (EnumeratedVariable*)((cjk->getVar(0) == xj) ? cjk->getVar(1) : cjk->getVar(0));
                        unsigned int ida = xk->toIndex(stra);
                        Value xka = xk->toValue(ida);
                        unsigned int idb = xk->toIndex(strb);
                        Value xkb = xk->toValue(idb);
                        for (EnumeratedVariable::iterator iterk = xk->begin(); iterk != xk->end(); ++iterk) {
                            if (*iterk != xka && *iterk != xkb && cjk->getCost(xj, xk, xja, *iterk) != cjk->getCost(xj, xk, xjb, *iterk)) {
                                return false;
                            }
                        }
                    } else if (!ctr->isAllDiff()) {
                        // Any other kind of constraint breaks the symmetry
                        return false;
                    }
                }
            }
        }
        return true;
    }

    // Returns the cost of the tuple s (0 if it is a valid alldifferent assignment, else the upper bound or a relaxed cost)
    Cost eval(const Tuple& s) override
    {
        // returns the cost of the corresponding assignment s

        if (isSquare) {
            Cost res = MIN_COST;
            Cost nbsame = 0; // number of violations (forbidden values or repeated values)

            vector<uint8_t> alreadyUsed(arity_, 0);
            for (int varIndex = 0; varIndex < arity_; varIndex++) {
                int valIndex = unionIndex(scope[varIndex]->getValueName(s[varIndex]));
                if (valIndex < 0) { // forbidden value (removed from the union)
                    nbsame++;
                    continue;
                }

                if (alreadyUsed[valIndex]) {
                    nbsame += 1;

                } else {
                    alreadyUsed[valIndex] = 1;
                }
            }
            if (nbsame > 0) {

                if (nbsame > 0 && Original_ub < wcsp->getUb() && 1.0L * Original_ub * nbsame < wcsp->getUb()) {
                    res = Original_ub * nbsame; // VNS-like methods may exploit a relaxation of the constraint
                } else {
                    res = wcsp->getUb();
                }
            }
            assert(res <= wcsp->getUb());
            assert(res >= MIN_COST);
            return res;
        }

        // Non-square case: the cost also includes the extended unary costs of the unused values
        // Cost res = -lb + assigneddeltas;
        Cost res = 0;
        Cost nbsame = 0; // number of violations
        vector<int> alreadyUsed(NbValues, 0); // number of variables taking each value
        for (int varIndex = 0; varIndex < arity_; varIndex++) {
            int valIndex = unionIndex(scope[varIndex]->getValueName(s[varIndex]));
            if (valIndex < 0) { // forbidden value (removed from the union)
                nbsame++;
                continue;
            }

            if (alreadyUsed[valIndex]) {
                int it = excepted == 0 ? 0 : isExceptedVal[valIndex];
                if (it == 0) {
                    nbsame += 1;
                }
            }
            alreadyUsed[valIndex] += 1;
        }
        if (nbsame == 0) {
            // Valid assignment: add the extended unary cost of each unused value (per remaining capacity if excepted values exist)
            if (excepted) {
                for (int valIndex = 0; valIndex < NbValues; valIndex++) {
                    res += (ValuesCapacity[valIndex] - alreadyUsed[valIndex]) * deltaCosts[valIndex];
                }
            } else {
                for (int valIndex = 0; valIndex < NbValues; valIndex++) {
                    if (alreadyUsed[valIndex] == 0)
                        res += deltaCosts[valIndex];
                }
            }
        }

        if (nbsame > 0 || res > wcsp->getUb()) {

            if (nbsame > 0 && Original_ub < wcsp->getUb() && 1.0L * Original_ub * nbsame < wcsp->getUb()) {
                res = Original_ub * nbsame; // VNS-like methods may exploit a relaxation of the constraint
            } else {
                res = wcsp->getUb();
            }
        }
        assert(res <= wcsp->getUb());
        assert(res >= MIN_COST);
        return res;
    }

    Cost getMaxFiniteCost() override ///< \brief returns the maximum finite cost for any valid tuple less than wcsp->getUb()
    {
        if (isSquare && !excepted)
            return MIN_COST;
        Cost sumdelta = -lb + assigneddeltas;

        Cost m = *max_element(deltaCosts.begin(), deltaCosts.end());
        if (m > MIN_COST)
            sumdelta += m * (NbValues - arity_);

        if (CUT(sumdelta, wcsp->getUb()))
            return MAX_COST;
        else
            return sumdelta;
    }

    // Returns the extended unary cost of value val for the variable at position index (always 0 in the square case)
    Cost getCost(int index, Value val)
    {
        assert(index >= 0 && index < arity_);
        if (isSquare)
            return MIN_COST;
        int valIndex = unionIndex(scope[index]->getValueName(scope[index]->toIndex(val)));
        return (valIndex < 0) ? MIN_COST : (Cost)deltaCosts[valIndex];
    }

    double computeTightness() override { return 0; } // TODO: compute factorial(n)/(n**n)?

    // TODO: needed for dominance test by DEE
    // pair<pair<Cost, Cost>, pair<Cost, Cost>> getMaxCost(int index, Value a, Value b)

    // void setInfiniteCost(Cost ub)
    // Updates the initial upper bound (used to compute the cost of a violation)
    void setInfiniteCost(Cost ub) override
    {
        Original_ub = min(ub, Original_ub);
    }

    // Called when the variable at position varIndex is assigned: disconnects it and either projects the
    // constraint into an n-ary cost function (few variables left) or propagates it
    void assign(int varIndex) override
    {
        if (ToulBar2::verbose >= 7) {
            cout << "assign " << scope[varIndex]->getName() << " in " << *this << endl;
        }
        if (connected(varIndex)) {
            deconnect(varIndex);
            assert(getNonAssigned() >= 0);
            if (getNonAssigned() <= NARYPROJECTIONSIZE && (getNonAssigned() <= 1 || prodInitDomSize <= NARYPROJECTIONPRODDOMSIZE || maxInitDomSize <= NARYPROJECTION3MAXDOMSIZE || (getNonAssigned() == 2 && maxInitDomSize <= NARYPROJECTION2MAXDOMSIZE))) {
                deconnect();
                projectNary();
            } else {
                // TODO: incremental bound propagation
                // remove the assigned value in other domains
                Value value = scope[varIndex]->getValue();
                int removeval = true;
                if(excepted) {
                    auto* variable = scope[varIndex];
                    int valIndex = unionIndex(variable->getValueName(variable->toIndex(value)));
                    if(isExceptedVal[valIndex]){
                        // don't remove an excepted value
                        removeval = false;
                    }
                }
                if(removeval){
                    for (int i = 0; i < arity_; i++) {
                        if (i != varIndex && scope[i]->unassigned() && scope[i]->canbe(value)) {
                            scope[i]->remove(value);
                        }
                    }
                }
                propagate();
                if (ToulBar2::FullEAC)
                    reviseEACGreedySolution();
            }
        }
    }

    /**
     * @brief Remove assigned values from unassigned
     *        variable domains, detecting contradictions, and triggering n-ary
     *        projection when thresholds are met.
     * @return true if propagation succeeded without n-ary projection, false otherwise.
     * @throws THROWCONTRADICTION if two variables share the same assigned value.
     */
    bool RemoveAssignVar()
    {
        // Reset tracking vectors for this propagation pass
        isAssignedValue = vector<bool>(NbValues, false);
        varAlreadyProcessed = vector<bool>(arity_, false);

        AssignedVar.clear();
        AssignedVar.reserve(arity_);
        AssignedVal.clear();
        AssignedVal.reserve(arity_);
        NoAssignedVar.clear();
        NoAssignedVar.reserve(arity_);
        unordered_map<int, int> conflictvar; // Maps value index -> variable index for contradiction detection

        // Partition variables into assigned / unassigned
        for (int varIndex = 0; varIndex < arity_; ++varIndex) {
            auto* variable = scope[varIndex];

            if (variable->unassigned()) {
                NoAssignedVar.push_back(varIndex);
            } else {
                // Variable assigned: find the index of its assigned value
                int valIndex = unionIndex(variable->getValueName(variable->toIndex(variable->getValue())));
                if (valIndex < 0) { // assigned to a forbidden value
                    conflictWeights[varIndex]++;
                    THROWCONTRADICTION;
                }
                int it = excepted == 0 ? 0 : isExceptedVal[valIndex];
                if (!isAssignedValue[valIndex]) {
                    AssignedVar.push_back(varIndex);
                    if (it == 0) {
                        isAssignedValue[valIndex] = true; // Mark value as taken
                        AssignedVal.push_back(valIndex);
                        conflictvar[valIndex] = varIndex; // Record owner for conflict detection
                    }
                } else {
                    // Two variables share the same non-excepted value: contradiction
                    if (it == 0) {
                        int varIndex2 = conflictvar[valIndex];
                        conflictWeights[varIndex2]++;
                        conflictWeights[varIndex]++;

                        THROWCONTRADICTION;
                    }
                }
            }
        }

        NbAssigned = AssignedVar.size();
        bool VarAssigned;
        bool NaryPro = false;

        // Iteratively remove assigned values from unassigned variable domains
        do {
            VarAssigned = false;
            vector<int> newlyAssignedVars; // Variables that become assigned during this iteration

            for (int valIndex = 0; valIndex < NbValues; ++valIndex) {
                if (!isAssignedValue[valIndex])
                    continue; // Only process taken values

                for (int varIndex : NoAssignedVar) {
                    if (varAlreadyProcessed[varIndex])
                        continue;

                    auto* variable = scope[varIndex];
                    Value value = toVarValue(varIndex, valIndex);

                    // Remove value if present in domain and not excepted
                    if (value != WRONG_VAL && variable->canbe(value) && (find(exceptedValIndex.begin(), exceptedValIndex.end(), valIndex) == exceptedValIndex.end())) {
                        variable->remove(value);

                        if (variable->assigned()) {
                            // Variable became singleton: schedule for propagation
                            newlyAssignedVars.push_back(varIndex);
                            varAlreadyProcessed[varIndex] = true;
                        }
                    }
                }
            }

            // Integrate newly assigned variables into the constraint state
            for (int varIndex : newlyAssignedVars) {
                if (connected(varIndex)) {
                    deconnect(varIndex);
                    NbNoAssigned--;
                    NbAssigned++;

                    // Switch to full n-ary projection if size thresholds are reached
                    if (getNonAssigned() <= NARYPROJECTIONSIZE && (getNonAssigned() <= 1 || prodInitDomSize <= NARYPROJECTIONPRODDOMSIZE || maxInitDomSize <= NARYPROJECTION3MAXDOMSIZE || (getNonAssigned() == 2 && maxInitDomSize <= NARYPROJECTION2MAXDOMSIZE))) {
                        deconnect();
                        projectNary();
                        NaryPro = true;
                        break;
                    } else {
                        VarAssigned = true; // Signal another propagation iteration
                        auto* variable = scope[varIndex];
                        // Variable assigned: find the index of its assigned value
                        int valIndex = unionIndex(variable->getValueName(variable->toIndex(variable->getValue())));
                        if (valIndex < 0) { // assigned to a forbidden value
                            conflictWeights[varIndex]++;
                            THROWCONTRADICTION;
                        }
                        int it = excepted == 0 ? 0 : isExceptedVal[valIndex];
                        if (it == 0) {
                            isAssignedValue[valIndex] = true;
                            AssignedVal.push_back(valIndex);
                        }
                        AssignedVar.push_back(varIndex);
                    }
                } else {
                    return false; // Variable disconnected externally: abort
                }
            }

            // Rebuild unassigned variable list after each iteration
            if (VarAssigned && !NaryPro) {
                vector<int> updatedNoAssigned;
                for (int varIndex : NoAssignedVar) {
                    if (!scope[varIndex]->assigned()) {
                        updatedNoAssigned.push_back(varIndex);
                    }
                }
                NoAssignedVar.swap(updatedNoAssigned);
            }

        } while (VarAssigned);

        return (!NaryPro);
    }

    /**
     * @brief Bimodal Dijkstra on the forward residual graph.
     *        Computes shortest distances from a source variable to all other
     *        variable nodes. At each settled node, chooses between iterating
     *        over the (sparse) adjacency list or the (small) unsettled set,
     *        whichever is cheaper.
     * @param source    Index of the source variable node.
     * @param VarList   VarList[val] = list of variable indices that have val as a non-matching arc.
     * @param ValList   ValList[var][val] = 1 if arc (var,val) exists in the residual graph.
     * @param SlackAssignedVals Values matched to the slack variable (slack formulation only).
     *
     *      Source : Bimodal Depth-First Search for Scalable GAC for AllDifferent.
     *      Sulian Le Bozec Chiffoleau; Nicolas Beldiceanu; Charles Prud'homme;
     *      Gilles Simonin; and Xavier Lorca, roceedings of the Thirty-Fourth
     *      International Joint Conference on Artificial Intelligence, IJCAI 2025,
     *      Montreal, Canada, August 16-22, 2025.
     */
    void BimodalDijkstra(int source, vector<vector<int>>& VarList, vector<vector<uint8_t>>& ValList, const vector<int>& SlackAssignedVals = {})
    {
        const int SlackVarInd = NbNoAssigned; // index of the dummy slack variable
        const int Nbvars = NbNoAssigned + (isSlack ? 1 : 0);

        if (dijkstraHandles.size() < (size_t)Nbvars)
            dijkstraHandles.resize(Nbvars);

        auto& handles = dijkstraHandles;
        PQ pq;

        distanceToVar.assign(Nbvars, MAX_COST);
        visited.assign(Nbvars, 0);

        // trackList = dense list of the nodes not yet settled, trackPos = position of each node in trackList (-1 if settled)
        trackList.resize(Nbvars);
        trackPos.resize(Nbvars);
        iota(trackList.begin(), trackList.end(), 0);
        iota(trackPos.begin(), trackPos.end(), 0);

        // Removes a node from trackList in O(1) (swap with the last element)
        auto trackErase = [&](int var) {
            int p = trackPos[var];
            if (p < 0)
                return;
            int last = trackList.back();
            trackList[p] = last;
            trackPos[last] = p;
            trackList.pop_back();
            trackPos[var] = -1;
        };

        inHeap.assign(Nbvars, 0);
        distanceToVar[source] = 0;
        handles[source] = pq.push({ 0, source });
        inHeap[source] = 1;

        // Relaxes the arc (val -> nextVar) from a node at distance d. Returns true if the distance of nextVar improved.
        auto relax = [&](int nextVar, Cost d, int val) -> bool {
            Cost alt = d + ReduceCostMatrix[nextVar * NbNoAssignedVal + val];
            if (alt < distanceToVar[nextVar]) {
                distanceToVar[nextVar] = alt;
                if (!inHeap[nextVar]) {
                    handles[nextVar] = pq.push({ alt, nextVar });
                    inHeap[nextVar] = 1;
                } else {
                    pq.update(handles[nextVar], { alt, nextVar });
                }
                return true;
            }
            return false;
        };

        // Expands the non-matching arcs of value val from a settled node at distance d.
        // Bimodal: iterates over the sparse adjacency list VarList[val] if it is smaller than the set of
        // unsettled nodes, otherwise over the unsettled nodes checking ValList.
        // If settleZeroCost is true, a node reached through a zero reduced-cost arc is settled immediately.
        auto expandValue = [&](Cost d, int val, bool settleZeroCost) {
            auto& adj = VarList[val];

            if (adj.size() < trackList.size()) {
                // Sparse mode
                for (int nextVar : adj) {
                    if (visited[nextVar])
                        continue;
                    if (relax(nextVar, d, val) && settleZeroCost && ReduceCostMatrix[nextVar * NbNoAssignedVal + val] == 0) {
                        trackErase(nextVar);
                        visited[nextVar] = 1;
                    }
                }
            } else {
                // Dense mode
                size_t k = 0;
                while (k < trackList.size()) {
                    int nextVar = trackList[k];
                    if (ValList[nextVar][val]) {
                        if (relax(nextVar, d, val) && settleZeroCost && ReduceCostMatrix[nextVar * NbNoAssignedVal + val] == 0) {
                            trackErase(nextVar);
                            visited[nextVar] = 1;
                            continue; // trackList[k] has been replaced by another node
                        }
                    }
                    ++k;
                }
            }
        };

        while (!pq.empty()) {
            auto [d, var] = pq.top();
            pq.pop();

            visited[var] = 1;
            trackErase(var);

            if (isSlack && var == SlackVarInd) {
                // The slack variable is matched to several values
                for (int val : SlackAssignedVals)
                    expandValue(d, val, false);
            } else {
                // Matching arc: traverse to the matched value
                expandValue(d, rowSol[var], true);
            }
        }
    }

    /**
     * @brief Bimodal Dijkstra on the reverse residual graph.
     *        Computes shortest distances from the value matched to the source
     *        variable back to all other value nodes, used by the Régin landmark
     *        upper-bound test.
     * @param source  Index of the source variable (its matched value is the start node).
     * @param ValList ValList[var][val] = 1 if arc (var,val) exists in the residual graph.
     * @param SlackAssignedVals Values matched to the slack variable (slack formulation only).
     */
    void BimodalDijkstraReverseGraph(int source, vector<vector<uint8_t>>& ValList, const vector<int>& SlackAssignedVals = {})
    {
        const int SlackVarInd = NbNoAssigned;
        const int Nbvars = NbNoAssigned + (isSlack ? 1 : 0);

        if (dijkstraHandles.size() < (size_t)max(Nbvars, NbNoAssignedVal))
            dijkstraHandles.resize(max(Nbvars, NbNoAssignedVal));

        auto& handles = dijkstraHandles;
        PQ pq;

        distanceToVal.assign(NbNoAssignedVal, MAX_COST);
        visited.assign(NbNoAssignedVal, 0);
        inHeap.assign(NbNoAssignedVal, 0);
        int visitSlack = -1; // first value matched to the slack variable which has been settled

        // Build inverse matching
        if (!isSlack)
            for (int var = 0; var < NbNoAssigned; ++var)
                colSol[rowSol[var]] = var;

        // Initialize tracking structures
        trackList.resize(NbNoAssignedVal);
        trackPos.resize(NbNoAssignedVal);
        iota(trackList.begin(), trackList.end(), 0);
        iota(trackPos.begin(), trackPos.end(), 0);

        // Removes a node from trackList in O(1) (swap with the last element)
        auto trackErase = [&](int val) {
            int p = trackPos[val];
            if (p < 0)
                return;
            int last = trackList.back();
            trackList[p] = last;
            trackPos[last] = p;
            trackList.pop_back();
            trackPos[val] = -1;
        };

        // Initialize source
        if (isSlack && source == SlackVarInd) {
            for (int val : SlackAssignedVals) {
                distanceToVal[val] = 0;
                handles[val] = pq.push({ 0, val });
                inHeap[val] = 1;
            }
        } else {
            int val = rowSol[source];
            distanceToVal[val] = 0;
            handles[val] = pq.push({ 0, val });
            inHeap[val] = 1;
        }

        while (!pq.empty()) {
            auto [d, val] = pq.top();
            pq.pop();

            // All the values matched to the slack variable share the same distance: only expand the first one
            if (colSol[val] == SlackVarInd) {
                if (visitSlack > -1) {
                    distanceToVal[val] = distanceToVal[visitSlack];
                    continue;
                } else {
                    visitSlack = val;
                }
            }

            visited[val] = 1;
            trackErase(val);

            int var = colSol[val];

            // Dense traversal over unsettled values
            size_t k = 0;
            while (k < trackList.size()) {
                int nextVal = trackList[k];

                if (ValList[var][nextVal]) {
                    Cost alt = d + ReduceCostMatrix[var * NbNoAssignedVal + nextVal];

                    if (alt < distanceToVal[nextVal]) {
                        distanceToVal[nextVal] = alt;

                        if (!inHeap[nextVal]) {
                            handles[nextVal] = pq.push({ alt, nextVal });
                            inHeap[nextVal] = 1;
                        } else {
                            pq.update(handles[nextVal], { alt, nextVal });
                        }

                        if (alt == d) { // Zero-cost arc: settle immediately
                            trackErase(nextVal);
                            visited[nextVal] = 1;
                            continue;
                        }
                    }
                }
                ++k;
            }
        }
    }

    /**
     * @brief Standard Dijkstra on the forward residual graph (sparse variant).
     *        Used when the bimodal strategy is not needed (e.g. dense graphs
     *        where VarList is always smaller than trackingList).
     * @param source  Index of the source variable node.
     * @param VarList VarList[val] = list of variable indices reachable via value val.
     */
    void Dijkstra(int source, vector<vector<int>>& VarList)
    {
        vector<PQ::handle_type> handles(NbNoAssigned);
        PQ pq;

        distanceToVar.assign(NbNoAssigned, MAX_COST);
        visited.assign(NbNoAssigned, 0);
        inHeap.assign(NbNoAssigned, 0);
        distanceToVar[source] = 0;
        handles[source] = pq.push({ 0, source });
        inHeap[source] = 1;

        while (!pq.empty()) {
            auto [d, var] = pq.top();
            pq.pop();
            visited[var] = 1;

            int val = rowSol[var]; // Matching arc: traverse to matched value
            auto& varlist = VarList[val]; // Variables reachable from that value (non-matching arcs)

            for (int nextVar : varlist) {
                if (visited[nextVar])
                    continue;
                Cost alt = d + ReduceCostMatrix[nextVar * NbNoAssignedVal + val];
                if (alt < distanceToVar[nextVar]) {
                    distanceToVar[nextVar] = alt;
                    if (!inHeap[nextVar]) {
                        handles[nextVar] = pq.push({ alt, nextVar });
                        inHeap[nextVar] = 1;
                    } else {
                        pq.update(handles[nextVar], { alt, nextVar });
                    }
                    if (alt == d) // Zero-cost arc: settle immediately without re-queuing
                        visited[nextVar] = 1;
                }
            }
        }
    }

    /**
     * @brief Régin landmark upper-bound test.
     *        Runs Bimodal Dijkstra in both directions from the source variable,
     *        then checks whether the sum of the two maximum distances is within
     *        the remaining cost budget H. If so, no arc in the residual graph
     *        can be pruned and filtering can stop early.
     * @param source  Source variable index used as the landmark.
     * @param VarList VarList[val] = variables reachable via value val (forward graph).
     * @param ValList ValList[var][val] = 1 if arc (var,val) exists in residual graph.
     * @param H       Remaining cost budget (wcsp upper bound minus current lower bound).
     * @param SlackAssignedVals Values matched to the slack variable (slack formulation only).
     * @return true if all values are consistent (no pruning possible), false otherwise.
     */
    bool ReginLandmarkUpperBound(int source, vector<vector<int>>& VarList, vector<vector<uint8_t>>& ValList, Cost H, const vector<int>& SlackAssignedVals = {})
    {
        // Forward pass: distances from source to all variable nodes
        BimodalDijkstra(source, VarList, ValList, SlackAssignedVals);
        // Backward pass: distances from source's matched value to all value nodes
        BimodalDijkstraReverseGraph(source, ValList, SlackAssignedVals);

        // An unreachable node (infinite distance) means we cannot conclude: avoid MAX_COST + MAX_COST overflow
        Cost maxDistanceFromVar = 0;
        for (Cost d : distanceToVar) {
            if (d >= MAX_COST)
                return false;
            maxDistanceFromVar = max(maxDistanceFromVar, d);
        }
        Cost maxDistanceToVar = 0;
        for (Cost d : distanceToVal) {
            if (d >= MAX_COST)
                return false;
            maxDistanceToVar = max(maxDistanceToVar, d);
        }

        // If max_forward + max_backward <= H - MaxReducedCost, no arc exceeds the budget
        bool AllValConsistent = maxDistanceFromVar + maxDistanceToVar <= H - MaxReducedCost;
        return AllValConsistent;
    }

    // Make the stored value the EAC support of the variable
    void updateSupport(int varIndex)
    {
        auto* variable = scope[varIndex];
        Value optimalValue = storeLastAssignment[varIndex];
        if (variable->getSupport() != optimalValue) {
            if (ToulBar2::verbose > 0)
                cout << "CHANGE ALLDIFF SUPPORT " << variable->getName() << " from " << variable->getSupport() << " to " << optimalValue << endl;
#ifndef NDEBUG
            variable->queueEAC1(); // EAC support may have been lost
#endif
            variable->setSupport(optimalValue);
            assert(variable->getCost(variable->getSupport()) == MIN_COST);
        }
    }


    // Updates the EAC support of every currently unassigned variable (NoAssignedVar) from the stored assignment
    void updateSupportsOfUnassigned()
    {
        for (int varInd = 0; varInd < NbNoAssigned; ++varInd)
            updateSupport(NoAssignedVar[varInd]);
    }

    // Updates the EAC support of every variable of the scope from the stored assignment
    void updateAllSupports()
    {
        for (int varIndex = 0; varIndex < arity_; ++varIndex)
            updateSupport(varIndex);
    }

    // Returns true if the previously stored optimal assignment is still valid
    // (all its values are still in the domains and have zero unary cost), so that propagation can be skipped
    bool canSkipPropagation()
    {
        if (!storeAssignment)
            return false;
        for (int varIndex = 0; varIndex < arity_; varIndex++) {
            if (scope[varIndex]->cannotbe(storeLastAssignment[varIndex]) || scope[varIndex]->getCost(storeLastAssignment[varIndex]) > MIN_COST)
                return false;
        }
        return true;
    }

    // Looks for a connected variable which is now assigned and calls assign() on it.
    // Returns true if one was found (assign() already propagates the constraint).
    bool assignFirstAssignedVariable()
    {
        for (int varIndex = 0; connected() && varIndex < arity_; varIndex++) {
            if (connected(varIndex) && scope[varIndex]->assigned()) {
                assign(varIndex);
                return true;
            }
        }
        return false;
    }

    // Stores the LAP solution (rowSol) of the unassigned variables, plus the values of the assigned variables
    // NoAssignedVal[c] = value index (in UnionVarDomain) of the column c of the reduced matrix
    void storePartialAssignment(const vector<int>& NoAssignedVal)
    {
        for (int varIndex = 0; varIndex < NbNoAssigned; ++varIndex) {
            storeLastAssignment[NoAssignedVar[varIndex]] = toVarValue(NoAssignedVar[varIndex], NoAssignedVal[rowSol[varIndex]]);
        }
        for (int varIndex = 0; varIndex < NbAssigned; ++varIndex) {
            storeLastAssignment[AssignedVar[varIndex]] = toVarValue(AssignedVar[varIndex], AssignedVal[varIndex]);
        }
        storeAssignment = true;
    }

    // Zero cost matrix with some assigned variables: any injective assignment of the unassigned variables is optimal
    void storePartialZeroCost(const vector<int>& NoAssignedVal)
    {
        for (int k = 0; k < NbNoAssigned; ++k)
            storeLastAssignment[NoAssignedVar[k]] = toVarValue(NoAssignedVar[k], NoAssignedVal[k]);
        // Assigned variables keep their value
        for (int k = 0; k < NbAssigned; ++k)
            storeLastAssignment[AssignedVar[k]] = toVarValue(AssignedVar[k], AssignedVal[k]);
        storeAssignment = true;
        updateSupportsOfUnassigned();
    }

    // Stores the LAP solution (rowSol) when all variables are unassigned
    void storeFullAssignment()
    {
        for (int varIndex = 0; varIndex < arity_; varIndex++) {
            storeLastAssignment[varIndex] = toVarValue(varIndex, rowSol[varIndex]);
        }
        storeAssignment = true;
    }

    // Zero cost matrix with no assigned variable: every variable can take every value at no cost
    void storeZeroCostFullAssignment()
    {
        for (int varIndex = 0; varIndex < arity_; ++varIndex) {
            int valIndex = (varIndex < NbValues) ? varIndex : exceptedValIndex[0];
            storeLastAssignment[varIndex] = toVarValue(varIndex, valIndex);
        }
        storeAssignment = true;
        updateAllSupports();
    }

    // Called when the LAP solver proves the constraint infeasible:
    // increases the conflict weight of the conflicting variables, revises the constraint and throws a contradiction.
    // sol: solution array returned by the LAP solver (rowSol or colSol), nbConflict: findConflict,
    // nbVars: entries of sol >= nbVars are dummy (not real variables), varMap: optional translation to scope indices.
    void onLAPFailure(const int* sol, int nbConflict, int nbVars, const vector<int>* varMap)
    {
        for (int i = 0; i < nbConflict; i++) {
            int v = sol[i];
            if (v >= nbVars)
                continue;
            conflictWeights[varMap ? (*varMap)[v] : v]++;
        }
        wcsp->revise(this);
        THROWCONTRADICTION;
    }

    /**
     * @brief Writes the dummy ("slack") row/column of the cost matrix from deltaCosts.
     *        deltaCosts are saturated at ub. If all the values that can be used are saturated
     *        (nb == nbVars) every other value is forbidden; if there are more, a contradiction is raised.
     * @param valIndexes Value index (in UnionVarDomain) of each matrix column.
     * @param nbVars     Number of real variables (the dummy variable has position nbVars).
     * @param ub         Current upper bound (wcsp->getUb() - wcsp->getLb()).
     * @param idx        idx(var, col) = position in costMatrix (encodes the layout of the matrix).
     */
    void setDeltaCostsInMatrix(const vector<int>& valIndexes, int nbVars, Cost ub, const function<size_t(int, int)>& idx)
    {
        int nbSaturated = 0;
        for (int col = 0; col < (int)valIndexes.size(); ++col) {
            int valIndex = valIndexes[col];
            if (deltaCosts[valIndex] >= ub) {
                deltaCosts[valIndex] = ub;
                nbSaturated++;
            }
            costMatrix[idx(nbVars, col)] = deltaCosts[valIndex];
        }
        if (nbSaturated > nbVars) {
            THROWCONTRADICTION;
        } else if (nbSaturated == nbVars) {
            for (int col = 0; col < (int)valIndexes.size(); ++col) {
                if (deltaCosts[valIndexes[col]] < ub) {
                    for (int varInd = 0; varInd < nbVars; ++varInd)
                        costMatrix[idx(varInd, col)] = ub;
                }
            }
        }
    }

    // Fills the cost matrix of the full problem (all variables x all values) with the unary costs, capped at ub
    void fillFullCostMatrix(Cost ub, const function<size_t(int, int)>& idx)
    {
        for (int varIndex = 0; varIndex < arity_; ++varIndex) {
            auto* variable = scope[varIndex];
            for (int valIndex = 0; valIndex < NbValues; ++valIndex) {
                Value value = toVarValue(varIndex, valIndex);
                Cost valCost = (value != WRONG_VAL && variable->canbe(value)) ? variable->getCost(value) : ub;
                costMatrix[idx(varIndex, valIndex)] = valCost < ub ? valCost : ub;
            }
        }
    }

    // Projects/extends the reduced costs (row dual + column dual) onto the unary costs of the full problem
    void projectReducedCostsFull(Cost ub, const function<size_t(int, int)>& idx)
    {
        for (int varIndex = 0; varIndex < arity_; ++varIndex) {
            for (int valIndex = 0; valIndex < NbValues; ++valIndex) {
                Value value = toVarValue(varIndex, valIndex);
                if (value == WRONG_VAL)
                    continue;
                if (costMatrix[idx(varIndex, valIndex)] < ub) {
                    ExtOrProJ(varIndex, value, (ReduceCostRow[varIndex] + ReduceCostCol[valIndex]));
                }
            }
        }
    }

       /**
     * @brief Partial case: projects the reduced costs onto the unary costs of the unassigned variables
     * @param NoAssignedVal Value index of each column of the reduced matrix.
     * @param current_ub    wcsp->getUb() - wcsp->getLb() before the LAP cost projection.
     * @param newCurrentUb  current_ub minus the LAP optimal cost (remaining budget).
     * @param idx           idx(var, col) = position in costMatrix.
     */
    void projectReducedCostsPartial(const vector<int>& NoAssignedVal, Cost current_ub,  const function<size_t(int, int)>& idx)
    {

        for (int varInd = 0; varInd < NbNoAssigned; ++varInd) {
            int varIndex = NoAssignedVar[varInd];
            for (int valInd = 0; valInd < NbNoAssignedVal; ++valInd) {
                size_t pos = idx(varInd, valInd);
                if (costMatrix[pos] < current_ub) {
                    int valIndex = NoAssignedVal[valInd];
                    Value value = toVarValue(varIndex, valIndex);
                    ExtOrProJ(varIndex, value, (ReduceCostRow[varInd] + ReduceCostCol[valInd]));

                }
            }
        }
    } 


    /**
     * @brief Partial case: projects the reduced costs onto the unary costs of the unassigned variables
     *        and, if reduced-cost filtering is enabled, builds the residual graph (VarList, ValList, ReduceCostMatrix, MaxReducedCost).
     * @param NoAssignedVal Value index of each column of the reduced matrix.
     * @param current_ub    wcsp->getUb() - wcsp->getLb() before the LAP cost projection.
     * @param newCurrentUb  current_ub minus the LAP optimal cost (remaining budget).
     * @param idx           idx(var, col) = position in costMatrix.
     */
    void projectAndBuildResidualGraph(const vector<int>& NoAssignedVal, Cost current_ub, Cost newCurrentUb, const function<size_t(int, int)>& idx, vector<vector<int>>& VarList, vector<vector<uint8_t>>& ValList)
    {
        MaxReducedCost = 0;
        const bool buildGraph = residualGraphEnabled();
        for (int varInd = 0; varInd < NbNoAssigned; ++varInd) {
            int varIndex = NoAssignedVar[varInd];
            ValList[varInd] = vector<uint8_t>(NbNoAssignedVal, 0);
            for (int valInd = 0; valInd < NbNoAssignedVal; ++valInd) {
                size_t pos = idx(varInd, valInd);
                if (costMatrix[pos] < current_ub) {
                    int valIndex = NoAssignedVal[valInd];
                    Value value = toVarValue(varIndex, valIndex);
                    ExtOrProJ(varIndex, value, (ReduceCostRow[varInd] + ReduceCostCol[valInd]));

                    if (buildGraph) {
                        ReducedCost = costMatrix[pos] - (ReduceCostRow[varInd] + ReduceCostCol[valInd]);
                        if (ReducedCost < newCurrentUb) {
                            ReduceCostMatrix[varInd * NbNoAssignedVal + valInd] = ReducedCost;
                            if (ReducedCost > MaxReducedCost)
                                MaxReducedCost = ReducedCost;
                            if (rowSol[varInd] != valInd) { // non-matching arc
                                VarList[valInd].push_back(varInd);
                                ValList[varInd][valInd] = 1;
                            }
                        }
                    }
                }
            }
        }
    }

    // Slack formulation: adds the arcs of the dummy (slack) variable to the residual graph
    // and collects the values matched to the slack variable in SlackAssignedVals
    void addSlackArcsToResidualGraph(Cost Current_ub, Cost newCurrentUb, const function<size_t(int, int)>& idx, vector<vector<int>>& VarList, vector<vector<uint8_t>>& ValList, vector<int>& SlackAssignedVals)
    {
        ValList[NbNoAssigned] = vector<uint8_t>(NbNoAssignedVal, 0);
        for (int valInd = 0; valInd < NbNoAssignedVal; ++valInd) {
            if (colSol[valInd] == NbNoAssigned) {
                SlackAssignedVals.push_back(valInd);
            } else {
                if(costMatrix[idx(NbNoAssigned, valInd)] < Current_ub){
                    ReducedCost = costMatrix[idx(NbNoAssigned, valInd)] - (ReduceCostRow[NbNoAssigned] + ReduceCostCol[valInd]);
                    if (ReducedCost < newCurrentUb) {
                        ReduceCostMatrix[NbNoAssigned * NbNoAssignedVal + valInd] = ReducedCost;
                        MaxReducedCost = max(ReducedCost, MaxReducedCost);
                        VarList[valInd].push_back(NbNoAssigned);
                        ValList[NbNoAssigned][valInd] = 1;
                    }
                }
            }
        }
    }


    // Recursive Tarjan algorithm for Strongly Connected Components – Wikipedia Version
    //source : https://courses.cms.caltech.edu/cs130/lectures-2024wi/CS130-Wi2024-Lec19.pdf

    void strongconnect(int v, const vector<vector<int>>& g, vector<int>& index, vector<int>& lowlink, vector<int>& comp, vector<uint8_t>& onStack, vector<int>& stk,int& counter, int& nbComp)
    {
        // Set index, initial lowlinklink for v
        index[v]  = counter;
        lowlink[v]  = counter;
        counter++;

        // Record v on the stack
        stk.push_back(v);
        onStack[v] = 1;

        // Consider neighbors of v to compute v.lowlink
        for (int w : g[v]) {

            if (index[w] == -1) {

                // Successor w has not yet been visited
                strongconnect(
                    w, g, index, lowlink, comp,
                    onStack, stk, counter, nbComp
                );

                // Update lowlink-link after returning from recursion
                lowlink[v] = min(lowlink[v], lowlink[w]);
            }
            else if (onStack[w]) {

                // Successor w has been visited, and is
                // also on stack S and is therefore in
                // the current SCC.
                lowlink[v] = min(lowlink[v], index[w]);
            }
        }

        // If v is a root of an SCC, pop vertices off
        // the stack to generate/record the SCC.
        if (lowlink[v] == index[v]) {

            int w;

            do {
                w = stk.back();
                stk.pop_back();

                onStack[w] = 0;
                comp[w] = nbComp;

            } while (w != v);
            //store or the current SCC 
            nbComp++;
        }
    }


    /**
     * @brief GAC by strongly connected components (Régin) on the residual graph of the LAP, standard or slack formulation.
     *        Step 1 of reducedCostFiltering(): structural criterion (alternating cycles) instead of reduced-cost shortest paths.
     *
     *        Nodes: real variables 0..NbNoAssigned-1, the slack variable NbNoAssigned (only if isSlack), then the columns.
     *        Arcs:  - matched:      var -> rowSol[var]               (real variables)
     *                               slack -> val for val in SlackAssignedVals
     *               - non-matching: val -> var for var in VarList[val] (the slack variable has index NbNoAssigned)
     *
     *        The slack variable stands for NbNoAssignedVal - NbNoAssigned copies, each matched to one column and adjacent to
     *        every column. Merging the copies into one node requires an arc val -> slack from each slack-matched column
     *        when there are at least two copies. Every column is matched (to a variable or to the slack), so there is no
     *        free column: a non-matching arc (x, v) is consistent iff x and v belong to the same SCC. Other arcs are pruned.
     *        Pruning is sound: only arcs of reduced cost < newCurrentUb are in the graph, and any solution within the
     *        budget uses only such arcs.
     *
     * @param NoAssignedVal     Value index (in UnionVarDomain) of each column.
     * @param VarList           VarList[val] = variables having a non-matching arc to val (pruned arcs are erased).
     * @param ValList           ValList[var][val] = 1 if the arc is in the residual graph (pruned arcs are reset to 0).
     * @param SlackAssignedVals Columns matched to the slack variable (empty when isSlack == 0).
     * @param newCurrentUb      Remaining cost budget (infinite cost projected on removed values).
     * @return true if at least one arc has been pruned.
     */
    bool gacFilteringSCC(const vector<int>& NoAssignedVal, vector<vector<int>>& VarList, vector<vector<uint8_t>>& ValList, const vector<int>& SlackAssignedVals, Cost newCurrentUb)
    {
        if (!residualGraphEnabled())
            return false;

        const int nVar = NbNoAssigned + (isSlack ? 1 : 0); // variable nodes (real + slack)
        const int slackNode = NbNoAssigned;
        const int m = NbNoAssignedVal; // column nodes: nVar .. nVar + m - 1
        const int N = nVar + m;

        // --- Directed graph
        vector<vector<int>> g(N);
        for (int var = 0; var < NbNoAssigned; ++var)
            g[var].push_back(nVar + rowSol[var]); // matched arc of a real variable
        if (isSlack) {
            for (int val : SlackAssignedVals)
                g[slackNode].push_back(nVar + val); // matched arcs of the slack variable
            if (SlackAssignedVals.size() >= 2) // at least two copies: they reach each other through the columns
                for (int val : SlackAssignedVals)
                    g[nVar + val].push_back(slackNode);
        }
        for (int val = 0; val < m; ++val)
            for (int var : VarList[val])
                g[nVar + val].push_back(var); // non-matching arcs (including those of the slack variable)


        // Strongly Connected Components (Recursive Tarjan)

        vector<int> index(N, -1);
        vector<int> lowlink(N, 0);
        vector<int> comp(N, -1);

        vector<uint8_t> onStack(N, 0);
        vector<int> stk;

        int counter = 0;
        int nbComp = 0;

        for (int root = 0; root < N; ++root) {

            if (index[root] == -1) {

                strongconnect(root, g, index, lowlink, comp, onStack, stk, counter, nbComp);
            }
        }

        // --- Pruning of the non-matching arcs between different components 
        bool removed = false;
        for (int valInd = 0; valInd < m; ++valInd) {
            auto& varlist = VarList[valInd];
            bool valremoved = false;
            for (int& varInd : varlist) {
                if (comp[varInd] == comp[nVar + valInd])
                    continue; // same SCC: the arc belongs to an alternating cycle
                int valIndex = NoAssignedVal[valInd];
                ValList[varInd][valInd] = 0;
                valremoved = true;
                if (varInd >= NbNoAssigned) {
                    // arc of the slack variable: increase the extended unary cost of the value
                    deltaCosts[valIndex] += newCurrentUb;
                } else {
                    int varIndex = NoAssignedVar[varInd];
                    auto* variable = scope[varIndex];
                    Value value = toVarValue(varIndex, valIndex);
                    if (value != WRONG_VAL && variable->canbe(value)) {
                        if (ToulBar2::verbose > 0)
                            cout << "GAC REMOVE VALUE " << value << " from " << variable->getName() << endl;
                        ExtOrProJ(varIndex, value, -newCurrentUb); // project an infinite cost on this value
                    }
                }
                varInd = -1; // mark arc for removal
            }
            if (valremoved) {
                varlist.erase(std::remove(varlist.begin(), varlist.end(), -1), varlist.end());
                removed = true;
            }
        }
        return removed;
    }


    /**
     * @brief Moves the optimal dual solution (ReduceCostRow = u per variable, ReduceCostCol = v per value) of the LAP
     *        into the relative interior of the set of dual optimizers (Algorithm 1 of Dlask & Savchynskyy,
     *        "Relative-interior solution for the (incomplete) linear assignment problem with applications to the
     *        quadratic assignment problem"). Works for the square case and for the slack formulation.
     *
     *        In the relative-interior solution, an edge has a zero reduced cost if and only if it belongs to
     *        at least one optimal assignment; every other edge gets a strictly positive reduced cost.
     *        The dual objective (hence the lower bound) and the optimal matching are unchanged.
     *
     *        Slack formulation (isSlack and more values than variables): the slack variable is ONE node (index n),
     *        with capacity m - n, owning ALL the values not taken by a real variable (the columns j with
     *        colSol[j] == n). Its copies are interchangeable (same costs, same dual u_slack), so they belong to the
     *        same strongly connected component; contracting them into a single node is therefore exact.
     *        When the node is shifted by alpha, u_slack decreases by alpha while its m - n owned values increase by
     *        alpha each: the dual objective sum_i cap_i * u_i + sum_j v_j is unchanged because the capacity equals
     *        the number of owned values.
     *
     *        Steps:
     *        1. Build the directed graph G on the nodes (real variables + slack): arc k -> i if the edge (i, c) has a
     *           zero reduced cost for a column c owned by k (matched edges are contracted).
     *        2. Compute the strongly connected components (Tarjan). Edges inside a component belong to an optimal
     *           assignment and must stay tight.
     *        3. Build the condensation and a topological order of it.
     *        4. For each component V having ingoing arcs in the condensation, taken in reverse topological order
     *           (sinks first), shift its duals: u_r -= alpha for r in V and v_c += alpha for every column c owned by V.
     *           This keeps all edges inside V (and the matched edges) tight, keeps the dual objective, and makes
     *           the ingoing tight arcs strictly positive. alpha is half of the smallest current reduced cost of the edges
     *           (node outside V, column owned by V), which are the only ones decreased by the shift.
     *
     *        Costs are integers: alpha = floor(min/2), so an arc of reduced cost 1 cannot be made strictly positive.
     *        The result is still an optimal and feasible dual solution.
     *
     *        Convention (same as everywhere else in this class):
     *        reducedCost(row, col) = costMatrix[idx(row, col)] - (ReduceCostRow[row] + ReduceCostCol[col]) >= 0,
     *        zero on the matching. The slack row has index n. Call it after rowSol (and colSol in the slack case) has
     *        been filled and BEFORE deltaCosts are recomputed from the duals of the slack variable.
     *
     * @param n   Number of real variables (rows of the LAP).
     * @param m   Number of values (columns of the LAP); m == n unless the slack formulation is used.
     * @param ub  Cost at which an edge is considered absent (costMatrix >= ub: forbidden edge, no dual constraint).
     * @param idx idx(row, col) = position in costMatrix (the slack row is n).
     */
    void makeDualRelativeInterior(int n, int m, Cost ub, const function<size_t(int, int)>& idx)
    {
        const bool hasSlack = isSlack && m > n;
        if (!hasSlack && m != n)
            return; // non-square problem without slack: not handled
        const int N = n + (hasSlack ? 1 : 0); // nodes: real variables, then the slack variable
        const int slackNode = n;
        auto present = [&](int row, int col) { return costMatrix[idx(row, col)] < ub; };
        auto reducedCosts = [&](int row, int col) { return costMatrix[idx(row, col)] - (ReduceCostRow[row] + ReduceCostCol[col]); };

        // Columns owned by each node: the matched column of a real variable, all the values not taken by a real variable for the slack
        vector<vector<int>> owned(N);
        for (int k = 0; k < n; ++k)
            owned[k].push_back(rowSol[k]);
        if (hasSlack) {
            for (int col = 0; col < m; ++col)
                if (colSol[col] == slackNode)
                    owned[slackNode].push_back(col);
        }

        // 1. Graph on the nodes: arc k -> i if the edge (i, c) is tight for a column c owned by k
        vector<vector<int>> g(N);
        for (int k = 0; k < N; ++k) {
            for (int col : owned[k]) {
                for (int i = 0; i < N; ++i) {
                    if (i != k && present(i, col) && reducedCosts(i, col) == MIN_COST)
                        g[k].push_back(i);
                }
            }
        }

        // 2. Strongly connected components
        vector<int> index(N, -1);
        vector<int> lowlink(N, 0);
        vector<int> comp(N, -1);
        vector<uint8_t> onStack(N, 0);
        vector<int> stk;
        int counter = 0;
        int nbComp = 0;
        for (int root = 0; root < N; ++root) {
            if (index[root] == -1) {
                strongconnect(root, g, index, lowlink, comp, onStack, stk, counter, nbComp);
            }
        }

        // 3. Build the list of nodes belonging to each SCC and the condensation graph.
        vector<vector<int>> components(nbComp);
        for (int i = 0; i < N; ++i) {
            components[comp[i]].push_back(i);
        }

        vector<vector<int>> dag(nbComp);
        vector<int> indegree(nbComp, 0);
        for (int k = 0; k < N; ++k) {
            for (int i : g[k]) {
                int from = comp[k];
                int to = comp[i];
                if (from != to) {
                    dag[from].push_back(to);
                    ++indegree[to];
                }
            }
        }
        // A component is shifted only if it has ingoing arcs in the condensation
        vector<uint8_t> hasIn(nbComp, 0);
        for (int c = 0; c < nbComp; ++c)
            hasIn[c] = (indegree[c] > 0);

        // Topological ordering of the condensation graph (Kahn's algorithm: sources first).
        queue<int> q;
        for (int c = 0; c < nbComp; ++c) {
            if (indegree[c] == 0) {
                q.push(c);
            }
        }
        vector<int> topo;
        while (!q.empty()) {
            int c = q.front();
            q.pop();
            topo.push_back(c);
            for (int next : dag[c]) {
                if (--indegree[next] == 0) {
                    q.push(next);
                }
            }
        }

        // 4. Process SCCs in reverse topological order (sinks first).
        for (auto it = topo.rbegin(); it != topo.rend(); ++it) {
            int c = *it;
            if (!hasIn[c])
                continue;

            // Smallest current reduced cost among the edges decreased by the shift:
            // (node outside the component, column owned by the component).
            Cost minRC = MAX_COST;
            for (int i : components[c]) {
                for (int col : owned[i]) {
                    for (int k = 0; k < N; ++k) {
                        if (comp[k] != c && present(k, col)) {
                            minRC = min(minRC, reducedCosts(k, col));
                        }
                    }
                }
            }

            // Maximum admissible shift.
            Cost alpha = (minRC >= MAX_COST) ? 1 : minRC / 2;
            if (alpha <= MIN_COST) {
                continue;
            }

            // Shift the dual variables of this SCC (for the slack node: u_slack and all its owned values).
            for (int i : components[c]) {
                ReduceCostRow[i] -= alpha;
                for (int col : owned[i])
                    ReduceCostCol[col] += alpha;
            }
        }
    }


    /**
     * @brief Filtering of variable domains with the Sellmann or Cambazard method (reduced-cost based filtering).
     *        Works for both the standard and the slack formulation (SlackAssignedVals is empty when isSlack == 0).
     *        Works with assigned variables (partial case) as well as with none (full case, see projectAndFilterFullProblem()).
     *        Step 1: GAC by strongly connected components (gacFilteringSCC), always run when the residual graph exists.
     *        Step 2: reduced-cost filtering with Dijkstra on the remaining arcs, only if some arc can still be pruned.
     *        Does nothing if filtering is disabled (see residualGraphEnabled()).
     *
     *        Source : Claus, G.; Cambazard, H.; and Jost, V. 2020. Analysis of Reduced Costs Filtering
     *        for Alldifferent and Minimum Weight Alldifferent Global Constraints. ECAI 2020, vol. 325, 323-330.
     *
     * @param NoAssignedVal     Value index of each column of the reduced matrix.
     * @param VarList           VarList[val] = variables having a non-matching arc to val (updated when arcs are pruned).
     * @param ValList           ValList[var][val] = 1 if arc (var,val) is in the residual graph (updated when arcs are pruned).
     * @param SlackAssignedVals Columns matched to the slack variable.
     * @param newCurrentUb      Remaining cost budget.
     */
    void reducedCostFiltering(const vector<int>& NoAssignedVal, vector<vector<int>>& VarList, vector<vector<uint8_t>>& ValList, const vector<int>& SlackAssignedVals, Cost newCurrentUb)
    {
        if (!residualGraphEnabled())
            return;

        // Step 1: GAC by strongly connected components
        gacFilteringSCC(NoAssignedVal, VarList, ValList, SlackAssignedVals, newCurrentUb);

        // Step 2: reduced-cost filtering, only if the largest reduced cost can exceed the budget
        if (!(MaxReducedCost >= newCurrentUb / NbNoAssigned))
            return;

        vector<int> VariableList(NbNoAssigned);
        iota(VariableList.begin(), VariableList.end(), 0);

        // Number of Dijkstra sources: all variables, or a random subset
        if (FiltLevel >= 1. - (double)ToulBar2::epsilon) {
            Q = NbNoAssigned;
        } else {
            Q = min(NbNoAssigned, 1 + static_cast<int>(NbNoAssigned * FiltLevel));
            myrearrange(VariableList);
        }

        int comp = 0; // number of Dijkstra runs already done
        int numzeroremoved = 0; // number of consecutive runs without any pruning
        int stop =  (Q / 2) + 1;
        bool AllValConsistent = false;

        // Prunes the arcs (varInd, valInd) of value valInd whose reduced cost exceeds the budget,
        // given the distances computed from the current source. 'var' is the node whose matched value is valInd.
        // Returns true if at least one arc was removed.
        auto pruneArcsOfValue = [&](int valInd, int var) -> bool {
            bool valremoved = false;
            int position = -1;
            auto& varlist = VarList[valInd];
            for (int varInd : varlist) {
                position++;
                if (distanceToVar[varInd] >= MAX_COST)
                    continue;
                Cost reducedCost = ReduceCostMatrix[varInd * NbNoAssignedVal + valInd] - distanceToVar[varInd] + distanceToVar[var];
                if (reducedCost >= newCurrentUb) {
                    int valIndex = NoAssignedVal[valInd];
                    varlist[position] = -1; // mark arc for removal
                    ValList[varInd][valInd] = 0;
                    valremoved = true;
                    if (varInd >= NbNoAssigned) {
                        // arc of the slack variable: increase the extended unary cost of the value
                        deltaCosts[valIndex] += newCurrentUb;
                    } else {
                        int varIndex = NoAssignedVar[varInd];
                        auto* variable = scope[varIndex];
                        Value value = toVarValue(varIndex, valIndex);
                        if (value != WRONG_VAL && variable->canbe(value)) {
                            if (ToulBar2::verbose > 0)
                                cout << "REMOVE VALUE " << value << " from " << variable->getName() << endl;
                            ExtOrProJ(varIndex, value, -newCurrentUb); // SdG: project infinite cost on this value and avoid to skip and reenter the AllDiff constraint without finishing the current filtering
                        }
                    }
                }
            }
            if (valremoved)
                varlist.erase(std::remove(varlist.begin(), varlist.end(), -1), varlist.end());
            return valremoved;
        };

        for (int ind = 0; ind < Q; ++ind) {
            if (numzeroremoved >= stop)
                break;
            int varInde = VariableList[ind];
            numzeroremoved++;

            if (comp == 0) {
                // Apply Régin Upper Bounds of Shortest Paths algorithm with landmark
                AllValConsistent = ReginLandmarkUpperBound(varInde, VarList, ValList, newCurrentUb, SlackAssignedVals);
            } else {
                // Bimodal Dijkstra’s shortest path algorithm
                BimodalDijkstra(varInde, VarList, ValList, SlackAssignedVals);
            }
            if (AllValConsistent)
                break;
            comp++;

            const int nbNodes = NbNoAssigned + (isSlack ? 1 : 0);
            for (int var = 0; var < nbNodes; ++var) {
                if (distanceToVar[var] >= MAX_COST)
                    continue;
                if (isSlack && var == NbNoAssigned) {
                    // the slack variable is matched to several values
                    for (int matchedVal : SlackAssignedVals) {
                        if (pruneArcsOfValue(matchedVal, var))
                            numzeroremoved = 0;
                    }
                } else if (pruneArcsOfValue(rowSol[var], var)) {
                    numzeroremoved = 0;
                }
            }
        }
    }

    // True if the residual graph (VarList, ValList, ReduceCostMatrix) is built and used for filtering:
    // reduced-cost filtering must be enabled
    bool residualGraphEnabled() const
    {
        return FiltLevel > 0;
    }


    /**
     * @brief Full case (no variable considered as assigned): projects the reduced costs onto the unary costs,
     *        builds the residual graph over ALL the variables and ALL the values, and runs reducedCostFiltering().
     *        The whole problem is viewed as a "reduced" problem: NbNoAssigned = arity_, NbNoAssignedVal = NbValues,
     *        NoAssignedVar and NoAssignedVal are the identity. Works for the standard (square) and the slack formulation.
     * @param current_ub   wcsp->getUb() - wcsp->getLb() before the LAP cost projection.
     * @param newCurrentUb current_ub minus the LAP optimal cost (remaining budget).
     * @param idx          idx(var, col) = position in costMatrix.
     */
    void projectAndFilterFullProblem(Cost current_ub, Cost newCurrentUb, const function<size_t(int, int)>& idx)
    {

        NbNoAssigned = arity_;
        NbNoAssignedVal = NbValues;
        NoAssignedVar.resize(arity_);
        iota(NoAssignedVar.begin(), NoAssignedVar.end(), 0);
        vector<int> NoAssignedVal(NbValues);
        iota(NoAssignedVal.begin(), NoAssignedVal.end(), 0);

        vector<vector<int>> VarList(NbNoAssignedVal);
        vector<vector<uint8_t>> ValList(NbNoAssigned + (isSlack ? 1 : 0));
        vector<int> SlackAssignedVals;
        projectAndBuildResidualGraph(NoAssignedVal, current_ub, newCurrentUb, idx, VarList, ValList);
        if (isSlack && residualGraphEnabled())
            addSlackArcsToResidualGraph(current_ub, newCurrentUb, idx, VarList, ValList, SlackAssignedVals);
        reducedCostFiltering(NoAssignedVal, VarList, ValList, SlackAssignedVals, newCurrentUb);

    }

    /**
     * @brief Standard case, some variables are assigned (same domains, no excepted values):
     *        reduces the problem to the unassigned variables and the free values,
     *        solves the LAP, projects the costs and filters the domains.
     */
    void propagatePartialAssignment()
    {
        // Collect unassigned values
        vector<int> NoAssignedVal;
        for (int valIndex = 0; valIndex < NbValues; ++valIndex) {
            if (!isAssignedValue[valIndex]) {
                NoAssignedVal.push_back(valIndex);
            }
        }

        NbNoAssignedVal = NoAssignedVal.size();
        if (NbNoAssignedVal < NbNoAssigned) {
            THROWCONTRADICTION;
        }

        // Initialize cost matrix (rows = unassigned variables, columns = free values) for the Jonker algorithm
        int removevals = 0;
        Cost current_ub = wcsp->getUb() - wcsp->getLb();

        for (int valInd = 0; valInd < NbNoAssignedVal; ++valInd) {
            int valIndex = NoAssignedVal[valInd];
            int comptremoveval = 0; // number of variables which cannot take this value
            for (int varInd = 0; varInd < NbNoAssigned; ++varInd) {
                int varIndex = NoAssignedVar[varInd];
                auto* variable = scope[varIndex];
                Value val = toVarValue(varIndex, valIndex);
                if (val != WRONG_VAL && variable->canbe(val)) {
                    auto valCost = variable->getCost(val);
                    costMatrix[varInd * NbNoAssignedVal + valInd] = valCost < current_ub ? valCost : current_ub;
                    if (valCost >= current_ub) {
                        comptremoveval++;
                    }
                } else {
                    comptremoveval++;
                    costMatrix[varInd * NbNoAssignedVal + valInd] = current_ub;
                }
            }

            // Nobody can take this value anymore: its extended unary cost is projected to the lower bound
            if (comptremoveval == NbNoAssigned) {
                if (isSquare) {
                    THROWCONTRADICTION;
                }
                Cost Delta = deltaCosts[valIndex];
                if (Delta >= current_ub) {
                    THROWCONTRADICTION;
                }
                deltaCosts[valIndex] = 0;
                projectLB(Delta);
                removevals++;
            }
        }

        int NbRestVal = NbNoAssignedVal - removevals;
        if (NbRestVal < NbNoAssigned) {
            THROWCONTRADICTION;
        }

        auto idx = [this](int var, int col) { return (size_t)var * NbNoAssignedVal + col; };

        if (any_of(costMatrix.begin(), costMatrix.begin() + NbNoAssigned * NbNoAssignedVal, [](auto x) { return x != MIN_COST; })) {
            // Solve assignment problem using Jonker algorithm
            Cost TotalCost;
            if (!isSquare && NbRestVal > NbNoAssigned) {
                // more values than variables: add the dummy row with the extended unary costs
                setDeltaCostsInMatrix(NoAssignedVal, NbNoAssigned, current_ub, idx);
                TotalCost = lapjv(NbNoAssigned + 1, NbNoAssignedVal, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub, findConflict);
            } else {
                TotalCost = lapjv(NbNoAssigned, NbNoAssignedVal, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub, findConflict);
            }

            if (TotalCost >= current_ub) {
                onLAPFailure(rowSol, findConflict, NbNoAssigned, &NoAssignedVar);
            } else if (TotalCost >= 0) {
                storePartialAssignment(NoAssignedVal);

                if (TotalCost > 0 || (any_of(ReduceCostRow, ReduceCostRow + NbNoAssigned, [](auto x) { return x > MIN_COST; }))) {
                    if (!isSquare && NbRestVal > NbNoAssigned) {
                        Cost mindelta = ReduceCostRow[NbNoAssigned];
                        for (int valInd = 0; valInd < NbNoAssignedVal; ++valInd) {
                            deltaCosts[NoAssignedVal[valInd]] -= (ReduceCostCol[valInd] + mindelta);
                        }
                        if (mindelta > 0) {
                            TotalCost += mindelta * (NbValues - arity_ - 1);
                            if (TotalCost >= current_ub) {
                                wcsp->revise(this);
                                THROWCONTRADICTION;
                            }
                        }
                    }
                    // Square case: move the duals to the relative interior of the dual optimizers
                    if (relativeInterior && isSquare)
                        makeDualRelativeInterior(NbNoAssigned, NbNoAssignedVal, current_ub, idx);
                    // Update the lower bound with the Jonker algorithm's total cost
                    projectLB(TotalCost);
                    Cost newCurrentUb = current_ub - TotalCost;
                    
                      if (residualGraphEnabled()){
                        // Adjust unary costs for unassigned variables based on reduced costs and build the residual graph
                        vector<vector<int>> VarList(NbNoAssignedVal);
                        vector<vector<uint8_t>> ValList(NbNoAssigned);
                        projectAndBuildResidualGraph(NoAssignedVal, current_ub, newCurrentUb, idx, VarList, ValList);

                        // Filtering of variables domains with Sellmann or Cambazard method
                        reducedCostFiltering(NoAssignedVal, VarList, ValList, vector<int>(), newCurrentUb);
                      } else {
                        projectReducedCostsPartial(NoAssignedVal, current_ub, idx);
                      }
                }

                // Update support values if needed for unassigned variables
                updateSupportsOfUnassigned();
            }
        } else {
            // Zero cost matrix: every unassigned variable can take every free value at no cost,
            // so any injective assignment is optimal (no cost to project)
            storePartialZeroCost(NoAssignedVal);
        }
    }

    /**
     * @brief Standard case, all variables unassigned (or excepted values, or different domains):
     *        solves the LAP on the whole problem, projects the costs and updates the supports.
     */
    void propagateFullAssignment()
    {
        Cost current_ub = wcsp->getUb() - wcsp->getLb();
        auto idx = [this](int var, int col) { return (size_t)var * NbValues + col; };

        // Initialize the cost matrix for all variables and their domain values
        fillFullCostMatrix(current_ub, idx);

        if (any_of(costMatrix.begin(), costMatrix.begin() + arity_ * NbValues, [](auto x) { return x != MIN_COST; })) {
            // Solve the Linear Assignment Problem (LAP) using the Jonker algorithm
            Cost TotalCost;

            if (!isSquare && !slack) {
                // more values than variables: add the dummy row with the extended unary costs
                vector<int> allValues(NbValues);
                iota(allValues.begin(), allValues.end(), 0);
                setDeltaCostsInMatrix(allValues, arity_, current_ub, idx);

                if (excepted) {
                    TotalCost = lapjv(arity_ + 1, NbValues, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub, exceptedValIndex, findConflict);
                } else {
                    TotalCost = lapjv(arity_ + 1, NbValues, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub, findConflict);
                }
            } else {
                TotalCost = lapjv(arity_, NbValues, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub, findConflict);
            }

            if (TotalCost >= current_ub) {
                onLAPFailure(rowSol, findConflict, arity_, nullptr);
            } else if (TotalCost >= 0) {
                storeFullAssignment();

                if (TotalCost > 0 || (any_of(ReduceCostRow, ReduceCostRow + arity_, [](auto x) { return x > MIN_COST; }))) {
                    if (!isSquare && !slack) {
                        Cost mindelta = ReduceCostRow[arity_];
                        for (int valInd = 0; valInd < NbValues; ++valInd) {
                            deltaCosts[valInd] -= (ReduceCostCol[valInd] + mindelta);
                        }

                        if (mindelta > 0) {
                            int nbexcep = exceptedValIndex.size();
                            if (excepted) {
                                TotalCost += (mindelta * (NbValues - arity_ - nbexcep + (nbexcep * arity_)) - mindelta);
                            } else {
                                TotalCost += mindelta * (NbValues - arity_ - 1);
                            }
                            if (TotalCost >= current_ub) {
                                wcsp->revise(this);
                                THROWCONTRADICTION;
                            }
                        }
                    }

                    // Square case: move the duals to the relative interior of the dual optimizers
                    if (relativeInterior && isSquare )
                        makeDualRelativeInterior(arity_, NbValues, current_ub, idx);
                    // Update lower bound with the total cost found
                    projectLB(TotalCost);
                    
                    // Adjust unary costs using reduced row and column costs
                    // (the non-square case without slack never uses filtering: FiltLevel > 0 implies slack)
                    if (residualGraphEnabled()){
                         projectAndFilterFullProblem(current_ub, current_ub - TotalCost, idx); // + GAC and reduced-cost filtering
                    }    
                    else {
                        projectReducedCostsFull(current_ub, idx);
                    }
                }

                // Update support if needed for all variables
                updateAllSupports();
            }
        } else {
            // Zero cost matrix: every variable can take every value at no cost.
            storeZeroCostFullAssignment();
        }
    }


    /**
     * @brief Slack formulation, some variables are assigned (same domains):
     *        same as propagatePartialAssignment() but with a dummy slack variable absorbing the unused values.
     *        The cost matrix is transposed (row = value, column = variable + 1 dummy column).
     */
    void propagatePartialAssignmentSlack()
    {
        // Collect unassigned values and initialize the cost matrix for the Jonker algorithm
        vector<int> NoAssignedVal;
        int valInd = 0;
        Cost current_ub = wcsp->getUb() - wcsp->getLb();

        for (int valIndex = 0; valIndex < NbValues; ++valIndex) {
            if (isAssignedValue[valIndex])
                continue;
            int comptremoveval = 0; // number of variables which cannot take this value
            for (int varInd = 0; varInd < NbNoAssigned; ++varInd) {
                int varIndex = NoAssignedVar[varInd];
                auto* variable = scope[varIndex];
                Value val = toVarValue(varIndex, valIndex);
                if (val != WRONG_VAL && variable->canbe(val)) {
                    auto valCost = variable->getCost(val);
                    costMatrix[valInd * (NbNoAssigned + 1) + varInd] = valCost < current_ub ? valCost : current_ub;
                    if (valCost >= current_ub) {
                        comptremoveval++;
                    }
                } else {
                    comptremoveval++;
                    costMatrix[valInd * (NbNoAssigned + 1) + varInd] = current_ub;
                }
            }

            // Nobody can take this value anymore: its extended unary cost is projected to the lower bound
            if (comptremoveval == NbNoAssigned) {
                Cost Delta = deltaCosts[valIndex];
                if (Delta >= current_ub) {
                    THROWCONTRADICTION;
                }
                deltaCosts[valIndex] = 0;
                projectLB(Delta);
            } else {
                NoAssignedVal.push_back(valIndex);
                valInd++;
            }
        }

        NbNoAssignedVal = NoAssignedVal.size();
        if (NbNoAssignedVal < NbNoAssigned) {
            THROWCONTRADICTION;
        }

        // The matrix stride depends on whether the dummy column is still present
        auto idx = [this](int var, int col) { return (size_t)col * (NbNoAssigned + (isSlack ? 1 : 0)) + var; };

        if (NbNoAssignedVal == NbNoAssigned) {
            // As many free values as unassigned variables: the problem is square, no more slack variable
            isSlack = 0;
            // Compact the matrix (remove the dummy column)
            for (int col = 1; col < NbNoAssignedVal; ++col) {
                for (int varInd = 0; varInd < NbNoAssigned; ++varInd) {
                    costMatrix[col * NbNoAssigned + varInd] = costMatrix[col * (NbNoAssigned + 1) + varInd];
                }
            }
        } else {
            // Dummy column with the extended unary costs
            setDeltaCostsInMatrix(NoAssignedVal, NbNoAssigned, current_ub, idx);
            capacity.assign(NbNoAssigned + 1, 1);
            capacity[NbNoAssigned] = NbNoAssignedVal - NbNoAssigned;
        }

        if (any_of(costMatrix.begin(), costMatrix.begin() + (NbNoAssigned + (isSlack ? 1 : 0)) * NbNoAssignedVal, [](auto x) { return x != MIN_COST; })) {
            // Solve assignment problem using Jonker algorithm
            Cost TotalCost;
            if (isSlack == 0) {
                TotalCost = lapjv(NbNoAssignedVal, NbNoAssigned, costMatrix, colSol, ReduceCostCol, ReduceCostRow, current_ub, findConflict, true);
            } else {
                TotalCost = lapjv_ub(NbNoAssignedVal, NbNoAssigned + 1, costMatrix, colSol, ReduceCostCol, ReduceCostRow, current_ub, capacity, findConflict, true);
            }

            if (TotalCost >= current_ub) {
                onLAPFailure(colSol, findConflict, NbNoAssigned, &NoAssignedVar);
            } else if (TotalCost >= 0) {

                // The LAP solver returns colSol (value -> variable): derive rowSol (variable -> value)
                for (int valIndex = 0; valIndex < NbNoAssignedVal; ++valIndex) {
                    int var = colSol[valIndex];
                    if (var < NbNoAssigned) {
                        rowSol[var] = valIndex;
                    }
                }

                storePartialAssignment(NoAssignedVal);

                if (TotalCost > 0 || (any_of(ReduceCostCol, ReduceCostCol + NbNoAssignedVal, [](auto x) { return x > MIN_COST; }))) {
                    // Slack (or locally square) formulation: move the duals to the relative interior of the dual
                    // optimizers, BEFORE the extended unary costs are derived from the duals of the slack variable
                    if (relativeInterior)
                        makeDualRelativeInterior(NbNoAssigned, NbNoAssignedVal, current_ub, idx);

                    if (isSlack == 1) {
                        // New extended unary costs from the duals of the dummy variable
                        for (int valIndex = 0; valIndex < NbNoAssignedVal; ++valIndex) {
                            deltaCosts[NoAssignedVal[valIndex]] = costMatrix[idx(NbNoAssigned, valIndex)] - (ReduceCostRow[NbNoAssigned] + ReduceCostCol[valIndex]);
                        }
                    }

                    // Update the lower bound with the Jonker algorithm's total cost
                    projectLB(TotalCost);
                    Cost newCurrentUb = current_ub - TotalCost;

                    if (residualGraphEnabled()) {
                        vector<vector<int>> VarList(NbNoAssignedVal);
                        vector<vector<uint8_t>> ValList(NbNoAssigned + (isSlack ? 1 : 0));
                        // Adjust unary costs for unassigned variables based on reduced costs and build the residual graph
                        projectAndBuildResidualGraph(NoAssignedVal, current_ub, newCurrentUb, idx, VarList, ValList);
                        if(isSlack){
                            // Values matched to the slack variable
                            vector<int> SlackAssignedVals;
                            addSlackArcsToResidualGraph(current_ub, newCurrentUb, idx, VarList, ValList, SlackAssignedVals);
                            reducedCostFiltering(NoAssignedVal, VarList, ValList, SlackAssignedVals, newCurrentUb);
                        } else {
                        // Filtering of variables domains with Sellmann or Cambazard method
                        reducedCostFiltering(NoAssignedVal, VarList, ValList, {}, newCurrentUb);
                        }
                    } else{
                        projectReducedCostsPartial(NoAssignedVal, current_ub, idx);
                    }     
                }

                // Update support values if needed for unassigned variables
                updateSupportsOfUnassigned();
            }
        } else {
            // Zero cost matrix: every unassigned variable can take every free value at no cost
            storePartialZeroCost(NoAssignedVal);
        }
    }

    /**
     * @brief Slack formulation, all variables unassigned (or different domains):
     *        solves the LAP with a dummy slack variable on the whole problem.
     *        The cost matrix is transposed (row = value, column = variable + 1 dummy column).
     */
    void propagateFullAssignmentSlack()
    {
        Cost current_ub = wcsp->getUb() - wcsp->getLb();
        auto idx = [this](int var, int col) { return (size_t)col * (arity_ + 1) + var; };

        // Initialize the cost matrix for all variables and their domain values + dummy column
        fillFullCostMatrix(current_ub, idx);
        vector<int> allValues(NbValues);
        iota(allValues.begin(), allValues.end(), 0);
        setDeltaCostsInMatrix(allValues, arity_, current_ub, idx);
        capacity.assign(arity_ + 1, 1);
        capacity[arity_] = NbValues - arity_;

        if (any_of(costMatrix.begin(), costMatrix.begin() + (arity_ + 1) * NbValues, [](auto x) { return x != MIN_COST; })) {
            // Solve the Linear Assignment Problem (LAP) using the Jonker algorithm
            Cost TotalCost = lapjv_ub(NbValues, arity_ + 1, costMatrix, colSol, ReduceCostCol, ReduceCostRow, current_ub, capacity, findConflict, true);

            if (TotalCost >= current_ub) {
                onLAPFailure(colSol, findConflict, arity_, nullptr);
            } else if (TotalCost >= 0) {
                // The LAP solver returns colSol (value -> variable): derive rowSol (variable -> value)
                for (int valIndex = 0; valIndex < NbValues; ++valIndex) {
                    int var = colSol[valIndex];
                    if (var < arity_) {
                        rowSol[var] = valIndex;
                    }
                }
                storeFullAssignment();

                // (in the non-square case the duals have the opposite sign convention: do not test their sign)
                if (TotalCost > 0 || (any_of(ReduceCostCol, ReduceCostCol + NbValues, [](auto x) { return x > MIN_COST; }))) {
                    // Move the duals to the relative interior of the dual optimizers,
                    // BEFORE the extended unary costs are derived from the duals of the slack variable
                    if (relativeInterior)
                        makeDualRelativeInterior(arity_, NbValues, current_ub, idx);

                    // New extended unary costs from the duals of the dummy variable
                    for (int valIndex = 0; valIndex < NbValues; ++valIndex) {
                        deltaCosts[valIndex] = costMatrix[idx(arity_, valIndex)] - (ReduceCostRow[arity_] + ReduceCostCol[valIndex]);
                    }
                    // Update lower bound with the total cost found
                    projectLB(TotalCost);
                    Cost newCurrentUb = current_ub - TotalCost;
                    // Adjust unary costs using reduced row and column costs
                    if (FiltLevel > 0)
                        projectAndFilterFullProblem(current_ub, newCurrentUb, idx); // + GAC and reduced-cost filtering
                    else
                        projectReducedCostsFull(current_ub, idx);
                }

                // Update support if needed for all variables
                updateAllSupports();
            }
        } else {
            // Zero cost matrix: every variable can take every value at no cost.
            storeZeroCostFullAssignment();
        }
    }


    /**
     * @brief Standard formulation (isSlack == 0): dispatches between the partial and the full case.
     */
    void propagateWithoutSlack()
    {
        // Initialize number of unassigned variables
        NbNoAssigned = getNonAssigned();
        bool filtreExcepted = false;

        // If there are assigned variables, handle accordingly
        if (NbNoAssigned < arity_) {
            if (excepted) {
                filtreExcepted = RemoveAssignVar();
            } else if (SameDomain && RemoveAssignVar()) {
                propagatePartialAssignment();
            }
        }

        // Case when all variables are unassigned or filtering exception or domains differ
        if (connected() && (NbNoAssigned == arity_ || filtreExcepted || !SameDomain)) {
            propagateFullAssignment();
        }
    }


    /**
     * @brief Slack formulation (isSlack == 1): dispatches between the partial and the full case.
     */
    void propagateWithSlack()
    {
        // Initialize number of unassigned variables
        NbNoAssigned = getNonAssigned();

        if (NbNoAssigned < arity_) {
            // If there are assigned variables, handle accordingly
            if (SameDomain && RemoveAssignVar()) {
                propagatePartialAssignmentSlack();
            }
        } else if (connected() && (NbNoAssigned == arity_ || !SameDomain)) {
            // Case when all variables are unassigned or domains differ
            propagateFullAssignmentSlack();
        }
    }


    /**
     * @brief Propagate the AllDifferent constraint to enforce consistency and improve problem lower bound.
     *
     * This method performs constraint propagation by checking assigned variables,
     * removing inconsistent values, computing lower bounds using the Jonker algorithm,
     * and updating unary costs and supports accordingly.
     *
     * It handles different cases depending on whether variables are assigned or not,
     * if the domains are the same, and whether "excepted" values are considered.
     * The work is delegated to propagateWithoutSlack() or propagateWithSlack().
     *
     * @throws TimeOut if the propagation is interrupted.
     * @throws ContradictionException if the constraint becomes unsatisfiable.
     */
    void propagate() override
    {
        if (ToulBar2::dumpWCSP % 2) // skip propagation if the problem is dumped before preprocessing
            return;

        if (ToulBar2::interrupted && !ToulBar2::isZ) {
            throw TimeOut();
        }

        // Propagate only if the constraint is connected
        if (!connected())
            return;

        if (ToulBar2::verbose >= 3) {
            cout << "propagate " << *this << endl;
        }

        // not enough distinct values for the variables (without excepted values)
        if (!excepted && NbValues < arity_)
            THROWCONTRADICTION;

        // Attempt to assign variables already connected and assigned (assign() propagates by itself)
        if (assignFirstAssignedVariable() || !connected())
            return;

        // Skip propagation if the stored optimal assignment is still valid
        if (canSkipPropagation())
            return;

        if (isSlack == 0)
            propagateWithoutSlack();
        else
            propagateWithSlack();
    }

    // checks that the constraint is still satisfiable (called by WCSP::verify in Debug mode at each search node)
    bool verify() override
    {
        if (!storeAssignment) {
            return false;
        }
        vector<bool> alreadyUsed(NbValues, false);
        for (int i = 0; i < arity_; i++) {
            int valIndex = unionIndex(scope[i]->getValueName(scope[i]->toIndex(storeLastAssignment[i])));
            if (valIndex < 0) {
                cout << "variable " << scope[i]->getName() << " value " << storeLastAssignment[i] << " is forbidden!" << endl;
                return false;
            }
            if ((alreadyUsed[valIndex] && (!excepted || find(exceptedValIndex.begin(), exceptedValIndex.end(), valIndex) == exceptedValIndex.end())) || scope[i]->cannotbe(storeLastAssignment[i]) || (ToulBar2::LcLevel >= LC_DAC && scope[i]->getCost(storeLastAssignment[i]) > MIN_COST)) {
                if (alreadyUsed[valIndex]) {
                    cout << "variable " << scope[i]->getName() << " value " << storeLastAssignment[i] << " used twice!" << endl;
                } else if (scope[i]->cannotbe(storeLastAssignment[i])) {
                    cout << "variable " << scope[i]->getName() << " value " << storeLastAssignment[i] << " has been removed!" << endl;
                } else if (scope[i]->getCost(storeLastAssignment[i]) > MIN_COST) {
                    cout << "variable " << scope[i]->getName() << " value " << storeLastAssignment[i] << " has nonzero cost!" << endl;
                }
                return false;
            } else {
                alreadyUsed[valIndex] = true;
            }
        }
        return true;
    }

    // Events on variable domains: a bound change is handled as a value removal
    void increase(int index) override
    {
        remove(index);
    }
    void decrease(int index) override
    {
        remove(index);
    }
    // A value has been removed from the variable at position index: propagate only if the stored assignment is broken
    void remove(int index) override
    {
        if (scope[index]->unassigned()) {
            bool revise = ToulBar2::FullEAC && getVar(index)->cannotbe(getVar(index)->getSupport());
            if (!storeAssignment || scope[index]->cannotbe(storeLastAssignment[index])) {
                propagate();
            }
            if (revise)
                reviseEACGreedySolution();
        } else {
            assign(index);
        }
    }
    // The unary cost of a value of the variable at position index increased from zero: propagate only if the stored assignment is broken
    void projectFromZero(int index) override
    {
        // TODO: incremental cost propagation
        bool revise = ToulBar2::FullEAC && (getVar(index)->cannotbe(getVar(index)->getSupport()) || getVar(index)->getCost(getVar(index)->getSupport()) > MIN_COST);
        if (!storeAssignment || scope[index]->cannotbe(storeLastAssignment[index]) || scope[index]->getCost(storeLastAssignment[index]) > MIN_COST) {
            propagate();
        }
        if (revise)
            reviseEACGreedySolution();
    }

    // bool checkEACGreedySolution(int index = -1, Value supportValue = 0) FINAL; // TODO

    // Checks the EAC greedy solution and unsets the full-EAC flag of the concerned variable(s) if it is broken
    bool reviseEACGreedySolution(int index = -1, Value supportValue = 0) FINAL
    {
        bool result = checkEACGreedySolution(index, supportValue);
        if (!result) {
            if (index >= 0) {
                getVar(index)->unsetFullEAC();
            } else {
                int a = arity();
                for (int i = 0; i < a; i++) {
                    getVar(i)->unsetFullEAC();
                }
            }
        }
        return result;
    }

    // returns the non-zero extended unary costs indexed by value (a value is reported through the first variable that can still take it)
    map<Value, Cost> nonZeroDeltaCosts() const
    {
        map<Value, Cost> result;
        if (isSquare)
            return result;
        for (int valIndex = 0; valIndex < NbValues; valIndex++) {
            Cost delta = deltaCosts[valIndex];
            if (delta == MIN_COST)
                continue;
            for (int varIndex = 0; varIndex < arity_; varIndex++) {
                Value value = toVarValue(varIndex, valIndex);
                if (value != WRONG_VAL && scope[varIndex]->canbe(value)) {
                    result[value] = delta;
                    break;
                }
            }
        }
        return result;
    }

    // List of the values removed from the domains of ALL the variables of the scope
    vector<Value> removedValues()
    {
        map<Value, bool> stillThere; // value -> true if at least one variable can still take it
        for (int i = 0; i < arity_; i++) {
            auto* var = scope[i];
            int initSize = (int)var->getDomainInitSize();
            for (int k = 0; k < initSize; k++) {
                Value v = var->toValue(k);
                bool can = var->canbe(v);
                auto res = stillThere.emplace(v, can);
                if (!res.second && can)
                    res.first->second = true;
            }
        }
        vector<Value> removed;
        for (const auto& [v, ok] : stillThere) {
            if (!ok)
                removed.push_back(v);
        }
        return removed;
    }

    // Prints the constraint (scope, extended unary costs, tightness, weights, excepted values, arity) for debugging
    void print(ostream& os) override
    {
        os << this << " alldifferent(";

        int unassigned_ = 0;
        for (int i = 0; i < arity_; i++) {
            if (scope[i]->unassigned())
                unassigned_++;
            os << wcsp->getName(scope[i]->wcspIndex);
            if (i < arity_ - 1)
                os << ",";
        }
        if (!isSquare) {
            os << ") "
               << " cost: " << -lb << " + " << assigneddeltas << " + (";
            // for (int i = 0; i < arity_; i++) {
            for (unsigned int j = 0; j < deltaCosts.size(); j++) {
                os << deltaCosts[j];
                if (j < deltaCosts.size() - 1)
                    os << "|";
            }
            // if (i < arity_ - 1)
            // os << ",";
            //}
            os << ") ";
        }
        os << "/" << getTightness();
        if (ToulBar2::weightedDegree) {
            os << "/" << getConflictWeight();
            for (int i = 0; i < arity_; i++) {
                os << "," << conflictWeights[i];
            }
        }
        os << " exceptedvalues:";
        for (Value v : exceptedValues) {
            os << " " << v;
        }
        os << " arity: " << arity_;
        os << " unassigned: " << getNonAssigned() << "/" << unassigned_ << endl;
    }

    // Dumps the constraint in WCSP format (original or current problem)
    void dump(ostream& os, bool original = true) override
    {
        map<Value, Cost> mapValuesDeltaCosts;
        if (!isSquare)
            mapValuesDeltaCosts = nonZeroDeltaCosts();
        if (original) {
            os << arity_;
            for (int i = 0; i < arity_; i++)
                os << " " << scope[i]->wcspIndex;
            os << " -1 alldiff ";
        } else {
            os << getNonAssigned();
            for (int i = 0; i < arity_; i++)
                if (scope[i]->unassigned())
                    os << " " << scope[i]->getCurrentVarId();
            os << " -1 alldiff ";
        }

        vector<Value> removedVals = removedValues();
        // the counters are positional in read(): write the excepted counter as soon as a later section exists
        // excepted values: number followed by the list
        if (!exceptedValues.empty()) {
            os << exceptedValues.size();
            for (Value v : exceptedValues)
                os << " " << v;
        } else if (!removedVals.empty() || !mapValuesDeltaCosts.empty())
            os << 0;

        // removed values: number followed by the list
        if (!removedVals.empty()) {
            os << " " << removedVals.size();
            for (Value v : removedVals)
                os << " " << v;
        } else if (!mapValuesDeltaCosts.empty())
            os << " " << 0;

        // delta costs (non-square case only)
        if (!isSquare) {
            if (!mapValuesDeltaCosts.empty()) {

                os << " " << mapValuesDeltaCosts.size();
                for (auto& [key, delta] : mapValuesDeltaCosts) {
                    os << " " << key;
                    os << " " << delta;
                }
            }
        }
        os << endl;
    }

    // Dumps the constraint in CFN (JSON) format (original or current problem)
    void dump_CFN(ostream& os, bool original = true) override
    {
        bool printed = false;
        os << "\"F_";
        if (original) {
            for (int i = 0; i < arity_; i++) {
                if (printed)
                    os << "_";
                os << scope[i]->wcspIndex;
                printed = true;
            }

            os << "\":{\"scope\":[";
            printed = false;
            for (int i = 0; i < arity_; i++) {
                if (printed)
                    os << ",";
                os << "\"" << name2cfn(scope[i]->getName()) << "\"";
                printed = true;
            }
            os << "],\"type\":\"alldiff\",\"params\":{";
        } else {
            for (int i = 0; i < arity_; i++)
                if (scope[i]->unassigned()) {
                    if (printed)
                        os << "_";
                    os << scope[i]->getCurrentVarId();
                    printed = true;
                }
            os << "\":{\"scope\":[";
            printed = false;
            for (int i = 0; i < arity_; i++)
                if (scope[i]->unassigned()) {
                    if (printed)
                        os << ",";
                    os << "\"" << name2cfn(scope[i]->getName()) << "\"";
                    printed = true;
                }
            os << "],\"type\":\"alldiff\",\"params\":{";
        }

        // Separator shared by the three optional parameters,
        // so that commas are always valid JSON.
        bool firstParam = true;
        auto sep = [&]() {
            if (!firstParam)
                os << ",";
            firstParam = false;
        };

        // excepted values
        if (!exceptedValues.empty()) {
            sep();
            os << "\"exceptedvalues\":[";
            printed = false;
            for (Value v : exceptedValues) {
                if (printed)
                    os << ",";
                os << v;
                printed = true;
            }
            os << "]";
        }

        // removed values
        vector<Value> removedVals = removedValues();
        if (!removedVals.empty()) {
            sep();
            os << "\"removedvalues\":[";
            for (size_t i = 0; i < removedVals.size(); i++) {
                if (i > 0)
                    os << ",";
                os << removedVals[i];
            }
            os << "]";
        }

        // delta costs (non-square case only)
        if (!isSquare) {
            map<Value, Cost> mapValuesDeltaCosts = nonZeroDeltaCosts();
            if (!mapValuesDeltaCosts.empty()) {
                sep();
                os << "\"deltacosts\":[";
                printed = false;
                for (auto& [key, delta] : mapValuesDeltaCosts) {
                    if (printed)
                        os << ",";
                    os << "[" << key << "," << delta << "]";
                    printed = true;
                }
                os << "]";
            }
        }
        os << "}},\n";
    }
};

#endif /*TB2ALLDIFFERENT_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
