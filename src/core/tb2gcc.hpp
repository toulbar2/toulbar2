/** \file tb2gcc.hpp
 *  \brief Propagates a Global Cardinality Constraint.
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

#ifndef TB2GCC_HPP_
#define TB2GCC_HPP_

#pragma once
#include <utility>
#include <variant>
#include <cstdint>
#include "tb2abstractconstr.hpp"
#include "tb2ternaryconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2wcsp.hpp"
#include "search/tb2clusters.hpp"
#include "utils/tb2lapjv.hpp"
#include <unordered_set>
#include <unordered_map>
#include <boost/heap/pairing_heap.hpp>
#include <set>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <numeric>
#include <cassert>
#include <random>
#include <algorithm>

using namespace std;
using PQ = boost::heap::pairing_heap<pair<Cost, int>, boost::heap::compare<greater<>>>;

class GlobalCardinalityConstraint : public AbstractNaryConstraint {
    Cost Original_ub; // Initial upper bound when creating the constraint
    StoreCost lb; // Projected cost to problem lower bound (if zero, all deltaCosts must be zero)
    StoreCost assigneddeltas; // Accumulated deltas from assigned values (used in cost propagation)
    vector<Long> conflictWeights; // Used by weighted degree heuristics to prioritize variables
    vector<StoreCost> deltaCosts; // Extended unary costs to all values
    vector<StoreValue> storeLastAssignment; // Stores the last optimal assignment of variables
    StoreInt storeAssignment; // True if an optimal assignment is currently stored
    int* rowSol = nullptr; // Row solution array: rowSol[var] = value index assigned to variable var by the LAP
    int* colSol = nullptr; // Column solution array: colSol[val] = variable index assigned to value val by the LAP
    Cost* ReduceCostRow = nullptr; // Reduced costs per row (per variable) after solving the LAP
    Cost* ReduceCostCol = nullptr; // Reduced costs per column (per value) after solving the LAP
    vector<Cost> costMatrix; // Flattened row-major cost matrix passed to the Jonker LAP
    int NbValues = 0; // Total number of distinct values across all variable domains (size of UnionVarDomain)
    vector<int> AssignedVar; // Indices (in scope) of currently assigned variables
    vector<int> AssignedVal; // Value indices (in UnionVarDomain) assigned to each variable in AssignedVar
    vector<int> NoAssignedVar; // Indices (in scope) of currently unassigned variables
    int NbAssigned; // Number of currently assigned variables
    int NbNoAssigned; // Number of currently unassigned variables
    int NbNoAssignedVal = 0; // Number of values not yet taken by any assigned variable (columns in the reduced cost matrix)
    vector<bool> isAssignedValue; // isAssignedValue[v] = true if value v (index in UnionVarDomain) is already taken by an assigned variable
    vector<bool> varAlreadyProcessed; // varAlreadyProcessed[i] = true if variable i has already been processed during domain filtering
    bool SameDomain = true; // True if all variables have the same domain
    map<Value, pair<int, int>> bounds; // lower and upper bound capacities for every value
    vector<int> capacity; // upper bound capacities for every value
    vector<int> CurrentCapacity; // remaining upper bound capacities for every value (updated during propagation)
    vector<int> demand; // lower bound capacities for every value
    int sumub = 0, sumlb = 0, Gcc_NbValues = 0; // sum of upper bound capacities, sum of lower bound capacities and number of values given in the gcc.
    bool isSquare = true; // True if the cost matrix is square (number of variables == number of values)
    vector<string> UnionVarDomain; // Sorted union of all variable domains (unique value names)
    vector<int> VarDomainSize; // Number of kept values (present in the union) for each variable (indexed by variable position in scope)
    unordered_map<string, int> mapDomainValToIndex; // Maps each value name to its index in UnionVarDomain (index >= NbValues means the value has been removed)
    vector<Cost> ReduceCostMatrix; // Reduced cost matrix after LAP solving: ReduceCostMatrix[var * NbNoAssignedVal + val] = reduced cost of assigning value val to variable var
    vector<uint8_t> visited; // visited[v] = 1 if variable node v has been finalized by Dijkstra
    vector<int> visitVal; // visitVal[v] = 1 if value has been used by other variable by Dijkstra
    vector<uint8_t> AlreadyUse; // AlreadyUse[val] = 1 if value has already been used by other variable
    vector<uint8_t> inHeap; // inHeap[v] = 1 if variable/value node v is currently in the priority queue
    vector<Cost> distanceToVar; // Shortest distances from the Dijkstra source variable to each other variable node in the residual graph
    vector<Cost> distanceToVal; // Shortest distances from the Dijkstra source variable to each value node in the residual graph (reverse graph)
    int Q = 0; // Number of Dijkstra runs to perform during reduced-cost filtering (derived from FiltLevel and NbNoAssigned)
    double FiltLevel; // Fraction of variables used as Dijkstra sources during filtering (0 = disabled, 1 = all variables)
    vector<int> trackList; // Variables not yet finalized by Dijkstra (dense array, used by the bimodal traversal strategy)
    vector<int> trackPos; // trackPos[var] = position of var in trackList, or -1 if finalized
    vector<PQ::handle_type> dijkstraHandles; // Handles of the priority queue, reused between two Dijkstra runs
    Cost MaxReducedCost = 0; // Maximum reduced cost observed over all arcs in the current residual graph (used as filtering threshold)
    Cost ReducedCost = 0; // Temporary reduced cost of a single arc computed during cost projection
    Cost VarMaxReducedCost = 0; // Maximum reduced cost among all arcs incident to a given variable (used to guide source selection)
    vector<int> varValIdx; // varValIdx[var * NbValues + val] = index of value val (in UnionVarDomain) inside the initial domain of variable var, or -1 if var cannot take it
    vector<int> lastConflictVars; // Indices of the variables involved in the last detected conflict (used for conflict-driven weight updates)
    int findConflict = 0; // Flag returned by the LAP solver: non-zero if a conflict (infeasible sub-problem) was detected, value encodes the number of conflicting variables

    void projectLB(Cost c)
    {
        if (c > MIN_COST) {
            if (!isSquare)
                lb += c;
            Constraint::projectLB(c);
        }
    }

    // returns true if the constraint can be projected to small local cost function in extension
    // bool canbeProjectedInExtension();

    // Depending of the value and the cost, extend or project the cost on the index value of the variable var
    void ExtOrProJ(int var, Value value, Cost C)
    {

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
    GlobalCardinalityConstraint(const GlobalCardinalityConstraint&) = delete; // owns raw arrays
    GlobalCardinalityConstraint& operator=(const GlobalCardinalityConstraint&) = delete;

    GlobalCardinalityConstraint(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in)
        : AbstractNaryConstraint(wcsp, scope_in, arity_in)
        , Original_ub(wcsp->getUb())
        , lb(0)
        , assigneddeltas(0)
        , storeAssignment(false)
        , NbValues(0)
        , NbAssigned(0)
        , NbNoAssigned(arity_in)
        , FiltLevel(0.0)
    {
        if (arity_in > 0) {
            SameDomain = true;
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

            int rate = ToulBar2::ReducedCostsFiltering;
            if (rate > 0)
                FiltLevel = rate / 100.0;

            // Initialize
            storeLastAssignment = vector<StoreValue>(arity_in, StoreValue(WRONG_VAL));
            NoAssignedVar = vector<int>(arity_in, -1);
            capacity = vector<int>(NbValues, arity_in);
            CurrentCapacity = vector<int>(NbValues, 0);
            demand = vector<int>(NbValues, 0);
            AssignedVar = vector<int>(arity_in, -1);
            AssignedVal = vector<int>(arity_in, -1);
            costMatrix = vector<Cost>((arity_in + 1) * NbValues, MAX_COST);
            ReduceCostMatrix = vector<Cost>(arity_in * NbValues, MAX_COST);

            // Allocate memory
            colSol = new int[NbValues];
            rowSol = new int[arity_ + 1];
            ReduceCostRow = new Cost[arity_ + 1];
            ReduceCostCol = new Cost[NbValues];

        } else {
            deconnect();
        }
    }

    virtual ~GlobalCardinalityConstraint()
    {
        delete[] colSol;
        delete[] rowSol;
        delete[] ReduceCostRow;
        delete[] ReduceCostCol;
    }
    void read(istream& files) // TODO: add a parameter for controlling the level of propagation if necessary
    {
        string line;
        getline(files, line);
        stringstream file(line);

        Gcc_NbValues = 0;
        sumlb = 0;
        sumub = 0;
        int lower = 0;
        int upper = 0;
        file >> Gcc_NbValues;
        bool findval;
        Cost Max_cost = wcsp->getUb();
        for (int v = 0; v < Gcc_NbValues; v++) {

            Value value;
            file >> value;
            file >> lower;
            sumlb += lower;
            file >> upper;
            if (lower > upper)
                THROWCONTRADICTION;
            sumub += upper;
            findval = false;
            bounds[value] = { lower, upper };
            for (int varIndex = 0; varIndex < arity_; ++varIndex) {
                auto* variable = scope[varIndex];
                if (variable->canbe(value) && (variable->getCost(value) < Max_cost)) {
                    int idx = unionIndex(variable->getValueName(variable->toIndex(value)));
                    if (idx < 0)
                        continue;
                    capacity[idx] = upper;
                    demand[idx] = lower;
                    findval = true;
                    break;
                }
            }
            if (!findval && lower > 0)
                THROWCONTRADICTION;
        }
        if (sumlb > arity_)
            THROWCONTRADICTION;
        if (sumub < arity_ && Gcc_NbValues == NbValues)
            THROWCONTRADICTION;

        // Values upper bound capacity update.
        int zeroCapacity = 0;
        if (sumlb > 0) {
            int newCapacity;
            int newsumub = 0;
            for (int valIndex = 0; valIndex < NbValues; valIndex++) {
                newCapacity = min((arity_ - (sumlb - demand[valIndex])), capacity[valIndex]);
                if (newCapacity <= 0) {
                    capacity[valIndex] = 0;
                    zeroCapacity++;
                    Value value;
                    for (int varIndex = 0; varIndex < arity_; varIndex++) {
                        auto* variable = scope[varIndex];
                        value = toVarValue(varIndex, valIndex);
                        if (value != WRONG_VAL && variable->canbe(value)) {
                            variable->project(value, wcsp->getUb(), true);
                        }
                    }
                } else {
                    capacity[valIndex] = capacity[valIndex] > newCapacity ? newCapacity : capacity[valIndex];
                    newsumub += capacity[valIndex];
                }
            }
            sumub = newsumub;
        } else {
            sumub = 0;
            for (int valIndex = 0; valIndex < NbValues; valIndex++) {
                sumub += capacity[valIndex]; // also counts values not listed in the gcc (capacity = arity)
                if (capacity[valIndex] <= 0) {
                    zeroCapacity++;
                    Value value;
                    for (int varIndex = 0; varIndex < arity_; varIndex++) {
                        auto* variable = scope[varIndex];
                        value = toVarValue(varIndex, valIndex);
                        if (value != WRONG_VAL && variable->canbe(value)) {
                            variable->project(value, wcsp->getUb(), true);
                        }
                    }
                }
            }
        }

        if (zeroCapacity > 0) {
            // Move the values of null capacity at the end of the union (indexes >= NbRestValues): they are then ignored
            int NbRestValues = NbValues - zeroCapacity;
            int restCapacity = 0;
            int zeroIndex = NbRestValues;

            vector<int> tempCapacity(capacity.begin(), capacity.begin() + NbValues);
            vector<int> tempDemand(demand.begin(), demand.begin() + NbValues);
            vector<string> tempUnionVarDomain(UnionVarDomain.begin(), UnionVarDomain.begin() + NbValues);

            for (int valIndex = 0; valIndex < NbValues; ++valIndex) {
                int dst = (tempCapacity[valIndex] == 0) ? zeroIndex++ : restCapacity++;
                capacity[dst] = tempCapacity[valIndex];
                demand[dst] = tempDemand[valIndex];
                UnionVarDomain[dst] = tempUnionVarDomain[valIndex];
                mapDomainValToIndex[tempUnionVarDomain[valIndex]] = dst;
            }

            NbValues = NbRestValues;
            rebuildVarMapping();
        }

        // If more values than variables, cost matrix is not square
        if ((Gcc_NbValues < NbValues && sumlb == 0) || sumub > arity_) {
            isSquare = false;
            // Initialize extended cost vector
            deltaCosts.assign(NbValues, MIN_COST);
        }
        if (!isSquare) {
            int nbDelta = 0;
            file >> nbDelta;
            if (nbDelta) {
                for (int v = 0; v < nbDelta; v++) {
                    Value value;
                    Cost delta;
                    file >> value;
                    file >> delta;
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
    bool isGCC() const FINAL { return true; }
    bool isGCCSquare() const FINAL { return isSquare; }

    const map<Value, pair<int, int>>& getBounds() const { return bounds; }

    Cost getDefCost() FINAL { return MAX_COST; }

    Long getConflictWeight() const override { return Constraint::getConflictWeight(); }
    Long getConflictWeight(int varIndex) const override
    {
        assert(varIndex >= 0);
        assert(varIndex < arity_);
        return conflictWeights[varIndex] + Constraint::getConflictWeight();
    }
    void incConflictWeight(Constraint* from) override
    {
        assert(from != NULL);
        if (from == this) {
            if (getNonAssigned() == arity_ || deconnected()) {
                if (arity_ < (int)wcsp->numberOfVariables()) { // SdG: no need to increase weightedDegree of all problem variables
                    Constraint::incConflictWeight(1);
                }
            } else {
                for (int i = 0; i < arity_; i++) {
                    // if (connected(i)) {
                    auto* variable = scope[i];

                    if (!variable->unassigned()) {

                        conflictWeights[i]++;
                    }
                }
            }
        } else if (deconnected()) {
            for (int i = 0; i < from->arity(); i++) {
                int index = getIndex(from->getVar(i));
                if (index >= 0) { // the last conflict constraint may be derived from two binary constraints (boosting search), each one derived from an n-ary constraint with a scope which does not include parameter constraint from
                    assert(index < arity_);
                    conflictWeights[index]++;
                }
            }
        }
    }
    void resetConflictWeight() override
    {
        conflictWeights.assign(conflictWeights.size(), 0);
        Constraint::resetConflictWeight();
    }

    /// \brief returns true if constraint always satisfied and has (less than) zero cost only
    // bool universal(Cost zero = MIN_COST) override;

    bool implies(Constraint* ctr) FINAL
    {
        if (scopeIncluded(ctr) && ctr->isBinary()) {
            BinaryConstraint* bctr = (BinaryConstraint*)ctr;
            EnumeratedVariable* x = (EnumeratedVariable*)(bctr->getVar(0));
            EnumeratedVariable* y = (EnumeratedVariable*)(bctr->getVar(1));
            if (x->isValueNames() && y->isValueNames()) {
                for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
                    for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
                        if ((x->getValueName(x->toIndex(*iterx)) != y->getValueName(y->toIndex(*itery)) || bounds.find(*iterx) == bounds.end() || bounds[*iterx].second >= 2) && bctr->getCost(*iterx, *itery) > MIN_COST) {
                            return false;
                        }
                    }
                }
            } else {
                for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
                    for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
                        if ((*iterx != *itery || bounds.find(*iterx) == bounds.end() || bounds[*iterx].second >= 2) && bctr->getCost(*iterx, *itery) > MIN_COST) {
                            return false;
                        }
                    }
                }
            }
            return true;
        }
        return false;
    }

    void projects(Constraint* ctr) FINAL
    {
        if (scopeIncluded(ctr) && ctr->isBinary()) {
            BinaryConstraint* bctr = (BinaryConstraint*)ctr;
            EnumeratedVariable* x = (EnumeratedVariable*)(bctr->getVar(0));
            EnumeratedVariable* y = (EnumeratedVariable*)(bctr->getVar(1));
            Cost mult_ub = (wcsp->getUb() < (MAX_COST / MEDIUM_COST)) ? (max(LARGE_COST, wcsp->getUb() * MEDIUM_COST)) : wcsp->getUb();
            if (x->isValueNames() && y->isValueNames()) {
                for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
                    string s = x->getValueName(x->toIndex(*iterx));
                    unsigned int yindex = y->toIndex(s);
                    Value yval = y->toValue(yindex);
                    if (bounds.find(*iterx) != bounds.end() && bounds[*iterx].second == 1 && y->canbe(yval) && !CUT(bctr->getCost(*iterx, yval), wcsp->getUb())) {
                        bctr->addcost(*iterx, yval, mult_ub - bctr->getCost(*iterx, yval));
                    }
                }
            } else {
                for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
                    if (bounds.find(*iterx) != bounds.end() && bounds[*iterx].second == 1 && y->canbe(*iterx) && !CUT(bctr->getCost(*iterx, *iterx), wcsp->getUb())) {
                        bctr->addcost(*iterx, *iterx, mult_ub - bctr->getCost(*iterx, *iterx));
                    }
                }
            }
        }
    }

    Cost eval(const Tuple& s) override
    {
        // returns the cost of the corresponding assignment s
        if (isSquare) {
            Cost res = MIN_COST;
            Cost nbsame = 0;
            vector<int> countUsed(NbValues, 0);
            for (int varIndex = 0; varIndex < arity_; varIndex++) {
                int valIndex = unionIndex(scope[varIndex]->getValueName(s[varIndex]));
                if (valIndex < 0) { // forbidden value (removed from the union)
                    nbsame++;
                    continue;
                }
                countUsed[valIndex]++;
                if (countUsed[valIndex] > capacity[valIndex]) {
                    nbsame++;
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

        Cost res = 0;

        Cost nbsame = 0;
        vector<int> countUsed(NbValues, 0);
        for (int varIndex = 0; varIndex < arity_; varIndex++) {
            int valIndex = unionIndex(scope[varIndex]->getValueName(s[varIndex]));
            if (valIndex < 0) { // forbidden value (removed from the union)
                nbsame++;
                continue;
            }
            countUsed[valIndex]++;
            if (countUsed[valIndex] > capacity[valIndex]) {
                nbsame++;
            }
        }
        if (sumlb > 0) {
            for (int valIndex = 0; valIndex < NbValues; valIndex++) {
                if (countUsed[valIndex] < demand[valIndex]) {
                    nbsame += (demand[valIndex] - countUsed[valIndex]);
                }
            }
        }

        if (nbsame == 0) {
            if (sumlb == 0) {
                for (int valIndex = 0; valIndex < NbValues; valIndex++) {
                    res += (capacity[valIndex] - countUsed[valIndex]) * deltaCosts[valIndex];
                }
            }
            if (sumlb > 0) {
                res = -lb + assigneddeltas;
                for (int valIndex = 0; valIndex < NbValues; valIndex++) {
                    res += countUsed[valIndex] * (deltaCosts[valIndex]);
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
        if (isSquare)
            return MIN_COST;
        Cost sumdelta = abs(-lb + assigneddeltas);
        Cost m = *max_element(deltaCosts.begin(), deltaCosts.end());
        if (m > MIN_COST)
            sumdelta += m;
        if (CUT(sumdelta, wcsp->getUb()))
            return MAX_COST;
        else
            return sumdelta;
    }

    Cost getCost(int index, Value val)
    {
        assert(index >= 0 && index < arity_);
        if (isSquare)
            return MIN_COST;
        // deltaCosts is indexed by the position of the value name in UnionVarDomain
        int idx = unionIndex(scope[index]->getValueName(scope[index]->toIndex(val)));
        if (idx < 0)
            return MIN_COST;
        return deltaCosts[idx];
    }

    double computeTightness() override { return 0; } // TODO: compute factorial(n)/(n**n)?

    // TODO: needed for dominance test by DEE
    // pair<pair<Cost, Cost>, pair<Cost, Cost>> getMaxCost(int index, Value a, Value b)

    // void setInfiniteCost(Cost ub)
    void setInfiniteCost(Cost ub) override
    {
        Original_ub = min(ub, Original_ub);
    }

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
                propagate();
                if (ToulBar2::FullEAC)
                    reviseEACGreedySolution();
            }
        }
    }

    /**
     * @brief Filters and removes assigned values from unassigned variables in the GCC constraint.
     *
     * This function enforces the GCC constraint by:
     * - Verifying that assigned variables do not share the same value.
     * - Removing values from unassigned variable domains if already assigned elsewhere.
     * - Propagating newly assigned variables and updating internal state.
     * - Triggering full n-ary projection if thresholds on domain sizes or variable count are met.
     *
     * @return true if constraint propagation was successful without needing full projection.
     * @return false if the constraint switched to full n-ary projection.
     *
     * @throws Throws a contradiction exception (THROWCONTRADICTION) if two variables are assigned the same value.
     */
    bool RemoveAssignVar()
    {
        varAlreadyProcessed.assign(arity_, false);

        // Clear and reserve space for assigned and unassigned variables containers
        AssignedVar.clear();
        AssignedVar.reserve(arity_);
        NoAssignedVar.clear();
        NoAssignedVar.reserve(arity_);
        AssignedVal.clear();
        AssignedVal.reserve(arity_);

        // Capacity of a value exceeded: increase the weight of the variable and of all the variables already using this value
        auto conflictOn = [&](int varIndex, int valIndex) {
            conflictWeights[varIndex]++;
            for (size_t k = 0; k < AssignedVar.size(); ++k) {
                if (AssignedVal[k] == valIndex)
                    conflictWeights[AssignedVar[k]]++;
            }
        };

        // Values without remaining capacity that still have to be removed from the unassigned variables
        vector<int> saturated;
        // Identify assigned and unassigned variables
        for (int varIndex = 0; varIndex < arity_; ++varIndex) {
            auto* variable = scope[varIndex];

            if (variable->unassigned()) {
                // Variable not assigned yet: add to unassigned list
                NoAssignedVar.push_back(varIndex);
            } else {
                // Variable assigned: find the index of its assigned value
                int valIndex = unionIndex(variable->getValueName(variable->toIndex(variable->getValue())));
                if (valIndex < 0) { // assigned to a forbidden value
                    conflictWeights[varIndex]++;
                    THROWCONTRADICTION;
                }
                CurrentCapacity[valIndex] -= 1;
                if (CurrentCapacity[valIndex] > -1) {
                    AssignedVar.push_back(varIndex);
                    AssignedVal.push_back(valIndex);
                    if (CurrentCapacity[valIndex] == 0)
                        saturated.push_back(valIndex); // this value is now saturated
                } else {
                    conflictOn(varIndex, valIndex);
                    THROWCONTRADICTION;
                }
            }
        }

        NbAssigned = AssignedVar.size();
        bool VarAssigned;
        bool NaryPro = false;

        // Propagation loop: remove saturated values from unassigned variables' domains
        // (a value only needs to be processed once: it is processed again only if it becomes saturated later)
        do {
            VarAssigned = false;
            vector<int> newlyAssignedVars;

            for (int valIndex : saturated) {
                for (int varIndex : NoAssignedVar) {
                    if (varAlreadyProcessed[varIndex])
                        continue; // Skip already processed vars

                    auto* variable = scope[varIndex];
                    Value value = toVarValue(varIndex, valIndex);

                    if (value != WRONG_VAL && variable->canbe(value)) {
                        // Remove the saturated value from the domain of the unassigned variable
                        variable->remove(value);
                        if (variable->assigned()) {
                            // If variable becomes assigned after removal, mark for propagation
                            newlyAssignedVars.push_back(varIndex);
                            varAlreadyProcessed[varIndex] = true;
                        }
                    }
                }
            }
            saturated.clear();

            // Process newly assigned variables after domain filtering
            for (int varIndex : newlyAssignedVars) {
                if (connected(varIndex)) {
                    // Disconnect variable from constraint and update counts
                    deconnect(varIndex);
                    NbNoAssigned--;
                    NbAssigned++;

                    // Check whether to switch to full n-ary projection based on thresholds
                    if (getNonAssigned() <= NARYPROJECTIONSIZE && (getNonAssigned() <= 1 || prodInitDomSize <= NARYPROJECTIONPRODDOMSIZE || maxInitDomSize <= NARYPROJECTION3MAXDOMSIZE || (getNonAssigned() == 2 && maxInitDomSize <= NARYPROJECTION2MAXDOMSIZE))) {
                        deconnect();
                        projectNary();
                        NaryPro = true;
                        break;
                    } else {
                        // Update assigned values and variables tracking
                        VarAssigned = true;
                        auto* variable = scope[varIndex];
                        int valIndex = unionIndex(variable->getValueName(variable->toIndex(variable->getValue())));
                        if (valIndex < 0) { // assigned to a forbidden value
                            conflictWeights[varIndex]++;
                            THROWCONTRADICTION;
                        }
                        CurrentCapacity[valIndex] -= 1;
                        if (CurrentCapacity[valIndex] > -1) {
                            AssignedVar.push_back(varIndex);
                            AssignedVal.push_back(valIndex);
                            if (CurrentCapacity[valIndex] == 0)
                                saturated.push_back(valIndex); // this value is now saturated
                        } else {
                            conflictOn(varIndex, valIndex);
                            THROWCONTRADICTION;
                        }
                    }
                } else {
                    return false;
                }
            }

            // The constraint has been fully projected: stop here
            if (NaryPro)
                return false;

            // Update list of unassigned variables after propagation
            if (VarAssigned) {
                vector<int> updatedNoAssigned;
                for (int varIndex : NoAssignedVar) {
                    if (!scope[varIndex]->assigned()) {
                        updatedNoAssigned.push_back(varIndex);
                    }
                }
                NoAssignedVar.swap(updatedNoAssigned);
            }

        } while (VarAssigned);

        return true;
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
     *
     *      Source : Bimodal Depth-First Search for Scalable GAC for AllDifferent.
     *      Sulian Le Bozec Chiffoleau; Nicolas Beldiceanu; Charles Prud'homme;
     *      Gilles Simonin; and Xavier Lorca, Proceedings of the Thirty-Fourth
     *      International Joint Conference on Artificial Intelligence, IJCAI 2025,
     *      Montreal, Canada, August 16-22, 2025.
     *
     */
    void BimodalDijkstra(int source, vector<vector<int>>& VarList, vector<vector<uint8_t>>& ValList)
    {
        if (dijkstraHandles.size() < (size_t)NbNoAssigned)
            dijkstraHandles.resize(NbNoAssigned);
        auto& handles = dijkstraHandles;
        PQ pq;

        distanceToVar.assign(NbNoAssigned, MAX_COST);
        visited.assign(NbNoAssigned, 0);
        visitVal.assign(NbNoAssignedVal, -1);

        // Variables not yet finalized: dense array with O(1) insertion/removal (swap with the last element)
        trackList.resize(NbNoAssigned);
        trackPos.resize(NbNoAssigned);
        iota(trackList.begin(), trackList.end(), 0);
        iota(trackPos.begin(), trackPos.end(), 0);
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

        inHeap.assign(NbNoAssigned, 0);
        distanceToVar[source] = 0;
        handles[source] = pq.push({ distanceToVar[source], source });
        inHeap[source] = 1;

        while (!pq.empty()) {
            auto [d, var] = pq.top();
            pq.pop();
            visited[var] = 1;
            trackErase(var);
            int val = rowSol[var];
            auto& varlist = VarList[val];

            if (visitVal[val] == -1) {
                visitVal[val] = var;
            } else {
                distanceToVar[var] = distanceToVar[visitVal[val]];
                continue;
            }
            if (varlist.size() < trackList.size()) {
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
                        if (alt == d) {
                            trackErase(nextVar);
                            visited[nextVar] = 1;
                        }
                    }
                }
            } else {
                size_t k = 0;
                while (k < trackList.size()) {
                    int nextVar = trackList[k];
                    if (ValList[nextVar][val]) {
                        Cost alt = d + ReduceCostMatrix[nextVar * NbNoAssignedVal + val];
                        if (alt < distanceToVar[nextVar]) {
                            distanceToVar[nextVar] = alt;
                            if (!inHeap[nextVar]) {
                                handles[nextVar] = pq.push({ alt, nextVar });
                                inHeap[nextVar] = 1;
                            } else {
                                pq.update(handles[nextVar], { alt, nextVar });
                            }
                            if (alt == d) {
                                trackErase(nextVar); // the last element now sits at position k
                                visited[nextVar] = 1;
                                continue;
                            }
                        }
                    }
                    ++k;
                }
            }
        }
    }

    /**
     * @brief Propagate the Global Cardinality constraint to enforce consistency and improve problem lower bound.
     *
     * This method performs constraint propagation by checking assigned variables,
     * removing inconsistent values, computing lower bounds using the Jonker algorithm,
     * and updating unary costs and supports accordingly.
     *
     * It handles different cases depending on whether variables are assigned or not,
     * if the domains are the same.
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
        if (connected()) {
            if (ToulBar2::verbose >= 3) {
                cout << "propagate " << *this << endl;
            }

            bool b = false;

            // Attempt to assign variables already connected and assigned
            for (int varIndex = 0; !b && connected() && varIndex < arity_; varIndex++) {
                if (connected(varIndex) && scope[varIndex]->assigned()) {
                    assign(varIndex);
                    b = true;
                }
            }

            if (!b && connected()) {
                // Check whether propagation should be skipped based on stored results
                bool skipPropagation = false;
                if (storeAssignment) {
                    skipPropagation = true;
                    for (int varIndex = 0; varIndex < arity_; varIndex++) {
                        if (scope[varIndex]->cannotbe(storeLastAssignment[varIndex]) || scope[varIndex]->getCost(storeLastAssignment[varIndex]) > MIN_COST) {
                            skipPropagation = false;
                            break;
                        }
                    }
                }

                if (!skipPropagation) {
                    // Initialize number of unassigned variables
                    NbNoAssigned = getNonAssigned();
                    // If there are assigned variables, handle accordingly
                    if (NbNoAssigned < arity_) {
                        for (int valIndex = 0; valIndex < NbValues; ++valIndex) {
                            CurrentCapacity[valIndex] = capacity[valIndex];
                        }
                        if (SameDomain && RemoveAssignVar()) {
                            // Collect unassigned values
                            vector<int> NoAssignedVal;
                            vector<int> NoAsscapacity;
                            vector<int> NoAssdemand;
                            int lbcap;
                            for (int valIndex = 0; valIndex < NbValues; ++valIndex) {
                                if (CurrentCapacity[valIndex] > 0) {
                                    NoAssignedVal.push_back(valIndex);
                                    NoAsscapacity.push_back(CurrentCapacity[valIndex]);
                                    if (sumlb > 0) {
                                        lbcap = demand[valIndex] - (capacity[valIndex] - CurrentCapacity[valIndex]);
                                        NoAssdemand.push_back(lbcap > 0 ? lbcap : 0);
                                    } else {
                                        NoAssdemand.push_back(0);
                                    }
                                }
                            }

                            NbNoAssignedVal = NoAssignedVal.size();
                            if (NbNoAssignedVal == 0)
                                THROWCONTRADICTION;

                            // Initialize cost matrix for the Jonker algorithm
                            Cost current_ub = wcsp->getUb() - wcsp->getLb();
                            const bool doFilter = FiltLevel > 0 && isSquare;

                            for (int varInd = 0; varInd < NbNoAssigned; ++varInd) {
                                int varIndex = NoAssignedVar[varInd];
                                auto* variable = scope[varIndex];

                                for (int valInd = 0; valInd < NbNoAssignedVal; ++valInd) {
                                    int valIndex = NoAssignedVal[valInd];
                                    Value val = toVarValue(varIndex, valIndex);
                                    if (val != WRONG_VAL && variable->canbe(val)) {
                                        auto valCost = variable->getCost(val);
                                        costMatrix[varInd * NbNoAssignedVal + valInd] = valCost < current_ub ? valCost : current_ub;
                                    } else {
                                        costMatrix[varInd * NbNoAssignedVal + valInd] = current_ub;
                                    }
                                }
                            }
                            if(any_of(costMatrix.begin(),costMatrix.begin() + arity_ * NbNoAssignedVal,  [](auto x) { return x > MIN_COST; })){
                                // Solve assignment problem using Jonker algorithm
                                Cost TotalCost;
                                if (!isSquare) {

                                    if (sumlb == 0) {
                                        for (int valInd = 0; valInd < NbNoAssignedVal; ++valInd) {
                                            costMatrix[NbNoAssigned * NbNoAssignedVal + valInd] = deltaCosts[NoAssignedVal[valInd]];
                                        }
                                        TotalCost = lapjv_ub(NbNoAssigned + 1, NbNoAssignedVal, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub, NoAsscapacity, findConflict);

                                    } else {
                                        TotalCost = lapjv_gcc(NbNoAssigned, NbNoAssignedVal, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub, NoAsscapacity, NoAssdemand, findConflict);
                                    }
                                } else {
                                    TotalCost = lapjv_ub(NbNoAssigned, NbNoAssignedVal, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub, NoAsscapacity, findConflict);
                                }

                                if (TotalCost >= current_ub) {
                                    if (findConflict) {
                                        int varIndex;
                                        for (int var = 0; var < findConflict; var++) {
                                            if (rowSol[var] >= NbNoAssigned)
                                                continue;
                                            varIndex = NoAssignedVar[rowSol[var]];
                                            conflictWeights[varIndex]++;
                                        }
                                    }
                                    wcsp->revise(this);
                                    THROWCONTRADICTION;

                                } else if (TotalCost >= 0) {
                                    Cost jonker = TotalCost;

                                    // Store results from Jonker algorithm
                                    for (int varIndex = 0; varIndex < NbNoAssigned; ++varIndex) {
                                        storeLastAssignment[NoAssignedVar[varIndex]] = toVarValue(NoAssignedVar[varIndex], NoAssignedVal[rowSol[varIndex]]);
                                    }

                                    for (int varIndex = 0; varIndex < NbAssigned; ++varIndex) {
                                        storeLastAssignment[AssignedVar[varIndex]] = toVarValue(AssignedVar[varIndex], AssignedVal[varIndex]);
                                    }

                                    storeAssignment = true;

                                    if (TotalCost > 0 || (any_of(ReduceCostRow, ReduceCostRow + NbNoAssigned,  [](auto x) { return x > MIN_COST; }))) {

                                        if (!isSquare) {
                                            int valIndex;
                                            if (sumlb == 0) {
                                                Cost mindelta = ReduceCostRow[NbNoAssigned];

                                                for (int valInd = 0; valInd < NbNoAssignedVal; ++valInd) {
                                                    valIndex = NoAssignedVal[valInd];
                                                    if (capacity[valIndex] > 0)
                                                        deltaCosts[valIndex] -= (ReduceCostCol[valInd] + mindelta);
                                                }
                                                if (mindelta > 0) {

                                                    jonker += mindelta * (sumub - arity_ - 1);

                                                    if (jonker >= current_ub) {
                                                        wcsp->revise(this);
                                                        THROWCONTRADICTION;
                                                    }
                                                }
                                            } else {
                                                Cost assignedPart = MIN_COST;
                                                for (int valInd = 0; valInd < NbNoAssignedVal; ++valInd) {
                                                    valIndex = NoAssignedVal[valInd];
                                                    if (capacity[valIndex] > 0) {
                                                        deltaCosts[valIndex] += ReduceCostCol[valInd];
                                                        assignedPart += ReduceCostCol[valInd] * (capacity[valIndex] - CurrentCapacity[valIndex]);
                                                    }
                                                }
                                                assigneddeltas -= assignedPart;
                                            }
                                        }
                                        // Update the lower bound with the Jonker algorithm's total cost
                                        projectLB(jonker);
                                        Cost lastCurrent_ub = current_ub;
                                        current_ub -= jonker;

                                        // Adjust unary costs for unassigned variables based on reduced costs

                                        MaxReducedCost = 0;
                                        vector<vector<int>> VarList(NbNoAssignedVal);
                                        vector<vector<uint8_t>> ValList(NbNoAssigned);
                                        vector<int> usedValList;
                                        vector<char> usedFlag;
                                        for (int varInd = 0; varInd < NbNoAssigned; ++varInd) {
                                            int varIndex = NoAssignedVar[varInd];
                                            assigneddeltas += ReduceCostRow[varInd];
                                            VarMaxReducedCost = 0;
                                            if (doFilter)
                                                ValList[varInd].assign(NbNoAssignedVal, 0);
                                            for (int valInd = 0; valInd < NbNoAssignedVal; ++valInd) {
                                                int valIndex = NoAssignedVal[valInd];
                                                Value value = toVarValue(varIndex, valIndex);
                                                if (value != WRONG_VAL && costMatrix[varInd * NbNoAssignedVal + valInd] < lastCurrent_ub) {
                                                    ExtOrProJ(varIndex, value, (ReduceCostRow[varInd] + ReduceCostCol[valInd]));
                                                    if (doFilter) {
                                                        ReducedCost = costMatrix[varInd * NbNoAssignedVal + valInd] - (ReduceCostRow[varInd] + ReduceCostCol[valInd]);
                                                        ReduceCostMatrix[varInd * NbNoAssignedVal + valInd] = ReducedCost;
                                                        if (ReducedCost > VarMaxReducedCost)
                                                            VarMaxReducedCost = ReducedCost;
                                                        if (rowSol[varInd] != valInd) {
                                                            VarList[valInd].push_back(varInd);
                                                            ValList[varInd][valInd] = 1;
                                                        }
                                                    }
                                                }
                                            }
                                            if (MaxReducedCost < current_ub)
                                                MaxReducedCost += VarMaxReducedCost;
                                        }

                                        if (doFilter && MaxReducedCost >= current_ub) {
                                            usedFlag.assign(NbNoAssignedVal, 0);
                                            for (int varInd = 0; varInd < NbNoAssigned; ++varInd) {
                                                if (!usedFlag[rowSol[varInd]]) {
                                                    usedFlag[rowSol[varInd]] = 1;
                                                    usedValList.push_back(rowSol[varInd]);
                                                    colSol[rowSol[varInd]] = varInd;
                                                }
                                            }
                                        }

                                        /* (BEGIN) : Filtering of variables domains with Sellmann or Cambazard method

                                        Source : Claus, G.; Cambazard, H.; and Jost, V. 2020. Analysis of Re-
                                                duced Costs Filtering for Alldifferent and Minimum Weight
                                                Alldifferent Global Constraints. In Giacomo, G. D.; Catal´a,
                                                A.; Dilkina, B.; Milano, M.; Barro, S.; Bugar´ın, A.; and
                                                Lang, J., eds., ECAI 2020 - 24th European Conference on
                                                Artificial Intelligence, volume 325, 323–330. Santiago de
                                                Compostela, Spain: IOS Press.  */

                                        if (doFilter && MaxReducedCost >= current_ub) {

                                            int position;
                                            vector<int> VariableList(NbNoAssigned);
                                            iota(VariableList.begin(), VariableList.end(), 0);
                                            if (FiltLevel >= 1. - (double)ToulBar2::epsilon) {
                                                Q = NbNoAssigned;
                                            } else {
                                                Q = min(NbNoAssigned, 1 + static_cast<int>(NbNoAssigned * FiltLevel));
                                                myrearrange(VariableList);
                                            }
                                            AlreadyUse.assign(NbNoAssignedVal, 0);
                                            int varInde;
                                            for (int ind = 0; ind < Q; ++ind) {
                                                varInde = VariableList[ind];
                                                if (AlreadyUse[rowSol[varInde]] == 0) {
                                                    AlreadyUse[rowSol[varInde]] = 1;
                                                } else {
                                                    continue;
                                                }
                                                // Dijkstra(varInde, VarList);
                                                // Bimodal Dijkstra’s shortest path algorithm
                                                BimodalDijkstra(varInde, VarList, ValList);

                                                for (int valInd : usedValList) {

                                                    if (distanceToVar[colSol[valInd]] >= MAX_COST)
                                                        continue;
                                                    position = -1;
                                                    auto& varlist = VarList[valInd];
                                                    for (int varInd : varlist) {
                                                        position++;
                                                        if (distanceToVar[varInd] == MAX_COST || costMatrix[varInd * NbNoAssignedVal + valInd] >= current_ub)
                                                            continue;
                                                        Cost reducedCost = ReduceCostMatrix[varInd * NbNoAssignedVal + valInd] - distanceToVar[varInd] + distanceToVar[colSol[valInd]];

                                                        if (reducedCost >= current_ub) {
                                                            int varIndex = NoAssignedVar[varInd];
                                                            int valIndex = NoAssignedVal[valInd];
                                                            auto* variable = scope[varIndex];
                                                            varlist[position] = -1;
                                                            ValList[varInd][valInd] = 0;

                                                            Value value = toVarValue(varIndex, valIndex);

                                                            if (value != WRONG_VAL && variable->canbe(value)) {
                                                                if (ToulBar2::verbose > 0)
                                                                    cout << "REMOVE VALUE " << value << " from " << variable->getName() << endl;
                                                                ExtOrProJ(varIndex, value, -current_ub); // SdG: project infinite cost on this value and avoid to skip and reenter the GCC constraint without finishing the current filtering
                                                            }
                                                        }
                                                    }
                                                    varlist.erase(std::remove(varlist.begin(), varlist.end(), -1), varlist.end());
                                                }
                                            }
                                        }
                                        /* (END) : Filtering of variables domains with Sellmann or Cambazard method */
                                    }
                                    // Update support values if needed for unassigned variables
                                    for (int varInd = 0; varInd < NbNoAssigned; ++varInd) {
                                        int varIndex = NoAssignedVar[varInd];
                                        auto* variable = scope[varIndex];
                                        Value optimalValue = storeLastAssignment[varIndex];
                                        if (variable->getSupport() != optimalValue) {
                                            if (ToulBar2::verbose > 0)
                                                cout << "CHANGE GCC SUPPORT " << variable->getName() << " from " << variable->getSupport() << " to " << optimalValue << endl;
    #ifndef NDEBUG
                                            variable->queueEAC1(); // EAC support may have been lost
    #endif
                                            variable->setSupport(optimalValue);
                                            assert(variable->getCost(variable->getSupport()) == MIN_COST);
                                        }
                                    }
                                }
                            } else if(!isSquare){
                               // cout << "GCC SASTISFY\n";
                                bool sastisfy = true;
                                for(int valIndex = 0; valIndex < NbNoAssignedVal; valIndex ++){
                                    if(NoAsscapacity[valIndex] < NbNoAssigned){
                                        sastisfy = false;
                                        break;
                                    }
                                }
                                if(sastisfy){
                                    if (ToulBar2::verbose > 0)
                                        cout << "GCC SASTISFY\n";
                                    deconnected();
                                }
                            }
                            
                        }
                    }

                    if (connected() && (NbNoAssigned == arity_ || (!SameDomain))) {

                        Cost current_ub = wcsp->getUb() - wcsp->getLb();
                        // Values of the union that do not belong to a variable domain are forbidden
                        fill(costMatrix.begin(), costMatrix.begin() + (size_t)arity_ * NbValues, current_ub);
                        // Initialize the cost matrix for all variables and their domain values
                        for (int varIndex = 0; varIndex < arity_; ++varIndex) {
                            auto* variable = scope[varIndex];

                            for (int valIndex = 0; valIndex < NbValues; ++valIndex) {
                                Value value = toVarValue(varIndex, valIndex);
                                if (value == WRONG_VAL)
                                    continue; // the variable cannot take this value: stays at current_ub (see fill above)
                                Cost valCost = variable->canbe(value) ? variable->getCost(value) : current_ub;
                                costMatrix[varIndex * NbValues + valIndex] = valCost < current_ub ? valCost : current_ub;
                            }
                        }

                        if(any_of(costMatrix.begin(),costMatrix.begin() + arity_ * NbValues,  [](auto x) { return x > MIN_COST; })){
                            // Solve the Linear Assignment Problem (LAP) using the Jonker algorithm
                            Cost TotalCost;
                            if (!isSquare) {

                                if (sumlb == 0) {
                                    for (int valInd = 0; valInd < NbValues; ++valInd) {
                                        costMatrix[arity_ * NbValues + valInd] = deltaCosts[valInd];
                                    }
                                    TotalCost = lapjv_ub(arity_ + 1, NbValues, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub, capacity, findConflict);

                                } else {
                                    TotalCost = lapjv_gcc(arity_, NbValues, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub, capacity, demand, findConflict);
                                }
                            } else {
                                TotalCost = lapjv_ub(arity_, NbValues, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub, capacity, findConflict);
                            }

                            if (TotalCost >= current_ub) {
                                if (findConflict) {
                                    int varIndex;
                                    for (int var = 0; var < findConflict; var++) {
                                        varIndex = rowSol[var];
                                        if (varIndex >= arity_)
                                            continue;
                                        conflictWeights[varIndex]++;
                                    }
                                }
                                wcsp->revise(this);
                                THROWCONTRADICTION;

                            } else if (TotalCost > 0 || (any_of(ReduceCostRow, ReduceCostRow + arity_,  [](auto x) { return x > MIN_COST; }))) {
                                // Store results from Jonker algorithm
                                for (int varIndex = 0; varIndex < arity_; varIndex++) {
                                    storeLastAssignment[varIndex] = toVarValue(varIndex, rowSol[varIndex]);
                                }

                                storeAssignment = true;

                                if (TotalCost > 0 || (any_of(ReduceCostRow, ReduceCostRow + arity_,  [](auto x) { return x > MIN_COST; }))) {
                                    if (!isSquare) {
                                        if (sumlb == 0) {
                                            Cost mindelta = ReduceCostRow[arity_];
                                            for (int valInd = 0; valInd < NbValues; ++valInd) {
                                                if (capacity[valInd] > 0)
                                                    deltaCosts[valInd] -= (ReduceCostCol[valInd] + mindelta);
                                            }

                                            if (mindelta > 0) {
                                                TotalCost += mindelta * (sumub - arity_ - 1);
                                                if (TotalCost >= current_ub) {
                                                    wcsp->revise(this);
                                                    THROWCONTRADICTION;
                                                }
                                            }

                                        } else {
                                            for (int valInd = 0; valInd < NbValues; ++valInd) {
                                                if (capacity[valInd] > 0)
                                                    deltaCosts[valInd] += ReduceCostCol[valInd];
                                            }
                                        }
                                    }

                                    // Update lower bound with the total cost found
                                    projectLB(TotalCost);

                                    // Adjust unary costs using reduced row and column costs
                                    for (int varIndex = 0; varIndex < arity_; ++varIndex) {
                                        assigneddeltas += ReduceCostRow[varIndex];
                                        for (int valIndex = 0; valIndex < NbValues; ++valIndex) {
                                            Value value = toVarValue(varIndex, valIndex);
                                            if (value == WRONG_VAL)
                                                continue;
                                            if (costMatrix[varIndex * NbValues + valIndex] < current_ub) {
                                                ExtOrProJ(varIndex, value, (ReduceCostRow[varIndex] + ReduceCostCol[valIndex]));
                                            }
                                        }
                                    }

                                    // Update support if needed for all variables
                                    for (int varIndex = 0; varIndex < arity_; ++varIndex) {
                                        auto* variable = scope[varIndex];
                                        Value optimalValue = storeLastAssignment[varIndex];
                                        if (variable->getSupport() != optimalValue) {
                                            if (ToulBar2::verbose > 0)
                                                cout << "CHANGE GCC SUPPORT " << variable->getName() << " from " << variable->getSupport() << " to " << optimalValue << endl;
        #ifndef NDEBUG
                                            variable->queueEAC1(); // EAC support may have been lost
        #endif
                                            variable->setSupport(optimalValue);
                                            assert(variable->getCost(variable->getSupport()) == MIN_COST);
                                        }
                                    }
                                }
                            }
                        }  else if(!isSquare){
                                //cout << "GCC SASTISFY\n";
                                bool sastisfy = true;
                                for(int valIndex = 0; valIndex < NbValues; valIndex ++){
                                    if(capacity[valIndex] < arity_){
                                        sastisfy = false;
                                        break;
                                    }
                                }
                                if(sastisfy){
                                    if (ToulBar2::verbose > 0)
                                        cout << "GCC SASTISFY\n";
                                    deconnected();
                                }
                            }
                    }
                }
            }
        }
    }

    // checks that the constraint is still satisfiable (called by WCSP::verify in Debug mode at each search node)
    bool verify() override
    {
        if (!storeAssignment) {
            return false;
        }
        map<Value, int> count;
        for (int i = 0; i < arity_; i++) {
            if (count.find(storeLastAssignment[i]) != count.end()) {
                count[storeLastAssignment[i]]++;
            } else {
                count[storeLastAssignment[i]] = 1;
            }
            if (scope[i]->cannotbe(storeLastAssignment[i]) || scope[i]->getCost(storeLastAssignment[i]) > MIN_COST) {
                if (scope[i]->cannotbe(storeLastAssignment[i])) {
                    cout << "variable " << scope[i]->getName() << " value " << storeLastAssignment[i] << " has been removed!" << endl;
                } else if (scope[i]->getCost(storeLastAssignment[i]) > MIN_COST) {
                    cout << "variable " << scope[i]->getName() << " value " << storeLastAssignment[i] << " has nonzero cost!" << endl;
                }
                return false;
            }
        }
        for (auto e : bounds) {
            int mycount = 0;
            if (count.find(e.first) != count.end()) {
                mycount = count[e.first];
            }
            if (mycount < e.second.first || mycount > e.second.second) {
                if (mycount < e.second.first) {
                    cout << " value " << e.first << " has underpassed the cardinality constraint of " << e.second.first << endl;
                } else if (mycount > e.second.second) {
                    cout << " value " << e.first << " has overpassed the cardinality constraint of " << e.second.second << endl;
                }
                return false;
            }
        }
        return true;
    }

    void increase(int index) override
    {
        remove(index);
    }
    void decrease(int index) override
    {
        remove(index);
    }
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

    void print(ostream& os) override
    {
        os << this << " gcc(";

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
            for (unsigned int j = 0; j < deltaCosts.size(); j++) {
                os << deltaCosts[j];
                if (j < deltaCosts.size() - 1)
                    os << "|";
            }
            os << ") ";
        }
        os << "/" << getTightness();
        if (ToulBar2::weightedDegree) {
            os << "/" << getConflictWeight();
            for (int i = 0; i < arity_; i++) {
                os << "," << conflictWeights[i];
            }
        }
        os << " " << bounds.size() << ":[";
        for (const auto& [key, bound] : bounds) {
            os << " (" << key;
            os << "," << bound.first;
            os << "," << bound.second;
            os << " )";
        }
        os << "] arity: " << arity_;
        os << " unassigned: " << getNonAssigned() << "/" << unassigned_ << endl;
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

    void dump(ostream& os, bool original = true) override
    {
        if (original) {
            os << arity_;
            for (int i = 0; i < arity_; i++)
                os << " " << scope[i]->wcspIndex;
            os << " -1 gcc ";

            os << bounds.size();
            for (const auto& [key, bound] : bounds) {
                os << " " << key;
                os << " " << bound.first;
                os << " " << bound.second;
            }
        } else {
            os << getNonAssigned();
            for (int i = 0; i < arity_; i++)
                if (scope[i]->unassigned())
                    os << " " << scope[i]->getCurrentVarId();
            os << " -1 gcc ";

            int compt = 0;
            for (const auto& [key, bound] : bounds) {
                string name = scope[0]->getValueName(scope[0]->toIndex(key));
                int idx = unionIndex(name);
                if (idx < 0)
                    continue;
                if (capacity[idx] > 0)
                    compt++;
            }

            os << compt;
            for (const auto& [key, bound] : bounds) {
                string name = scope[0]->getValueName(scope[0]->toIndex(key));
                int idx = unionIndex(name);
                if (idx < 0)
                    continue;
                if (capacity[idx] == 0)
                    continue;
                os << " " << key;
                os << " " << bound.first;
                os << " " << bound.second;
            }
        }

        if (!isSquare) {
            map<Value, Cost> mapValuesdeltaCosts = nonZeroDeltaCosts();
            if (!mapValuesdeltaCosts.empty()) {
                int current_val = (int)mapValuesdeltaCosts.size();
                os << " " << current_val;
                for (auto& [key, delta] : mapValuesdeltaCosts) {
                    os << " " << key;
                    os << " " << delta;
                }
            } else {
                os << " " << 0;
            }
        } else {
            os << " " << 0;
        }
        os << endl;
    }

    void dump_CFN(ostream& os, bool original = true) override
    {
        bool printed = false;
        os << "\"F_";

        if (original) {
            printed = false;
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
            os << "],\"type\":\"gcc\",\"params\":{\"bounds\":[";
            printed = false;
            for (const auto& [key, bound] : bounds) {
                if (printed)
                    os << ",";
                os << "[" << key;
                os << "," << bound.first;
                os << "," << bound.second;
                os << "]";
                printed = true;
            }
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
            os << "],\"type\":\"gcc\",\"params\":{\"bounds\":[";
            printed = false;
            for (const auto& [key, bound] : bounds) {
                string name = scope[0]->getValueName(scope[0]->toIndex(key));
                int idx = unionIndex(name);
                if (idx < 0)
                    continue;
                if (capacity[idx] == 0)
                    continue;
                if (printed)
                    os << ",";
                os << "[" << key;
                os << "," << bound.first;
                os << "," << bound.second;
                os << "]";
                printed = true;
            }
        }

        os << "],\"deltaCosts\":[";
        if (isSquare) {
            os << "]}},\n";
        } else {
            map<Value, Cost> mapValuesdeltaCosts = nonZeroDeltaCosts();
            if (!mapValuesdeltaCosts.empty()) {
                printed = false;
                for (auto& [key, delta] : mapValuesdeltaCosts) {
                    if (printed)
                        os << ",";
                    os << "[" << key;
                    os << "," << delta;
                    os << "]";
                    printed = true;
                }
            }
            os << "]}},\n";
        }
    }
};

#endif /*TB2GCC_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
