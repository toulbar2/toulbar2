/** \file tb2gcc.hpp
 *  \brief Propagates a Global Cardinality Constraint.
 *
 *  It maintains different levels of consistency:
 *  - tests if the constraint is always satisfied or violated when enough variables are assigned (NC level)
 *  - projects the optimal cost of the underlying linear assignment problem (LAP) onto the lower bound and onto the unary costs
 *  - GAC by strongly connected components and reduced-cost filtering on the residual graph (only in the "square" case)
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
#include <functional>
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
#include <queue>
#include <cassert>
#include <random>
#include <algorithm>

using namespace std;
using PQ = boost::heap::pairing_heap<pair<Cost, int>, boost::heap::compare<greater<>>>;

class GlobalCardinalityConstraint : public AbstractNaryConstraint {
    Cost Original_ub; // Initial upper bound when creating the constraint
    StoreCost lb; // Projected cost to problem lower bound (if zero, all deltaCosts must be zero)
    StoreCost assigneddeltas; // Accumulated deltas from assigned values (used in cost propagation)
    StoreCost sumcostsprojectonvar; // Accumulated costs projected on variables (case: non square and sumlb > 0)
    vector<Long> conflictWeights; // Used by weighted degree heuristics to prioritize variables
    vector<StoreCost> deltaCosts; // Extended unary costs to all values
    vector<StoreValue> storeLastAssignment; // Stores the last optimal assignment of variables
    StoreInt storeAssignment; // True if an optimal assignment is currently stored
    int* rowSol = nullptr; // Row solution array: rowSol[var] = value index assigned to variable var by the LAP solver
    Cost* ReduceCostRow = nullptr; // Reduced costs per row (per variable) after solving the LAP
    Cost* ReduceCostCol = nullptr; // Reduced costs per column (per value) after solving the LAP
    vector<Cost> costMatrix; // Flattened row-major cost matrix passed to the Jonker LAP solver
    int NbValues = 0; // Total number of distinct values across all variable domains (size of UnionVarDomain)
    vector<int> AssignedVar; // Indices (in scope) of currently assigned variables
    vector<int> AssignedVal; // Value indices (in UnionVarDomain) assigned to each variable in AssignedVar
    vector<int> NoAssignedVar; // Indices (in scope) of currently unassigned variables
    int NbAssigned; // Number of currently assigned variables
    int NbNoAssigned; // Number of currently unassigned variables
    int NbNoAssignedVal = 0; // Number of values with remaining capacity (columns in the reduced cost matrix)
    vector<bool> varAlreadyProcessed; // varAlreadyProcessed[i] = true if variable i has already been processed during domain filtering
    bool SameDomain = true; // True if all variables have the same domain
    map<Value, pair<int, int>> bounds; // lower and upper bound capacities for every value
    map<Value, pair<int, int>> boundsTemp; // lower and upper bound capacities for every value (as read in the file)
    vector<int> capacity; // upper bound capacities for every value
    vector<int> CurrentCapacity; // remaining upper bound capacities for every value (updated during propagation)
    vector<int> demand; // lower bound capacities for every value
    int sumub = 0, sumlb = 0, Gcc_NbValues = 0; // sum of upper bound capacities, sum of lower bound capacities and number of values given in the gcc.
    bool isSquare = true; // True if the cost matrix is square (sum of upper capacities == number of variables)
    vector<string> UnionVarDomain; // Sorted union of all variable domains (unique value names)
    vector<int> VarDomainSize; // Number of kept values (present in the union) for each variable (indexed by variable position in scope)
    unordered_map<string, int> mapDomainValToIndex; // Maps each value name to its index in UnionVarDomain (index >= NbValues means the value has been removed)
    vector<Cost> ReduceCostMatrix; // Reduced cost matrix after LAP solving: ReduceCostMatrix[var * NbNoAssignedVal + val] = reduced cost of assigning value val to variable var
    vector<uint8_t> visited; // visited[v] = 1 if variable node v has been finalized by Dijkstra
    vector<int> visitVal; // visitVal[val] = first variable matched to val which has been settled by Dijkstra (-1 if none)
    vector<uint8_t> inHeap; // inHeap[v] = 1 if variable node v is currently in the priority queue
    vector<Cost> distanceToVar; // Shortest distances from the Dijkstra source variable to each other variable node in the residual graph
    vector<Cost> distanceToVal; // Shortest distances from each value node to the value matched to the Dijkstra source variable in the residual graph (reverse graph)
    int Q = 0; // Number of Dijkstra runs to perform during reduced-cost filtering (derived from FiltLevel and NbNoAssigned)
    double FiltLevel; // Fraction of variables used as Dijkstra sources during filtering (0 = disabled, 1 = all variables)
    vector<int> trackList; // Variables not yet finalized by Dijkstra (dense array, used by the bimodal traversal strategy)
    vector<int> trackPos; // trackPos[var] = position of var in trackList, or -1 if finalized
    vector<PQ::handle_type> dijkstraHandles; // Handles of the priority queue, reused between two Dijkstra runs
    Cost MaxReducedCost = 0; // Maximum reduced cost observed over all arcs in the current residual graph (used as filtering threshold)
    Cost ReducedCost = 0; // Temporary reduced cost of a single arc computed during cost projection
    vector<int> varValIdx; // varValIdx[var * NbValues + val] = index of value val (in UnionVarDomain) inside the initial domain of variable var, or -1 if var cannot take it
    bool relativeInterior = false; // If true, in the square case the dual solution of the LAP is moved to the relative interior of the dual optimizers before projecting costs
    int findConflict = 0; // Flag returned by the LAP solver: non-zero if a conflict (infeasible sub-problem) was detected, value encodes the number of conflicting variables

    // Projects cost c to the problem lower bound (only if it is strictly positive)
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
            if (count == 0) {
                THROWCONTRADICTION;
            }
        }
    }

    // Greedy feasible assignment for a zero cost matrix: demands (lower bounds) first, then capacities (upper bounds).
    // varIdx / valIdx map local indices to scope / UnionVarDomain indices (nullptr = identity).
    // Feasibility (sum of demands <= nbVars <= sum of capacities) must already be checked.
    void fillFeasibleSupport(int nbVars, const int* varIdx, int nbVals, const int* valIdx,
        const vector<int>& cap, const vector<int>& dem)
    {
        int x = 0; // next variable to place
        auto place = [&](int j) {
            int v = varIdx ? varIdx[x] : x;
            storeLastAssignment[v] = toVarValue(v, valIdx ? valIdx[j] : j);
            ++x;
        };
        // Phase 1: satisfy demands
        if (sumlb > 0) {
            for (int j = 0; j < nbVals; ++j)
                for (int c = 0; c < dem[j] && x < nbVars; ++c)
                    place(j);
        }
        // Phase 2: fill remaining capacities
        for (int j = 0; j < nbVals && x < nbVars; ++j)
            for (int c = dem[j]; c < cap[j] && x < nbVars; ++c)
                place(j);
        assert(x == nbVars);
        storeAssignment = true;
    }

    // Make the stored value the EAC support of the variable
    void updateSupport(int varIndex)
    {
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

    // Stores the LAP solution (rowSol) when all variables are unassigned
    void storeFullAssignment()
    {
        for (int varIndex = 0; varIndex < arity_; varIndex++) {
            storeLastAssignment[varIndex] = toVarValue(varIndex, rowSol[varIndex]);
        }
        storeAssignment = true;
    }

    // If the problem is not square and every value can still host all the (nbVars) variables, the cardinality
    // constraint cannot be violated anymore: the constraint is disconnected.
    void deconnectIfAlwaysSatisfied(const vector<int>& cap, int nbVals, int nbVars)
    {
        if (isSquare)
            return;
        for (int valIndex = 0; valIndex < nbVals; valIndex++) {
            if (cap[valIndex] < nbVars)
                return;
        }
        if (ToulBar2::verbose > 0)
            cout << "GCC SATISFIED\n";
        deconnect();
    }

    // Called when the LAP solver proves the constraint infeasible:
    // increases the conflict weight of the conflicting variables, revises the constraint and throws a contradiction.
    // sol: solution array returned by the LAP solver (rowSol), nbConflict: findConflict,
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
     *        (used when the residual graph is not built).
     * @param NoAssignedVal Value index of each column of the reduced matrix.
     * @param ub            wcsp->getUb() - wcsp->getLb() before the LAP cost projection.
     * @param idx           idx(var, col) = position in costMatrix.
     */
    void projectReducedCostsPartial(const vector<int>& NoAssignedVal, Cost ub, const function<size_t(int, int)>& idx)
    {
        for (int varInd = 0; varInd < NbNoAssigned; ++varInd) {
            int varIndex = NoAssignedVar[varInd];
            for (int valInd = 0; valInd < NbNoAssignedVal; ++valInd) {
                if (costMatrix[idx(varInd, valInd)] < ub) {
                    Value value = toVarValue(varIndex, NoAssignedVal[valInd]);
                    if (value == WRONG_VAL)
                        continue;
                    ExtOrProJ(varIndex, value, (ReduceCostRow[varInd] + ReduceCostCol[valInd]));
                }
            }
        }
    }

    // True if the residual graph (VarList, ValList, ReduceCostMatrix) is built and used for filtering:
    // reduced-cost filtering must be enabled, and the problem must be square (all the value capacities are saturated)
    bool residualGraphEnabled() const
    {
        return FiltLevel > 0 && isSquare;
    }

    /**
     * @brief Projects the reduced costs onto the unary costs of the unassigned variables
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
            if (buildGraph)
                ValList[varInd] = vector<uint8_t>(NbNoAssignedVal, 0);
            for (int valInd = 0; valInd < NbNoAssignedVal; ++valInd) {
                size_t pos = idx(varInd, valInd);
                if (costMatrix[pos] < current_ub) {
                    int valIndex = NoAssignedVal[valInd];
                    Value value = toVarValue(varIndex, valIndex);
                    if (value == WRONG_VAL)
                        continue;
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
     * @brief Moves the optimal dual solution (ReduceCostRow = u per variable, ReduceCostCol = v per value) of the
     *        SQUARE capacitated assignment problem of the gcc into the relative interior of the set of dual optimizers
     *        (Algorithm 1 of Dlask & Savchynskyy, "Relative-interior solution for the (incomplete) linear assignment
     *        problem with applications to the quadratic assignment problem").
     *
     *        In the relative-interior solution, an edge has a zero reduced cost if and only if it belongs to
     *        at least one optimal assignment; every other edge gets a strictly positive reduced cost.
     *        The dual objective sum_i u_i + sum_j cap_j * v_j (hence the lower bound) and the optimal matching are unchanged.
     *
     *        Square case: every value j is used by exactly cap[j] variables. A value of capacity cap[j] stands for cap[j]
     *        interchangeable copies sharing the same dual v_j. The nodes of the graph are the variables; a variable k
     *        "owns" its matched column rowSol[k]. All the variables matched to the same column are mutually connected
     *        (the matched edge (k', j) is tight), so they always belong to the same strongly connected component.
     *        A component V therefore contains whole groups of variables sharing a column, and
     *        #variables(V) = sum of cap[j] over the distinct columns j owned by V.
     *        Shifting V by alpha (u_r -= alpha for r in V, v_j += alpha ONCE per distinct owned column j) keeps the dual
     *        objective unchanged: -alpha * #variables(V) + alpha * sum cap[j] = 0.
     *
     *        Steps:
     *        1. Build the directed graph G on the variables: arc k -> i if the edge (i, rowSol[k]) has a zero reduced cost.
     *        2. Compute the strongly connected components (Tarjan). Edges inside a component belong to an optimal
     *           assignment and must stay tight.
     *        3. Build the condensation and a topological order of it.
     *        4. For each component V having ingoing arcs in the condensation, taken in reverse topological order
     *           (sinks first), shift its duals as above. alpha is half of the smallest current reduced cost of the edges
     *           (variable outside V, column owned by V), which are the only ones decreased by the shift.
     *
     *        Costs are integers: alpha = floor(min/2), so an arc of reduced cost 1 cannot be made strictly positive.
     *        The result is still an optimal and feasible dual solution.
     *
     *        Convention (same as projectReducedCostsFull):
     *        reducedCost(row, col) = costMatrix[idx(row, col)] - (ReduceCostRow[row] + ReduceCostCol[col]) >= 0,
     *        zero on the matching. Call it after rowSol has been filled and only when isSquare.
     *
     * @param n   Number of variables (rows of the LAP).
     * @param m   Number of values (columns of the LAP).
     * @param cap cap[j] = (remaining) capacity of column j; the matching must use column j exactly cap[j] times.
     * @param ub  Cost at which an edge is considered absent (costMatrix >= ub: forbidden edge, no dual constraint).
     * @param idx idx(row, col) = position in costMatrix.
     */
    void makeDualRelativeInterior(int n, int m, const vector<int>& cap, Cost ub, const function<size_t(int, int)>& idx)
    {
        auto present = [&](int row, int col) { return costMatrix[idx(row, col)] < ub; };
        auto reducedCosts = [&](int row, int col) { return costMatrix[idx(row, col)] - (ReduceCostRow[row] + ReduceCostCol[col]); };

        // Variables matched to each column; in the square case every column is saturated
        vector<vector<int>> byCol(m);
        for (int k = 0; k < n; ++k)
            byCol[rowSol[k]].push_back(k);
        for (int j = 0; j < m; ++j)
            if ((int)byCol[j].size() != cap[j])
                return; // matching not saturating the capacities (should not happen): leave the duals unchanged

        // 1. Graph on the variables: arc k -> i if the edge (i, rowSol[k]) is tight (computed once per column)
        vector<vector<int>> g(n);
        for (int j = 0; j < m; ++j) {
            if (byCol[j].empty())
                continue;
            vector<int> tightRows;
            for (int i = 0; i < n; ++i)
                if (present(i, j) && reducedCosts(i, j) == MIN_COST)
                    tightRows.push_back(i);
            for (int k : byCol[j])
                for (int i : tightRows)
                    if (i != k)
                        g[k].push_back(i);
        }

        // 2. Strongly connected components
        vector<int> index(n, -1);
        vector<int> lowlink(n, 0);
        vector<int> comp(n, -1);
        vector<uint8_t> onStack(n, 0);
        vector<int> stk;
        int counter = 0;
        int nbComp = 0;
        for (int root = 0; root < n; ++root) {
            if (index[root] == -1) {
                strongconnect(root, g, index, lowlink, comp, onStack, stk, counter, nbComp);
            }
        }

        // 3. Build the list of variables belonging to each SCC and the condensation graph.
        vector<vector<int>> components(nbComp);
        for (int i = 0; i < n; ++i) {
            components[comp[i]].push_back(i);
        }

        vector<vector<int>> dag(nbComp);
        vector<int> indegree(nbComp, 0);
        for (int k = 0; k < n; ++k) {
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
        vector<int> colMark(m, -1); // colMark[j] == c if column j has already been collected for component c
        for (auto it = topo.rbegin(); it != topo.rend(); ++it) {
            int c = *it;
            if (!hasIn[c])
                continue;

            // Distinct columns owned by the component
            vector<int> cols;
            for (int i : components[c]) {
                int j = rowSol[i];
                if (colMark[j] != c) {
                    colMark[j] = c;
                    cols.push_back(j);
                }
            }

            // Smallest current reduced cost among the edges decreased by the shift:
            // (variable outside the component, column owned by the component).
            Cost minRC = MAX_COST;
            for (int j : cols) {
                for (int k = 0; k < n; ++k) {
                    if (comp[k] != c && present(k, j)) {
                        minRC = min(minRC, reducedCosts(k, j));
                    }
                }
            }

            // Maximum admissible shift.
            Cost alpha = (minRC >= MAX_COST) ? 1 : minRC / 2;
            if (alpha <= MIN_COST) {
                continue;
            }

            // Shift the dual variables of this SCC: every variable once, every owned column ONCE.
            for (int i : components[c])
                ReduceCostRow[i] -= alpha;
            for (int j : cols)
                ReduceCostCol[j] += alpha;
        }
    }

    /**
     * @brief GAC by strongly connected components (Régin) on the residual graph of the LAP.
     *        Step 1 of reducedCostFiltering(): structural criterion (alternating cycles) instead of reduced-cost shortest paths.
     *
     *        Nodes: variables 0..NbNoAssigned-1, then the columns (values).
     *        Arcs:  - matched:      var -> rowSol[var]               (several variables may be matched to the same value)
     *               - non-matching: val -> var for var in VarList[val]
     *
     *        The residual graph is only built when the problem is square: every value capacity is saturated, so
     *        there is no free capacity and a non-matching arc (x, v) is consistent iff x and v belong to the same SCC.
     *        Other arcs are pruned. Pruning is sound: only arcs of reduced cost < newCurrentUb are in the graph,
     *        and any solution within the budget uses only such arcs.
     *
     * @param NoAssignedVal Value index (in UnionVarDomain) of each column.
     * @param VarList       VarList[val] = variables having a non-matching arc to val (pruned arcs are erased).
     * @param ValList       ValList[var][val] = 1 if the arc is in the residual graph (pruned arcs are reset to 0).
     * @param newCurrentUb  Remaining cost budget (infinite cost projected on removed values).
     * @return true if at least one arc has been pruned.
     */
    bool gacFilteringSCC(const vector<int>& NoAssignedVal, vector<vector<int>>& VarList, vector<vector<uint8_t>>& ValList, Cost newCurrentUb)
    {
        if (!residualGraphEnabled())
            return false;

        const int nVar = NbNoAssigned; // variable nodes
        const int m = NbNoAssignedVal; // column nodes: nVar .. nVar + m - 1
        const int N = nVar + m;

        // --- Directed graph
        vector<vector<int>> g(N);
        for (int var = 0; var < nVar; ++var)
            g[var].push_back(nVar + rowSol[var]); // matched arc
        for (int val = 0; val < m; ++val)
            for (int var : VarList[val])
                g[nVar + val].push_back(var); // non-matching arcs

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
                int varIndex = NoAssignedVar[varInd];
                auto* variable = scope[varIndex];
                Value value = toVarValue(varIndex, valIndex);
                if (value != WRONG_VAL && variable->canbe(value)) {
                    if (ToulBar2::verbose > 0)
                        cout << "GAC REMOVE VALUE " << value << " from " << variable->getName() << endl;
                    ExtOrProJ(varIndex, value, -newCurrentUb); // project an infinite cost on this value
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
     * @brief Bimodal Dijkstra on the forward residual graph.
     *        Computes shortest distances from a source variable to all other
     *        variable nodes. At each settled node, chooses between iterating
     *        over the (sparse) adjacency list or the (small) unsettled set,
     *        whichever is cheaper.
     *        Several variables can be matched to the same value (capacity > 1): only the first settled one
     *        expands the value, the others inherit its distance.
     * @param source    Index of the source variable node.
     * @param VarList   VarList[val] = list of variable indices that have val as a non-matching arc.
     * @param ValList   ValList[var][val] = 1 if arc (var,val) exists in the residual graph.
     *
     *      Source : Bimodal Depth-First Search for Scalable GAC for AllDifferent.
     *      Sulian Le Bozec Chiffoleau; Nicolas Beldiceanu; Charles Prud'homme;
     *      Gilles Simonin; and Xavier Lorca, Proceedings of the Thirty-Fourth
     *      International Joint Conference on Artificial Intelligence, IJCAI 2025,
     *      Montreal, Canada, August 16-22, 2025.
     */
    void BimodalDijkstra(int source, vector<vector<int>>& VarList, vector<vector<uint8_t>>& ValList)
    {
        const int Nbvars = NbNoAssigned;

        if (dijkstraHandles.size() < (size_t)Nbvars)
            dijkstraHandles.resize(Nbvars);

        auto& handles = dijkstraHandles;
        PQ pq;

        distanceToVar.assign(Nbvars, MAX_COST);
        visited.assign(Nbvars, 0);
        visitVal.assign(NbNoAssignedVal, -1);

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
        // A node reached through a zero reduced-cost arc is settled immediately.
        auto expandValue = [&](Cost d, int val) {
            auto& adj = VarList[val];

            if (adj.size() < trackList.size()) {
                // Sparse mode
                for (int nextVar : adj) {
                    if (visited[nextVar])
                        continue;
                    if (relax(nextVar, d, val) && ReduceCostMatrix[nextVar * NbNoAssignedVal + val] == 0) {
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
                        if (relax(nextVar, d, val) && ReduceCostMatrix[nextVar * NbNoAssignedVal + val] == 0) {
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

            int val = rowSol[var]; // Matching arc: traverse to the matched value
            if (visitVal[val] == -1) {
                visitVal[val] = var;
            } else {
                // The value has already been expanded by another variable matched to it
                distanceToVar[var] = distanceToVar[visitVal[val]];
                continue;
            }
            expandValue(d, val);
        }
    }

    /**
     * @brief Bimodal-style Dijkstra on the reverse residual graph (square case).
     *        Computes, for every value node, the shortest distance from this value to the value matched to the
     *        source variable (the "landmark"), used by the Régin landmark upper-bound test.
     *        Reverse expansion of a settled value u: every variable y matched to u (there may be several since
     *        capacities can be > 1) gives access to every value w having a non-matching arc w -> y, at cost
     *        ReduceCostMatrix[y][w].
     * @param source  Index of the source variable (its matched value is the start node).
     * @param ValList ValList[var][val] = 1 if the non-matching arc (val -> var) exists in the residual graph.
     */
    void BimodalDijkstraReverseGraph(int source, vector<vector<uint8_t>>& ValList)
    {
        const int Nbvals = NbNoAssignedVal;

        if (dijkstraHandles.size() < (size_t)max(NbNoAssigned, Nbvals))
            dijkstraHandles.resize(max(NbNoAssigned, Nbvals));

        auto& handles = dijkstraHandles;
        PQ pq;

        distanceToVal.assign(Nbvals, MAX_COST);
        visited.assign(Nbvals, 0);
        inHeap.assign(Nbvals, 0);

        // Variables matched to each value (several variables per value when the capacity is > 1)
        vector<vector<int>> matchedVars(Nbvals);
        for (int var = 0; var < NbNoAssigned; ++var)
            matchedVars[rowSol[var]].push_back(var);

        // trackList = dense list of the values not yet settled, trackPos = position of each value in trackList (-1 if settled)
        trackList.resize(Nbvals);
        trackPos.resize(Nbvals);
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

        // Initialize source: the value matched to the source variable
        int sourceVal = rowSol[source];
        distanceToVal[sourceVal] = 0;
        handles[sourceVal] = pq.push({ 0, sourceVal });
        inHeap[sourceVal] = 1;

        while (!pq.empty()) {
            auto [d, val] = pq.top();
            pq.pop();

            visited[val] = 1;
            trackErase(val);

            for (int var : matchedVars[val]) {
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
                                continue; // trackList[k] has been replaced by another node
                            }
                        }
                    }
                    ++k;
                }
            }
        }
    }

    /**
     * @brief Régin landmark upper-bound test (square case).
     *        Runs Dijkstra in both directions from the source variable (forward on the variables, backward on the
     *        values), then checks whether the sum of the two maximum distances is within the remaining cost budget H.
     *        If so, no arc of the residual graph can have a reduced cost >= H and the filtering can stop early.
     *        With capacities > 1, the forward distance of a value is the minimum over the variables matched to it.
     *        The forward distances (distanceToVar) are left valid for the pruning that follows.
     * @param source  Source variable index used as the landmark.
     * @param VarList VarList[val] = variables reachable via value val (forward graph).
     * @param ValList ValList[var][val] = 1 if arc (var,val) exists in residual graph.
     * @param H       Remaining cost budget (wcsp upper bound minus current lower bound).
     * @return true if all values are consistent (no pruning possible), false otherwise.
     */
    bool ReginLandmarkUpperBound(int source, vector<vector<int>>& VarList, vector<vector<uint8_t>>& ValList, Cost H)
    {
        // Forward pass: distances from source to all variable nodes
        BimodalDijkstra(source, VarList, ValList);
        // Backward pass: distances from all value nodes to the value matched to source
        BimodalDijkstraReverseGraph(source, ValList);

        // Forward distance of each value = minimum over the variables matched to it
        vector<Cost> distanceFromVal(NbNoAssignedVal, MAX_COST);
        for (int var = 0; var < NbNoAssigned; ++var)
            distanceFromVal[rowSol[var]] = min(distanceFromVal[rowSol[var]], distanceToVar[var]);

        // An unreachable node (infinite distance) means we cannot conclude: avoid MAX_COST + MAX_COST overflow
        Cost maxDistanceFromVar = 0;
        for (Cost d : distanceFromVal) {
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

        // Pruning needs a reduced cost >= H: if max_forward + max_backward + MaxReducedCost < H, no arc can be pruned
        return maxDistanceFromVar + maxDistanceToVar < H - MaxReducedCost;
    }

    /**
     * @brief Filtering of variable domains with the Sellmann or Cambazard method (reduced-cost based filtering).
     *        Works with assigned variables (partial case) as well as with none (full case, see projectAndFilterFullProblem()).
     *        Step 1: GAC by strongly connected components (gacFilteringSCC), always run when the residual graph exists.
     *        Step 2: reduced-cost filtering with Dijkstra on the remaining arcs, only if some arc can still be pruned.
     *        Does nothing if filtering is disabled (see residualGraphEnabled()).
     *
     *        Source : Claus, G.; Cambazard, H.; and Jost, V. 2020. Analysis of Reduced Costs Filtering
     *        for Alldifferent and Minimum Weight Alldifferent Global Constraints. ECAI 2020, vol. 325, 323-330.
     *
     * @param NoAssignedVal Value index of each column of the reduced matrix.
     * @param VarList       VarList[val] = variables having a non-matching arc to val (updated when arcs are pruned).
     * @param ValList       ValList[var][val] = 1 if arc (var,val) is in the residual graph (updated when arcs are pruned).
     * @param newCurrentUb  Remaining cost budget.
     */
    void reducedCostFiltering(const vector<int>& NoAssignedVal, vector<vector<int>>& VarList, vector<vector<uint8_t>>& ValList, Cost newCurrentUb)
    {
        if (!residualGraphEnabled())
            return;

        // Step 1: GAC by strongly connected components
        gacFilteringSCC(NoAssignedVal, VarList, ValList, newCurrentUb);

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
        int stop = (Q / 2) + 1;
        bool AllValConsistent = false;
        vector<uint8_t> valUsedAsSource(NbNoAssignedVal, 0); // a value is used as Dijkstra source only once
        vector<uint8_t> valDone; // values already pruned for the current Dijkstra run

        // Prunes the arcs (varInd, valInd) of value valInd whose reduced cost exceeds the budget,
        // given the distances computed from the current source. 'var' is a node matched to valInd.
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
                    int varIndex = NoAssignedVar[varInd];
                    auto* variable = scope[varIndex];
                    Value value = toVarValue(varIndex, valIndex);
                    if (value != WRONG_VAL && variable->canbe(value)) {
                        if (ToulBar2::verbose > 0)
                            cout << "REMOVE VALUE " << value << " from " << variable->getName() << endl;
                        ExtOrProJ(varIndex, value, -newCurrentUb); // SdG: project infinite cost on this value and avoid to skip and reenter the GCC constraint without finishing the current filtering
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
            if (valUsedAsSource[rowSol[varInde]])
                continue;
            valUsedAsSource[rowSol[varInde]] = 1;
            numzeroremoved++;

            if (comp == 0) {
                // Apply Régin Upper Bounds of Shortest Paths algorithm with landmark
                AllValConsistent = ReginLandmarkUpperBound(varInde, VarList, ValList, newCurrentUb);
            } else {
                // Bimodal Dijkstra's shortest path algorithm
                BimodalDijkstra(varInde, VarList, ValList);
            }
            if (AllValConsistent)
                break;
            comp++;

            valDone.assign(NbNoAssignedVal, 0);
            for (int var = 0; var < NbNoAssigned; ++var) {
                if (distanceToVar[var] >= MAX_COST)
                    continue;
                int valInd = rowSol[var];
                if (valDone[valInd])
                    continue; // all the variables matched to a value share the same distance
                valDone[valInd] = 1;
                if (pruneArcsOfValue(valInd, var))
                    numzeroremoved = 0;
            }
        }
    }

    /**
     * @brief Full case (no variable considered as assigned): projects the reduced costs onto the unary costs,
     *        builds the residual graph over ALL the variables and ALL the values, and runs reducedCostFiltering().
     *        The whole problem is viewed as a "reduced" problem: NbNoAssigned = arity_, NbNoAssignedVal = NbValues,
     *        NoAssignedVar and NoAssignedVal are the identity.
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
        vector<vector<uint8_t>> ValList(NbNoAssigned);
        projectAndBuildResidualGraph(NoAssignedVal, current_ub, newCurrentUb, idx, VarList, ValList);
        reducedCostFiltering(NoAssignedVal, VarList, ValList, newCurrentUb);
    }

    // Non-square case, partial assignment: updates the extended unary costs (deltaCosts) from the LAP duals
    // and adjusts TotalCost (the lower bound that can be projected). Throws a contradiction if it exceeds current_ub.
    void updateDeltaCostsPartial(const vector<int>& NoAssignedVal, Cost& TotalCost, Cost current_ub)
    {
        if (sumlb == 0) {
            Cost mindelta = ReduceCostRow[NbNoAssigned];
            for (int valInd = 0; valInd < NbNoAssignedVal; ++valInd) {
                int valIndex = NoAssignedVal[valInd];
                if (capacity[valIndex] > 0)
                    deltaCosts[valIndex] -= (ReduceCostCol[valInd] + mindelta);
            }
            if (mindelta > 0) {
                TotalCost += mindelta * (sumub - arity_ - 1);
                if (TotalCost >= current_ub) {
                    wcsp->revise(this);
                    THROWCONTRADICTION;
                }
            }
        } else {
            Cost assignedPart = MIN_COST;
            for (int valInd = 0; valInd < NbNoAssignedVal; ++valInd) {
                int valIndex = NoAssignedVal[valInd];
                if (capacity[valIndex] > 0) {
                    deltaCosts[valIndex] += ReduceCostCol[valInd];
                    assignedPart += ReduceCostCol[valInd] * (capacity[valIndex] - CurrentCapacity[valIndex]);
                }
            }
            sumcostsprojectonvar -= assignedPart;
        }
        for (int varInd = 0; varInd < NbNoAssigned; ++varInd)
            sumcostsprojectonvar += ReduceCostRow[varInd];
    }

    // Non-square case, all variables: updates the extended unary costs (deltaCosts) from the LAP duals
    // and adjusts TotalCost (the lower bound that can be projected). Throws a contradiction if it exceeds current_ub.
    void updateDeltaCostsFull(Cost& TotalCost, Cost current_ub)
    {
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

    /**
     * @brief Some variables are assigned (same domains):
     *        reduces the problem to the unassigned variables and the values with remaining capacity,
     *        solves the LAP, projects the costs and filters the domains.
     */
    void propagatePartialAssignment()
    {
        // Collect the values with a remaining capacity, with their remaining capacity and demand
        vector<int> NoAssignedVal;
        vector<int> NoAsscapacity;
        vector<int> NoAssdemand;
        int sumlbtemp = 0;
        int sumubtemp = 0;
        for (int valIndex = 0; valIndex < NbValues; ++valIndex) {
            if (CurrentCapacity[valIndex] > 0) {
                NoAssignedVal.push_back(valIndex);
                NoAsscapacity.push_back(CurrentCapacity[valIndex]);
                sumubtemp += CurrentCapacity[valIndex];
                int lbcap = 0;
                if (sumlb > 0) {
                    lbcap = max(0, demand[valIndex] - (capacity[valIndex] - CurrentCapacity[valIndex]));
                    sumlbtemp += lbcap;
                }
                NoAssdemand.push_back(lbcap);
            }
        }

        NbNoAssignedVal = NoAssignedVal.size();
        if (NbNoAssignedVal == 0 || sumlbtemp > NbNoAssigned || sumubtemp < NbNoAssigned)
            THROWCONTRADICTION;

        // Initialize cost matrix (rows = unassigned variables, columns = values with remaining capacity) for the Jonker algorithm
        Cost current_ub = wcsp->getUb() - wcsp->getLb();
        auto idx = [this](int var, int col) { return (size_t)var * NbNoAssignedVal + col; };

        for (int varInd = 0; varInd < NbNoAssigned; ++varInd) {
            int varIndex = NoAssignedVar[varInd];
            auto* variable = scope[varIndex];
            for (int valInd = 0; valInd < NbNoAssignedVal; ++valInd) {
                Value val = toVarValue(varIndex, NoAssignedVal[valInd]);
                if (val != WRONG_VAL && variable->canbe(val)) {
                    auto valCost = variable->getCost(val);
                    costMatrix[idx(varInd, valInd)] = valCost < current_ub ? valCost : current_ub;
                } else {
                    costMatrix[idx(varInd, valInd)] = current_ub;
                }
            }
        }

        if (any_of(costMatrix.begin(), costMatrix.begin() + NbNoAssigned * NbNoAssignedVal, [](auto x) { return x > MIN_COST; })) {
            // Solve assignment problem using Jonker algorithm
            Cost TotalCost;
            if (!isSquare) {
                if (sumlb == 0) {
                    // dummy row with the extended unary costs
                    for (int valInd = 0; valInd < NbNoAssignedVal; ++valInd)
                        costMatrix[idx(NbNoAssigned, valInd)] = deltaCosts[NoAssignedVal[valInd]];
                    TotalCost = lapjv_ub(NbNoAssigned + 1, NbNoAssignedVal, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub, NoAsscapacity, findConflict);
                } else {
                    TotalCost = lapjv_gcc(NbNoAssigned, NbNoAssignedVal, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub, NoAsscapacity, NoAssdemand, findConflict);
                }
            } else {
                TotalCost = lapjv_ub(NbNoAssigned, NbNoAssignedVal, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub, NoAsscapacity, findConflict);
            }

            if (TotalCost >= current_ub) {
                onLAPFailure(rowSol, findConflict, NbNoAssigned, &NoAssignedVar);
            } else if (TotalCost >= 0) {
                storePartialAssignment(NoAssignedVal);

                if (TotalCost > 0 || (any_of(ReduceCostRow, ReduceCostRow + NbNoAssigned, [](auto x) { return x > MIN_COST; }))) {
                    if (!isSquare)
                        updateDeltaCostsPartial(NoAssignedVal, TotalCost, current_ub);

                    // Square case: move the duals to the relative interior of the dual optimizers
                    if (relativeInterior && isSquare)
                        makeDualRelativeInterior(NbNoAssigned, NbNoAssignedVal, NoAsscapacity, current_ub, idx);

                    // Update the lower bound with the Jonker algorithm's total cost
                    projectLB(TotalCost);
                    Cost newCurrentUb = current_ub - TotalCost;

                    if (residualGraphEnabled()) {
                        // Adjust unary costs for unassigned variables based on reduced costs and build the residual graph
                        vector<vector<int>> VarList(NbNoAssignedVal);
                        vector<vector<uint8_t>> ValList(NbNoAssigned);
                        projectAndBuildResidualGraph(NoAssignedVal, current_ub, newCurrentUb, idx, VarList, ValList);

                        // GAC and filtering of variables domains with Sellmann or Cambazard method
                        reducedCostFiltering(NoAssignedVal, VarList, ValList, newCurrentUb);
                    } else {
                        projectReducedCostsPartial(NoAssignedVal, current_ub, idx);
                    }
                }

                // Update support values if needed for unassigned variables
                updateSupportsOfUnassigned();
            }
        } else {
            // Zero cost matrix: build a feasible assignment quickly
            fillFeasibleSupport(NbNoAssigned, NoAssignedVar.data(), NbNoAssignedVal, NoAssignedVal.data(), NoAsscapacity, NoAssdemand);
            // Assigned variables keep their value
            for (int k = 0; k < NbAssigned; ++k)
                storeLastAssignment[AssignedVar[k]] = toVarValue(AssignedVar[k], AssignedVal[k]);
            updateSupportsOfUnassigned();
            deconnectIfAlwaysSatisfied(NoAsscapacity, NbNoAssignedVal, NbNoAssigned);
        }
    }

    /**
     * @brief All variables unassigned (or different domains):
     *        solves the LAP on the whole problem, projects the costs, filters the domains and updates the supports.
     */
    void propagateFullAssignment()
    {
        Cost current_ub = wcsp->getUb() - wcsp->getLb();
        auto idx = [this](int var, int col) { return (size_t)var * NbValues + col; };

        // Initialize the cost matrix for all variables and their domain values
        fillFullCostMatrix(current_ub, idx);

        if (any_of(costMatrix.begin(), costMatrix.begin() + arity_ * NbValues, [](auto x) { return x > MIN_COST; })) {
            // Solve the Linear Assignment Problem (LAP) using the Jonker algorithm
            Cost TotalCost;
            if (!isSquare) {
                if (sumlb == 0) {
                    // dummy row with the extended unary costs
                    for (int valInd = 0; valInd < NbValues; ++valInd)
                        costMatrix[idx(arity_, valInd)] = deltaCosts[valInd];
                    TotalCost = lapjv_ub(arity_ + 1, NbValues, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub, capacity, findConflict);
                } else {
                    TotalCost = lapjv_gcc(arity_, NbValues, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub, capacity, demand, findConflict);
                }
            } else {
                TotalCost = lapjv_ub(arity_, NbValues, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub, capacity, findConflict);
            }

            if (TotalCost >= current_ub) {
                onLAPFailure(rowSol, findConflict, arity_, nullptr);
            } else if (TotalCost >= 0) {
                storeFullAssignment();

                if (TotalCost > 0 || (any_of(ReduceCostRow, ReduceCostRow + arity_, [](auto x) { return x > MIN_COST; }))) {
                    if (!isSquare)
                        updateDeltaCostsFull(TotalCost, current_ub);

                    // Square case: move the duals to the relative interior of the dual optimizers
                    if (relativeInterior && isSquare)
                        makeDualRelativeInterior(arity_, NbValues, capacity, current_ub, idx);

                    // Update lower bound with the total cost found
                    projectLB(TotalCost);

                    if (!isSquare) {
                        for (int varIndex = 0; varIndex < arity_; ++varIndex)
                            sumcostsprojectonvar += ReduceCostRow[varIndex];
                    }

                    // Adjust unary costs using reduced row and column costs
                    if (residualGraphEnabled())
                        projectAndFilterFullProblem(current_ub, current_ub - TotalCost, idx); // + GAC and reduced-cost filtering
                    else
                        projectReducedCostsFull(current_ub, idx);
                }

                // Update support for all variables (also when the optimal LAP cost is zero)
                updateAllSupports();
            }
        } else {
            // Zero cost matrix: build a feasible assignment quickly
            fillFeasibleSupport(arity_, nullptr, NbValues, nullptr, capacity, demand);
            updateAllSupports();
            deconnectIfAlwaysSatisfied(capacity, NbValues, arity_);
        }
    }

public:
    GlobalCardinalityConstraint(const GlobalCardinalityConstraint&) = delete; // owns raw arrays
    GlobalCardinalityConstraint& operator=(const GlobalCardinalityConstraint&) = delete;

    // Constructor for the GlobalCardinalityConstraint class
    GlobalCardinalityConstraint(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in)
        : AbstractNaryConstraint(wcsp, scope_in, arity_in)
        , Original_ub(wcsp->getUb())
        , lb(0)
        , assigneddeltas(0)
        , sumcostsprojectonvar(0)
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

            // Percentage of variables used as Dijkstra sources in reduced-cost filtering (0 = disabled)
            int rate = ToulBar2::ReducedCostsFiltering;
            if (rate > 0)
                FiltLevel = rate / 100.0;

            // Relative-interior dual solution (applied only in the square case, see makeDualRelativeInterior)
            if (ToulBar2::RelativeInteriorSolution > 0)
                relativeInterior = true;

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
            rowSol = new int[arity_ + 1];
            ReduceCostRow = new Cost[arity_ + 1];
            ReduceCostCol = new Cost[NbValues];

        } else {
            deconnect();
        }
    }

    // Destructor: frees the arrays used by the LAP solver
    virtual ~GlobalCardinalityConstraint()
    {
        delete[] rowSol;
        delete[] ReduceCostRow;
        delete[] ReduceCostCol;
    }

    // Reads the parameters of the constraint from a WCSP file: bounds of each value, removed values and (non-square case) extended unary costs
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
        Cost Max_cost = wcsp->getUb();
        unordered_map<Value, Value> mapValuetempValue;

        vector<Value> ListValTemp;
        ListValTemp.reserve(Gcc_NbValues);
        for (int v = 0; v < Gcc_NbValues; v++) {
            Value valueTemp;
            file >> valueTemp;
            file >> lower;
            sumlb += lower;
            file >> upper;
            if (lower > upper)
                THROWCONTRADICTION;
            sumub += upper;
            boundsTemp[valueTemp] = { lower, upper };
            ListValTemp.push_back(valueTemp);
        }

        // Values removed from the domains
        int Nb_removeValues = 0;
        file >> Nb_removeValues;
        vector<Value> ListremoveVal;
        ListremoveVal.reserve(Nb_removeValues);
        for (int val = 0; val < Nb_removeValues; val++) {
            Value rmval;
            file >> rmval;
            ListremoveVal.push_back(rmval);
        }

        // Shift each value by the number of removed values smaller than it
        vector<Value> ListVal;
        ListVal.reserve(Gcc_NbValues);
        for (int v = 0; v < Gcc_NbValues; v++) {
            int distance = 0;
            for (int r = 0; r < Nb_removeValues; r++) {
                if (ListremoveVal[r] < ListValTemp[v])
                    distance++;
            }
            ListVal.push_back(ListValTemp[v] - distance);
        }

        // Set capacity and demand of each listed value
        for (int v = 0; v < Gcc_NbValues; v++) {
            Value valueTemp = ListValTemp[v];
            Value value = ListVal[v];
            bounds[value] = boundsTemp[valueTemp];
            mapValuetempValue[valueTemp] = value;
            lower = bounds[value].first;
            upper = bounds[value].second;

            bool findval = false;
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

            if (!findval) {
                if (lower > 0) {
                    THROWCONTRADICTION;
                } else {
                    sumub -= upper;
                }
            }
        }
        if (sumlb > arity_) {
            THROWCONTRADICTION;
        }
        if (sumub < arity_ && Gcc_NbValues == NbValues) {
            THROWCONTRADICTION;
        }

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
                sumub += capacity[valIndex];
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
            // Move the values of null capacity at the end of the union: they are then ignored
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

        for (int valIndex = 0; valIndex < NbValues; ++valIndex) {
            CurrentCapacity[valIndex] = capacity[valIndex];
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
            if (nbDelta > 0) {
                for (int v = 0; v < nbDelta; v++) {
                    Value value;
                    Cost delta;
                    file >> value;
                    file >> delta;
                    for (int varIndex = 0; varIndex < arity_; ++varIndex) {
                        auto* variable = scope[varIndex];
                        if (variable->canbe(mapValuetempValue[value])) {
                            int idx = unionIndex(variable->getValueName(variable->toIndex(mapValuetempValue[value])));
                            if (idx < 0)
                                continue;
                            deltaCosts[idx] = delta;
                            break;
                        }
                    }
                }
                Cost SumCostsProjectVar = 0, Lb = 0;
                file >> Lb;
                file >> SumCostsProjectVar;
                lb += Lb;
                sumcostsprojectonvar += SumCostsProjectVar;
            }
        }
    }

    bool extension() const FINAL { return false; } // TODO: allows functional variable elimination but no other preprocessing
    bool isGCC() const FINAL { return true; } // this constraint is a gcc
    bool isGCCSquare() const FINAL { return isSquare; } // true if the sum of the upper capacities == #variables

    // Returns the lower and upper capacities of every value
    const map<Value, pair<int, int>>& getBounds() const { return bounds; }

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
                    auto* variable = scope[i];
                    if (!variable->unassigned()) {
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

    // Returns true if the binary constraint ctr is implied by this gcc
    // (i.e. ctr has zero cost on every pair of different values, or when the value has an upper capacity >= 2)
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
                        if ((x->getValueName(x->toIndex(*iterx)) != y->getValueName(y->toIndex(*itery)) || bounds.find(*iterx) == bounds.end() || bounds[*iterx].second >= 2) && bctr->getCost(*iterx, *itery) > MIN_COST) {
                            return false;
                        }
                    }
                }
            } else {
                // Values are compared directly
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

    // Projects the gcc (forbidden equal values for values with upper capacity 1) onto the binary constraint ctr included in its scope
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
                    if (bounds.find(*iterx) != bounds.end() && bounds[*iterx].second == 1 && y->canbe(yval) && !CUT(bctr->getCost(*iterx, yval), wcsp->getUb())) {
                        bctr->addcost(*iterx, yval, mult_ub - bctr->getCost(*iterx, yval));
                    }
                }
            } else {
                // Values are matched directly
                for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
                    if (bounds.find(*iterx) != bounds.end() && bounds[*iterx].second == 1 && y->canbe(*iterx) && !CUT(bctr->getCost(*iterx, *iterx), wcsp->getUb())) {
                        bctr->addcost(*iterx, *iterx, mult_ub - bctr->getCost(*iterx, *iterx));
                    }
                }
            }
        }
    }

    // Returns the cost of the tuple s (0 if it satisfies the cardinalities, else the upper bound or a relaxed cost)
    Cost eval(const Tuple& s) override
    {
        // returns the cost of the corresponding assignment s
        if (isSquare) {
            Cost res = MIN_COST;
            Cost nbsame = 0; // number of violations (forbidden values or capacity exceeded)
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

        // Non-square case: the cost also includes the extended unary costs
        Cost res = 0;
        Cost nbsame = 0; // number of violations
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
                res = -lb + sumcostsprojectonvar;
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
        Cost sumdelta = abs(-lb + sumcostsprojectonvar);
        Cost m = *max_element(deltaCosts.begin(), deltaCosts.end());
        if (m > MIN_COST)
            sumdelta += m;
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
        // deltaCosts is indexed by the position of the value name in UnionVarDomain
        int idx = unionIndex(scope[index]->getValueName(scope[index]->toIndex(val)));
        if (idx < 0)
            return MIN_COST;
        return deltaCosts[idx];
    }

    double computeTightness() override { return 0; } // TODO: compute factorial(n)/(n**n)?

    // TODO: needed for dominance test by DEE
    // pair<pair<Cost, Cost>, pair<Cost, Cost>> getMaxCost(int index, Value a, Value b)

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
                auto* variable = scope[varIndex];
                int valIndex = unionIndex(variable->getValueName(variable->toIndex(value)));
                if(CurrentCapacity[valIndex] > 0){
                    CurrentCapacity[valIndex] = CurrentCapacity[valIndex] - 1;
                    if(CurrentCapacity[valIndex] == 0){
                        for (int i = 0; i < arity_; i++) {
                            if (i != varIndex && scope[i]->unassigned() && scope[i]->canbe(value)) {
                                scope[i]->remove(value);
                            }
                        }
                    }   
                } else {
                    THROWCONTRADICTION;
                }
                propagate();
                if (ToulBar2::FullEAC)
                    reviseEACGreedySolution();
            }
        }
    }

    /**
     * @brief Filters and removes saturated values from unassigned variables in the GCC constraint.
     *
     * This function enforces the GCC constraint by:
     * - Verifying that no value is used more often than its upper capacity.
     * - Removing values from unassigned variable domains once their capacity is saturated.
     * - Propagating newly assigned variables and updating internal state.
     * - Triggering full n-ary projection if thresholds on domain sizes or variable count are met.
     *
     * @return true if constraint propagation was successful without needing full projection.
     * @return false if the constraint switched to full n-ary projection.
     *
     * @throws Throws a contradiction exception (THROWCONTRADICTION) if a value capacity is exceeded.
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
                    return false; // Variable disconnected externally: abort
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
     * @brief Propagate the Global Cardinality constraint to enforce consistency and improve problem lower bound.
     *
     * This method performs constraint propagation by checking assigned variables,
     * removing inconsistent values, computing lower bounds using the Jonker algorithm,
     * and updating unary costs and supports accordingly.
     *
     * It handles different cases depending on whether variables are assigned or not,
     * and if the domains are the same. The work is delegated to propagatePartialAssignment()
     * or propagateFullAssignment().
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

        // Attempt to assign variables already connected and assigned (assign() propagates by itself)
        if (assignFirstAssignedVariable() || !connected())
            return;

        // Skip propagation if the stored optimal assignment is still valid
        if (canSkipPropagation())
            return;

        // Initialize number of unassigned variables
        NbNoAssigned = getNonAssigned();

        // If there are assigned variables (and same domains), handle accordingly
        if (NbNoAssigned < arity_) {
            for (int valIndex = 0; valIndex < NbValues; ++valIndex) {
                CurrentCapacity[valIndex] = capacity[valIndex];
            }
            if (SameDomain && RemoveAssignVar()) {
                propagatePartialAssignment();
            }
        }

        // Case when all variables are unassigned or domains differ
        if (connected() && (NbNoAssigned == arity_ || !SameDomain)) {
            propagateFullAssignment();
        }
        for (int valIndex = 0; valIndex < NbValues; ++valIndex) {
            CurrentCapacity[valIndex] = capacity[valIndex];
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

    // Prints the constraint (scope, extended unary costs, tightness, weights, bounds, arity) for debugging
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
               << " cost: " << -lb << " + " << sumcostsprojectonvar << " + (";
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
        os << " " << boundsTemp.size() << ":[";
        for (const auto& [key, bound] : boundsTemp) {
            os << " (" << key;
            os << "," << bound.first;
            os << "," << bound.second;
            os << " )";
        }
        os << "] arity: " << arity_;
        os << " unassigned: " << getNonAssigned() << "/" << unassigned_ << endl;
    }

    // returns the non-zero extended unary costs indexed by value (a value is reported through the first variable that can still take it)
    map<Value, Cost> nonZeroDeltaCosts()
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

    // Dumps the constraint in WCSP format (original or current problem)
    void dump(ostream& os, bool original = true) override
    {
        // a bound is kept in the current problem only if its value still has a capacity
        auto keepBound = [&](Value key) {
            string name = scope[0]->getValueName(scope[0]->toIndex(key));
            int idx = unionIndex(name);
            return idx >= 0 && capacity[idx] > 0;
        };

        if (original) {
            // original problem: bounds (+ delta costs if any); removed values are omitted
            os << arity_;
            for (int i = 0; i < arity_; i++)
                os << " " << scope[i]->wcspIndex;
            os << " -1 gcc ";

            os << boundsTemp.size();
            for (const auto& [key, bound] : boundsTemp) {
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
            for (const auto& [key, bound] : bounds)
                if (keepBound(key))
                    compt++;

            os << compt;
            for (const auto& [key, bound] : bounds) {
                if (!keepBound(key))
                    continue;
                os << " " << key;
                os << " " << bound.first;
                os << " " << bound.second;
            }
        }

        map<Value, Cost> mapValuesdeltaCosts;
        if (!isSquare)
            mapValuesdeltaCosts = nonZeroDeltaCosts();

        // removed values: written only if non-empty, or if delta costs follow
        // (the count is then needed, even 0, so that the reader stays in sync)
        if (!original) {
            vector<Value> removedVals = removedValues();
            if (!removedVals.empty() || !mapValuesdeltaCosts.empty()) {
                os << " " << removedVals.size();
                for (Value v : removedVals)
                    os << " " << v;
            }
        }

        // delta costs, lb and sumcostsprojectonvar: only if there are non-zero delta costs
        if (!mapValuesdeltaCosts.empty()) {
            os << " " << mapValuesdeltaCosts.size();
            for (auto& [key, delta] : mapValuesdeltaCosts) {
                os << " " << key;
                os << " " << delta;
            }
            os << " " << lb;
            os << " " << sumcostsprojectonvar;
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
        }
        os << "],\"type\":\"gcc\",\"params\":{";

        // Separator shared by all optional parameters: commas are always valid JSON.
        bool firstParam = true;
        auto sep = [&]() {
            if (!firstParam)
                os << ",";
            firstParam = false;
        };

        // bounds (always written, even if empty)
        ostringstream bs;
        bool anyBound = false;
        auto addBound = [&](Value key, const auto& bound) {
            if (anyBound)
                bs << ",";
            bs << "[" << key << "," << bound.first << "," << bound.second << "]";
            anyBound = true;
        };
        if (original) {
            for (const auto& [key, bound] : boundsTemp)
                addBound(key, bound);
        } else {
            for (const auto& [key, bound] : bounds) {
                string name = scope[0]->getValueName(scope[0]->toIndex(key));
                int idx = unionIndex(name);
                if (idx < 0 || capacity[idx] <= 0)
                    continue;
                addBound(key, bound);
            }
        }
        sep();
        os << "\"bounds\":[" << bs.str() << "]";

        // removed values (current problem only, written only if non-empty)
        if (!original) {
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
        }

        // delta costs (non-square problems), then lb / sumcostsprojectonvar if non-zero
        map<Value, Cost> mapValuesdeltaCosts;
        if (!isSquare)
            mapValuesdeltaCosts = nonZeroDeltaCosts();
        if (!mapValuesdeltaCosts.empty()) {
            sep();
            os << "\"deltacosts\":[";
            printed = false;
            for (auto& [key, delta] : mapValuesdeltaCosts) {
                if (printed)
                    os << ",";
                os << "[" << key << "," << delta << "]";
                printed = true;
            }
            os << "]";
            if (lb != 0) {
                sep();
                os << "\"lb\":" << lb;
            }
            if (sumcostsprojectonvar != 0) {
                sep();
                os << "\"sumcostsprojectonvar\":" << sumcostsprojectonvar;
            }
        }

        os << "}},\n";
    }
};

#endif /*TB2GCC_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
