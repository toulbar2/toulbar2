/** \file tb2lapjv.hpp
 *  \brief Solves linear assignment problems using Jonker and Volgenant's algorithm.
 *
 */

#ifndef TB2LAPJV_HPP_
#define TB2LAPJV_HPP_

#include <bits/stdc++.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <iostream>
#include <cstdlib> 
#include "core/tb2enumvar.hpp"
#include <boost/heap/pairing_heap.hpp>

using namespace std;
using PQ = boost::heap::pairing_heap<pair<Cost, int>, boost::heap::compare<greater<>>>;

/**
 * @brief Implements the Jonker Algorithm (LAPJV) to solve the assignment problem.
 *
 * The Jonker Algorithm provides an optimal solution to the linear assignment problem,
 * which consists of assigning tasks to agents in such a way that the total cost is minimized.
 * through a well-defined series of steps.
 * Original source code from https://github.com/scipy/scipy/tree/main/scipy/optimize/rectangular_lsap
 * 
 * Paper : 
 * 
 * 
 */


 static intptr_t
augmenting_path(intptr_t dim_val, const vector<Cost>& cost, vector<Cost>& u,
                vector<Cost>& v, vector<intptr_t>& path,
                vector<intptr_t>& row4col,
                vector<Cost>& shortestPathCosts, intptr_t i,
                vector<bool>& SR, vector<bool>& SC,
                vector<intptr_t>& remaining, Cost* p_minVal,Cost MAX_COST)
{
    Cost minVal = 0;

    // Crouse's pseudocode uses set complements to keep track of remaining
    // nodes.  Here we use a vector, as it is more efficient in C++.
    intptr_t num_remaining = dim_val;
    for (intptr_t it = 0; it < dim_val; it++) {
        // Filling this up in reverse order ensures that the solution of a
        // constant cost matrix is the identity matrix .
        remaining[it] = dim_val - it - 1;
    }

    fill(SR.begin(), SR.end(), false);
    fill(SC.begin(), SC.end(), false);
    fill(shortestPathCosts.begin(), shortestPathCosts.end(), MAX_COST);

    // find shortest augmenting path
    intptr_t sink = -1;
    while (sink == -1) {

        intptr_t index = -1;
        Cost lowest = MAX_COST;
        SR[i] = true;

        for (intptr_t it = 0; it < num_remaining; it++) {
            intptr_t j = remaining[it];

            Cost r = minVal + cost[i * dim_val + j] - u[i] - v[j];
            if (r < shortestPathCosts[j]) {
                path[j] = i;
                shortestPathCosts[j] = r;
            }

            // When multiple nodes have the minimum cost, we select one which
            // gives us a new sink node. This is particularly important for
            // integer cost matrices with small co-efficients.
            if (shortestPathCosts[j] < lowest ||
                (shortestPathCosts[j] == lowest && row4col[j] == -1)) {
                lowest = shortestPathCosts[j];
                index = it;
            }
        }

        minVal = lowest;
        if (minVal == MAX_COST) { // MAX_COSTeasible cost matrix
            return -1;
        }

        intptr_t j = remaining[index];
        if (row4col[j] == -1) {
            sink = j;
        } else {
            i = row4col[j];
        }

        SC[j] = true;
        remaining[index] = remaining[--num_remaining];
    }

    *p_minVal = minVal;
    return sink;
}

static Cost lapjv(intptr_t dim_var, intptr_t dim_val, const vector<Cost>& cost,  int* b,  Cost* usol, Cost* vsol, Cost MAX_COST, int& findConflict)
{

    // initialize variables
    findConflict = 0;
    vector<Cost> u(dim_var, 0);
    vector<Cost> v(dim_val, 0);
    vector<Cost> shortestPathCosts(dim_val);
    vector<intptr_t> path(dim_val, -1);
    vector<intptr_t> col4row(dim_var, -1);
    vector<intptr_t> row4col(dim_val, -1);
    vector<bool> SR(dim_var);
    vector<bool> SC(dim_val);
    vector<intptr_t> remaining(dim_val);

    // iteratively build the solution
    for (intptr_t curRow = 0; curRow < dim_var; curRow++) {

        Cost minVal;
        intptr_t sink = augmenting_path(dim_val, cost, u, v, path, row4col,
                                        shortestPathCosts, curRow, SR, SC,
                                        remaining, &minVal, MAX_COST);
                                        
        if (sink < 0) {
            findConflict = 1;
            b[0] = curRow;
            for(int col=0; col< dim_val; col++){
                if(cost[curRow * dim_val + col] < MAX_COST){
                    b[findConflict] = row4col[col];
                    findConflict++;
                }
            }     
            return MAX_COST;
        }

        // update dual variables
        u[curRow] += minVal;
        for (intptr_t i = 0; i < dim_var; i++) {
            if (SR[i] && i != curRow) {
                u[i] += minVal - shortestPathCosts[col4row[i]];
            }
        }

        for (intptr_t j = 0; j < dim_val; j++) {
            if (SC[j]) {
                v[j] -= minVal - shortestPathCosts[j];
            }
        }

        // augment previous solution
        intptr_t j = sink;
        while (1) {
            intptr_t i = path[j];
            row4col[j] = i;
            swap(col4row[i], j);
            if (i == curRow) {
                break;
            }
        }
    }

    Cost Total_cost = 0;
    for (intptr_t i = 0; i < dim_var; i++) {

            b[i] = col4row[i];
            Total_cost+= cost[i * dim_val + col4row[i]];
            usol[i] = u[i];
    }

    for (intptr_t j = 0; j < dim_val; j++) {

        vsol[j] = v[j];
    }
    return Total_cost;
}




static intptr_t
augmenting_path_gcc(intptr_t dim_val,
const vector<Cost>& cost,
vector<Cost>& u,
vector<Cost>& v,
intptr_t i,
vector<intptr_t>& path,
vector<Cost>& shortestPathCost,
vector<bool>& SC,
vector<bool>& SR,
const vector<int>& col4row,
const vector<vector<int>>& row4col,
                     Cost MAX_COST,
const vector<int>& capacity)
{
    /* --- Initialisation --------------------------------------------------- */

    // Clear visited sets and reset path costs for this augmentation round.
    fill(SR.begin(), SR.end(), false);
    fill(SC.begin(), SC.end(), false);
    fill(shortestPathCost.begin(), shortestPathCost.end(), MAX_COST);

    // Mark the source row as visited so it is never relaxed through again.
    SR[i] = true;

    // Seed shortest-path costs with the reduced cost from source row i to
    // every column:  rc(i,j) = cost(i,j) - u[i] - v[j].
    for (intptr_t j = 0; j < dim_val; ++j) {
        if(capacity[j] > 0){
            path[j] = -1;
            Cost r = cost[i * dim_val + j] - u[i] - v[j];
            if (r < shortestPathCost[j]) {
                shortestPathCost[j] = r;
                path[j] = i;
            }
        }
    }

    /* --- Main label-setting (Dijkstra) loop ------------------------------- */
    intptr_t sink = -1;
    while (true) {

        // Extract-min: pick the unsettled column with the smallest
        Cost lowest = MAX_COST;
        intptr_t index = -1;
        for (intptr_t j = 0; j < dim_val; ++j) {
            if (!SC[j] && shortestPathCost[j] < lowest){
                lowest = shortestPathCost[j];
                index = j;
            }
        }

        // No reachable unsettled column remains → problem is infeasible
        // from this source row.
        if (index == -1) {
            return -1;
        }

        // Settle column `index`.
        SC[index] = true;

        // If column `index` still has free capacity it can absorb one more
        // assignment, so it is a valid augmenting-path sink.
        if ((int)row4col[index].size() < capacity[index]) {
            sink = index;
            break;
        }

        for (int r : row4col[index]) {
            if (SR[r]) continue;   // row r already settled, skip
            SR[r] = true;

            // Reduced cost of the backward arc (undoing assignment r→index).
            Cost backward = cost[r * dim_val + index] - u[r] - v[index];

            for (intptr_t k = 0; k < dim_val; ++k) {
                if (SC[k] || capacity[k] <= 0) continue;   // column k already settled, skip

                // Reduced cost of the forward arc r→k.
                Cost forward = cost[r * dim_val + k] - u[r] - v[k];

                // Candidate shortest-path cost to reach column k by
                // re-routing row r from column `index` to column k.
                Cost cand = shortestPathCost[index] - backward + forward;

                if (cand < shortestPathCost[k]) {
                    shortestPathCost[k] = cand;
                    path[k] = r;   // row r is the predecessor of column k
                }
            }
        }
    }
    return sink;
}


static Cost lapjv_ub(intptr_t dim_var, intptr_t dim_val,
                     const vector<Cost>& cost,
                     int* b,
                     Cost* usol, Cost* vsol,
                     Cost MAX_COST,
                     vector<int>& capacity,
                     int& findConflict)
{
    // Fast path: all capacities are 1 → standard (unit-capacity) LAP.
    int maxcap = *max_element(capacity.begin(), capacity.end());
    int mincap = *min_element(capacity.begin(), capacity.end());
    if (maxcap == 1 && mincap > 0) {
        return lapjv(dim_var, dim_val, cost, b, usol, vsol, MAX_COST, findConflict);
    }

    findConflict = 0;

    // col4row[i]   : column currently assigned to row i  (-1 = unassigned).
    // row4col[j]   : list of rows currently assigned to column j.
    // count_col[j] : current load of column j  (== row4col[j].size()).
    // u[i], v[j]   : dual potentials used to compute reduced costs.
    // shortestPathCost[j] : tentative shortest-path cost to column j.
    // path[j]             : predecessor row of column j on the shortest path.
    // SC[j]               : true if column j has been settled.
    // SR[i]               : true if row    i has been settled.

    vector<int>          col4row(dim_var, -1);
    vector<vector<int>>  row4col(dim_val);
    vector<int>          count_col(dim_val, 0);
    vector<Cost>         u(dim_var, 0), v(dim_val, 0);

    vector<intptr_t> path(dim_val, -1);
    vector<Cost>     shortestPathCost(dim_val);
    vector<bool>     SC(dim_val, false);
    vector<bool>     SR(dim_var, false);

    /* --- Main loop: assign each row --------------------------------------- */
    for (intptr_t curRow = 0; curRow < dim_var; ++curRow) {

        // Row already assigned (can happen if a previous augmentation
        if (col4row[curRow] != -1) continue;

        // Find the shortest augmenting path from curRow to some column
        // that still has remaining capacity.

        intptr_t sink = augmenting_path_gcc(dim_val, cost, u, v, curRow,
                                            path, shortestPathCost, SC, SR,
                                            col4row, row4col, MAX_COST, capacity);

        // No augmenting path found → problem is infeasible.
        // Report curRow and every row that competes with it for the same
        // columns so the caller can diagnose the conflict.
        if (sink < 0) {
            findConflict = 1;
            b[0] = curRow;
            for (int col = 0; col < dim_val; col++) {
                if (cost[curRow * dim_val + col] < MAX_COST) {
                    for (int var : row4col[col]) {
                        b[findConflict] = var;
                        findConflict++;
                    }
                }
            }
            return MAX_COST;
        }

        /* --- Dual variable update (complementary slackness) --------------- */

        // shortestPathCost[sink] is the length of the augmenting path; all
        // dual corrections below use this as the reference distance.
        Cost minVal = shortestPathCost[sink];

        // Shift the dual of the source row: u[curRow] increases by minVal
        // so that the reduced cost of the newly added arc becomes 0.
        u[curRow] += minVal;


        for (intptr_t r = 0; r < dim_var; ++r) {
            if (!SR[r] || r == curRow) continue;
            int colAssigned = col4row[r];   // pre-augmentation assignment
            if (colAssigned >= 0) {
                u[r] += (minVal - shortestPathCost[colAssigned]);
            }
        }

        // For every settled column j, shift v[j] so that the reduced cost
        // of its incoming arc (from the path) remains 0 after the update.
        for (intptr_t j = 0; j < dim_val; ++j) {
            if (!SC[j] || capacity[j] <= 0 ) continue;
            v[j] -= (minVal - shortestPathCost[j]);
        }

        /* --- Path augmentation -------------------------------------------- */

        // Walk the predecessor chain from `sink` back to `curRow`,
        // re-assigning each row along the path to the next column.
        intptr_t j = sink;
        while (true) {
            int r       = path[j];       // predecessor row of column j
            int old_col = col4row[r];    // column that row r was assigned to

            // Assign row r to column j.
            col4row[r] = j;
            row4col[j].push_back(r);
            count_col[j] = (int)row4col[j].size();

            // Remove row r from its previous column (if any).
            if (old_col != -1) {
                auto& vec = row4col[old_col];
                auto  it  = find(vec.begin(), vec.end(), r);
                if (it != vec.end()) vec.erase(it);
                count_col[old_col] = (int)vec.size();
            }

            // Stop once we reach the free row (no previous column).
            if (old_col == -1) break;
            j = old_col;
        }
    }

    /* --- Build output arrays ---------------------------------------------- */

    Cost total_cost = 0;
    for (intptr_t i = 0; i < dim_var; ++i) {
        int col      = col4row[i];
        b[i]         = col;
        total_cost  += cost[i * dim_val + col];
        usol[i]      = u[i];
    }

    for (intptr_t j = 0; j < dim_val; ++j) { 
        vsol[j] = capacity[j] > 0 ? v[j] : MAX_COST;

    }
    return total_cost;
}


//Linear assignment problem with excepted values

static Cost lapjv(intptr_t dim_var, intptr_t dim_val, const vector<Cost>& cost,  int* b,  Cost* usol, Cost* vsol, Cost MAX_COST, vector<int>& exceptedValIndex, int& findConflict)
{
    
    vector<int> capacity(dim_val, 1);
    for(int val: exceptedValIndex){
        capacity[val] = dim_var;
    } 
    return lapjv_ub(dim_var, dim_val, cost, b, usol, vsol, MAX_COST, capacity, findConflict);
}

/* -----------------------------------------------------------------------------
 * shortestPaths
 *
 * Bellman-based search that finds, from a "source" column `source` (which
 * has an unsatisfied lower-bound demand), the cheapest row `dest` that can be
 * re-routed to `source` without breaking any upper-bound constraint.
 *
 * Parameters
 * ----------
 * dim_var      Number of rows.
 * dim_val      Number of columns.
 * cost         Flat row-major cost matrix.
 * VarList      VarList[j] = list of rows that have a finite cost to column j
 *              (pre-computed adjacency list used for efficiency).
 * col4row      Current assignment row → column.
 * path         Output: predecessor array along the shortest path, size dim_var.
 *              path[v] = predecessor row of v, or -2 if v is directly adjacent
 *              to the source column.
 * MAX_COST     Sentinel for infinity.
 * source       The infeasible column whose demand is not yet met.
 *
 * Returns
 * -------
 * Index of the destination row `dest` whose rerouting resolves the demand
 * deficit at `source`, or -1 if no such row exists (infeasible).
 * -----------------------------------------------------------------------------
 */
static intptr_t BellmanShortestPaths(intptr_t dim_var, intptr_t dim_val,
                                    const vector<Cost>& cost,
                                    const vector<vector<int>>& VarList,
                                    int* col4row,
                                    vector<int>& path,
                                    vector<int>& demand,
                                    vector<int>& count_col,
                                    Cost MAX_COST,
                                    int source)
{
    // dist[v] : best known distance from `source` to row v.
    // Initialised to MAX_COST (= +∞) for all rows.
    vector<Cost> dist(dim_var, MAX_COST);
    path.assign(dim_var, -1);

    /* --- Initialisation: seed rows directly adjacent to `source` --------- */

    // A row is "directly reachable" from `source` if it has a finite-cost arc
    // to `source` AND is not already assigned to `source` (assigning it again
    // would not help satisfy source's demand).
    for (int v : VarList[source]) {
        if (col4row[v] == source) continue;
        Cost d = cost[v * dim_val + source];
        if (d < dist[v]) {
            dist[v] = d;
            path[v] = -2;   // -2 : predecessor is the virtual source node
        }
    }

    /* --- Bellman-Ford relaxation loop ------------------------------------ */

    // Standard Bellman-Ford runs (|V| - 1) relaxation passes over all edges.
    // Here |V| = dim_var (row nodes only; columns are implicit via col4row).
    // Each pass iterates over every row u that has been reached (dist[u] < MAX_COST),
    // then relaxes outgoing edges through u's currently assigned column.
    //
    // Convergence is guaranteed in at most (dim_var - 1) passes because the
    // longest simple path in the exchange graph visits at most dim_var nodes.
    for (int iter = 0; iter < dim_var - 1; ++iter) {

        bool updated = false;   // early-exit flag: stop if no relaxation occurred

        for (int u = 0; u < dim_var; ++u) {
            if (dist[u] == MAX_COST) continue;   

            int val = col4row[u];      
            if (val == source) continue; 

            Cost base = dist[u] - cost[u * dim_val + val];

            // Relax each row v that is a neighbour of column `val`.
            for (int v : VarList[val]) {
                if (col4row[v] == val) continue;   // v already assigned to val — skip

                // Guard against overflow when base is already large.
                if (base == MAX_COST) continue;

                Cost alt = base + cost[v * dim_val + val];

                if (alt < dist[v]) {
                    dist[v]  = alt;
                    path[v]  = u;   // u is the predecessor of v on this path
                    updated  = true;
                }
            }
        }

        if (!updated) break;   // no improvement in this pass → converged early
    }

    /* --- Identify best destination --------------------------------------- */

    // Scan all rows; pick the one that:
    //   1. Has been reached (dist[v] < MAX_COST).
    //   2. Its current column is NOT `source` (we need to pull flow FROM elsewhere).
    //   3. Its current column has a surplus over its lower-bound demand
    //      (count_col[col] > demand[col]), so removing it won't create a new
    //      infeasibility.
    //   4. The net distance to that column (dist[v] - cost[v][col4row[v]])
    //      is non-negative (ensures the path cost accounting is consistent)
    //      and is the minimum found so far.
    Cost distmin = MAX_COST;
    int  dest    = -1;
    for (int var = 0; var < dim_var; var++) {
        if (dist[var] == MAX_COST) continue;

        int   col      = col4row[var];
        Cost  distToVal = dist[var] - cost[dim_val * var + col];

        if (   col != source                            // not the source column
            && demand[col] < count_col[col]             // column has surplus
            && distToVal >= 0                           // consistent path cost
            && distToVal < distmin)
        {
            distmin = distToVal;
            dest    = var;
        }
    }
    return dest;
}



/* -----------------------------------------------------------------------------
 * sendFlow
 *
 * Applies the augmenting path found by shortestPaths() to the current
 * assignment, effectively rerouting one unit of flow so that the infeasible
 * column `notFeasVal` gains one more assigned row.
 *
 * The path is stored in the `path` array produced by shortestPaths():
 *   path[v] = predecessor of v on the shortest path (-2 means "directly
 *   adjacent to the source column").
 *
 * The operation performed is:
 *   1. Remove `dest` from its current column's assignment.
 *   2. Walk the path from `dest` backward, shifting each row one step along
 *      the path (each row takes its predecessor's old column).
 *   3. Assign the path's tail row to `notFeasVal`.
 *
 * Parameters
 * ----------
 * path         Predecessor array from shortestPaths(), size dim_var.
 * col4row      Current row → column assignment.  Modified in-place.
 * row4col      Current column → list of rows assignment.  Modified in-place.
 * count_col    Current load per column.  Modified in-place.
 * notFeasVal   The infeasible column that needs one more row.
 * (dest)       The destination row identified by shortestPaths() — assumed
 *              to be in scope as an outer variable in the original code.
 * -----------------------------------------------------------------------------
 */
static void sendVarFlow(vector<int>& path, int* col4row, vector<vector<int>>& row4col,
         vector<int>& count_col, int notFeasVal, int dest)
{
    // Remove `dest` from its current column (it will be re-assigned below).
    int oldVal = col4row[dest];
    int vtx = dest;
    auto &vec = row4col[oldVal];
    auto it = find(vec.begin(), vec.end(), dest);
    if (it != vec.end()) vec.erase(it);

    // Walk the path: each row along the path inherits the column of its
    // predecessor, effectively shifting assignments one step toward `source`.
    while (path[vtx] != -2) {
        int prev = path[vtx];
        int col  = col4row[prev];
        col4row[vtx] = col;
        auto &vec = row4col[col];
        auto it = std::find(vec.begin(), vec.end(), prev);
        if (it != vec.end()) {
            *it = vtx;   // replace prev with vtx in-place
        }
        vtx = prev;
    }

    // The tail of the path is assigned to the infeasible column.
    col4row[vtx] = notFeasVal;
    row4col[notFeasVal].push_back(vtx);
    count_col[notFeasVal]++;
    count_col[oldVal]--;
}


/* -----------------------------------------------------------------------------
 * checkFlow
 *
 * Scans the current assignment to find the first column whose load is below
 * its lower-bound demand.
 *
 * Parameters
 * ----------
 * dim_val      Number of columns.
 * demand       Lower-bound demand per column.
 * count_col    Current number of rows assigned to each column.
 *
 * Returns
 * -------
 * Index of the first infeasible column (count_col[j] < demand[j]), or -1 if
 * all demands are satisfied.
 * -----------------------------------------------------------------------------
 */
static intptr_t checkFlow(intptr_t dim_val,
                           vector<int>& demand,
                           vector<int>& count_col)
{
    for (intptr_t j = 0; j < dim_val; ++j) {
        if (count_col[j] < demand[j]) {
            return j;   // first column with unmet demand
        }
    }
    return -1;   // all demands are satisfied
}


/* =============================================================================
 * LAP-JV with Global Cardinality Constraints (GCC)
 *
 * Solves a min-cost assignment solver that handles:
 *   - upper-bound capacities per column  (each column j can serve at most
 *     capacity[j] rows)
 *   - lower-bound demands per column     (each column j must serve at least
 *     demand[j] rows)
 *
 * The base optimality structure follows:
 *   J-C. Régin, "Cost-Based Arc Consistency for Global Cardinality
 *   Constraints", Constraints 7, 2002, pp. 387-405.
 *
 * Pipeline
 * -------------------
 *  1. lapjv_gcc()
 *       |
 *       +-- lapjv_ub()          find a min-cost assignment respecting
 *       |                       upper-bound capacities only
 *       |
 *       +-- checkFlow()         verify every lower-bound demand is met
 *       |
 *       +-- (repair loop)
 *       |     shortestPaths()   Dijkstra from an infeasible column to find
 *       |                       the cheapest rerouting path
 *       |     sendFlow()        apply the rerouting along that path
 *       |__   checkFlow()       repeat until all demands are satisfied
 
 *
 *
 *
 * Algorithm outline
 * -----------------
 *  Step 1. Call lapjv_ub() to find a min-cost assignment satisfying only the
 *          upper-bound capacities.  This gives an optimal (but potentially
 *          demand-infeasible) solution with duals u, v.
 *
 *  Step 2. Build the row4col / count_col structures from the assignment b[].
 *
 *  Step 3. Call checkFlow() to detect the first column whose lower-bound
 *          demand is not met.
 *
 *  Step 4. While there exists an infeasible column (notFeasVal ≥ 0):
 *            a. Run shortestPaths() (Dijkstra on the exchange graph) from
 *               notFeasVal to find the cheapest row that can be rerouted to it
 *               without creating a new infeasibility elsewhere.
 *            b. Apply sendFlow() to reroute that flow.
 *            c. Re-check feasibility with checkFlow().
 *
 * Parameters
 * ----------
 * dim_var      Number of rows / variables (n).
 * dim_val      Number of columns / values (m).
 * cost         Flat row-major cost matrix of size (n × m).
 * b            Output array size n; b[i] = column assigned to row i.
 * usol         Output row dual potentials, size n.
 * vsol         Output column dual potentials, size m.
 * MAX_COST     Sentinel for infeasible / forbidden assignments.
 * capacity     Upper-bound capacities per column, size m.
 * demand       Lower-bound demands per column, size m.
 * findConflict Output: 0 = feasible solution found; otherwise the number of
 *              conflicting variables (stored in b[0..findConflict-1]).
 *
 * Returns
 * -------
 * Total assignment cost, or MAX_COST if no feasible solution exists.
 * =============================================================================
 */
static Cost lapjv_gcc(intptr_t dim_var, intptr_t dim_val,
                         const vector<Cost>& cost,
                         int* b,
                         Cost* usol, Cost* vsol,
                         Cost MAX_COST,
                         vector<int>& capacity,
                         vector<int>&  demand,
                         int& findConflict)
{
    vector<vector<int>> row4col(dim_val);      // column → list of assigned rows
    vector<int> count_col(dim_val, 0);         // current load per column
    Cost total_cost = 0;
    int  notFeasVal = -1;
    int  val;

    /* Step 1: Find a min-cost assignment under upper-bound constraints only */
    total_cost = lapjv_ub(dim_var, dim_val, cost, b, usol, vsol,
                          MAX_COST, capacity, findConflict);

    if(total_cost >= MAX_COST) return total_cost;   // already infeasible

    /* Step 2: Reconstruct row4col and count_col from output b[] */
    for(int var = 0; var < dim_var; var++){  
        val = b[var],                          
        row4col[val].push_back(var);
        count_col[val] += 1;
    }

    /* Step 3: Check whether all lower-bound demands are satisfied */
    notFeasVal = checkFlow(dim_val, demand, count_col);

    if (notFeasVal == -1) return total_cost;   // already feasible — done

    /* Build adjacency list VarList[j] = rows with finite cost to column j.
     * This avoids re-scanning the full cost matrix inside the repair loop. */
    vector<int> path;
    vector<vector<int>> VarList(dim_val);

    for(int val = 0; val < dim_val; val++){
        for(int var = 0; var < dim_var; var++){
            if (cost[dim_val * var + val] < MAX_COST){
                VarList[val].push_back(var);
            }
        }
    }

    /* Step 4: Repair loop — restore feasibility one column at a time */
    while(notFeasVal > -1){

        // (a) Find the cheapest rerouting from the infeasible column.
        int dest = BellmanShortestPaths(dim_var, dim_val, cost, VarList,
                                        b, path, demand, count_col, MAX_COST, notFeasVal);

        if(dest == -1) return MAX_COST;   // no rerouting possible → infeasible

        // (b) Apply the rerouting: shift flow along the discovered path.
        sendVarFlow(path, b, row4col, count_col, notFeasVal, dest);

        // (c) Re-check feasibility; continue if another infeasible column exists.
        notFeasVal = checkFlow(dim_val, demand, count_col);
    }

    return lapjv_ub(dim_var, dim_val, cost, b, usol, vsol,
                    MAX_COST, count_col, findConflict);


    return total_cost;
}

//Linear assignment problem with excepted values


#endif // LAPJV_HPP_

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
