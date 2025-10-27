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
using namespace std;

/**
 * @brief Implements the Jonker Algorithm (LAPJV) to solve the assignment problem.
 *
 * The Jonker Algorithm provides an optimal solution to the assignment problem,
 * which consists of assigning tasks to agents in such a way that the total cost is minimized.
 * This class handles the transformation of a square cost matrix and applies the algorithm
 * through a well-defined series of steps.
 * Original source code from https://github.com/scipy/scipy/tree/main/scipy/optimize/rectangular_lsap
 */
 static intptr_t
augmenting_path(intptr_t dim_val, vector<Cost>& cost, vector<Cost>& u,
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
        // constant cost matrix is the identity matrix (c.f. #11602).
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

Cost lapjv(intptr_t dim_var, intptr_t dim_val, vector<Cost>& cost,  int* b,  Cost* usol, Cost* vsol, Cost MAX_COST)
{

    // initialize variables
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

//Linear assignment problem with excepted values
static intptr_t
augmenting_path(intptr_t dim_val, vector<Cost>& cost, vector<Cost>& u,
                vector<Cost>& v, vector<intptr_t>& path,
                vector<intptr_t>& row4col,
                vector<Cost>& shortestPathCosts, intptr_t i,
                vector<bool>& SR, vector<bool>& SC,
                vector<intptr_t>& remaining, Cost* p_minVal,Cost MAX_COST, vector<int>& exceptedValIndex)
{
    Cost minVal = 0;

    // Crouse's pseudocode uses set complements to keep track of remaining
    // nodes.  Here we use a vector, as it is more efficient in C++.
    intptr_t num_remaining = dim_val;
    for (intptr_t it = 0; it < dim_val; it++) {
        // Filling this up in reverse order ensures that the solution of a
        // constant cost matrix is the identity matrix (c.f. #11602).
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
           auto val = find(exceptedValIndex.begin(), exceptedValIndex.end(), j);

            // When multiple nodes have the minimum cost, we select one which
            // gives us a new sink node. This is particularly important for
            // integer cost matrices with small co-efficients.
            if (shortestPathCosts[j] < lowest ||
                (shortestPathCosts[j] == lowest && (row4col[j] == -1 || (val != exceptedValIndex.end())))) {
                lowest = shortestPathCosts[j];
                index = it;
            }
        }

        minVal = lowest;
        if (minVal == MAX_COST) { // MAX_COSTeasible cost matrix
            return -1;
        }

        intptr_t j = remaining[index];
		 auto val = find(exceptedValIndex.begin(), exceptedValIndex.end(), j);
        if (row4col[j] == -1 || (val != exceptedValIndex.end())) {
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

Cost lapjv(intptr_t dim_var, intptr_t dim_val, vector<Cost>& cost,  int* b,  Cost* usol, Cost* vsol, Cost MAX_COST, vector<int>& exceptedValIndex)
{

    // initialize variables
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
                                        remaining, &minVal, MAX_COST, exceptedValIndex);
        if (sink < 0) {
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



/*
  Solve assignment with column capacities.

  Parameters:
    dim_var   : number of rows to assign (n)
    dim_val   : number of original columns (m)
    cost      : vector of size (dim_var * dim_val) storing row-major cost matrix
    b         : output array of size dim_var; b[i] = assigned original column for row i
    usol      : output potentials for rows (size dim_var)
    vsol      : output potentials for original columns (size dim_val)
    MAX_COST  : large sentinel cost returned on impossible cases
    capacity  : vector<int> size dim_val, capacity[j] >= 0

  Returns:
    total cost (sum of assigned costs) or MAX_COST if not feasible (sum capacities < dim_var)
*/
static intptr_t
augmenting_path_gcc(intptr_t dim_val,                 // number of columns (m)
                     const vector<Cost>& cost,       // row-major (rows * cols) costs
                     vector<Cost>& u,                // row potentials, size dim_var
                     vector<Cost>& v,                // column potentials, size dim_val
                     intptr_t i,               // free row to start from
                     vector<intptr_t>& path,        // predecessor row for each column (size dim_val, -1 if none)
                     vector<Cost>& shortestPathCost,     // shortest reduced cost to each column
                     vector<bool>& SC,               // visited columns in this search
                     vector<bool>& SR,               // visited rows in this search
                     const vector<int>& col4row, // current assignment row -> col (-1 if free)
                     const vector<vector<int>>& row4col, // current rows assigned to each column
                     Cost MAX_COST,
                     const vector<int>& capacity)    // capacity per column
{
    // initialize
    fill(SR.begin(), SR.end(), false);
    fill(SC.begin(), SC.end(), false);
    fill(shortestPathCost.begin(), shortestPathCost.end(), MAX_COST);

    SR[i] = true;
    for (intptr_t j = 0; j < dim_val; ++j) {
        path[j] = -1;
        Cost r = cost[i * dim_val + j] - u[i] - v[j];
        if (r < shortestPathCost[j]) {
            shortestPathCost[j] = r;
            path[j] = i;
        }
    }

    // find shortest augmenting path
    intptr_t sink = -1;
    while (true) {
        Cost lowest = MAX_COST;
        intptr_t index = -1;
        for (intptr_t j = 0; j < dim_val; ++j) {
            if (!SC[j] && shortestPathCost[j] < lowest) {
                lowest = shortestPathCost[j];
                index = j;
            }
        }
        if (index == -1) {
            return -1;
        }

        SC[index] = true;

        if ((int)row4col[index].size() < capacity[index]) {
            sink = index;
            break;
        }

        for (int r : row4col[index]) {
            if (SR[r]) continue;
            SR[r] = true;
            for (intptr_t k = 0; k < dim_val; ++k) {
                if (SC[k]) continue;
                Cost cand = shortestPathCost[index] + (cost[r * dim_val + k] - u[r] - v[k]);
                if (cand < shortestPathCost[k]) {
                    shortestPathCost[k] = cand;
                    path[k] = r;
                }
            }
        }
    }

    return sink;
}


Cost lapjv_gcc(intptr_t dim_var, intptr_t dim_val,
                         const vector<Cost>& cost,
                         int* b,           // size dim_var, output: column assigned to row i
                         Cost* usol, Cost* vsol,      // output duals (u size n, v size m)
                         Cost MAX_COST,                // Top
                         vector<int>& capacity)         // capacities size m
{
    if (dim_var <= 0) return 0;
    if (dim_val <= 0) return MAX_COST;

    // feasibility check:
    long long total_cap = 0;
    for (intptr_t j = 0; j < dim_val; ++j) total_cap += capacity[j];
    if (total_cap < dim_var) return MAX_COST;

    vector<int> col4row(dim_var, -1);       // row -> col (or -1)
    vector<vector<int>> row4col(dim_val);      // list of rows assigned to column j
    vector<int> count_col(dim_val, 0);
    vector<Cost> u(dim_var, 0), v(dim_val, 0);    // duals  solutions
    vector<intptr_t> path(dim_val, -1);           // predecessor row for each column during search
    vector<Cost> shortestPathCost(dim_val);      // shortest reduced cost to each column
    vector<bool> SC(dim_val, false);                  // visited columns
    vector<bool> SR(dim_var, false);                  // visited rows

    // main loop: assign each row
    for (intptr_t curRow = 0; curRow < dim_var; ++curRow) {
        if (col4row[curRow] != -1) continue;

        // find shortest augmenting path from curRow to some column with free capacity
        intptr_t sink = augmenting_path_gcc(dim_val, cost, u, v, curRow,
                                          path, shortestPathCost, SC, SR,
                                          col4row, row4col,MAX_COST, capacity);
        if (sink < 0) {
            return MAX_COST;
        }

        Cost minVal = shortestPathCost[sink];
        u[curRow] += minVal;
        for (intptr_t i = 0; i < dim_var; ++i) {
            if (!SR[i] || i == curRow) continue;
            int colAssigned = col4row[i];
            if (colAssigned >= 0) {
                u[i] += (minVal - shortestPathCost[colAssigned]);
            }
        }
        for (intptr_t j = 0; j < dim_val; ++j) {
            if (!SC[j]) continue;
            v[j] -= (minVal - shortestPathCost[j]);
        }


        intptr_t j = sink;
        while (true) {
            int i = path[j];
            int old_col = col4row[i];
            col4row[i] = j;
            row4col[j].push_back(i);
            count_col[j] = (int)row4col[j].size();
            if (old_col != -1) {
                auto &vec = row4col[old_col];
                auto it = find(vec.begin(), vec.end(), i);
                if (it != vec.end()) vec.erase(it);
                count_col[old_col] = (int)vec.size();
            }
            if (old_col == -1) break;
            j = old_col;
        }
    }
    // Build outputs and compute total cost
    Cost total_cost = 0;
    for (intptr_t i = 0; i < dim_var; ++i) {
        int col = col4row[i];
        b[i] = col;
        total_cost += cost[i * dim_val + col];
        usol[i] = u[i];
    }
    for (intptr_t j = 0; j < dim_val; ++j) {
        vsol[j] = v[j];
    }

    return total_cost;
}

#endif // LAPJV_HPP_
