#ifndef TB2LAPJV_HPP_
#define TB2LAPJV_HPP_

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
        if (minVal == MAX_COST) { // infeasible cost matrix
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
        if (minVal == MAX_COST) { // infeasible cost matrix
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



static intptr_t
augmenting_path_gcc(intptr_t dim_val, vector<Cost>& cost, vector<Cost>& u,
                vector<Cost>& v, vector<intptr_t>& path,
                vector<intptr_t>& row4col,
                vector<Cost>& shortestPathCosts, intptr_t i,
                vector<bool>& SR, vector<bool>& SC,
                vector<intptr_t>& remaining, Cost* p_minVal, Cost MAX_COST, vector<int>& capacity)
{
    Cost minVal = 0;
    intptr_t num_remaining = dim_val;
    for (intptr_t it = 0; it < dim_val; it++) {
        remaining[it] = dim_val - it - 1;
    }

    fill(SR.begin(), SR.end(), false);
    fill(SC.begin(), SC.end(), false);
    fill(shortestPathCosts.begin(), shortestPathCosts.end(), MAX_COST);

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

            if (shortestPathCosts[j] < lowest ||
                (shortestPathCosts[j] == lowest && (row4col[j] == -1 || capacity[j] != 0))) {
                lowest = shortestPathCosts[j];
                index = it;
            }
        }

        minVal = lowest;
        if (minVal >= MAX_COST) {
            return -1;
        }

        intptr_t j = remaining[index];
        if (row4col[j] == -1 || (capacity[j] != 0)) {
            sink = j;
            capacity[j]--;
        } else {
            i = row4col[j];
        }

        SC[j] = true;
        remaining[index] = remaining[--num_remaining];
    }

    *p_minVal = minVal;
    return sink;
}

Cost lapjv_gcc(intptr_t dim_var, intptr_t dim_val, vector<Cost>& cost,
           int* b, Cost* usol, Cost* vsol, Cost MAX_COST, vector<int>& capacity)
{
    vector<Cost> u(dim_var, 0);
    vector<Cost> v(dim_val, 0);
    vector<Cost> shortestPathCosts(dim_val);
    vector<intptr_t> path(dim_val, -1);
    vector<intptr_t> col4row(dim_var, -1);
    vector<intptr_t> row4col(dim_val, -1);
    vector<bool> SR(dim_var);
    vector<bool> SC(dim_val);
    vector<intptr_t> remaining(dim_val);

    for (intptr_t curRow = 0; curRow < dim_var; curRow++) {
        Cost minVal;
        intptr_t sink = augmenting_path_gcc(dim_val, cost, u, v, path, row4col,
                                        shortestPathCosts, curRow, SR, SC,
                                        remaining, &minVal, MAX_COST, capacity);
        if (sink < 0) {
            return MAX_COST;
        }

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
        Total_cost += cost[i * dim_val + col4row[i]];
        usol[i] = u[i];
    }
    for (intptr_t j = 0; j < dim_val; j++) {
        vsol[j] = v[j];
    }      
    return Total_cost;
 }
#endif // LAPJV_HPP_
