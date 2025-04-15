#include "tb2hungarian.hpp"

// Constructor initializing the Hungarian algorithm with a disallowed cost.
// The disallowed cost (MAX_COST) is used to mark invalid or forbidden assignments.
Hungarian::Hungarian(Cost DISALLOWED) : n(0), Z0_r(0), Z0_c(0), MAX_COST(DISALLOWED) {}

// Executes the Hungarian algorithm to find the optimal assignment with minimal cost.
// Returns the total cost of the assignment or MAX_COST if no valid solution exists.
Cost Hungarian::compute(vector<vector<Cost>> &cost_matrix, vector<int> &supports) {
    C = cost_matrix;
    n = C.size();

    row_covered.assign(n, false);
    col_covered.assign(n, false);
    marked.assign(n, vector<int>(n, 0));
    path.assign(2 * n, vector<int>(2, 0));
    reduce_cost_row.assign(n, 0);
    reduce_cost_col.assign(n, 0);
    assignment.assign(n, 0);

    int step = 2;
    bool done = false;

    // Main step loop of the algorithm
    while (!done) {
        switch (step) {
            case 1: step = step1(); break;            // Step 1: Reduce each row
            case 2: step = step2(supports); break;    // Step 2: Star initial zeros
            case 3: step = step3(); break;            // Step 3: Cover columns with stars
            case 4: step = step4(); break;            // Step 4: Prime uncovered zeros
            case 5: step = step5(); break;            // Step 5: Build and augment paths
            case 6: step = step6(); break;            // Step 6: Adjust matrix
            case 7: done = true; break;               // Solution found
            case 8: done = true; break;               // No solution possible
        }
    }

    if (step == 8) {
        return MAX_COST; // No feasible assignment found
    }

    // Compute total cost from reduced costs
    Cost total_cost = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (marked[i][j] == 1) {
                assignment[i] = j;
            }
        }
        total_cost += -(reduce_cost_row[i] + reduce_cost_col[i]);
    }

    return total_cost;
}



// Returns the row reduction values (row dual solution u).
vector<Cost> Hungarian::getReduceCostRow() const {
    return reduce_cost_row;
}

// Returns the column reduction values (column dual solution v).
vector<Cost> Hungarian::getReduceCostCol() const {
    return reduce_cost_col;
}

// Returns the optimal column assignment for each row(primal solution x).
vector<int> Hungarian::getAssignment() const {
    return assignment;
}


// Step 1: Row reduction.
// For each row, subtract the smallest value from all elements in that row.
int Hungarian::step1() {
    for (int i = 0; i < n; i++) {
        Cost minval = *min_element(C[i].begin(), C[i].end());
        reduce_cost_row[i] -= minval;
        for (int j = 0; j < n; j++)
            if (C[i][j] < MAX_COST)
                C[i][j] -= minval;
    }
    return 2;
}

// Step 2: Star independent zeros in the matrix.
// First use the preferred supports vector, then find additional zeros if needed.
int Hungarian::step2(vector<int> &supports) {
    int count = 0;
    for (int i = 0; i < n; i++) {
        if (col_covered[supports[i]]) continue;
        count++;
        marked[i][supports[i]] = 1;
        row_covered[i] = true;
        col_covered[supports[i]] = true;
    }

    if (count == n)
        return 7;

    for (int i = 0; i < n; i++) {
        if (row_covered[i]) continue;
        for (int j = 0; j < n; j++) {
            if (!col_covered[j] && C[i][j] == 0) {
                marked[i][j] = 1;
                row_covered[i] = true;
                col_covered[j] = true;
                break;
            }
        }
    }

    clear_covers();
    return 3;
}

// Step 3: Cover columns containing starred zeros.
// If all columns are covered, the assignment is complete.
int Hungarian::step3() {
    int count = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (marked[i][j] == 1 && !col_covered[j]) {
                col_covered[j] = true;
                count++;
            }

    return (count >= n) ? 7 : 4;
}

// Step 4: Find and prime uncovered zeros.
// If a primed zero has no starred zero in its row, go to Step 5.
int Hungarian::step4() {
    while (true) {
        auto [row, col] = find_a_zero();
        if (row == -1)
            return 6;

        marked[row][col] = 2; // Prime the zero
        int star_col = find_star_in_row(row);
        if (star_col >= 0) {
            row_covered[row] = true;
            col_covered[star_col] = false;
        } else {
            Z0_r = row;
            Z0_c = col;
            return 5;
        }
    }
}

// Step 5: Build an alternating path of primed and starred zeros, then update the marks.
int Hungarian::step5() {
    int count = 0;
    path[count][0] = Z0_r;
    path[count][1] = Z0_c;

    bool done = false;
    while (!done) {
        int row = find_star_in_col(path[count][1]);
        if (row >= 0) {
            count++;
            path[count][0] = row;
            path[count][1] = path[count - 1][1];
        } else {
            done = true;
        }

        if (!done) {
            int col = find_prime_in_row(path[count][0]);
            count++;
            path[count][0] = path[count - 1][0];
            path[count][1] = col;
        }
    }

    convert_path(path, count);
    clear_covers();
    erase_primes();
    return 3;
}

// Step 6: Modify the matrix to create more zeros.
// Add the smallest uncovered value to covered rows, subtract from uncovered columns.
int Hungarian::step6() {
    Cost minval = find_smallest();
    if (minval == MAX_COST)
        return 8;

    vector<bool> col_modified(n, false);
    int events = 0;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (C[i][j] == MAX_COST) continue;
            if (row_covered[i]) {
                C[i][j] += minval;
                events++;
            }
            if (!col_covered[j]) {
                C[i][j] -= minval;
                events++;
                if (!col_modified[j]) {
                    reduce_cost_col[j] -= minval;
                    col_modified[j] = true;
                }
            }
            if (row_covered[i] && !col_covered[j]) events -= 2;
        }
        if (row_covered[i]) {
            reduce_cost_row[i] += minval;
        }
    }

    return (events == 0) ? 8 : 4;
}

// Finds the smallest uncovered value in the cost matrix.
Cost Hungarian::find_smallest() {
    Cost minval = MAX_COST;
    for (int i = 0; i < n; i++)
        if (!row_covered[i])
            for (int j = 0; j < n; j++)
                if (!col_covered[j] && C[i][j] < MAX_COST)
                    minval = min(minval, C[i][j]);
    return minval;
}

// Finds the position of the next uncovered zero in the matrix.
pair<int, int> Hungarian::find_a_zero() {
    for (int i = 0; i < n; i++)
        if (!row_covered[i])
            for (int j = 0; j < n; j++)
                if (!col_covered[j] && C[i][j] == 0)
                    return {i, j};
    return {-1, -1};
}

// Finds a starred zero in a given row, or returns -1.
int Hungarian::find_star_in_row(int row) {
    for (int j = 0; j < n; j++)
        if (marked[row][j] == 1)
            return j;
    return -1;
}

// Finds a starred zero in a given column, or returns -1.
int Hungarian::find_star_in_col(int col) {
    for (int i = 0; i < n; i++)
        if (marked[i][col] == 1)
            return i;
    return -1;
}

// Finds a primed zero in a given row, or returns -1.
int Hungarian::find_prime_in_row(int row) {
    for (int j = 0; j < n; j++)
        if (marked[row][j] == 2)
            return j;
    return -1;
}

// Converts a path of alternating primes and stars by flipping their status.
void Hungarian::convert_path(vector<vector<int>>& path, int count) {
    for (int i = 0; i <= count; i++)
        marked[path[i][0]][path[i][1]] = (marked[path[i][0]][path[i][1]] == 1) ? 0 : 1;
}

// Resets all row and column cover indicators.
void Hungarian::clear_covers() {
    fill(row_covered.begin(), row_covered.end(), false);
    fill(col_covered.begin(), col_covered.end(), false);
}

// Removes all prime marks (value 2) from the marked matrix.
void Hungarian::erase_primes() {
    for (auto& row : marked)
        for (auto& cell : row)
            if (cell == 2) cell = 0;
}
