#include "tb2hungarian.hpp"

Hungarian::Hungarian(Cost limit_val) : n(0), Z0_r(0), Z0_c(0), total_cost(0), MAX_COST(limit_val) {}

// Computes the optimal assignment for a given cost matrix.
Cost Hungarian::compute(EnumeratedVariable** scope, int arity) {

    n = arity;
    
    row_covered.assign(n, false);
    col_covered.assign(n, false);
    marked.assign(n, vector<int>(n, 0));
    path.assign(2 * n, vector<int>(2, 0));
    reduce_cost_row.assign(n, 0); // Ligne : réduction
    reduce_cost_col.assign(n, 0); // Colonne : réduction
    results.assign(n, 0);

    int step = 2;
    bool done = false;

    while (!done) {
        switch (step) {
            case 2: step = step2(scope); break;
            case 3: step = step3(); break;
            case 4: step = step4(scope); break;
            case 5: step = step5(); break;
            case 6: step = step6(scope); break;
            case 7: done = true; break;
            case 8: done = true; break;
        }
    }

    if (step == 8) return MAX_COST;

    total_cost = 0;
	//Cost total_cost1 = 0;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (marked[i][j] == 1) {
                results[i] = j;

            }
        }
		 total_cost += -(reduce_cost_row[i] + reduce_cost_col[i]);
    }


    return total_cost;
}



// Gets the row reduction values.
vector<Cost> Hungarian::getReduceCostRow() const {
    return reduce_cost_row;
}
// Gets the column reduction values.
vector<Cost> Hungarian::getReduceCostCol() const {
    return reduce_cost_col;
}
// Get the optimal assignment
vector<int> Hungarian::getAssignment() const {
    return results;
}
// Gets the reduced cost matrix after applying row and column reductions.
/*vector<vector<int>> Hungarian::getReduceMatrix() const {
	return C;
}
*/
// Step 1.
// For each row of the matrix, find the smallest element and subtract it
// from every element in its row.  Go to Step 2
/*int Hungarian::step1() {
    for (int i = 0; i < n; i++) {
        Cost minval = *min_element(C[i].begin(), C[i].end());
        reduce_cost_row[i] -= minval;
        for (int j = 0; j < n; j++)
			if(C[i][j] < MAX_COST)
				C[i][j] -= minval;
    }
    return 2;
}*/

// Step 2.
// Find a zero (Z) in the matrix.  If there is no starred zero in its row
// or column, star Z.  Repeat for every element in the matrix.  Go to step 3.
int Hungarian::step2(EnumeratedVariable** scope) {
	
	int count = 0;
    for (int i = 0; i < n; i++) {
		int index_support = scope[i]->toIndex(scope[i]->getSupport());
		if (col_covered[index_support]) continue;
		count++;
		marked[i][index_support] = 1;
		row_covered[i] = true;
		col_covered[index_support] = true;
    }	
	if (count == n)
		return 7;
	
	
    for (int i = 0; i < n; i++) {
        if (row_covered[i]) continue;
        for (int j = 0; j < n; j++) {
            if (col_covered[j] || scope[i]->cannotbe(scope[i]->toValue(j))) continue;
            if (scope[i]->getCost(scope[i]->toValue(j)) == 0) {
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

// Step 3.
// Cover each column containing a starred zero.  If all columns are
// covered, the starred zeros describe a complete set of unique assignments.
// In this case, terminate the algorithm.  Otherwise, go to step 4.
int Hungarian::step3() {
    int count = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (marked[i][j] == 1 ) {
				if(!col_covered[j]){
					col_covered[j] = true;
				}
                count++;
            }

    return (count >= n) ? 7 : 4;
}

// Step 4.
// Find a noncovered zero and prime it.  If there is no starred zero in the
// row containing this primed zero, Go to Step 5.  Otherwise, cover this row
// and uncover the column containing the starred zero. Continue in this manner
// until there are no uncovered zeros left, then go to Step 6.
int Hungarian::step4(EnumeratedVariable** scope) {
	while (true) {
		auto [row, col] = find_a_zero(scope);
		if (row == -1)
			return 6;

		marked[row][col] = 2;
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

// Step 5.
// Construct a series of alternating primed and starred zeros
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
			path[count][1] = path[count-1][1];
		} else {
			done = true;
		}

		if (!done) {
			int col = find_prime_in_row(path[count][0]);
			count++;
			path[count][0] = path[count-1][0];
			path[count][1] = col;
		}
	}

	convert_path(path, count);
	clear_covers();
	erase_primes();
	return 3;
}

// Step 6
// Add the smallest uncovered value in the matrix to every element of each
// covered row, and subtract it from every element of each uncovered column.
// Return to Step 4 without altering any stars, primes, or covered lines.
int Hungarian::step6(EnumeratedVariable** scope) {
    Cost minval = find_smallest(scope);
    if (minval == MAX_COST)
        return 8;

    int events = 0;
    vector<bool> col_modified(n, false);

    for (int i = 0; i < n; ++i) {
        bool rowCovered = row_covered[i];
        if (rowCovered) {
            reduce_cost_row[i] += minval;
        }

        for (int j = 0; j < n; ++j) {
            if (scope[i]->cannotbe(scope[i]->toValue(j))) continue;

            bool colCovered = col_covered[j];

            if (rowCovered && !colCovered) {
                // Compensation entre +minval (ligne) et -minval (colonne)
                // donc aucune modification effective
                continue;
            }

            // Simulation d'événements via les décalages
            if (rowCovered) {
                ++events;
            }

            if (!colCovered) {
                ++events;
                if (!col_modified[j]) {
                    reduce_cost_col[j] -= minval;
                    col_modified[j] = true;
                }
            }
        }
    }

    return (events == 0) ? 8 : 4;
}


// Find the smallest uncovered cell in the matrix.
Cost Hungarian::find_smallest(EnumeratedVariable** scope) {
    Cost minval = MAX_COST;
    for (int i = 0; i < n; i++) {
        if (!row_covered[i]) {
            for (int j = 0; j < n; j++) {
                if (!col_covered[j] && scope[i]->canbe(scope[i]->toValue(j))) {
                    Cost reduced = scope[i]->getCost(scope[i]->toValue(j)) + reduce_cost_row[i] + reduce_cost_col[j];
                    if (reduced < minval)
                        minval = reduced;
                }
            }
        }
    }
    return minval;
}



// Find an uncovered zero and store its co-ordinates in (zeroRow, zeroCol)
pair<int, int> Hungarian::find_a_zero(EnumeratedVariable** scope) {
        for (int i = 0; i < n; i++)
            if (!row_covered[i])
                for (int j = 0; j < n; j++)
					if(!col_covered[j] && scope[i]->canbe(scope[i]->toValue(j))){
						if (scope[i]->getCost(scope[i]->toValue(j)) + reduce_cost_row[i] + reduce_cost_col[j] == 0)
							return {i, j};
					}
        return {-1, -1};
}

// Find a column in row 'row' containing a star, or return
int Hungarian::find_star_in_row(int row) {
	for (int j = 0; j < n; j++)
		if (marked[row][j] == 1)
			return j;
	return -1;
}

// Find a row in column 'col' containing a star, or return
int Hungarian::find_star_in_col(int col) {
	for (int i = 0; i < n; i++)
		if (marked[i][col] == 1)
			return i;
	return -1;
}

// Find a column in row containing a prime, or return
int Hungarian::find_prime_in_row(int row) {
	for (int j = 0; j < n; j++)
		if (marked[row][j] == 2)
			return j;
	return -1;
}

void Hungarian::convert_path(vector<vector<int>>& path, int count) {
	for (int i = 0; i <= count; i++)
		marked[path[i][0]][path[i][1]] = (marked[path[i][0]][path[i][1]] == 1) ? 0 : 1;
}

// Uncovery ever row and column in the matrix.
void Hungarian::clear_covers() {
	fill(row_covered.begin(), row_covered.end(), false);
	fill(col_covered.begin(), col_covered.end(), false);
}

// Remove the prime marks from every cell in the matrix.
void Hungarian::erase_primes() {
	for (auto& row : marked)
		for (auto& cell : row)
			if (cell == 2) cell = 0;
}
