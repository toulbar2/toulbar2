#ifndef TB2HUNGARIAN_HPP_
#define TB2HUNGARIAN_HPP_

#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>
#include <fstream>

using namespace std;

/**
 * @class Hungarian
 * @brief Implements the Hungarian Algorithm for solving the assignment problem.
 *
 * The Hungarian Algorithm is used to find the optimal assignment that minimizes cost.
 * It follows a series of steps to transform the cost matrix and find the best matching.
 * source : https://github.com/bmc/munkres
 */
class Hungarian {
private:
    vector<vector<int>> C; ///< Cost matrix.
    vector<bool> row_covered; ///< Boolean vector for row coverage.
    vector<bool> col_covered; ///< Boolean vector for column coverage.
    vector<vector<int>> marked; ///< Matrix for marking starred and primed zeros.
    vector<vector<int>> path; ///< Path matrix used in step 5.
    vector<int> reduce_cost_row; ///< Row reduction values.
    vector<int> reduce_cost_col; ///< Column reduction values.
    int n; ///< Matrix size.
    int Z0_r, Z0_c; ///< Coordinates for the zero used in step 5.
    int total_cost; ///< Total minimal cost.

public:
    /**
     * @brief Default constructor.
     * Initializes an empty Hungarian solver instance.
     */
    Hungarian();

    /**
     * @brief Computes the optimal assignment for a given cost matrix.
     * @param cost_matrix The matrix representing the cost of assignments.
     * @return A vector of pairs (row, column) representing the optimal assignments.
     */
    vector<int> compute(vector<vector<int>> cost_matrix);

    /**
     * @brief Gets the total minimal cost of the assignment.
     * @return The minimal total cost.
     */
    int getTotalCost() const;

    /**
     * @brief Gets the row reduction values.
     * @return A vector representing the reduction applied to each row.
     */
    vector<int> getReduceCostRow() const;

    /**
     * @brief Gets the column reduction values.
     * @return A vector representing the reduction applied to each column.
     */
    vector<int> getReduceCostCol() const;

    /**
     * @brief Gets the reduced cost matrix after applying row and column reductions.
     * @return The reduced cost matrix.
     */
    vector<vector<int>> getReduceMatrix() const;

private:
    /**
     * @brief Step 1: Row reduction.
     * Subtracts the minimum value in each row from all elements in that row.
     * @return The next step number.
     */
    int step1();

    /**
     * @brief Step 2: Star zeros.
     * Finds and marks zeros following the algorithm's rules.
     * @return The next step number.
     */
    int step2();

    /**
     * @brief Step 3: Cover columns with starred zeros.
     * Determines if an optimal assignment has been found.
     * @return The next step number or termination.
     */
    int step3();

    /**
     * @brief Step 4: Prime uncovered zeros.
     * Searches for uncovered zeros and primes them if needed.
     * @return The next step number.
     */
    int step4();

    /**
     * @brief Step 5: Augment paths and update assignments.
     * Builds an alternating sequence of primed and starred zeros.
     * @return The next step number.
     */
    int step5();

    /**
     * @brief Step 6: Adjust the matrix values.
     * Adds the smallest uncovered value to covered rows and subtracts it from uncovered columns.
     * @return The next step number.
     */
    int step6();

    /**
     * @brief Finds the smallest uncovered value in the matrix.
     * @return The smallest uncovered value.
     */
    int find_smallest();

    /**
     * @brief Finds an uncovered zero in the matrix.
     * @return A pair (row, column) of the uncovered zero.
     */
    pair<int, int> find_a_zero();

    /**
     * @brief Finds a starred zero in a given row.
     * @param row The row index.
     * @return The column index of the starred zero or -1 if none.
     */
    int find_star_in_row(int row);

    /**
     * @brief Finds a starred zero in a given column.
     * @param col The column index.
     * @return The row index of the starred zero or -1 if none.
     */
    int find_star_in_col(int col);

    /**
     * @brief Finds a primed zero in a given row.
     * @param row The row index.
     * @return The column index of the primed zero or -1 if none.
     */
    int find_prime_in_row(int row);

    /**
     * @brief Converts the alternating sequence of starred and primed zeros.
     * @param path The path matrix.
     * @param count Number of elements in the path.
     */
    void convert_path(vector<vector<int>>& path, int count);

    /**
     * @brief Clears all row and column covers.
     */
    void clear_covers();

    /**
     * @brief Erases all primed zeros in the matrix.
     */
    void erase_primes();
};

#endif // HUNGARIAN_HPP_
