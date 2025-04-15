#ifndef TB2HUNGARIAN_HPP_
#define TB2HUNGARIAN_HPP_

#include <vector>
#include <algorithm>
#include "core/tb2enumvar.hpp"

/**
 * @class Hungarian
 * @brief Implements the Hungarian Algorithm (Munkres Algorithm) to solve the assignment problem.
 *
 * The Hungarian Algorithm provides an optimal solution to the assignment problem,
 * which consists of assigning tasks to agents in such a way that the total cost is minimized.
 * This class handles the transformation of a square cost matrix and applies the algorithm
 * through a well-defined series of steps.
 * Source code in python: https://github.com/bmc/munkres
 */
class Hungarian {
private:
    vector<vector<Cost>> C; ///< Input cost matrix.
    vector<bool> row_covered; ///< Tracks whether each row is currently covered.
    vector<bool> col_covered; ///< Tracks whether each column is currently covered.
    vector<vector<int>> marked; ///< Matrix used to mark starred and primed zeros.
    vector<int> assignment; ///< Stores the column assigned to each row in the optimal solution.
    vector<vector<int>> path; ///< Used to store the alternating path of starred and primed zeros (Step 5).
    vector<Cost> reduce_cost_row; ///< Reduction values applied to each row.
    vector<Cost> reduce_cost_col; ///< Reduction values applied to each column.
    int n; ///< Size of the square cost matrix.
    int Z0_r, Z0_c; ///< Coordinates of an uncovered zero used in Step 5.
    Cost MAX_COST; ///< Value representing disallowed or invalid assignments.

public:
    /**
     * @brief Constructs a Hungarian solver instance with a disallowed cost value.
     * @param DISALLOWED The cost value used to indicate an invalid or forbidden assignment.
     */
    Hungarian(Cost DISALLOWED);

    /**
     * @brief Computes the optimal assignment based on the given cost matrix.
     * @param cost_matrix A square matrix representing the cost of assigning rows to columns.
     * @param supports A preference or weighting vector that may guide the initial zero marking.
     * @return The total cost of the optimal assignment, or MAX_COST if no valid assignment is found.
     */
    Cost compute(vector<vector<Cost>> &cost_matrix, vector<int> &supports);

    /**
     * @brief Retrieves the row reduction values applied during the algorithm.
     * @return A vector containing the reduction value for each row.
     */
    vector<Cost> getReduceCostRow() const;

    /**
     * @brief Retrieves the column reduction values applied during the algorithm.
     * @return A vector containing the reduction value for each column.
     */
    vector<Cost> getReduceCostCol() const;

    /**
     * @brief Returns the optimal assignments as computed by the algorithm.
     * @return A vector where each index represents a row and the value is the assigned column.
     */
    vector<int> getAssignment() const;

private:
    /**
     * @brief Step 1: Row reduction.
     * Subtracts the minimum value from each row to ensure at least one zero per row.
     * @return The next step number.
     */
    int step1();

    /**
     * @brief Step 2: Star initial zeros.
     * Identifies and marks independent zeros as "starred" based on the algorithm's rules.
     * @param supports A vector of row preferences used to guide zero selection.
     * @return The next step number.
     */
    int step2(vector<int> &supports);

    /**
     * @brief Step 3: Cover all columns with starred zeros.
     * Determines whether the current assignment is complete.
     * @return The next step number, or a signal to terminate if assignment is optimal.
     */
    int step3();

    /**
     * @brief Step 4: Prime uncovered zeros.
     * Searches for uncovered zeros, primes them, and prepares for path augmentation.
     * @return The next step number.
     */
    int step4();

    /**
     * @brief Step 5: Construct an alternating path and update starred zeros.
     * Builds and processes a sequence of alternating primed and starred zeros.
     * @return The next step number.
     */
    int step5();

    /**
     * @brief Step 6: Adjust matrix values.
     * Modifies the cost matrix to create additional zeros by adjusting uncovered elements.
     * @return The next step number.
     */
    int step6();

    /**
     * @brief Finds the smallest uncovered value in the matrix.
     * @return The smallest uncovered cost value.
     */
    Cost find_smallest();

    /**
     * @brief Finds an uncovered zero in the cost matrix.
     * @return A pair (row, column) of an uncovered zero location.
     */
    pair<int, int> find_a_zero();

    /**
     * @brief Finds a starred zero in the given row.
     * @param row The row index.
     * @return The column index of the starred zero, or -1 if none exists.
     */
    int find_star_in_row(int row);

    /**
     * @brief Finds a starred zero in the given column.
     * @param col The column index.
     * @return The row index of the starred zero, or -1 if none exists.
     */
    int find_star_in_col(int col);

    /**
     * @brief Finds a primed zero in the given row.
     * @param row The row index.
     * @return The column index of the primed zero, or -1 if none exists.
     */
    int find_prime_in_row(int row);

    /**
     * @brief Converts the current path of primed and starred zeros into a new star configuration.
     * @param path The path matrix.
     * @param count The number of points in the path.
     */
    void convert_path(vector<vector<int>>& path, int count);

    /**
     * @brief Clears all row and column coverage.
     */
    void clear_covers();

    /**
     * @brief Erases all primed zeros from the marked matrix.
     */
    void erase_primes();
};


#endif // HUNGARIAN_HPP_
