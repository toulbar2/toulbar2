#ifndef BI_CRITERIA_HPP
#define BI_CRITERIA_HPP

#include "multicfn.hpp"

/*!
 * \class Bicriteria
 * \brief functions to compute supported points in bi-criteria problem
 */
class Bicriteria {

private:
    /*!
     * \brief default constructor
     */
    Bicriteria();

public: /* types and enum */
    static constexpr Double Delta = 1e-7; // constant for defining weights to compute the optimal points individually in the objectives (should be equal to 10**(-resolution))
    /*!
     * \brief type representing a point in the objective space
     */
    typedef pair<Double, Double> Point;

    /*!
     * \brief type representing a pair of weights for the scalarization
     */
    typedef pair<Double, Double> Weights;

    /*!
     * \brief represenation of the optimization direction: maximization or minimization
     */
    enum OptimDir { Optim_Max,
        Optim_Min };

private: /* static members to store solutions and points */
    static std::vector<Point> _points; // points in the objective space

    static std::vector<Weights> _weights; // the weights used to compute the points

    static std::vector<MultiCFN::Solution> _solutions; // solutions computed

    static std::vector<Double> _lower_bounds; // lower_bounds returned for each solution, equals sol_cost[0]*sol_weight[0]+sol_cost[1]*sol_weight[1]

    static unsigned int _first_cfn_index;

    static unsigned int _second_cfn_index;

private: /* static members to for parameters and global start time */
    static int _showSolutions; // save ToulBar2 showSolutions print intermediate solutions found

    static int _vac; // save ToulBar2 vac

    static int _seed; // save ToulBar2 seed

    static int _verbose; // save ToulBar2 verbosity

    static int _sol_timeout; // timeout in seconds for each solution

    static int _global_timeout; // timeout in seconds for the resolution

    static int _nsol_max; // limit on number of compute solutions

    static double _startGlobalCpuTime; // global CPU-time consumed by computing a Pareto front

private: /* static functions */
    /*!
     * \brief sort the solutions obtained
     * \param optim_dir the optimization direction
     */
    static void sortSolutions(pair<OptimDir, OptimDir> optim_dir);

    /*!
     * \brief tests if two points are equal or not
     * \param p1
     * \param p2
     * \return true if the two points are not equal to each other
     */
    static bool notEqual(Point p1, Point p2);

    /*!
     * \brief tests if two points are equal or not
     * \param p1
     * \param p2
     * \return true if the two points are equal to each other
     */
    static bool equal(Point p1, Point p2);

    /*!
     * \brief returns true if p1 dominates p2
     * \param p1
     * \param p2
     * \param optim_dir the optimization direction on the two criteria
     */
    static bool dominates(Point p1, Point p2, pair<OptimDir, OptimDir> optim_dir);

    /*!
     * \brief solve the scalarization of two cost function networks registered in a combiner with provided weights
     * \param multicfn the problem containing (at least) the two cost function networks
     * \param weights the weights applied to the two networks
     * \param solution the solution returned by the solver, if not null
     * \param point the point in the objective space computed by the solver, if not null
     * \param lb the lower bound returned by the solver, if not null
     * \return true if a solution has beend found, false otherwise
     */
    static bool solveScalarization(MultiCFN* multicfn, pair<Double, Double> weights, MultiCFN::Solution* solution = nullptr, Point* point = nullptr, Double* lb = nullptr);

public: /* static functions */
    /*!
     * \brief compute a list of supported points for a bi-objective cost function network
     * \param multicfn the problem containing (at least) the two cost function networks
     * \param first_cfn_index index of the first cfn to optimize in the multicfn object
     * \param second_cfn_index index of the second cfn to optimize in the multicfn object
     * \param optim_dir the optimization direction of the two objectives: Optim_Max or Optim_Dir
     * \param delta constant for defining weights to compute the optimal points individually in the objectives
     */
    static void computeSupportedPoints(MultiCFN* multicfn, unsigned int first_cfn_index, unsigned int second_cfn_index, pair<OptimDir, OptimDir> optim_dir, Double delta = Delta);

    /*!
     * \brief compute a list of supported points for a bi-objective cost function network when optimizing respectively the first and the second cfn in the multicfn object
     * \param multicfn the problem containing (at least) the two cost function networks
     * \param optim_dir the optimization direction of the two objectives: Optim_Max or Optim_Dir
     * \param delta constant for defining weights to compute the optimal points individually in the objectives
     */
    static void computeSupportedPoints(MultiCFN* multicfn, pair<OptimDir, OptimDir> optim_dir, Double delta = Delta);

    /*!
     * \brief compute additional (potentially non dominated) solutions via enumeration in a nondominated triangle
     * \param multicfn the bicriteria cost function network
     * \param optim_dir the direction of the optimization for the two criteria
     * \param solIndex the index of the solution from which searching new solutions
     * \param nbLimit maximum number of solutions to obtain
     * \param pct the percentage of the search space, 1.0 = complete triangle
     */
    static void computeAdditionalSolutions(MultiCFN* multicfn, pair<Bicriteria::OptimDir, Bicriteria::OptimDir> optim_dir, unsigned int solIndex, unsigned int nbLimit = 100, Double pct = 1.0);

    /*!
     * \brief compute all the non supported solutions
     * \param multicfn the bicriteria problem
     * \param optim_dir the optimization direction
     * \param nbLimit maximum number of points to enumerate
     */
    static void computeNonSupported(MultiCFN* multicfn, pair<OptimDir, OptimDir> optim_dir, unsigned int nbLimit = 100);

    /*!
     * \brief get the list of solutions computed
     * \return a vector of the solutions
     */
    static std::vector<MultiCFN::Solution> getSolutions();

    /*!
     * \brief get the list of points computed in the objective space
     * \return a vector of the points
     */
    static std::vector<Point> getPoints();

    /*!
     * \brief get the list of weights used to obtain the supported points
     * \return a vector of the pairs of weights
     */
    static std::vector<Weights> getWeights();

    /*!
     * \brief get the list of lower bounds returned for each solution (lb = weight_1 * cost_1 + weight_2*cost_2)
     * \return a vector of the lower bounds as Double
     */
    static std::vector<Double> getLowerBounds();

    /*!
     * \brief set a CPU timeout for the computation of each solution
     * \param timeout the timeout in seconds, 0 by default, meaning no timeout
     */
    static void setSolutionTimeout(int timeout) {
        _sol_timeout = timeout;
    }

    /*!
     * \brief set a global CPU timeout
     * \param timeout the timeout in seconds, 0 by default, meaning no timeout
     */
    static void setGlobalTimeout(int timeout) {
        _global_timeout = timeout;
    }

    /*!
     * \brief set the maximum number of solution to compute
     * \param solCount maximum number of solution, 0 by default, meaning no limit
     */
    static void setMaxSolutionCount(int solCount) {
        _nsol_max = solCount;
    }

    /*!
     * \brief set the _vac parameter
     * \param vac enable VAC in preprocessing and during search, 0 by default, meaning no VAC
     */
    static void setVAC(int vac) {
        _vac = vac;
    }

    /*!
     * \brief set the _seed parameter
     * \param random seed used in the solver, 1 by default
     */
    static void setSeed(int seed) {
        _seed = seed;
    }

    /*!
     * \brief set the _verbose parameter
     * \param verbose verbosity level, -1 by default, meaning no output messages
     */
    static void setVerbose(int verbose) {
        _verbose = verbose;
    }

    /*!
     * \brief set the _showSolutions parameter
     * \param showSolutions print intermediate solutions found during search, 0 by default, meaning solutions are not shown
     */
    static void setShowSolutions(int showSolutions) {
        _showSolutions = showSolutions;
    }

};

#endif // BI_CRITERIA_HPP
