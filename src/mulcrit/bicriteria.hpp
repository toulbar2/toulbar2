#ifndef BI_CRITERIA_HPP
#define BI_CRITERIA_HPP

#include "mulcrit/multiwcsp.hpp"

namespace mulcrit {

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

    static constexpr Double delta = 1e-3; // constant for defining weights to compute the optimal points individually in the objectives

  public: /* types and enum */

    /*!
    * \brief type representing a point in the objective space
    */
    typedef pair<Double, Double> Point;

    /*!
    * \brief represenation of the optimization direction: maximization or minimization
    */
    enum OptimDir {Optim_Max, Optim_Min};

  public: /* static functions */

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
    * \brief solve the scalarization of two cost function networks registered in a combiner with provided weights
    * \param multiwcsp the problem containing (at least) the two cost function networks
    * \param weights the weights applied to the two networks
    * \param solution the solution returned by the solver, if not null
    * \return point the point in the objective space of the solution
    */
    static Point solveScalarization(mulcrit::MultiWCSP* multiwcsp, pair<Double,Double> weights, mulcrit::Solution* solution = nullptr);

    /*!
    * \brief compute a list of supported points for a bi-objective cost function network
    * \param multiwcsp the problem containing (at least) the two cost function networks
    * \param optim_dir the optimization direction of the two objectives: Optim_Max or Optim_Dir
    * \param supported_points the list of points computed by the algorithm
    * \param solutions solutions associated to the supported points
    */
    static void computeSupportedPoints(mulcrit::MultiWCSP* multiwcsp, pair<OptimDir, OptimDir> optim_dir, std::vector<Point>* supported_points = nullptr, std::vector<mulcrit::Solution>* solutions = nullptr);

};

} // namespace mulcrit

#endif // BI_CRITERIA_HPP