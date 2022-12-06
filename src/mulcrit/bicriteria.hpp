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

    // static constexpr Double delta = 1e-3; // constant for defining weights to compute the optimal points individually in the objectives
    

  public: /* types and enum */

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
    enum OptimDir {Optim_Max, Optim_Min};


  private: /* static members to store solutions and points */

    static std::vector<Point> _points; // points in the objective space

    static std::vector<Weights> _weights; // the weights used to compute the points

    static std::vector<mulcrit::Solution> _solutions; // solutions computed

  private: /* static functions */

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
     * \param point the point in the objective space computed by the solver, if not null
     * \return true if a solution has beend found, false otherwise
     */
    static bool solveScalarization(mulcrit::MultiWCSP* multiwcsp, pair<Double,Double> weights, mulcrit::Solution* solution = nullptr, Point* point = nullptr);



  public: /* static functions */

    /*!
     * \brief compute a list of supported points for a bi-objective cost function network
     * \param multiwcsp the problem containing (at least) the two cost function networks
     * \param optim_dir the optimization direction of the two objectives: Optim_Max or Optim_Dir
     * \param delta constant for defining weights to compute the optimal points individually in the objectives
     */
    static void computeSupportedPoints(mulcrit::MultiWCSP* multiwcsp, pair<OptimDir, OptimDir> optim_dir, Double delta = 1e-3);

    /*!
     * \brief get the list of solutions computed
     * \return a vector of the solutions
     */
    static std::vector<mulcrit::Solution> getSolutions();

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

};

} // namespace mulcrit

#endif // BI_CRITERIA_HPP