/**
 * \file multiwcsp.hpp
 * \brief Data structure combining multiple wcsp objects in the same problem with weights
 */

#ifndef MULTI_WCSP_HPP
#define MULTI_WCSP_HPP

#include <string>

#include "toulbar2lib.hpp"

namespace mulcrit {

class MultiWCSP; // forward delaration

typedef std::map<std::string, std::string> Solution; // type representing a solution of a multi-cfn

/*!
 * \class Var
 * \brief store variable data: name, value names: the values can only be defined by their name (str)
 */
class Var {

  public:

    /*!
     * \brief constructor
     * \param multiwcsp the global problem
     */
    Var(MultiWCSP* multiwcsp);

    /*!
     * \brief print the variable data
     * \param os the stream to print to
     */
    void print(ostream& os);

    /*!
     * \brief number of values in the domaine of the variable
     * \return the number of values
     */
    unsigned int nbValues();


  public:

    MultiWCSP* multiwcsp; //!< pointer to the MultiWCSP instance

    std::string name;
    std::vector<std::string> domain_str;
    std::map<string, int> str_to_index;

};

/*!
 * \class CostFunction
 * \brief store cost function data: name, scope, costs
 */
class CostFunction {

  public:

    /*!
     * \brief constructor
     */
    CostFunction(MultiWCSP* multiwcsp);

    /*!
     * \brief print the cost function data
     * \brief os the stream to print to
     */
    void print(std::ostream& os);

    /*!
     * \brief return the cost of a given tuple
     * \param tuple the tuple
     */
    Double getCost(std::vector<unsigned int>& tuple);

    /*!
     * \brief return the arity of the cost function
     */
    unsigned int arity();

  public:

    MultiWCSP* multiwcsp;

    std::string name;
    std::vector<unsigned int> scope;


    // cost table
    Double default_cost;
    std::vector<Double> costs;
    std::vector<std::vector<unsigned int> > tuples; // value indexes of the variables

};


/*!
 * \class MultiWCSP
 * \brief store a combination of several cost function network
 */
class MultiWCSP {

  public:

    /*!
     * \brief default constructor
     */
    MultiWCSP();

    /*!
     * \brief constructor: build a multiwcsp combining wcsps given as input
     * \param wcsps a vector of wcsp objects to combine with weighted sums
     * \param weights a list of weights for each wcsp
     */
    MultiWCSP(std::vector<WCSP*>& wcsps, std::vector<Double>& weights);

    /*!
     * \brief add a wcsp to the network, create the variables if they do not exist
     * \param wcsp the wcsp to add
     * \param weight the weight of the wcsp in the objective function (sum of the cost functions)
     */
    void addWCSP(WCSP* wcsp, double weight);

    /*!
     * \brief set the weight of the cost functions of one of the network
     * \param wcsp_index the index of the network to modify
     * \param weight the new weight
     */
    void setWeight(unsigned int wcsp_index, double weight);

    /*!
     * \brief number of networks loaded in the combiner
     * \return the number of network
     */
    unsigned int nbNetworks();

    /*!
     * \brief number of variables in the problem
     * \return the number of variables
     */
    unsigned int nbVariables();

    /*!
     * \brief load the cfn from toulbar2
     * \param wcsp the original tb2 csp
     */
    void load_cfn(WCSP* wcsp);

    /*!
     * \brief print the cfn
     * \brief os the stream to print to
     */
    void print(std::ostream& os);

    /*!
     * \brief make a wcsp from the convex combination of all the wcsps
     */
    WeightedCSP* makeWeightedCSP();


    /*!
     * \brief convert a tuple to a cost index, rightmost value indexed first
     * \param variables the list of variables from the tuple
     * \param tuple the tuple: value indexes for each variable
     */
    unsigned int tupleToIndex(std::vector<Var*> variables, std::vector<unsigned int> tuple);

    /*!
     * \brief compute the solution of the separated problems from the combined tb2 wcsp
     * \param solver the tb2 solver
     * \param obj_value the values of the sub cost function networks, not provided if null
     * \param solution optional, the solution returned by the algorithm  
     */
    void getSolution(WeightedCSPSolver* solver, std::vector<Double>* obj_value = nullptr, Solution* solution = nullptr);

  private: /* private methods */

    /*!
     * \brief send the cfn to toulbar2
     * \param wcsp tb2 wcsp
     */
    void exportToWCSP(WCSP* wcsp);

  private:

    /*!
     * \brief add a cost function to the network
     * \param wcsp the tb2 wcsp
     * \param cstr the original tb2 cost function
     */
    void addCostFunction(WCSP* wcsp, Constraint* cstr);

    /*!
     * \brief compute a TOP (infinity) value for the internal representation of the cfns
     */
    Double computeTop();

  public: // public attributes

    // variables
    std::vector<Var> var; // variables
    std::map<std::string, int> var_index; // index of variables

    // cost functions
    std::vector<CostFunction> cost_function; // list of the cost functions
    std::map<std::string, unsigned int> cost_function_index; // map between cfn names and indices

  private: // private attributes

    std::vector<double> weights; // list of weights for all the loaded networks
    std::vector<std::string> network_names; // names of the networks
    std::vector<std::vector<unsigned int>> networks; // list of the cost function networks (function indexes for each network)
    std::vector<unsigned int> network_index; // index of the network for each cost function

    std::vector<double> _doriginal_lbs; // list of original lower bounds (as Double) for each network
    std::vector<Double> _original_costMultipliers; // list of cost multipliers of all the original wcsp

    unsigned int _tb2_decimalpoint; // precision of the wcsp

};

} // namespace mulcrit

#endif // MULTI_WCSP_HPP