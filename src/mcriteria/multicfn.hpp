/**
 * \file multicfn.hpp
 * \brief Data structure combining multiple wcsp objects in the same problem with weights
 */

#ifndef MULTI_CFN_HPP
#define MULTI_CFN_HPP

#include "toulbar2lib.hpp"

#ifdef ILOGCPLEX
#include <ilcplex/ilocplex.h>
#endif

class MultiCFN; // forward delaration

// namespace preventing collisions with tb2 Var and Cost functions classes
namespace mcriteria {

/*!
 * \class Var
 * \brief store variable data: name, value names: the values can only be defined by their name (str)
 */
class Var {

public:
    /*!
     * \brief constructor
     * \param multicfn the global problem
     */
    Var(MultiCFN* multicfn);

    /*!
     * \brief print the variable data
     * \param os the stream to print to
     */
    void print(ostream& os);

    /*!
     * \brief number of values in the domain of the variable
     * \return the number of values
     */
    unsigned int nbValues();

public:
    MultiCFN* multicfn; //!< pointer to the MultiCFN instance

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
    enum Type { Tuple,
        Linear };

public:
    /*!
     * \brief constructor
     */
    CostFunction(MultiCFN* multicfn, unsigned int net_index);

    /*!
     * \brief destructor
     */
    virtual ~CostFunction() {}

    /*!
     * \brief print the cost function data
     * \brief os the stream to print to
     */
    void print(std::ostream& os);

    /*!
     * \brief return the arity of the cost function
     */
    virtual unsigned int arity();

    /*!
     * \brief return the cost of a given tuple
     * \param tuple the tuple
     */
    virtual Double getCost(std::vector<unsigned int>& tuple) = 0;

    /*!
     * \brief return the type of the cost function
     */
    virtual Type getType() = 0;

public:
    MultiCFN* multicfn;
    std::string name;
    unsigned int net_index; /* index of the network this cost function belongs to */
    std::vector<unsigned int> scope; /* internal variable indexes */
};

/*!
 * \class TupleCostFunction
 * \brief store cost function data: name, scope, costs
 */
class TupleCostFunction : public CostFunction {

public:
    /*!
     * \brief constructor
     */
    TupleCostFunction(MultiCFN* multicfn, unsigned int net_index);

    /*!
     * \brief print the cost function data
     * \brief os the stream to print to
     */
    void print(std::ostream& os);

    /*!
     * \brief return the cost of a given tuple
     * \param tuple the tuple
     */
    virtual Double getCost(std::vector<unsigned int>& tuple);

    /*!
     * \brief get the max and min cost of the cost function
     * \param min_cost the minimum cost of the function returned by the method
     * \param min_cost the maximum cost of the function returned by the method
     */
    void getMinMaxCost(double& min_cost, double& max_cost);

    /*!
     * \brief compute the total number of tuples
     */
    size_t compute_n_tuples();

    /*!
     * \brief check in the costs if the cost function is a hard constraint
     */
    bool detectIfHard();

    /*!
     * \brief return the type of the cost function
     */
    virtual Type getType()
    {
        return CostFunction::Type::Tuple;
    }

public:
    // cost table
    Double default_cost;
    std::vector<Double> costs;
    std::vector<std::vector<unsigned int>> tuples; // value indexes of the variables
    size_t n_total_tuples; // total number of tuples (length of the domain cartesian product)
    bool all_tuples; // true if all tuples are stored, false otherwise (default_cost)
    bool hard; // true if all costs are either null or =infinity, the cost function is then used only as a constraint
};

/*!
 * \class LinearCostFunction
 * \brief store a linear cost function data: name, scope, coefficients
 */
class LinearCostFunction : public CostFunction {

public:
    /*!
     * \brief constructor
     */
    LinearCostFunction(MultiCFN* multicfn, unsigned int net_index);

    /*!
     * \brief print the cost function data
     * \brief os the stream to print to
     */
    void print(std::ostream& os);

    /*!
     * \brief return the cost of a given tuple
     * \param tuple the tuple
     */
    virtual Double getCost(std::vector<unsigned int>& tuple);

    /*!
     * \brief return the type of the cost function
     */
    virtual Type getType()
    {
        return CostFunction::Type::Linear;
    }

public:
    Double capacity;
    std::vector<std::vector<std::pair<unsigned int, Double>>> weights; // for each variable of the scope, a list of values and their weight
};

} // namespace mcriteria

/*!
 * \class MultiCFN
 * \brief store a combination of several cost function networks.
 * one new class: multiwcsp to create a wcsp as the linear combination of wcsp's given as input
 * makeWeightedCSPSolver: possibility to give a wcsp as input, which will be used as the instance solved by the solver (works only with the base Solver class, otherwise, a new weightedCSP is created withinin the solver object)
 * adding a function for the python API for weightedCSP to read an instance independently from the solver
 * adding the multiwcsp object and methods in the python API
 * modification of the python API in makeWeightedCSPSolver to take a weightedCSP as input
 * modification of the python API to create a weightedCSP object without solvers
 * modification of the python API to read an instance file directly from a weightedCSP object
 */
class MultiCFN {

public:
    // type representing a solution of a multicfn
    typedef std::map<std::string, std::string> Solution;

// ILP encoding type: direct or tuple
#ifdef ILOGCPLEX
    enum ILP_encoding { ILP_Direct,
        ILP_Tuple };
#endif

public:
    /*!
     * \brief default constructor
     */
    MultiCFN();

    /*!
     * \brief constructor: build a multicfn combining wcsps given as input
     * \param wcsps a vector of wcsp objects to combine with weighted sums
     * \param weights a list of weights for each wcsp
     */
    MultiCFN(std::vector<WCSP*>& wcsps, std::vector<Double>& weights);

    /*!
     * \brief destructor (delete all constraints)
     */
    ~MultiCFN();

    /*!
     * \brief add a wcsp to the network, create the variables if they do not exist, the wcsp is stored internally, the original wcsp will not be referenced
     * \param wcsp the wcsp to add
     * \param weight the weight of the wcsp in the objective function (sum of the cost functions)
     */
    void push_back(WCSP* wcsp, Double weight = 1.0);

    /*!
     * \brief set the weight of the cost functions of one of the network
     * \param wcsp_index the index of the network to modify
     * \param weight the new weight
     */
    void setWeight(unsigned int wcsp_index, Double weight);

    /*!
     * \brief get the wieght of a network
     * \param wcsp_index the index of the network (by adding order)
     * \return the weight assigned to the network
     */
    Double getWeight(unsigned int wcsp_index);

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
     * \brief number of cost functions in the problem
     * \return the number of cost functions
     */
    unsigned int nbCostFunctions();

    /*!
     * \brief get the name of one of the network added to the multiwcsp
     * \param index the index of the network
     * \return the name associated to the network
     */
    std::string getNetworkName(unsigned int index);

    /*!
     * \brief get the index of one of the variable added to the multiwcsp
     * \param name the name of the variable
     * \return the index associated to the variable or -1 if not found
     */
    int getVariableIndex(std::string name);

    /*!
     * \brief number of values in the domain of the variable
     * \param index the index of the variable (must be lower than nbVariables())
     * \return the number of values in the domain of the variable
     */
    unsigned int nbValues(unsigned int index);

    /*!
     * \brief scope of the cost function
     * \param index the index of the cost function (must be lower than nbCostFunctions())
     * \return the list of variable indexes in the scope of the cost function
     */
    std::vector<unsigned int> getScope(unsigned int index);

    /*!
     * \brief type of the cost function
     * \param index the index of the cost function (must be lower than nbCostFunctions())
     * \return the type of the cost function
     */
    mcriteria::CostFunction::Type getType(unsigned int index);

    /*!
     * \brief return the cost of a given tuple
     * \param index the index of the cost function (must be lower than nbCostFunctions())
     * \param tuple the tuple
     */
    Double getCost(unsigned int index, std::vector<unsigned int>& tuple);

    /*!
     * \brief return the precision used in the combined wcsp (max of the decimalPoint of the wcsp given as input)
     */
    unsigned int getDecimalPoint();

    /*!
     * \brief get the wcsp's unit cost precision as double
     */
    Double getUnitCost();

    /*!
     * \brief print the cfn
     * \brief os the stream to print to
     */
    void print(std::ostream& os);

    /*!
     * \brief make a wcsp from the convex combination of all the wcsps
     * \param vars the optional set of variable indexes to extract the induced graph (if missing or empty then no restriction)
     * \param scopes the optional list of allowed scopes to extract the partial graph (if missing or empty then no restriction)
     * \param constrs the optional list of allowed cost function indexes (same index as in cfn dump file) to extract the partial graph (if missing or empty then no restriction)
     */
    WeightedCSP* makeWeightedCSP(const set<unsigned int>& vars = {}, const vector<set<unsigned int>>& scopes = {}, const vector<unsigned int>& constrs = {});

    /*!
     * \brief fill a wcsp with the convex combination of all the wcsps already added
     * \param wcsp the weighted csp to be filled
     * \param vars the optional set of variable indexes to extract the induced graph (if missing or empty then no restriction)
     * \param scopes the optional list of allowed scopes to extract the partial graph (if missing or empty then no restriction)
     * \param constrs the optional list of allowed cost function indexes (same index as in cfn dump file) to extract the partial graph (if missing or empty then no restriction)
     */
    void makeWeightedCSP(WeightedCSP* wcsp, const set<unsigned int>& vars = {}, const vector<set<unsigned int>>& scopes = {}, const vector<unsigned int>& constrs = {});

#ifdef ILOGCPLEX

    /*!
     * \brief export the multicfn to a cplex model data structure
     * \param model the cplex model data structure
     * \param objectives list of networks to combine in the optimized criterion
     * \param constraints list of global constraints in the model: network index with two bounds, one bound may be infinity
     * \param domain_vars the list of cplex variables corresponding the multicfn variables
     */
    void makeIloModel(IloEnv& env, IloModel& model, ILP_encoding encoding, std::vector<size_t>& objectives, std::vector<std::pair<size_t, std::pair<Double, Double>>>& constraints, std::vector<IloNumVarArray>& domain_vars, std::vector<std::shared_ptr<IloNumVarArray>>& tuple_vars);

    /*!
     * \brief extract the cplex solution
     * \param cplex the cplex solver object
     * \param domain_vars the cplex variables
     * \param solution the MultiCFN solution
     */
    void getCplexSolution(IloCplex& cplex, std::vector<IloNumVarArray>& domain_vars, MultiCFN::Solution& solution);

    // debug
    Double computeCriteriaSol(IloCplex& cplex, size_t index, bool weighted, std::vector<IloNumVarArray>& domain_vars, std::vector<std::shared_ptr<IloNumVarArray>>& tuple_vars);

#endif

    /*!
     * \brief convert a tuple to a cost index, rightmost value indexed first
     * \param variables the list of variables from the tuple
     * \param tuple the tuple: value indexes for each variable
     * \return the index corresponding to the tuple
     */
    unsigned int tupleToIndex(const std::vector<mcriteria::Var*>& variables, const std::vector<unsigned int>& tuple);

    /*!
     * \brief get the solution of the created wcsp after being solved
     * \return the solution as a dictionary of variable names/value names
     * \pre the wcsp must have been solved and not been deleted
     */
    Solution getSolution();

    /*!
     * \brief get the objective values of the different cost function networks from the created wcsp after being solved
     * \return the objective values as a dictionary of variable names/value names
     * \pre the wcsp must have been solved and not been deleted
     */
    std::vector<Double> getSolutionValues();

    /*!
     * \brief compute the values of an existing solution
     * \param solution the solution given
     * \return the costs of the solution
     */
    std::vector<Double> computeSolutionValues(Solution& solution);

    void outputNetSolutionCosts(size_t index, Solution& solution);

    /*!
     * \brief convert a solution returned by ToulBar2 to a dictionary with variable names and values as labels
     * \param solution the solution given bu ToulBar2
     * \return the solution as a dictionary
     */
    Solution convertToSolution(std::vector<Value>& solution);

    /*!
     * \brief activate the noise parameter to add uniform random noise to each cost in every cost function
     * \param activation true activates the noise
     */
    void setNoiseActivation(bool activation);

    /*!
     * \brief set the min and max value for uniform noise perturbation of the costs
     */
    void setNoiseLevel(Double min_level, Double max_level);

private: /* private methods */
    /*!
     * \brief send the cfn to toulbar2
     * \param wcsp tb2 wcsp
     * \param vars the optional set of variable indexes to extract the induced graph (if missing or empty then no restriction)
     * \param scopes the optional list of allowed scopes to extract the partial graph (if missing or empty then no restriction)
     * \param constrs the optional list of allowed cost function indexes (same index as in cfn dump file) to extract the partial graph (if missing or empty then no restriction)
     * \warning export the problem lower bound only if no restrictions are given (or if scopes contains the emptyset)
     */
    void exportToWCSP(WCSP* wcsp, const set<unsigned int>& vars, const vector<set<unsigned int>>& scopes, const vector<unsigned int>& constrs);
    void exportToWCSP_(WCSP* wcsp, const set<unsigned int>& vars, const set<set<unsigned int>>& scopes, const set<unsigned int>& constrs);

    /*!
     * \brief export a tuple cost function to the wcsp
     */
    void exportTupleCostFunction(WCSP* wcsp, unsigned int func_ind, Double top, std::uniform_real_distribution<>& dis);

    /*!
     * \brief export a linear cost function to the wcsp
     */
    void exportLinearCostFunction(WCSP* wcsp, unsigned int func_ind);

    /*!
     * \brief  extract the solution and the objective values from the created wcsp
     */
    void extractSolution();

    /*!
     * \brief add a cost function to the network
     * \param wcsp the tb2 wcsp
     * \param cstr the original tb2 cost function
     */
    void addCostFunction(WCSP* wcsp, Constraint* cstr);

    /*!
     * \brief add a cost function to the network
     * \param wcsp the tb2 wcsp
     * \param cstr the original tb2 cost function
     */
    void addTupleCostFunction(WCSP* wcsp, Constraint* cstr);

    /*!
     * \brief compute a TOP (infinity) value and a minimum cost for the internal representation of the cfn costs
     * \return a pair containing the top value (first) and the min cost (second)
     */
    std::pair<Double, Double> computeTopMinCost();

    /*!
     * \brief check if a variable in a wcsp and a variable in the multicfn with the same name have the same domains, throw an exception otherwise
     * \param tb2_var the variable in the wcsp
     * \param multicfn_var the variable in the multicfn
     */
    void checkVariablesConsistency(EnumeratedVariable* tb2_var, mcriteria::Var& multicfn_var);

    /*!
     * \brief check if a linear constraint is verifed for a solution
     * \param sol the solution to check
     * \return true if the solution verifies the constraint
     */
    bool checkLinCostFuncConsistency(unsigned int func_ind, Solution& sol);

#ifdef ILOGCPLEX

    /*!
     * \brief fill ilog expression with terms representing the objective function of a network
     * \param expr the ilog expression to fill
     * \param index the index of the network to add in the expression
     * \param weighted true if the expression should be multiplied by the network weight, false otherwise
     * \parma negShift all the costs of every cost function to obtain only positive costs, and add the corresponding constant
     * \param domain_vars the ilog variables representing variable domains
     * \param tuple_vars the ilog variables representing cost function tuples
     */
    void addCriterion(IloExpr& expr, size_t index, bool weighted, bool negShift, std::vector<IloNumVarArray>& domain_vars, std::vector<std::shared_ptr<IloNumVarArray>>& tuple_vars);

#endif

public: // public attributes
    // variables
    std::vector<mcriteria::Var> var; // variables
    std::map<std::string, int> var_index; // index of variables

    // cost functions
    std::vector<mcriteria::CostFunction*> cost_function; // list of the cost functions
    std::map<std::string, unsigned int> cost_function_index; // map between cfn names and indices

private: // private attributes
    std::vector<Double> weights; // list of weights for all the loaded networks
    std::vector<std::string> network_names; // names of the networks
    std::vector<std::vector<unsigned int>> networks; // list of the cost function networks (function indexes for each network)
    std::vector<unsigned int> network_index; // index of the network for each cost function

    std::vector<Double> _doriginal_ubs; // list of original upper bounds (as Double) for each network
    std::vector<Double> _doriginal_lbs; // list of original lower bounds (as Double) for each network
    std::vector<Double> _original_costMultipliers; // list of cost multipliers of all the original wcsp

    unsigned int _tb2_decimalpoint; // precision of the wcsp
    Double _tb2_unit_cost; // toulbar2 cost precision as Double

    // parameter that add noise to the costs during exportation if activated
    bool add_noise;
    Double noise_min, noise_max;
    std::random_device rd;
    std::mt19937 gen;

    /* solution */
    WCSP* _wcsp; // pointer to the wcsp containing the combination of the input wcsps
    bool _sol_extraction; // indicates if the solution and objective value have already been extracted
    std::vector<Double> _obj_values;
    Solution _solution;
};

#endif // MULTI_CFN_HPP
