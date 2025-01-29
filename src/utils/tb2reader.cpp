/*
 * **************** Read wcsp format files **************************
 *
 */

#include "core/tb2wcsp.hpp"
#include "core/tb2enumvar.hpp"
#include "applis/tb2pedigree.hpp"
#include "applis/tb2haplotype.hpp"
#include "applis/tb2bep.hpp"
#include "core/tb2naryconstr.hpp"
#include "tb2randomgen.hpp"
#include "core/tb2globaldecomposable.hpp"
#include "core/tb2clqcover.hpp"
#include "core/tb2knapsack.hpp"
#include "lp-parser.hpp"
#ifdef BOOST
#define BOOST_IOSTREAMS_NO_LIB
#include <boost/version.hpp>
#include <boost/tokenizer.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

#if (BOOST_VERSION < 106500)
#define NO_LZMA
#endif

#ifndef NO_LZMA
#include <boost/iostreams/filter/lzma.hpp>
#include "lzma-cpp.inc"
#endif

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/compressed_pair.hpp>
#endif

typedef struct {
    EnumeratedVariable* var;
    vector<Cost> costs;
} TemporaryUnaryConstraint;

/**
 * \defgroup wcspformat Weighted Constraint Satisfaction Problem file format (wcsp)
 *
 * It is a text format composed of a list of numerical and string terms separated by spaces.
 * Instead of using names for making reference to variables, variable
 * indexes are employed. The same for domain values. All indexes start at
 * zero.
 *
 * Cost functions can be defined in intention (see below) or in extension, by their list of
 * tuples. A default cost value is defined per function in order to
 * reduce the size of the list. Only tuples with a different cost value
 * should be given (not mandatory). All the cost values must be positive. The arity of a cost function in extension may be
 * equal to zero. In this case, there is no tuples and the default cost
 * value is added to the cost of any solution. This can be used to represent
 * a global lower bound constant of the problem.
 *
 * The wcsp file format is composed of three parts: a problem header, the list of
 * variable domain sizes, and the list of cost functions.
 *
 * - Header definition for a given problem:
 * \verbatim
 <Problem name>
 <Number of variables (N)>
 <Maximum domain size>
 <Number of cost functions>
 <Initial global upper bound of the problem (UB)>
 \endverbatim
 * The goal is to find an assignment of all the variables with minimum total cost,
 * strictly lower than UB.
 * Tuples with a cost greater than or equal to UB are forbidden (hard constraint).
 *
 * - Definition of domain sizes
 * \verbatim
 <Domain size of variable with index 0>
 ...
 <Domain size of variable with index N - 1>
 \endverbatim
 *
 * Note : domain values range from zero to \e size-1
 *
 * Note : a negative domain size is interpreted as a variable with an interval domain in \f$[0,-size-1]\f$
 *
 * Warning : variables with interval domains are restricted to arithmetic and disjunctive cost functions in intention (see below)
 *
 * - General definition of cost functions
 *   - Definition of a cost function in extension
 * \verbatim
 <Arity of the cost function>
 <Index of the first variable in the scope of the cost function>
 ...
 <Index of the last variable in the scope of the cost function>
 <Default cost value>
 <Number of tuples with a cost different than the default cost>
 \endverbatim
 * followed by for every tuple with a cost different than the default cost:
 * \verbatim
 <Index of the value assigned to the first variable in the scope>
 ...
 <Index of the value assigned to the last variable in the scope>
 <Cost of the tuple>
 \endverbatim
 *
 * Note : Shared cost function: A cost function in extension can be shared by several cost functions with the same arity (and same domain sizes) but different scopes. In order to do that, the cost function to be shared must start by a negative scope size. Each shared cost function implicitly receives an occurrence number starting from 1 and incremented at each new shared definition. New cost functions in extension can reuse some previously defined shared cost functions in extension by using a negative number of tuples representing the occurrence number of the desired shared cost function. Note that default costs should be the same in the shared and new cost functions. Here is an example of 4 variables with domain size 4 and one AllDifferent hard constraint decomposed into 6 binary constraints.
 *
 *   - Shared CF used inside a small example in wcsp format:
 * \code
 AllDifferentDecomposedIntoBinaryConstraints 4 4 6 1
 4 4 4 4
 -2 0 1 0 4
 0 0 1
 1 1 1
 2 2 1
 3 3 1
 2 0 2 0 -1
 2 0 3 0 -1
 2 1 2 0 -1
 2 1 3 0 -1
 2 2 3 0 -1
 \endcode
 *   - Definition of a cost function in intension by replacing the default cost value by -1 and by giving its keyword name and its K parameters
 * \verbatim
 <Arity of the cost function>
 <Index of the first variable in the scope of the cost function>
 ...
 <Index of the last variable in the scope of the cost function>
 -1
 <keyword>
 <parameter1>
 ...
 <parameterK>
 \endverbatim
 *   .
 * .
 * Possible keywords of cost functions defined in intension followed by their specific parameters:
 * - >= \e cst \e delta to express soft binary constraint \f$x \geq y + cst\f$ with associated cost function \f$max( (y + cst - x \leq delta)?(y + cst - x):UB , 0 )\f$
 * - > \e cst \e delta to express soft binary constraint \f$x > y + cst\f$ with associated cost function  \f$max( (y + cst + 1 - x \leq delta)?(y + cst + 1 - x):UB , 0 )\f$
 * - <= \e cst \e delta to express soft binary constraint \f$x \leq y + cst\f$ with associated cost function  \f$max( (x - cst - y \leq delta)?(x - cst - y):UB , 0 )\f$
 * - < \e cst \e delta to express soft binary constraint \f$x < y + cst\f$ with associated cost function  \f$max( (x - cst + 1 - y \leq delta)?(x - cst + 1 - y):UB , 0 )\f$
 * - = \e cst \e delta to express soft binary constraint \f$x = y + cst\f$ with associated cost function  \f$(|y + cst - x| \leq delta)?|y + cst - x|:UB\f$
 * - disj \e cstx \e csty \e penalty to express soft binary disjunctive constraint \f$x \geq y + csty \vee y \geq x + cstx\f$ with associated cost function \f$(x \geq y + csty \vee y \geq x + cstx)?0:penalty\f$
 * - sdisj \e cstx \e csty \e xinfty \e yinfty \e costx \e costy to express a special disjunctive constraint with three implicit hard constraints \f$x \leq xinfty\f$ and \f$y \leq yinfty\f$ and \f$x < xinfty \wedge y < yinfty \Rightarrow (x \geq y + csty \vee y \geq x + cstx)\f$ and an additional cost function \f$((x = xinfty)?costx:0) + ((y= yinfty)?costy:0)\f$
 * - Global cost functions using a dedicated propagator:
 *     - clique \e 1 (\e nb_values (\e value)*)* to express a hard clique cut to restrict the number of variables taking their value into a given set of values (per variable) to at most \e 1 occurrence for all the variables (warning! it assumes also a clique of binary constraints already exists to forbid any two variables using both the restricted values)
 *     - knapsack \e capacity (\e weight)* to express a reverse knapsack constraint (i.e., a linear constraint on 0/1 variables with >= operator) with capacity and weights are positive or negative integer coefficients (use negative numbers to express a linear constraint with <= operator)
 *     - knapsackc \e capacity (\e weight)* \e nb_AMO (\e nb_variables (\e variable \e value)*)* to express a reverse knapsack constraint (i.e., a linear constraint on 0/1 variables with >= operator) combined with a list of non-overlapping at-most-one constraints
 *     - knapsackp \e capacity (\e nb_values (\e value \e weight)*)* to express a reverse knapsack constraint with for each variable the list of values to select the item in the knapsack with their corresponding weight
 *     - knapsackv \e capacity \e nb_triplets (\e variable \e value \e weight)* to express a reverse knapsack constraint with a list of triplets variable, value, and its corresponding weight
 *     - salldiffkp hard \e UB to express a hard alldifferent constraint (decomposes into knapsack cost functions)
 *     - wcsp \e lb \e ub \e duplicatehard \e strongduality \e wcsp to express a hard global constraint on the cost of an input weighted constraint satisfaction problem in wcsp format such that its valid solutions must have a cost value in [lb,ub[.
 *
 * - Global cost functions using a flow-based propagator:
 *     - salldiff var|dec|decbi \e cost to express a soft alldifferent constraint with either variable-based (\e var keyword) or decomposition-based (\e dec and \e decbi keywords) cost semantic with a given \e cost per violation (\e decbi decomposes into a binary cost function complete network)
 *     - sgcc var|dec|wdec \e cost \e nb_values (\e value \e lower_bound \e upper_bound (\e shortage_weight \e excess_weight)?)* to express a soft global cardinality constraint with either variable-based (\e var keyword) or decomposition-based (\e dec keyword) cost semantic with a given \e cost per violation and for each value its lower and upper bound (if \e wdec then violation cost depends on each value shortage or excess weights)
 *     - ssame \e cost \e list_size1 \e list_size2 (\e variable_index)* (\e variable_index)* to express a permutation constraint on two lists of variables of equal size (implicit variable-based cost semantic)
 *     - sregular var|edit \e cost \e nb_states \e nb_initial_states (\e state)* \e nb_final_states (\e state)* \e nb_transitions (\e start_state \e symbol_value \e end_state)* to express a soft regular constraint with either variable-based (\e var keyword) or edit distance-based (\e edit keyword) cost semantic with a given \e cost per violation followed by the definition of a deterministic finite automaton with number of states, list of initial and final states, and list of state transitions where symbols are domain values
 *     .
 * - Global cost functions using a dynamic programming DAG-based propagator:
 *     - sregulardp var \e cost \e nb_states \e nb_initial_states (\e state)* \e nb_final_states (\e state)* \e nb_transitions (\e start_state \e symbol_value \e end_state)* to express a soft regular constraint with a variable-based (\e var keyword) cost semantic with a given \e cost per violation followed by the definition of a deterministic finite automaton with number of states, list of initial and final states, and list of state transitions where symbols are domain values
 *     - sgrammar|sgrammardp var|weight \e cost \e nb_symbols \e nb_values \e start_symbol \e nb_rules ((0 \e terminal_symbol \e value)|(1 \e nonterminal_in \e nonterminal_out_left \e nonterminal_out_right)|(2 \e terminal_symbol \e value \e weight)|(3 \e nonterminal_in \e nonterminal_out_left \e nonterminal_out_right \e weight))* to express a soft/weighted grammar in Chomsky normal form
 *     - samong|samongdp var \e cost \e lower_bound \e upper_bound \e nb_values (\e value)* to express a soft among constraint to restrict the number of variables taking their value into a given set of values
 *     - salldiffdp var \e cost to express a soft alldifferent constraint with variable-based (\e var keyword) cost semantic with a given \e cost per violation (decomposes into samongdp cost functions)
 *     - sgccdp var \e cost \e nb_values (\e value \e lower_bound \e upper_bound)* to express a soft global cardinality constraint with variable-based (\e var keyword) cost semantic with a given \e cost per violation and for each value its lower and upper bound (decomposes into samongdp cost functions)
 *     - max|smaxdp \e defCost \e nbtuples (\e variable \e value \e cost)* to express a weighted max cost function to find the maximum cost over a set of unary cost functions associated to a set of variables (by default, \e defCost if unspecified)
 *     - MST|smstdp to express a spanning tree hard constraint where each variable is assigned to its parent variable index in order to build a spanning tree (the root being assigned to itself)
 *     .
 * - Global cost functions using a cost function network-based propagator:
 *     - wregular \e nb_states \e nb_initial_states (\e state and cost)* \e nb_final_states (\e state and cost)* \e nb_transitions (\e start_state \e symbol_value \e end_state \e cost)* to express a weighted regular constraint with weights on initial states, final states, and transitions, followed by the definition of a deterministic finite automaton with number of states, list of initial and final states with their costs, and list of weighted state transitions where symbols are domain values
 *     - walldiff hard|lin|quad \e cost to express a soft alldifferent constraint as a set of wamong hard constraint (\e hard keyword) or decomposition-based (\e lin and \e quad keywords) cost semantic with a given \e cost per violation
 *     - wgcc hard|lin|quad \e cost \e nb_values (\e value \e lower_bound \e upper_bound)* to express a soft global cardinality constraint as either a hard constraint (\e hard keyword) or with decomposition-based (\e lin and \e quad keyword) cost semantic with a given \e cost per violation and for each value its lower and upper bound
 *     - wsame hard|lin|quad \e cost to express a permutation constraint on two lists of variables of equal size (implicitly concatenated in the scope) using implicit decomposition-based cost semantic
 *     - wsamegcc hard|lin|quad \e cost \e nb_values (\e value \e lower_bound \e upper_bound)* to express the combination of a soft global cardinality constraint and a permutation constraint
 *     - wamong hard|lin|quad \e cost \e nb_values (\e value)* \e lower_bound \e upper_bound to express a soft among constraint to restrict the number of variables taking their value into a given set of values
 *     - wvaramong hard \e cost \e nb_values (\e value)* to express a hard among constraint to restrict the number of variables taking their value into a given set of values to be equal to the last variable in the scope
 *     - woverlap hard|lin|quad \e cost \e comparator \e righthandside  overlaps between two sequences of variables X, Y (i.e. set the fact that Xi and Yi take the same value (not equal to zero))
 *     - wsum hard|lin|quad \e cost \e comparator \e righthandside to express a soft sum constraint with unit coefficients to test if the sum of a set of variables matches with a given comparator and right-hand-side value
 *     - wvarsum hard \e cost \e comparator to express a hard sum constraint to restrict the sum to be \e comparator to the value of the last variable in the scope
 *     - wdiverse \e distance (\e value)* to express a hard diversity constraint using a dual encoding such that there is a given minimum Hamming distance to a given variable assignment
 *     - whdiverse \e distance (\e value)* to express a hard diversity constraint using a hidden encoding such that there is a given minimum Hamming distance to a given variable assignment
 *     - wtdiverse \e distance (\e value)* to express a hard diversity constraint using a ternary encoding such that there is a given minimum Hamming distance to a given variable assignment
 *
 *       Let us note <> the comparator, K the right-hand-side value associated to the comparator, and Sum the result of the sum over the variables. For each comparator, the gap is defined according to the distance as follows:
 *       -	if <> is == : gap = abs(K - Sum)
 *       -  if <> is <= : gap = max(0,Sum - K)
 *       -  if <> is < : gap = max(0,Sum - K - 1)
 *       -	if <> is != : gap = 1 if Sum != K and gap = 0 otherwise
 *       -  if <> is > : gap = max(0,K - Sum + 1);
 *       -	if <> is >= : gap = max(0,K - Sum);
 *       .
 *     .
 * .
 *
 * Warning : The decomposition of wsum and wvarsum may use an exponential size (sum of domain sizes).
 *
 * Warning :  \e list_size1 and \e list_size2 must be equal in \e ssame.
 *
 * Warning :  Cost functions defined in intention cannot be shared.
 *
 * Note More about network-based global cost functions can be found on ./misc/doc/DecomposableGlobalCostFunctions.html
 *
 * Examples:
 * - quadratic cost function \f$x0 * x1\f$ in extension with variable domains \f$\{0,1\}\f$ (equivalent to a soft clause \f$\neg x0 \vee \neg x1\f$): \code 2 0 1 0 1 1 1 1 \endcode
 * - simple arithmetic hard constraint \f$x1 < x2\f$: \code 2 1 2 -1 < 0 0 \endcode
 * - hard temporal disjunction\f$x1 \geq x2 + 2 \vee x2 \geq x1 + 1\f$: \code 2 1 2 -1 disj 1 2 UB \endcode
 * - clique cut ({x0,x1,x2,x3}) on Boolean variables such that value 1 is used at most once: \code 4 0 1 2 3 -1 clique 1 1 1 1 1 1 1 1 1 \endcode
 * - knapsack constraint (\f$2 * x0 + 3 * x1 + 4 * x2 + 5 * x3 >= 10\f$) on four Boolean 0/1 variables: \code 4 0 1 2 3 -1 knapsack 10 2 3 4 5 \endcode
 * - knapsackc constraint (\f$2 * x0 + 3 * x1 + 4 * x2 + 5 * x3 >= 10\f$, \f$x1 + x2 <= 1\f$) on four Boolean 0/1 variables: \code 4 0 1 2 3 -1 knapsackc 10 2 3 4 5 1 2 1 1 2 1\endcode
 * - knapsackp constraint (\f$2 * (x0=0) + 3 * (x1=1) + 4 * (x2=2) + 5 * (x3=0 \vee x3=1) >= 10\f$) on four {0,1,2}-domain variables: \code 4 0 1 2 3 -1 knapsackp 10 1 0 2 1 1 3 1 2 4 2 0 5 1 5\endcode
 * - knapsackv constraint (\f$2 * (x0=0) + 3 * (x1=1) + 4 * (x2=2) + 5 * (x3=0 \vee x3=1) >= 10\f$) on four {0,1,2}-domain variables: \code 4 0 1 2 3 -1 knapsackv 10 5 0 0 2 1 1 3 2 2 4 3 0 5 3 1 5\endcode
 * - wcsp constraint (\f$3 <= 2 * x1 * x2 + 3 * x1 * x4 + 4 * x2 * x4 < 5\f$) on three Boolean 0/1 variables: \code 3 1 2 4 -1 wcsp 3 5 0 0 name 3 2 3 1000 2 2 2 2 0 1 0 1 1 1 2 2 0 2 0 1 1 1 3 2 1 2 0 1 1 1 4\endcode
 * - soft_alldifferent({x0,x1,x2,x3}): \code 4 0 1 2 3 -1 salldiff var 1 \endcode
 * - soft_gcc({x1,x2,x3,x4}) with each value \e v from 1 to 4 only appearing at least v-1 and at most v+1 times: \code 4 1 2 3 4 -1 sgcc var 1 4 1 0 2 2 1 3 3 2 4 4 3 5 \endcode
 * - soft_same({x0,x1,x2,x3},{x4,x5,x6,x7}): \code 8 0 1 2 3 4 5 6 7 -1 ssame 1 4 4 0 1 2 3 4 5 6 7 \endcode
 * - soft_regular({x1,x2,x3,x4}) with DFA (3*)+(4*): \code 4 1 2 3 4 -1 sregular var 1 2 1 0 2 0 1 3 0 3 0 0 4 1 1 4 1 \endcode
 * - soft_grammar({x0,x1,x2,x3}) with hard cost (1000) producing well-formed parenthesis expressions: \code 4 0 1 2 3 -1 sgrammardp var 1000 4 2 0 6 1 0 0 0 1 0 1 2 1 0 1 3 1 2 0 3 0 1 0 0 3 1 \endcode
 * - soft_among({x1,x2,x3,x4}) with hard cost (1000) if \f$\sum_{i=1}^4(x_i \in \{1,2\}) < 1\f$ or \f$\sum_{i=1}^4(x_i \in \{1,2\}) > 3\f$: \code 4 1 2 3 4 -1 samongdp var 1000 1 3 2 1 2 \endcode
 * - soft max({x0,x1,x2,x3}) with cost equal to \f$\max_{i=0}^3((x_i!=i)?1000:(4-i))\f$: \code 4 0 1 2 3 -1 smaxdp 1000 4 0 0 4 1 1 3 2 2 2 3 3 1 \endcode
 * - wregular({x0,x1,x2,x3}) with DFA (0(10)*2*): \code 4 0 1 2 3 -1 wregular 3 1 0 0 1 2 0 9 0 0 1 0 0 1 1 1 0 2 1 1 1 1 0 0 1 0 0 1 1 2 0 1 1 2 2 0 1 0 2 1 1 1 2 1 \endcode
 * - wamong({x1,x2,x3,x4}) with hard cost (1000) if \f$\sum_{i=1}^4(x_i \in \{1,2\}) < 1\f$ or \f$\sum_{i=1}^4(x_i \in \{1,2\}) > 3\f$: \code 4 1 2 3 4 -1 wamong hard 1000 2 1 2 1 3 \endcode
 * - wvaramong({x1,x2,x3,x4}) with hard cost (1000) if \f$\sum_{i=1}^3(x_i \in \{1,2\}) \neq x_4\f$: \code 4 1 2 3 4 -1 wvaramong hard 1000 2 1 2 \endcode
 * - woverlap({x1,x2,x3,x4}) with hard cost (1000) if \f$\sum_{i=1}^2(x_i = x_{i+2}) \geq 1\f$: \code 4 1 2 3 4 -1 woverlap hard 1000 < 1\endcode
 * - wsum({x1,x2,x3,x4}) with hard cost (1000) if \f$\sum_{i=1}^4(x_i) \neq 4\f$: \code 4 1 2 3 4 -1 wsum hard 1000 == 4 \endcode
 * - wvarsum({x1,x2,x3,x4}) with hard cost (1000) if \f$\sum_{i=1}^3(x_i) \neq x_4\f$: \code 4 1 2 3 4 -1 wvarsum hard 1000 == \endcode
 * - wdiverse({x0,x1,x2,x3}) hard constraint on four variables with minimum Hamming distance of 2 to the value assignment (1,1,0,0): \code 4 0 1 2 3 -1 wdiverse 2 1 1 0 0 \endcode
 * .
 *
 * Latin Square 4 x 4 crisp CSP example in wcsp format:
 * \code
 latin4 16 4 8 1
 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
 4 0 1 2 3 -1 salldiff var 1
 4 4 5 6 7 -1 salldiff var 1
 4 8 9 10 11 -1 salldiff var 1
 4 12 13 14 15 -1 salldiff var 1
 4 0 4 8 12 -1 salldiff var 1
 4 1 5 9 13 -1 salldiff var 1
 4 2 6 10 14 -1 salldiff var 1
 4 3 7 11 15 -1 salldiff var 1
 \endcode
 *
 * 4-queens binary weighted CSP example with random unary costs in wcsp format:
 * \code
 4-WQUEENS 4 4 10 5
 4 4 4 4
 2 0 1 0 10
 0 0 5
 0 1 5
 1 0 5
 1 1 5
 1 2 5
 2 1 5
 2 2 5
 2 3 5
 3 2 5
 3 3 5
 2 0 2 0 8
 0 0 5
 0 2 5
 1 1 5
 1 3 5
 2 0 5
 2 2 5
 3 1 5
 3 3 5
 2 0 3 0 6
 0 0 5
 0 3 5
 1 1 5
 2 2 5
 3 0 5
 3 3 5
 2 1 2 0 10
 0 0 5
 0 1 5
 1 0 5
 1 1 5
 1 2 5
 2 1 5
 2 2 5
 2 3 5
 3 2 5
 3 3 5
 2 1 3 0 8
 0 0 5
 0 2 5
 1 1 5
 1 3 5
 2 0 5
 2 2 5
 3 1 5
 3 3 5
 2 2 3 0 10
 0 0 5
 0 1 5
 1 0 5
 1 1 5
 1 2 5
 2 1 5
 2 2 5
 2 3 5
 3 2 5
 3 3 5
 1 0 0 2
 1 1
 3 1
 1 1 0 2
 1 1
 2 1
 1 2 0 2
 1 1
 2 1
 1 3 0 2
 0 1
 2 1
 \endcode
 **/

#ifdef BOOST
class CFNStreamReader {

public:
    CFNStreamReader(istream& stream, WCSP* wcsp);
    ~CFNStreamReader();

    std::pair<int, string> getNextToken();
    void skipOBrace(); // checks if next token is an opening brace and spits an error otherwise.
    void skipCBrace(); // checks if next token is a  closing brace and spits an error otherwise.
    void testJSONTag(const std::pair<int, string>& token, const string& tag);
    void skipJSONTag(const string& tag);
    void testAndSkipFirstOBrace();
    bool isCost(const string& str);

    // WCSP2 reading methods
    Cost readHeader();
    pair<unsigned, unsigned> readVariables();
    unsigned readVariable(unsigned int varIndex);
    int readDomain(std::vector<string>& valueNames);
    int getValueIdx(int variableIdx, const string& token, int lineNumber);
    void readScope(vector<int>& scope);
    pair<unsigned, unsigned> readCostFunctions();
    void readZeroAryCostFunction(bool all, Cost defaultCost);
    void readNaryCostFunction(vector<int>& scope, bool all, Cost defaultCost);
    void readArithmeticCostFunction();
    void readGlobalCostFunction(vector<int>& scope, const std::string& globalCfnName, int line);

    void generateGCFStreamFromTemplate(vector<int>& scope, const string& funcName, string GCFTemplate, stringstream& stream);

    void generateGCFStreamSgrammar(vector<int>& scope, stringstream& stream);
    void generateGCFStreamSsame(vector<int>& scope, stringstream& stream);

    void readIntervalUnaryTable(int varIdx, vector<Value>& authorized);
    std::vector<Cost> readFunctionCostTable(vector<int> scope, bool all, Cost defaultCost, Cost& minCost);
    void enforceUB();

    std::map<std::string, int> varNameToIdx;
    std::vector<std::map<std::string, int>> varValNameToIdx;
    std::map<std::string, std::vector<pair<string, std::vector<int>>>> tableShares;
    vector<TemporaryUnaryConstraint> unaryCFs;

private:
    istream& iStream;
    WCSP* wcsp;
    bool getNextLine();
    unsigned int lineCount;
    string currentLine;
    boost::char_separator<char> sep;
    bool tok_init;
    boost::tokenizer<boost::char_separator<char>> tok;
    boost::tokenizer<boost::char_separator<char>>::iterator tok_iter;
    bool JSONMode;
    Cost upperBound;
};

CFNStreamReader::CFNStreamReader(istream& stream, WCSP* wcsp)
    : iStream(stream)
    , wcsp(wcsp)
    , tok(string())
{
    this->lineCount = 0;
    this->JSONMode = false;
    this->tok_init = false;
    this->sep = boost::char_separator<char>(" \n\f\r\t\":,", "{}[]");
    this->upperBound = readHeader();
    if (ToulBar2::costThresholdS.size())
        ToulBar2::costThreshold = wcsp->decimalToCost(ToulBar2::costThresholdS, 0);
    if (ToulBar2::costThresholdPreS.size())
        ToulBar2::costThresholdPre = wcsp->decimalToCost(ToulBar2::costThresholdPreS, 0);
    unsigned nvar, nval;
    tie(nvar, nval) = readVariables();
    unsigned ncf, maxarity;
    tie(ncf, maxarity) = readCostFunctions();

    // all negCosts are collected. We should be fine enforcing the UB
    enforceUB();
    if (ToulBar2::vnsOptimumS.size())
        ToulBar2::vnsOptimum = wcsp->decimalToCost(ToulBar2::vnsOptimumS, 0) + wcsp->getNegativeLb();

    // merge unary cost functions if they are on the same variable
    vector<int> seen(wcsp->numberOfVariables(), -1);
    vector<TemporaryUnaryConstraint> newunaryCFs;
    for (unsigned int u = 0; u < unaryCFs.size(); u++) {
        if (seen[unaryCFs[u].var->wcspIndex] == -1) {
            seen[unaryCFs[u].var->wcspIndex] = newunaryCFs.size();
            newunaryCFs.push_back(unaryCFs[u]);
        } else {
            for (unsigned int i = 0; i < unaryCFs[u].var->getDomainInitSize(); i++) {
                if (newunaryCFs[seen[unaryCFs[u].var->wcspIndex]].costs[i] < wcsp->getUb()) {
                    if (unaryCFs[u].costs[i] < wcsp->getUb())
                        newunaryCFs[seen[unaryCFs[u].var->wcspIndex]].costs[i] += unaryCFs[u].costs[i];
                    else
                        newunaryCFs[seen[unaryCFs[u].var->wcspIndex]].costs[i] = wcsp->getUb();
                }
            }
        }
    }

    unaryCFs = newunaryCFs;

    if (ToulBar2::sortDomains) {
        cout << "Error: cannot sort domains in preprocessing with CFN format (remove option -sortd)." << endl;
        throw BadConfiguration();
        //        if (maxarity > 2) {
        //            cout << "Error: cannot sort domains in preprocessing with non-binary cost functions." << endl;
        //            throw BadConfiguration();
        //        } else {
        //            ToulBar2::sortedDomains.clear();
        //            for (unsigned int u = 0; u < unaryCFs.size(); u++) {
        //                ToulBar2::sortedDomains[unaryCFs[u].var->wcspIndex] = unaryCFs[u].var->sortDomain(unaryCFs[u].costs);
        //            }
        //        }
    }

    for (auto& cf : unaryCFs) {

        bool bilevel2 = false;
        bool bilevel3 = false;
        Cost negcost23 = MIN_COST;
        Cost initlb23 = MIN_COST;
        if (ToulBar2::bilevel) {
            if (cf.var->wcspIndex < (int)wcsp->varsBLP[0].size()) { // cost function on leader variable(s) only
                if (ToulBar2::bilevel == 2) { // skip when reading Problem2
                    bilevel2 = true;
                    ToulBar2::bilevel = 5;
                    negcost23 = wcsp->getNegativeLb();
                    initlb23 = wcsp->getLb();
                } else if (ToulBar2::bilevel == 3) { // incorporate directly in Problem1
                    bilevel3 = true;
                    ToulBar2::bilevel = 1;
                    negcost23 = wcsp->getNegativeLb();
                    initlb23 = wcsp->getLb();
                }
            }
        }

        wcsp->postUnaryConstraint(cf.var->wcspIndex, cf.costs);

        if (bilevel2) { // we must ignore this constraint
            ToulBar2::bilevel = 2;
            Cost deltaneg = wcsp->getNegativeLb() - negcost23;
            if (deltaneg != MIN_COST) {
                wcsp->decreaseLb(-deltaneg);
            }
            Cost deltalb = wcsp->getLb() - initlb23;
            if (deltalb != MIN_COST) {
                wcsp->setLb(initlb23);
            }
        } else if (bilevel3) { // we move the corresponding negcost and lb to the restricted leader problem
            ToulBar2::bilevel = 3;
            Cost deltaneg = wcsp->getNegativeLb() - negcost23;
            if (deltaneg != MIN_COST) {
                ToulBar2::negCostBLP[0] += deltaneg;
                ToulBar2::initialUbBLP[0] += deltaneg;
                wcsp->decreaseLb(-deltaneg);
            }
            Cost deltalb = wcsp->getLb() - initlb23;
            if (deltalb != MIN_COST) {
                ToulBar2::initialLbBLP[0] += deltalb;
                wcsp->setLb(initlb23);
            }
        }
    }

    wcsp->sortConstraints();

    if (ToulBar2::verbose >= 0)
        cout << "Read " << nvar << " variables, with " << nval << " values at most, and " << ncf << " cost functions, with maximum arity " << maxarity << "." << endl;
}

// Reads a line. Skips comment lines starting with '#' and // too.
bool CFNStreamReader::getNextLine()
{
    string line;
    lineCount++;
    std::getline(iStream, line);

    while (line == "" || line.at(0) == '#') {
        if (std::getline(iStream, line))
            lineCount++;
        else
            return false;
    }
    size_t posJSONcomments = line.find("//");
    if (posJSONcomments == string::npos) {
        this->currentLine = line;
        return true;
    } else if (posJSONcomments != 0) {
        this->currentLine = line.substr(0, posJSONcomments);
        return true;
    } else
        return getNextLine(); // tail recurse
}

// Reads a token using lazily updated line by line reads
std::pair<int, string> CFNStreamReader::getNextToken()
{
    if (tok_init) {
        if (tok_iter != tok.end()) {
            string token = *tok_iter;
            tok_iter = std::next(tok_iter);
            return make_pair(lineCount, token);
        } else {
            tok_init = false;
            return getNextToken();
        }
    } else {
        if (this->getNextLine()) {
            tok = boost::tokenizer<boost::char_separator<char>>(currentLine, sep);
            tok_iter = tok.begin();
            tok_init = true;
            return getNextToken();
        } else {
            return make_pair(-1, "");
        }
    }
}

CFNStreamReader::~CFNStreamReader()
{
    // TODO clear vectors / maps
}

// Utilities: test opening and closing braces.
bool inline isOBrace(const string& token) { return ((token == "{") || (token == "[")); }
bool inline isCBrace(const string& token) { return ((token == "}") || (token == "]")); }

inline void yellOBrace(const string& token, const int& l)
{
    if (!isOBrace(token)) {
        cerr << "Error: expected a '{' or '[' instead of '" << token << "' at line " << l << endl;
        throw WrongFileFormat();
    }
}

inline void yellCBrace(const string& token, const int& l)
{
    if (!isCBrace(token)) {
        cerr << "Error: expected a ']' or ']' instead of '" << token << "' at line " << l << endl;
        throw WrongFileFormat();
    }
}

// checks if the next token is an opening brace
// and yells otherwise
void CFNStreamReader::skipOBrace()
{
    int l;
    string token;

    std::tie(l, token) = this->getNextToken();
    yellOBrace(token, l);
}
// checks if the next token is a closing brace
// and yells otherwise
void CFNStreamReader::skipCBrace()
{
    int l;
    string token;

    std::tie(l, token) = this->getNextToken();
    yellCBrace(token, l);
}

// Tests if a read token is the expected (JSON) tag and yells otherwise.
inline void CFNStreamReader::testJSONTag(const std::pair<int, string>& token, const string& tag)
{
    if (token.second != tag) {
        cerr << "Error: expected '" << tag << "' instead of '" << token.second << "' at line " << token.first << endl;
        throw WrongFileFormat();
    }
}

// In JSON mode, checks is the next token is the expected (JSON) tag and yells otherwise.
inline void CFNStreamReader::skipJSONTag(const string& tag)
{
    if (JSONMode) {
        testJSONTag(this->getNextToken(), tag);
    }
}

// Checks for the first internal opening brace.
// If it is preceded by a "problem" tag, activates JSON tag checking.
void CFNStreamReader::testAndSkipFirstOBrace()
{
    int l;
    string token;

    std::tie(l, token) = this->getNextToken();
    if (token == "problem") {
        JSONMode = true;
        std::tie(l, token) = this->getNextToken();
    }

    if (!isOBrace(token)) {
        cerr << "Error: expected a '{' or '[' instead of '" << token << "' at line " << l << endl;
        throw WrongFileFormat();
    }
}

// Tests if the token starts with a digit, '+, '-' or "." (is a Cost)
bool CFNStreamReader::isCost(const string& str)
{
    // Test if the first char of a string is a decimal digit, point or +/- sign.
    return (string("0123456789-+.").find(str[0]) != string::npos);
}

// Reads the problem header (problem name and global Bound) and returns the bound
// Starts: at the beginning of the stream
// Ends  : after the closing brace of the header
Cost CFNStreamReader::readHeader()
{
    int lineNumber;
    string token;

    skipOBrace();
    testAndSkipFirstOBrace(); // check if we are in JSON mode
    skipJSONTag("name");
    std::tie(lineNumber, token) = this->getNextToken();

    // TODO pas de méthode WCSP pour écrire le nom du problème
    if (ToulBar2::verbose >= 1)
        cout << "Read problem: " << token << endl;

    skipJSONTag("mustbe");
    Cost pbBound;

    std::tie(lineNumber, token) = this->getNextToken();
    if ((token[0] == '<') || token[0] == '>') {

        auto pos = token.find('.');
        string integerPart = token.substr(1, token.find('.'));
        string decimalPart;

        // precision option management :
        if (pos == string::npos) {
            ToulBar2::decimalPoint = 0;
        } else {
            decimalPart = token.substr(token.find('.') + 1);
            if (ToulBar2::verbose >= 0) {
                cout << "Initial cost precision of " << decimalPart.size() << " digits";
                // cout << " (primal bound: " << token << ")";
            }
            if (ToulBar2::resolution_Update) {
                if (ToulBar2::resolution >= 0) {
                    string str;
                    str = decimalPart.substr(0, ToulBar2::resolution);
                    decimalPart = str;
                    if (decimalPart.size() == 0) {
                        pos = string::npos;
                    }
                    if (ToulBar2::verbose >= 0) {
                        cout << " changed to " << decimalPart.size() << " digits";
                        // cout << " (new primal bound: " << integerPart << decimalPart << ")";
                    }
                } else {
                    cerr << "Sorry, cannot use a negative value for precision! (see option -precision)" << endl;
                    throw BadConfiguration();
                }
            }
            if (ToulBar2::verbose >= 0) {
                cout << "." << endl;
            }
            ToulBar2::decimalPoint = decimalPart.size();
        }

        try {
            if (pos != string::npos) {
                pbBound = (std::stoll(integerPart) * (Cost)powl(10, decimalPart.size()));
                pbBound += ((pbBound >= 0) ? std::stoll(decimalPart) : -std::stoll(decimalPart));
            } else {
                pbBound = std::stoll(integerPart);
            }
        } catch (const std::invalid_argument&) {
            cerr << "Error: invalid global bound '" << token << "' at line " << lineNumber << endl;
            throw WrongFileFormat();
        }
    } else {
        cerr << "Error: global bound '" << token << "' misses upper/lower bound comparator at line " << lineNumber << endl;
        throw WrongFileFormat();
    }

    if (token[0] == '>') {
        ToulBar2::setCostMultiplier(ToulBar2::costMultiplier * -1.0);
    }
    if (ToulBar2::bilevel == 3) {
        ToulBar2::setCostMultiplier(ToulBar2::costMultiplier * -1.0);
    }

    if (ToulBar2::verbose >= 1)
        cout << "Read bound: " << pbBound << " with precision " << ToulBar2::decimalPoint << endl;
    skipCBrace();

    long double fbound = pbBound * ToulBar2::costMultiplier;

    if (fbound > MAX_COST) {
        cerr << "Error: bound generates Cost overflow with -C multiplier = " << ToulBar2::costMultiplier << " ( " << this->upperBound << " )" << endl;
        throw WrongFileFormat();
    }
    return fbound;
}

// Reads the variables and domains and creates them.
// Starts: after the opening brace of the variables list.
// Ends:   after the closing brace of the variables list.
pair<unsigned, unsigned> CFNStreamReader::readVariables()
{
    skipJSONTag("variables");
    // Check first opening brace
    skipOBrace();

    unsigned domsize;
    unsigned maxdomsize = 0;
    unsigned nVar = 0;
    while ((domsize = readVariable(nVar)) != 0) {
        maxdomsize = max(maxdomsize, domsize);
        nVar++;
    }
    return make_pair(nVar, maxdomsize);
}

// Reads the description of the ith variable, creates it and returns the domain size (> 0 iff successful)
// Starts: after the opening brace of all variables or closing brace of previous variable.
// Ends:   after the closing brace of all variables or of the variable read otherwise.
unsigned CFNStreamReader::readVariable(unsigned i)
{
    string varName;
    int domainSize = 0;
    vector<string> valueNames;

    string token;
    int lineNumber;
    std::tie(lineNumber, token) = this->getNextToken();

    if (isCBrace(token)) { // End of variable list
        return 0;
    }

    // A domain or domain size is there: the variable has no name
    // we create an integer name that cannot clash with user names
    if (isOBrace(token) || isdigit(token[0])) {
        varName = to_string("x") + to_string(i);
    } else {
        varName = token;
        std::tie(lineNumber, token) = this->getNextToken();
    }
    // This is a list of symbols, read it
    if (isOBrace(token)) {
        domainSize = readDomain(valueNames);
    } else { // Just a domain size
        try {
            domainSize = stoi(token);
            if (domainSize >= 0)
                for (int ii = 0; ii < domainSize; ii++)
                    valueNames.push_back(to_string("v") + to_string(ii));
        } catch (std::invalid_argument&) {
            cerr << "Error: expected domain or domain size instead of '" << token << "' at line " << lineNumber << endl;
        }
    }

    unsigned int varIndex = wcsp->getVarIndex(varName);
    bool newvar = (varIndex >= wcsp->numberOfVariables());
    if (varIndex > wcsp->numberOfVariables()) {
        varName += "_neg";
    }
    if (ToulBar2::verbose >= 1)
        cout << "Variable " << varName << ((newvar) ? " new " : " known ") << "with domain size " << domainSize << " read";
    // Create the toulbar2 variable and store its name in the variable map.
    if (newvar) {
        varIndex = ((domainSize >= 0) ? this->wcsp->makeEnumeratedVariable(varName, 0, domainSize - 1) : this->wcsp->makeIntervalVariable(varName, 0, -domainSize - 1));
    } else if (ToulBar2::bilevel) {
        wcsp->varsBLP[ToulBar2::bilevel - 1].insert(varIndex);
    }
    if (ToulBar2::verbose >= 1)
        cout << " # " << varIndex << endl;
    if (not varNameToIdx.insert(std::pair<string, int>(varName, varIndex)).second) {
        cerr << "Error: variable name '" << varName << "' not unique at line " << lineNumber << endl;
        throw WrongFileFormat();
    }
    // set the value names (if any) in the Variable.values map
    varValNameToIdx.resize(max((size_t)varIndex + 1, varValNameToIdx.size() + 1));
    assert(varValNameToIdx.size() >= varIndex + 1);
    for (unsigned int ii = 0; ii < valueNames.size(); ++ii) {
        if (not varValNameToIdx[varIndex].insert(std::pair<string, int>(valueNames[ii], ii)).second) {
            cerr << "Error: duplicated value name '" << valueNames[ii] << "' for variable '" << wcsp->getName(varIndex) << "' at line " << lineNumber << endl;
            throw WrongFileFormat();
        }
    }
    if (newvar) {
        for (unsigned int ii = 0; ii < valueNames.size(); ++ii)
            wcsp->addValueName(varIndex, valueNames[ii]);
    } else {
        if (((EnumeratedVariable*)wcsp->getVar(varIndex))->getDomainInitSize() != (unsigned int)domainSize) {
            cerr << "Error: same variable has two different domain sizes " << ((EnumeratedVariable*)wcsp->getVar(varIndex))->getDomainInitSize() << ", " << domainSize << " for variable '" << wcsp->getName(varIndex) << "' at line " << lineNumber << endl;
            throw WrongFileFormat();
        }
        for (unsigned int ii = 0; ii < valueNames.size(); ++ii) {
            if (wcsp->getVar(varIndex)->getValueName(ii) != valueNames[ii]) {
                cerr << "Error: same variable has two different domains '" << valueNames[ii] << "' for variable '" << wcsp->getName(varIndex) << "' at line " << lineNumber << endl;
                throw WrongFileFormat();
            }
        }
    }

    return domainSize;
}

// Reads a domain defined as a set of symbolic values in the valueNames vector and returns domain size.
// Starts: after the domain opening brace has been read.
// Ends:   after the domain closing brace has been read.
int CFNStreamReader::readDomain(std::vector<string>& valueNames)
{
    int l;
    string token;
    std::tie(l, token) = this->getNextToken();

    while (!isCBrace(token)) {
        if (isdigit(token[0])) { // not a symbol !
            cerr << "Error: value name '" << token << "' starts with a digit at line " << l << endl;
            throw WrongFileFormat();
        } else {
            valueNames.push_back(token);
        }
        std::tie(l, token) = this->getNextToken();
    }
    return valueNames.size();
}
// Reads list of authorized values for a unary cost function over an interval variable
void CFNStreamReader::readIntervalUnaryTable(int varIdx, vector<Value>& authorized)
{
    int lineNumber;
    string token;

    tie(lineNumber, token) = this->getNextToken();
    while (!isCBrace(token)) {
        unsigned int valueIdx = getValueIdx(varIdx, token, lineNumber);
        authorized.push_back(valueIdx);

        tie(lineNumber, token) = this->getNextToken();
        Cost cost = wcsp->decimalToCost(token, lineNumber);
        if (cost != MIN_COST) {
            cerr << "Error: Unary cost function with non zero non default cost at line " << lineNumber << endl;
            throw WrongFileFormat();
        }
    }
}

// Reads a cost function table for the scope given.
// If all is set, a lexicographically ordered list of costs is expected.
// ELse, a list of tuples with their cost is expectd and all other tuples get the defaultCost
std::vector<Cost> CFNStreamReader::readFunctionCostTable(vector<int> scope, bool all, Cost defaultCost, Cost& minCost)
{
    int lineNumber;
    string token;
    minCost = MAX_COST;

    // Create a vector filled with defaultCost values
    std::vector<Cost> costVector;
    long unsigned int costVecSize = 1;
    for (auto i : scope) {
        costVecSize *= wcsp->getDomainInitSize(i);
    }
    costVector.resize(costVecSize);
    fill(costVector.begin(), costVector.end(), defaultCost);

    std::tie(lineNumber, token) = this->getNextToken();
    if (!all) { // will be a tuple:cost table
        int arity = scope.size();
        int scopeIdx = 0; // position in the scope
        int tableIdx = 0; // position in the cost table
        size_t nbCostInserted = 0;

        while (!isCBrace(token)) {
            // if we have read a full tuple and cost
            if (scopeIdx == arity) {
                Cost cost = wcsp->decimalToCost(token, lineNumber);
                // the same tuple has already been defined.
                if (costVector[tableIdx] != defaultCost) {
                    cerr << "Error: tuple on scope [ ";
                    for (int i : scope)
                        cerr << i << " ";
                    cerr << "] with cost " << cost << " redefined at line " << lineNumber << endl;
                    throw WrongFileFormat();
                } else {
                    costVector[tableIdx] = cost;
                }
                nbCostInserted++;
                minCost = min(cost, minCost);
            } else {
                // striding in the costs array
                if (scopeIdx != 0) {
                    tableIdx *= wcsp->getDomainInitSize(scope[scopeIdx]);
                }
                unsigned int valueIdx = getValueIdx(scope[scopeIdx], token, lineNumber);
                assert(valueIdx >= 0 && valueIdx < wcsp->getDomainInitSize(scope[scopeIdx]));
                tableIdx += valueIdx;
            }
            // If we just finished a tuple reset indices
            if (scopeIdx == arity) {
                scopeIdx = tableIdx = 0;
            } else {
                scopeIdx++;
            }
            std::tie(lineNumber, token) = this->getNextToken();
        }

        if (nbCostInserted < costVecSize) // there are some defaultCost remaining
            minCost = min(defaultCost, minCost);
    }
    // all is true: we expect a full costs list
    else {
        unsigned int tableIdx = 0;
        while (tableIdx < costVecSize && !isCBrace(token)) {
            Cost cost = wcsp->decimalToCost(token, lineNumber);

            minCost = min(cost, minCost);

            costVector[tableIdx] = cost;
            tableIdx++;

            std::tie(lineNumber, token) = this->getNextToken();
        }
        if (tableIdx != costVecSize) {
            cerr << "Error: incorrect number of costs in cost table ending at line " << lineNumber << endl;
            throw WrongFileFormat();
        }
        if (!isCBrace(token)) {
            cerr << "Error: too many costs in cost table ending at line " << lineNumber << endl;
            throw WrongFileFormat();
        }
    }

    // make all costs non negative and remember the shift
    for (Cost& c : costVector) {
        c -= minCost;
    }

    wcsp->negCost -= minCost;
    skipCBrace();

    //    for (size_t i = 0; i < costVector.size(); i++)
    //       cout << i << " " << costVector[i] << endl ;

    return costVector;
}

// this->upperBound is the raw bound from the header (unshifted, rescaled)
void CFNStreamReader::enforceUB()
{
    if (Add(this->upperBound, wcsp->negCost, &this->upperBound)) {
        cerr << "Error: bound generates Cost overflow with -C multiplier = " << ToulBar2::costMultiplier << " ( " << this->upperBound << " " << wcsp->negCost << " )" << endl;
        throw WrongFileFormat();
    }

    // if the shifted/scaled bound is less than zero, we equivalently set it to zero
    if (this->upperBound < MIN_COST)
        this->upperBound = MIN_COST;

    if (ToulBar2::bilevel == 0) {
        if (ToulBar2::externalUB.length() != 0) {
            this->upperBound = min(this->upperBound, wcsp->decimalToCost(ToulBar2::externalUB, 0) + wcsp->negCost);
        }
        if (ToulBar2::deltaUbS.length() != 0) {
            ToulBar2::deltaUbAbsolute = max(MIN_COST, wcsp->decimalToCost(ToulBar2::deltaUbS, 0));
            ToulBar2::deltaUb = max(ToulBar2::deltaUbAbsolute, (Cost)(ToulBar2::deltaUbRelativeGap * (Double)min(this->upperBound, wcsp->getUb())));
            if (ToulBar2::deltaUb > MIN_COST) {
                // as long as a true certificate as not been found we must compensate for the deltaUb in CUT
                this->upperBound += ToulBar2::deltaUb;
            }
        }
        wcsp->updateUb(this->upperBound);
    } else if (ToulBar2::bilevel >= 1 && ToulBar2::bilevel <= 2) {
        wcsp->setUb(this->upperBound);
    }
}

// Returns the index of the value name for the given variable
int CFNStreamReader::getValueIdx(int variableIdx, const string& token, int lineNumber)
{
    if (not isdigit(token[0])) {
        std::map<std::string, int>::iterator it;

        if ((it = varValNameToIdx[variableIdx].find(token)) != varValNameToIdx[variableIdx].end()) {
            return it->second;
        } else {
            cerr << "Error: value name '" << token << "' not in the domain of variable '" << wcsp->getName(variableIdx) << "' at line " << lineNumber << endl;
            throw WrongFileFormat();
        }
    } else {
        int valueIdx = -1;
        try {
            valueIdx = stoi(token);
        } catch (std::invalid_argument&) {
            cerr << "Error: value '" << token << "' is not a proper name/index for variable " << wcsp->getName(variableIdx) << " at line " << lineNumber << endl;
            throw WrongFileFormat();
        }
        if (valueIdx < 0 || (unsigned)valueIdx >= wcsp->getDomainInitSize(variableIdx)) {
            cerr << "Error: value '" << token << "' out of range of variable " << wcsp->getName(variableIdx) << " at line " << lineNumber << endl;
            throw WrongFileFormat();
        }
        return valueIdx;
    }
}

// Reads a scope.
// Starts: after the opening brace of the scope.
// Ends:   after the closing brace of the scope.
void CFNStreamReader::readScope(vector<int>& scope)
{
    int lineNumber;
    string token;

    std::tie(lineNumber, token) = this->getNextToken();
    while (!isCBrace(token)) {
        // It's a name, convert to index
        if (not isdigit(token[0])) {
            map<string, int>::iterator it;
            string tokenNeg = token + "_neg";
            if (ToulBar2::bilevel == 3 && (it = varNameToIdx.find(tokenNeg)) != varNameToIdx.end()) { // change variable names specific to NegProblem2 by adding "_neg" suffix to their names
                scope.push_back(it->second);
            } else if ((it = varNameToIdx.find(token)) != varNameToIdx.end()) {
                scope.push_back(it->second);
            } else {
                cerr << "Error: unknown variable with name '" << token << "' at line " << lineNumber << endl;
                throw WrongFileFormat();
            }
        } else {
            int varIdx = -1;
            try {
                varIdx = stoi(token);
            } catch (std::invalid_argument&) {
                cerr << "Error: not a variable name or index " << varIdx << " at line " << lineNumber << endl;
                throw WrongFileFormat();
            }
            if (varIdx < 0 || (unsigned)varIdx >= wcsp->numberOfVariables()) {
                cerr << "Error: unknown variable index " << varIdx << " at line " << lineNumber << endl;
                throw WrongFileFormat();
            } else {
                scope.push_back(varIdx);
            }
        }
        // prepare for next iteration (will ultimately read final CBrace)
        std::tie(lineNumber, token) = this->getNextToken();
    }
}
// Reads all cost functions.
// Starts: after the opening brace of all cost functions.
// Ends:   at the end of file or after the closing brace of the problem.
pair<unsigned, unsigned> CFNStreamReader::readCostFunctions()
{
    int lineNumber;
    string token;
    unsigned nbcf = 0;
    unsigned maxarity = 0;

    skipJSONTag("functions");
    skipOBrace();

    std::tie(lineNumber, token) = this->getNextToken(); // start the token pump!

    while ((lineNumber != -1) && !isCBrace(token)) {
        //  Read function name (if any) and move after next OBrace
        string funcName;
        if (!isOBrace(token)) {
            funcName = token;
            skipOBrace();
        }
        // Reads a scope
        skipJSONTag("scope");
        skipOBrace();

        //  Read variable names/indices until CBrace
        vector<int> scope;
        readScope(scope);
        maxarity = max(maxarity, static_cast<unsigned>(scope.size()));
        nbcf++;
        bool isUnaryInterval = (scope.size() == 1 && !wcsp->getVar(scope[0])->enumerated()); // Unary over interval variable are special
        bool isShared = false;

        // Test if function is shared
        if (funcName.size() != 0) {
            isShared = (tableShares.find(funcName) != tableShares.end());
        } else { // If no function name, generate one
            funcName = "f(";
            for (const auto& var : scope) {
                funcName += this->wcsp->getVar(var)->getName();
                funcName += ",";
            }
            funcName[funcName.size() - 1] = ')';
        }

        bool bilevel2 = false;
        bool bilevel3 = false;
        Cost negcost23 = MIN_COST;
        Cost initlb23 = MIN_COST;
        if (ToulBar2::bilevel) {
            bool inleader = true;
            for (int idx : scope) {
                if (idx >= (int)wcsp->varsBLP[0].size()) {
                    inleader = false;
                    break;
                }
            }
            if (inleader) { // cost function on leader variable(s) only
                if (ToulBar2::bilevel == 2) { // skip when reading Problem2
                    bilevel2 = true;
                    ToulBar2::bilevel = 5;
                    negcost23 = wcsp->getNegativeLb();
                    initlb23 = wcsp->getLb();
                } else if (ToulBar2::bilevel == 3) { // incorporate directly in Problem1
                    bilevel3 = true;
                    ToulBar2::bilevel = 1;
                    negcost23 = wcsp->getNegativeLb();
                    initlb23 = wcsp->getLb();
                }
            }
        }

        if (ToulBar2::verbose >= 1)
            cout << "Cost function header for " << funcName << " read" << endl;

        //  Test if a defaultCost is there (and tuples will be expected later)
        bool skipDefaultCost = false;
        Cost defaultCost = MIN_COST;
        std::tie(lineNumber, token) = this->getNextToken();

        if (JSONMode) {
            if (token == "defaultcost") {
                // read the defaultCost
                std::tie(lineNumber, token) = this->getNextToken();
            } else {
                skipDefaultCost = true;
            }
        } else {
            skipDefaultCost = !isCost(token);
        }

        if (!skipDefaultCost) { // Set default cost and skip to next token
            defaultCost = wcsp->decimalToCost(token, lineNumber);
            std::tie(lineNumber, token) = this->getNextToken();
        }

        if (isUnaryInterval && defaultCost == MIN_COST) {
            cerr << "Error: unary function " << funcName << " over an interval variable must have non zero default cost at line " << lineNumber << endl;
            throw WrongFileFormat();
        }

        // Discriminate between global/shared and table cost functions
        bool isGlobal = false;
        bool isReused = false;

        if (JSONMode) {
            if (token == "type") { // This is a global/arithmetic
                isGlobal = true;
                // read type
                std::tie(lineNumber, token) = this->getNextToken();
                skipJSONTag("params");
                skipOBrace();
                // ready to read params (after OBrace)
            } else if (token != "costs") {
                cerr << "Error: expected tag 'costs' instead of '" << token << "' at line " << lineNumber << endl;
                throw WrongFileFormat();
            } else { // cost table: can be reused or explicit
                std::tie(lineNumber, token) = this->getNextToken();
                isReused = !isOBrace(token); // no brace, so reused
                if (isReused) {
                    if (!skipDefaultCost) {
                        cerr << "Error: function " << funcName << " sharing cost tables with " << token << " cannot have default costs at line " << lineNumber << endl;
                        throw WrongFileFormat();
                    }
                    tableShares[token].push_back(make_pair(funcName, scope));
                    skipCBrace();
                }
            }
        } else if (!isOBrace(token)) { // Non JSON. No OBrace means reused or global (type)
            string token2;
            int lineNumber2;
            std::tie(lineNumber2, token2) = this->getNextToken();
            if (isOBrace(token2)) { // parameters are there, so global
                isGlobal = true;
            } else if (!isCBrace(token2)) { // reused, nothing expected after the shared fname
                cerr << "Error: expected closing brace after type at line " << lineNumber2 << endl;
            } else {
                if (!skipDefaultCost) {
                    cerr << "Error: function " << funcName << " sharing cost tables with " << token << " cannot have default costs at line " << lineNumber << endl;
                    throw WrongFileFormat();
                }
                isReused = true;
                tableShares[token].push_back(make_pair(funcName, scope));
            }
        }

        // Table cost function
        if (!isGlobal && !isReused) {
            if (scope.size() == 0) {
                this->readZeroAryCostFunction(skipDefaultCost, defaultCost);
            } else if (scope.size() > NARYPROJECTIONSIZE) {
                this->readNaryCostFunction(scope, skipDefaultCost, defaultCost);
            } else {
                Cost minCost;
                vector<Cost> costs;
                vector<Value> authorized; // list of authorized values for Interval Unary

                if (scope.size() != 1 || wcsp->getVar(scope[0])->enumerated())
                    costs = this->readFunctionCostTable(scope, skipDefaultCost, defaultCost, minCost);
                else
                    readIntervalUnaryTable(scope[0], authorized);

                switch (scope.size()) {
                case 1:
                    if (wcsp->getVar(scope[0])->enumerated()) {
                        TemporaryUnaryConstraint unarycf;
                        unarycf.var = (EnumeratedVariable*)wcsp->getVar(scope[0]);
                        assert(costs.size() == unarycf.var->getDomainInitSize());
                        unarycf.costs = costs;
                        unaryCFs.push_back(unarycf);
                        if (isShared) {
                            unsigned int domSize = wcsp->getDomainInitSize(scope[0]);
                            for (const auto& ns : tableShares[funcName]) {
                                if ((ns.second.size() == 1) && wcsp->getVar(ns.second[0])->enumerated() && wcsp->getDomainInitSize(ns.second[0]) == domSize) {
                                    TemporaryUnaryConstraint unarycf;
                                    unarycf.var = (EnumeratedVariable*)wcsp->getVar(ns.second[0]);
                                    assert(costs.size() == unarycf.var->getDomainInitSize());
                                    unarycf.costs = costs;
                                    unaryCFs.push_back(unarycf);
                                    // this->wcsp->postUnaryConstraint(s[0], costs);
                                    wcsp->negCost -= minCost;
                                    // TODO must remember name too
                                } else {
                                    cerr << "Error: cannot share cost function '" << funcName << "' with '" << ns.first << " on scope { ";
                                    for (auto v : ns.second)
                                        cerr << wcsp->getVar(v)->getName() << " ";
                                    cerr << "}" << endl;
                                    throw WrongFileFormat();
                                }
                            }
                        }
                    } else {
                        wcsp->postUnaryConstraint(scope[0], authorized.data(), authorized.size(), defaultCost);
                    }
                    break;
                case 2: {
                    int cfIdx = this->wcsp->postBinaryConstraint(scope[0], scope[1], costs);
                    if (cfIdx != INT_MAX)
                        this->wcsp->getCtr(cfIdx)->setName(funcName);
                    if (isShared) {
                        unsigned int domSize0 = wcsp->getDomainInitSize(scope[0]);
                        unsigned int domSize1 = wcsp->getDomainInitSize(scope[1]);
                        for (const auto& ns : tableShares[funcName]) {
                            if ((ns.second.size() == 2) && wcsp->getDomainInitSize(ns.second[0]) == domSize0 && wcsp->getDomainInitSize(ns.second[1]) == domSize1) {
                                cfIdx = this->wcsp->postBinaryConstraint(ns.second[0], ns.second[1], costs);
                                wcsp->negCost -= minCost;
                                if (cfIdx != INT_MAX)
                                    this->wcsp->getCtr(cfIdx)->setName(ns.first);
                            } else {
                                cerr << "Error: cannot share cost function '" << funcName << "' with '" << ns.first << "' on scope { ";
                                for (auto v : ns.second)
                                    cerr << wcsp->getVar(v)->getName() << " ";
                                cerr << "}" << endl;
                                throw WrongFileFormat();
                            }
                        }
                    }
                } break;
                case 3: {
                    int cfIdx = this->wcsp->postTernaryConstraint(scope[0], scope[1], scope[2], costs);
                    if (cfIdx != INT_MAX)
                        wcsp->getCtr(cfIdx)->setName(funcName);
                    if (isShared) {
                        unsigned int domSize0 = wcsp->getDomainInitSize(scope[0]);
                        unsigned int domSize1 = wcsp->getDomainInitSize(scope[1]);
                        unsigned int domSize2 = wcsp->getDomainInitSize(scope[2]);
                        for (const auto& ns : tableShares[funcName]) {
                            if ((ns.second.size() == 3) && wcsp->getDomainInitSize(ns.second[0]) == domSize0 && wcsp->getDomainInitSize(ns.second[1]) == domSize1 && wcsp->getDomainInitSize(ns.second[2]) == domSize2) {
                                cfIdx = this->wcsp->postTernaryConstraint(ns.second[0], ns.second[1], ns.second[2], costs);
                                wcsp->negCost -= minCost;
                                if (cfIdx != INT_MAX)
                                    wcsp->getCtr(cfIdx)->setName(ns.first);
                            } else {
                                cerr << "Error: cannot share cost function '" << funcName << "' on scope { ";
                                for (auto v : ns.second)
                                    cerr << wcsp->getVar(v)->getName() << " ";
                                cerr << "}" << endl;
                                throw WrongFileFormat();
                            }
                        }
                    }
                } break;
                }
            }
        } else if (isReused) {
            if ((scope.size() <= 1) || (scope.size() > NARYPROJECTIONSIZE) || isGlobal) {
                cerr << "Error: only unary, binary and ternary cost functions can share cost tables for '" << funcName << " at line " << lineNumber << endl;
                throw WrongFileFormat();
            }
        } else if (isGlobal) {
            this->readGlobalCostFunction(scope, token, lineNumber);
        }

        if (bilevel2) { // we must ignore this constraint
            ToulBar2::bilevel = 2;
            Cost deltaneg = wcsp->getNegativeLb() - negcost23;
            if (deltaneg != MIN_COST) {
                wcsp->decreaseLb(-deltaneg);
            }
            Cost deltalb = wcsp->getLb() - initlb23;
            if (deltalb != MIN_COST) {
                wcsp->setLb(initlb23);
            }
        } else if (bilevel3) { // we move the corresponding negcost and lb to the restricted leader problem
            ToulBar2::bilevel = 3;
            Cost deltaneg = wcsp->getNegativeLb() - negcost23;
            if (deltaneg != MIN_COST) {
                ToulBar2::negCostBLP[0] += deltaneg;
                ToulBar2::initialUbBLP[0] += deltaneg;
                wcsp->decreaseLb(-deltaneg);
            }
            Cost deltalb = wcsp->getLb() - initlb23;
            if (deltalb != MIN_COST) {
                ToulBar2::initialLbBLP[0] += deltalb;
                wcsp->setLb(initlb23);
            }
        }

        std::tie(lineNumber, token) = this->getNextToken();
    } // end of while (token != closing braces = EOF)

    return make_pair(nbcf, maxarity);
}

// Reads a 0ary function.
// Starts: after the cost table OBrace
// Ends:   after the function CBrace
void CFNStreamReader::readZeroAryCostFunction(bool all, Cost defaultCost)
{
    string token;
    int lineNumber;

    std::tie(lineNumber, token) = this->getNextToken();
    Cost zeroAryCost = 0;

    if (!isCBrace(token)) { // We have a cost
        zeroAryCost = wcsp->decimalToCost(token, lineNumber);
        skipCBrace();
    } else {
        if (all) { // We should have a cost
            cerr << "Error: no cost or default cost given for 0 arity function at line " << lineNumber << endl;
            throw WrongFileFormat();
        } else
            zeroAryCost = defaultCost;
    }
    if (zeroAryCost < 0) {
        wcsp->negCost -= zeroAryCost;
        zeroAryCost = 0;
    }
    wcsp->increaseLb(zeroAryCost);
    skipCBrace(); // read final function CBrace
}

// Reads a Nary cost function
// Starts:
// Ends  :
void CFNStreamReader::readNaryCostFunction(vector<int>& scope, bool all, Cost defaultCost)
{
    int lineNumber;
    string token;
    Cost minCost = MAX_COST;

    // Compute the cardinality of the cartesian product as unsigned long and floating point (log)
    // the size_t may overflow but all tuples cannot be available in this case
    long double logCard = 0.0;
    size_t card = 1;
    for (auto i : scope) {
        logCard += log(wcsp->getDomainInitSize(i));
        card *= (size_t)wcsp->getDomainInitSize(i);
    }

    unsigned int arity = scope.size();
    Tuple tup(arity);
    map<Tuple, Cost> costFunction;
    size_t nbTuples = 0;
    int scopeArray[arity];
    for (unsigned int i = 0; i < scope.size(); i++) {
        scopeArray[i] = scope[i];
    }

    // Start reading
    std::tie(lineNumber, token) = this->getNextToken();
    if (not all) {
        unsigned int scopeIdx = 0; // Index of the cost table tuple
        while (!isCBrace(token)) {
            // We have read a full tuple: finish the tuple
            if (scopeIdx == arity) {
                Cost cost = wcsp->decimalToCost(token, lineNumber);

                if (not costFunction.insert(pair<Tuple, Cost>(tup, cost)).second) {
                    cerr << "Error: tuple on scope [ ";
                    for (int i : scope)
                        cout << i << " ";
                    cout << "] with cost " << cost << " redefined at line " << lineNumber << endl;
                    throw WrongFileFormat();
                } else {
                    nbTuples++;
                    minCost = min(cost, minCost);
                }
            } else {
                unsigned int valueIdx = getValueIdx(scope[scopeIdx], token, lineNumber);
                assert(valueIdx >= 0 && valueIdx < wcsp->getDomainInitSize(scope[scopeIdx]));
                tup[scopeIdx] = valueIdx; // fill Tuple
            }
            scopeIdx = ((scopeIdx == arity) ? 0 : scopeIdx + 1);
            std::tie(lineNumber, token) = this->getNextToken();
        }
        // Is there any remaining default cost (either too many tuples or less than we need)
        if ((logCard > log(std::numeric_limits<size_t>::max())) || nbTuples < card) {
            minCost = min(minCost, defaultCost);
        }

        int naryIndex = this->wcsp->postNaryConstraintBegin(scopeArray, arity, defaultCost - minCost, nbTuples);
        for (auto it = costFunction.begin(); it != costFunction.end(); ++it) {
            this->wcsp->postNaryConstraintTupleInternal(naryIndex, it->first, it->second - minCost); // For each tuple
        }
        this->wcsp->postNaryConstraintEnd(naryIndex);
    }
    // all tuples in lexico order
    else {
        if (ToulBar2::verbose >= 3) {
            cout << "read nary cost function on ";
            for (unsigned int i = 0; i < arity; i++) {
                cout << scope[i] << " ";
            }
            cout << endl;
        }

        int cfIndex = this->wcsp->postNaryConstraintBegin(scopeArray, arity, MIN_COST, LONGLONG_MAX);
        NaryConstraint* nctr = (NaryConstraint*)this->wcsp->getCtr(cfIndex);
        Cost cost;
        vector<Cost> costs;

        // Read all costs
        while (!isCBrace(token)) {
            cost = wcsp->decimalToCost(token, lineNumber);
            costs.push_back(cost);
            minCost = min(minCost, cost);
            nbTuples++;
            std::tie(lineNumber, token) = this->getNextToken();
        }

        // Test if all tuples have been read
        if ((logCard > log(std::numeric_limits<size_t>::max())) || nbTuples < card) {
            cerr << "Error : incorrect number of tuples for scope : ";
            for (unsigned int i = 0; i < arity; i++) {
                cout << scope[i] << " ";
            }
            cout << endl;
            throw WrongFileFormat();
        }

        int j = 0;
        nctr->firstlex();
        while (nctr->nextlex(tup, cost)) {
            this->wcsp->postNaryConstraintTupleInternal(cfIndex, tup, costs[j]);
            j++;
        }
        if (ToulBar2::verbose >= 3)
            cout << "read arity " << arity << " table costs." << endl;
        this->wcsp->postNaryConstraintEnd(nctr->wcspIndex);
    }
    wcsp->negCost -= minCost;
    skipCBrace(); // Function final CBrace read.
}

// Reads a global/arithmetic cost function
// Starts: after the OBrace of parameters
// Ends:
void CFNStreamReader::readGlobalCostFunction(vector<int>& scope, const string& funcName, int line)
{
    unsigned int arity = scope.size();

    map<string, string> GCFTemplates = {
        { "clique", ":rhs:N:values:[v+]S" },
        { "knapsack", ":capacity:N:weights:[N]S" },
        { "knapsackv", ":capacity:N:weightedvalues:[VvN]+" },
        { "cfnconstraint", ":cfn:W:lb:c:ub:c:duplicatehard:N:strongduality:N" },
        { "salldiff", ":metric:K:cost:c" },
        { "sgcc", ":metric:K:cost:c:bounds:[vNN]+" }, // Read first keyword then special case processing
        { "ssame", "SPECIAL" }, // Special case processing
        { "sregular", ":metric:K:cost:c:nb_states:N:starts:[N]+:ends:[N]+:transitions:[NvN]+" },
        { "sregulardp", ":metric:K:cost:C:nb_states:N:starts:[N]+:ends:[N]+:transitions:[NvN]+" },
        { "sgrammar", "SPECIAL" }, // Special case processing
        { "sgrammardp", "SPECIAL" }, // Special case processing
        { "samong", ":metric:K:cost:c:min:N:max:N:values:[v]+" },
        { "samongdp", ":metric:K:cost:c:min:N:max:N:values:[v]+" },
        { "salldiffdp", ":metric:K:cost:c" },
        { "salldiffkp", ":metric:K:cost:c" },
        { "sgccdp", ":metric:K:cost:c:bounds:[vNN]+" },
        { "max", ":defaultcost:c:tuples:[Vvc]+" },
        { "smaxdp", ":defaultcost:c:tuples:[Vvc]+" },
        { "MST", "" },
        { "smstdp", "" },
        { "wregular", ":nb_states:N:starts:[NC]+:ends:[NC]+:transitions:[NvNC]+" },
        { "walldiff", ":metric:K:cost:c" },
        { "wgcc", ":metric:K:cost:c:bounds:[vNN]+" },
        { "wsame", ":metric:K:cost:c" },
        { "wsamegcc", ":metric:K:cost:c:bounds:[vNN]+" },
        { "wamong", ":metric:K:cost:c:values:[v]+:min:N:max:N" },
        { "wvaramong", ":metric:K:cost:c:values:[v]+" },
        { "woverlap", ":metric:K:cost:c:comparator:K:to:N" },
        { "wsum", ":metric:K:cost:c:comparator:K:to:N" },
        { "wvarsum", ":metric:K:cost:c:comparator:K" },
        { "wdiverse", ":distance:N:values:[v]S" },
        { "whdiverse", ":distance:N:values:[v]S" },
        { "wtdiverse", ":distance:N:values:[v]S" }
    };

    auto it = GCFTemplates.find(funcName);
    if (it != GCFTemplates.end()) {
        // Reads function using template and generates the corresponding stream
        stringstream paramsStream;

        this->generateGCFStreamFromTemplate(scope, funcName, GCFTemplates[funcName], paramsStream);

        int scopeArray[arity];
        for (unsigned int i = 0; i < arity; i++) {
            scopeArray[i] = scope[i];
        }

        if (funcName == "clique" && !(paramsStream.str().size() > 1 && paramsStream.str()[0] == '1' && paramsStream.str()[1] == ' ')) {
            cerr << "Error: the clique global constraint does not accept RHS different from 1 for now at line" << line << endl;
            throw WrongFileFormat();
        }

        int nbconstr = 0; // unused int for pointer ref
        this->wcsp->postGlobalFunction(scopeArray, (int)arity, funcName, paramsStream, &nbconstr, false);
    }
    // Arithmetic function
    else {
        const set<string> arithmeticFuncNames = { ">=", ">", "<=", "<", "=", "disj", "sdisj" };

        if (arithmeticFuncNames.find(funcName) == arithmeticFuncNames.end()) {
            cerr << "Error: unknown global cost function: " << funcName << " at line " << line << endl;
            throw WrongFileFormat();
        }

        if (arity != 2) {
            cerr << "Error : arithmetic function " << funcName << " has incorrect arity at line " << line << endl;
            throw WrongFileFormat();
        }

        if (ToulBar2::costMultiplier < 0.0 || ToulBar2::decimalPoint != 0) {
            cerr << "Error : arithmetic function " << funcName << " at line " << line << "cannot be used with decimal costs or in maximization mode." << endl;
            throw WrongFileFormat();
        }

        pair<int, string> token = this->getNextToken();
        vector<pair<int, string>> funcParams;

        while (!isCBrace(token.second)) {
            funcParams.push_back(token);
            token = this->getNextToken();
        }

        if (funcName == ">=") {
            if (funcParams.size() != 2) {
                cerr << "Error : arithmetic function " << funcName << " has incorrect number of parameters." << endl;
                throw WrongFileFormat();
            }
            try {
                wcsp->postSupxyc(scope[0], scope[1], stoi(funcParams[0].second), stoi(funcParams[1].second));
            } catch (std::invalid_argument&) {
                cerr << "Error: invalid parameters for '" << funcName << "' at line " << funcParams[0].first << endl;
                throw WrongFileFormat();
            }

        } else if (funcName == ">") {
            if (funcParams.size() != 2) {
                cerr << "Error : arithmetic function " << funcName << " has incorrect number of parameters." << endl;
                throw WrongFileFormat();
            }
            try {
                wcsp->postSupxyc(scope[0], scope[1], stoi(funcParams[0].second) + 1, stoi(funcParams[1].second));
            } catch (std::invalid_argument&) {
                cerr << "Error: invalid parameters for '" << funcName << "' at line " << funcParams[0].first << endl;
                throw WrongFileFormat();
            }

        } else if (funcName == "<=") {
            if (funcParams.size() != 2) {
                cerr << "Error : arithmetic function " << funcName << " has incorrect number of parameters." << endl;
                throw WrongFileFormat();
            }
            try {
                wcsp->postSupxyc(scope[0], scope[1], -stoi(funcParams[0].second), stoi(funcParams[1].second));
            } catch (std::invalid_argument&) {
                cerr << "Error: invalid parameters for '" << funcName << "' at line " << funcParams[0].first << endl;
                throw WrongFileFormat();
            }

        } else if (funcName == "<") {
            if (funcParams.size() != 2) {
                cerr << "Error : arithmetic function " << funcName << " has incorrect number of parameters." << endl;
                throw WrongFileFormat();
            }
            try {
                wcsp->postSupxyc(scope[0], scope[1], -stoi(funcParams[0].second) + 1, stoi(funcParams[1].second));
            } catch (std::invalid_argument&) {
                cerr << "Error: invalid parameters for '" << funcName << "' at line " << funcParams[0].first << endl;
                throw WrongFileFormat();
            }

        } else if (funcName == "=") {
            if (funcParams.size() != 2) {
                cerr << "Error : arithmetic function " << funcName << " has incorrect number of parameters." << endl;
                throw WrongFileFormat();
            }
            try {
                wcsp->postSupxyc(scope[0], scope[1], stoi(funcParams[0].second), stoi(funcParams[1].second));
                wcsp->postSupxyc(scope[1], scope[0], -stoi(funcParams[0].second), stoi(funcParams[1].second));
            } catch (std::invalid_argument&) {
                cerr << "Error: invalid parameters for '" << funcName << "' at line " << funcParams[0].first << endl;
                throw WrongFileFormat();
            }
        } else if (funcName == "disj") {
            if (funcParams.size() != 3) {
                cerr << "Error : arithmetic function " << funcName << " has incorrect number of parameters." << endl;
                throw WrongFileFormat();
            }
            Cost cost = wcsp->decimalToCost(funcParams[2].second, funcParams[2].first);
            try {
                wcsp->postDisjunction(scope[0], scope[1], stoi(funcParams[0].second), stoi(funcParams[1].second), cost);
            } catch (std::invalid_argument&) {
                cerr << "Error: invalid parameters for '" << funcName << "' at line " << funcParams[0].first << endl;
                throw WrongFileFormat();
            }

        } else if (funcName == "sdisj") {
            if (funcParams.size() != 6) {
                cerr << "Error : arithmetic function " << funcName << " has incorrect number of parameters." << endl;
                throw WrongFileFormat();
            }
            Cost cost1 = wcsp->decimalToCost(funcParams[4].second, funcParams[4].first);
            Cost cost2 = wcsp->decimalToCost(funcParams[5].second, funcParams[4].first);
            try {
                wcsp->postSpecialDisjunction(scope[0], scope[1], stoi(funcParams[0].second), stoi(funcParams[1].second),
                    stoi(funcParams[2].second), stoi(funcParams[3].second), cost1, cost2);
            } catch (std::invalid_argument&) {
                cerr << "Error: invalid parameters for '" << funcName << "' at line " << funcParams[0].first << endl;
                throw WrongFileFormat();
            }
        }
        skipCBrace();
    }
}

void CFNStreamReader::generateGCFStreamFromTemplate(vector<int>& scope, const string& funcType, string GCFTemplate, stringstream& stream)
{
    // -------------------- Special cases are treated separately
    if (funcType == "sgrammar" || funcType == "sgrammardp") {
        this->generateGCFStreamSgrammar(scope, stream);
        return;
    } else if (funcType == "ssame") {
        this->generateGCFStreamSsame(scope, stream);
        return;
    }

    // -------------------- Function reading using template
    int lineNumber = -1;
    string token;
    vector<char> repeatedSymbols;
    unsigned int numberOfTuplesRead = 0;
    bool isOpenedBrace = false;
    bool variableRepeat = false;
    vector<pair<char, string>> streamContentVec;

    // Main loop: read template string char by char, and read the CFN file accordingly to the pattern
    for (unsigned int i = 0; i < GCFTemplate.size(); i++) {

        if (isOpenedBrace) {
            if (GCFTemplate[i] == ']') {
                isOpenedBrace = false;
            } else if (GCFTemplate[i] == '+') {
                variableRepeat = true;
                repeatedSymbols.push_back(GCFTemplate[i]);
            } else {
                repeatedSymbols.push_back(GCFTemplate[i]);
            }
        }
        // ---------- Read keyword and add it to stream
        else if (GCFTemplate[i] == 'K') {

            std::tie(lineNumber, token) = this->getNextToken();
            streamContentVec.push_back(std::make_pair('K', token));

            // Special case of sgcc
            if (funcType == "sgcc") {
                if (token == "wdec") {
                    if (ToulBar2::verbose >= 2)
                        cout << "Updating template (wdec) : "
                             << ":metric:K:cost:c:bounds:[vNNcc]+" << endl;
                    GCFTemplate = ":metric:K:cost:c:bounds:[vNNcc]+";
                }
            }
        }
        // ---------- Read cost, transform it to cost and add it to stream
        else if (GCFTemplate[i] == 'C' || GCFTemplate[i] == 'c') {

            std::tie(lineNumber, token) = this->getNextToken();
            Cost cost = MIN_COST;
            if (funcType == "cfnconstraint") {
                int wcspind = stoi(streamContentVec.front().second);
                assert(WCSP::CollectionOfWCSP.find(wcspind) != WCSP::CollectionOfWCSP.end());
                cost = WCSP::CollectionOfWCSP[wcspind]->decimalToCost(token, lineNumber);
            } else {
                cost = wcsp->decimalToCost(token, lineNumber);
                if (GCFTemplate[i] == 'c' && cost < 0) {
                    cerr << "Error: the global cost function " << funcType << " cannot accept negative costs at line " << lineNumber << endl;
                    throw WrongFileFormat();
                }
            }
            streamContentVec.push_back(std::make_pair(GCFTemplate[i], to_string(cost)));
        }
        // ---------- Read variable and add it to stream
        else if (GCFTemplate[i] == 'V') {

            std::tie(lineNumber, token) = this->getNextToken();

            if (not isdigit(token[0])) {
                auto it = varNameToIdx.find(token);
                if (it != varNameToIdx.end()) {
                    token = it->second;
                } else {
                    cerr << "Error: unknown variable with name '" << token << "' at line " << lineNumber << endl;
                    throw WrongFileFormat();
                }
            }

            streamContentVec.push_back(std::make_pair('V', token));
        }
        // ---------- Read value and add it to stream
        else if (GCFTemplate[i] == 'v') {

            std::tie(lineNumber, token) = this->getNextToken();
            // V0 : value MUST be a number
            for (char c : token) {
                if (!isdigit(c)) {
                    cerr << "Error: value index required at line " << lineNumber << " but read " << token << endl;
                    throw WrongFileFormat();
                }
            }
            streamContentVec.push_back(std::make_pair('v', token));

        }
        // ---------- Read number and add it to stream
        else if (GCFTemplate[i] == 'N') {

            std::tie(lineNumber, token) = this->getNextToken();
            for (char c : token) {
                if (!isdigit(c) && c != '-') {
                    cerr << "Error: number required at line " << lineNumber << " but read " << token << endl;
                    throw WrongFileFormat();
                }
            }
            streamContentVec.push_back(std::make_pair('N', token));
        }
        // ---------- Read cost function network and add its index number to stream
        else if (GCFTemplate[i] == 'W') {

            WCSP* wcsp_ = dynamic_cast<WCSP*>(WeightedCSP::makeWeightedCSP(MAX_COST));
            CFNStreamReader fileReader(iStream, wcsp_); // read recursively a cost function network
            lineCount += fileReader.lineCount;
            std::tie(lineNumber, token) = this->getNextToken(); // skip last }

            assert(WCSP::CollectionOfWCSP.find(wcsp_->getIndex()) == WCSP::CollectionOfWCSP.end());
            WCSP::CollectionOfWCSP[wcsp_->getIndex()] = wcsp_; // memorize WCSP object
            streamContentVec.push_back(std::make_pair('W', to_string(wcsp_->getIndex())));
        }
        // ---------- Read JSON tag
        else if (GCFTemplate[i] == ':') {
            auto idx = GCFTemplate.find_first_of(':', i + 1);
            string jsonTag = GCFTemplate.substr(i + 1, idx - i - 1); // extract 'jsonTag'
            i += jsonTag.size() + 1; // Increase i according to length of token read
            skipJSONTag(jsonTag);
        }
        // ---------- Entering a repeated section
        else if (GCFTemplate[i] == '[') {
            isOpenedBrace = true;
        }
        // Read function repeated section
        else if ((GCFTemplate[i] == '+' || GCFTemplate[i] == 'S') && !isOpenedBrace) {
            vector<pair<char, string>> repeatedContentVec; // Function repeated params
            vector<pair<char, string>> variableRepeatVec; // Stack if internal tuples have unknown size
            // [ delimiting the start of the list
            skipOBrace();
            // Inside the list of parameter tuples

            std::tie(lineNumber, token) = this->getNextToken();
            while (token != "]") {
                // Each (non unary) tuple is inside []. Skip first [
                if (repeatedSymbols.size() > 1 || variableRepeat) {
                    if (!isOBrace(token)) {
                        cerr << "Error: expected '[/{' but read " << token << " at line " << lineNumber << endl;
                        throw WrongFileFormat();
                    } else
                        std::tie(lineNumber, token) = this->getNextToken();
                }

                size_t repeatIndex = 0;
                while ((repeatIndex < repeatedSymbols.size()) && !isCBrace(token)) {
                    if (repeatedSymbols[repeatIndex] == '+')
                        repeatIndex = 0;
                    char symbol = repeatedSymbols[repeatIndex];
                    if (symbol == 'N') {
                        for (char c : token) {
                            if (!isdigit(c) && c != '-') {
                                cerr << "Error: integer required at line " << lineNumber << " but read " << token << endl;
                                throw WrongFileFormat();
                            }
                        }
                        if (variableRepeat)
                            variableRepeatVec.push_back(std::make_pair('N', token));
                        else
                            repeatedContentVec.push_back(std::make_pair('N', token));
                    } else if (symbol == 'V') {
                        // If variable name (string)
                        if (not isdigit(token[0])) {
                            auto it = varNameToIdx.find(token);
                            if (it != varNameToIdx.end()) {
                                token = to_string(it->second);
                            } else {
                                cerr << "Error: unknown variable with name '" << token << "' at line " << lineNumber << endl;
                                throw WrongFileFormat();
                            }
                        }
                        if (variableRepeat)
                            variableRepeatVec.push_back(std::make_pair('V', token));
                        else
                            repeatedContentVec.push_back(std::make_pair('V', token));
                    } else if (symbol == 'v') {
                        // V0 : value MUST be a number
                        for (char c : token) {
                            if (!isdigit(c)) {
                                cerr << "Error: value index required at line " << lineNumber << " but read " << token << endl;
                                throw WrongFileFormat();
                            }
                        }
                        if (variableRepeat)
                            variableRepeatVec.push_back(std::make_pair('v', token));
                        else
                            repeatedContentVec.push_back(std::make_pair('v', token));
                    } else if ((symbol == 'C') || (symbol == 'c')) {
                        Cost c = wcsp->decimalToCost(token, lineNumber);
                        if (symbol == 'c' && c < 0) {
                            cerr << "Error: the global cost function " << funcType << " cannot accept negative costs at line " << lineNumber << endl;
                            throw WrongFileFormat();
                        }
                        if (variableRepeat)
                            variableRepeatVec.push_back(std::make_pair(symbol, to_string(c)));
                        else
                            repeatedContentVec.push_back(std::make_pair(symbol, to_string(c)));
                    }

                    repeatIndex++;
                    std::tie(lineNumber, token) = this->getNextToken();
                }

                if (repeatedSymbols.size() > 1 || variableRepeat) {
                    if (!isCBrace(token)) {
                        cerr << "Error: expected ']/}' but read " << token << " at line " << lineNumber << endl;
                        throw WrongFileFormat();
                    } else
                        std::tie(lineNumber, token) = this->getNextToken();
                }
                if (variableRepeat) { // we must push the size of the repeat and its contents
                    repeatedContentVec.push_back(std::make_pair('N', to_string(variableRepeatVec.size())));
                    repeatedContentVec.insert(repeatedContentVec.end(), variableRepeatVec.begin(), variableRepeatVec.end());
                    variableRepeatVec.clear();
                }
                numberOfTuplesRead++; // Number of tuples read
            }
            if (GCFTemplate[i] == 'S' && numberOfTuplesRead != scope.size()) {
                cerr << "Error: expected " << scope.size() << " tuples for '" << funcType << "' but read " << numberOfTuplesRead << " at line " << lineNumber << endl;
                throw WrongFileFormat();
            }
            // Add number of tuples before the list if the number of expected tuples is not known
            if (GCFTemplate[i] == '+')
                streamContentVec.push_back(std::make_pair('N', to_string(numberOfTuplesRead)));
            // Copy repeatedContentVec to streamContent
            for (pair<char, string> repContentPair : repeatedContentVec) {
                streamContentVec.push_back(std::make_pair(repContentPair.first, repContentPair.second));
            }

            // Reset vector for future usage
            numberOfTuplesRead = 0;
            repeatedSymbols.clear();

        } // end of if (GCFTemplate[i] == '+')
    } // end of for (unsigned int i=0; i < GCFTemplate.size(); i++)

    // End of params
    skipCBrace();
    // End of function
    skipCBrace();

    // -------------------- Data processing
    Cost minCost = MAX_COST;

    // DEBUG DISPLAY
    if (ToulBar2::verbose >= 2) {
        cout << "Output Data map :" << endl;
        for (unsigned int i = 0; i < streamContentVec.size(); i++) {
            cout << streamContentVec[i].first << "\t" << streamContentVec[i].second << endl;
        }
    }

    // FIND MIN COST
    bool minUpdated = false;
    for (pair<char, string> streamContentPair : streamContentVec) {
        if (streamContentPair.first == 'C') {
            Cost currentCost = (Cost)std::stoll(streamContentPair.second);
            minCost = MIN(currentCost, minCost);
            minUpdated = true;
        }
    }
    if (!minUpdated)
        minCost = 0;

    // WRITE ALL TO STREAM AND SUBSTRACT MIN COST TO ALL COSTS
    for (unsigned int i = 0; i < streamContentVec.size(); i++) {
        if (streamContentVec[i].first == 'C') {
            Cost currentCost = (Cost)std::stoll(streamContentVec[i].second);
            currentCost -= minCost;
            streamContentVec[i].second = to_string(currentCost);
        }
        stream << streamContentVec[i].second << " ";
    }

    // Correct for negative costs
    if (funcType == "wregular") { // regular: we can handle all costs. The number of transitions is known and we have one start and end state
        wcsp->negCost -= ((scope.size() + 2) * minCost); // TODO we could do better and compute different mins for initial/final/transitions.
    } else
        wcsp->negCost -= minCost;

    // STREAM DEBUG
    if (ToulBar2::verbose >= 1)
        cout << "Stream for " << funcType << ": '" << stream.str() << "'" << endl;
    return;
}

/*
 * Example :
 * metric : var
 * cost : 15
 * nb_symbols : 4
 * nb_values : 2
 * start_symbol : 0
 * terminals : [ [1 0] [3 1] ]
 * non_terminals : [ [0 0 0] [0 1 2] [0 1 3] [2 0 3] ]
 * return stream : [var|weight cost nb_symbols nb_values start_symbol nb_rules ((0 terminal_symbol value)|(1 nonterminal_in nonterminal_out_left nonterminal_out_right)|(2 terminal_symbol value weight)|(3 nonterminal_in nonterminal_out_left nonterminal_out_right weight))∗]
 */
void CFNStreamReader::generateGCFStreamSgrammar(vector<int>& scope, stringstream& stream)
{

    int lineNumber;
    string token;
    string metric;
    vector<string> terminal_rules;
    vector<string> non_terminal_rules;

    skipJSONTag("metric");
    std::tie(lineNumber, token) = this->getNextToken();
    metric = token;
    if (metric != "var" && metric != "weight") {
        cerr << "Error: sgrammar metric must be either 'var' or 'weight' at line " << lineNumber << endl;
        throw WrongFileFormat();
    }
    // Read cost
    skipJSONTag("cost");
    std::tie(lineNumber, token) = this->getNextToken();
    Cost cost = wcsp->decimalToCost(token, lineNumber);
    if (cost < 0) {
        cerr << "Error: sgrammar at line " << lineNumber << "uses a negative cost." << endl;
        throw WrongFileFormat();
    }
    // Read Nb Symbols
    skipJSONTag("nb_symbols");
    std::tie(lineNumber, token) = this->getNextToken();
    string nb_symbols = token;
    // Read Nb Values
    skipJSONTag("nb_values");
    std::tie(lineNumber, token) = this->getNextToken();
    string nb_values = token;
    // Read start symbol
    skipJSONTag("start");
    std::tie(lineNumber, token) = this->getNextToken();
    string start_symbol = token;

    skipJSONTag("terminals"); // 0 or 2
    skipOBrace();
    std::tie(lineNumber, token) = this->getNextToken(); // Second [ or ]
    while (token != "]") {

        string terminal_rule;
        yellOBrace(token, lineNumber);
        // Read terminal_symbol
        std::tie(lineNumber, token) = this->getNextToken();
        terminal_rule += token + " ";

        // Read value
        std::tie(lineNumber, token) = this->getNextToken();
        terminal_rule += token + " ";

        if (metric == "weight") {
            // Read weight
            std::tie(lineNumber, token) = this->getNextToken();
            Cost tcost = wcsp->decimalToCost(token, lineNumber);
            if (tcost < 0) {
                cerr << "Error: sgrammar at line " << lineNumber << "uses a negative cost." << endl;
                throw WrongFileFormat();
            }
            terminal_rule += to_string(tcost) + to_string(" ");
        }

        terminal_rules.push_back(terminal_rule);

        skipCBrace();
        std::tie(lineNumber, token) = this->getNextToken();
    }

    skipJSONTag("non_terminals"); // 1 or 3
    skipOBrace();
    std::tie(lineNumber, token) = this->getNextToken(); // Second [ or ]
    while (token != "]") {

        string non_terminal_rule;
        yellOBrace(token, lineNumber);

        // Read nonterminal_in
        std::tie(lineNumber, token) = this->getNextToken();
        non_terminal_rule += token + " ";

        // Read nonterminal_out_left
        std::tie(lineNumber, token) = this->getNextToken();
        non_terminal_rule += token + " ";

        // Read nonterminal_out_right
        std::tie(lineNumber, token) = this->getNextToken();
        non_terminal_rule += token + " ";

        if (metric == "weight") {
            // Read weight
            std::tie(lineNumber, token) = this->getNextToken();
            Cost tcost = wcsp->decimalToCost(token, lineNumber);
            if (tcost < 0) {
                cerr << "Error: sgrammar at line " << lineNumber << "uses a negative cost." << endl;
                throw WrongFileFormat();
            }
            non_terminal_rule += to_string(tcost) + to_string(" ");
        }

        non_terminal_rules.push_back(non_terminal_rule);

        skipCBrace();
        std::tie(lineNumber, token) = this->getNextToken();
    }

    // End of function, write to stream

    // End of params
    skipCBrace();
    // End of function
    skipCBrace();

    // Cost had no impact on negCost here
    stream << metric << " " << cost << " " << nb_symbols << " " << nb_values << " " << start_symbol
           << " " << to_string(terminal_rules.size() + non_terminal_rules.size()) << " ";
    if (metric == "var") {
        for (string terminal_rule : terminal_rules)
            stream << "0 " << terminal_rule;
        for (string non_terminal_rule : non_terminal_rules)
            stream << "1 " << non_terminal_rule;
    } else if (metric == "weight") {
        for (string terminal_rule : terminal_rules)
            stream << "2 " << terminal_rule;
        for (string non_terminal_rule : non_terminal_rules)
            stream << "3 " << non_terminal_rule;
    }

    if (ToulBar2::verbose >= 1)
        cout << "Stream for sgrammar : '" << stream.str() << "'" << endl;
}
/*
 * Example :
 * cost : 10.8
 * vars1 : [v1 v2 v3]
 * vars2 : [v4 v5 v6]
 * return stream : [cost list_size1 list_size2 (variable_index)∗ (variable_index)∗]
 */
void CFNStreamReader::generateGCFStreamSsame(vector<int>& scope, stringstream& stream)
{

    int lineNumber;
    string token;
    vector<string> variables1;
    vector<string> variables2;

    skipJSONTag("cost");
    std::tie(lineNumber, token) = this->getNextToken();
    Cost cost = wcsp->decimalToCost(token, lineNumber);
    // TODO Cost should be >= 0

    skipJSONTag("vars1");
    skipOBrace();
    std::tie(lineNumber, token) = this->getNextToken();

    while (token != "]") {

        if (not isdigit(token[0])) {
            map<string, int>::iterator it;
            if ((it = varNameToIdx.find(token)) != varNameToIdx.end()) {
                token = to_string(it->second);
            } else {
                cerr << "Error: unknown variable with name '" << token << "' at line " << lineNumber << endl;
                throw WrongFileFormat();
            }
        }
        variables1.push_back(token);

        std::tie(lineNumber, token) = this->getNextToken();
    }

    skipJSONTag("vars2");
    skipOBrace();
    std::tie(lineNumber, token) = this->getNextToken();

    while (token != "]") {

        if (not isdigit(token[0])) {
            map<string, int>::iterator it;
            if ((it = varNameToIdx.find(token)) != varNameToIdx.end()) {
                token = to_string(it->second);
            } else {
                cerr << "Error: unknown variable with name '" << token << "' at line " << lineNumber << endl;
                throw WrongFileFormat();
            }
        }
        variables2.push_back(token);

        std::tie(lineNumber, token) = this->getNextToken();
    }

    // End of params
    skipCBrace();
    // End of function
    skipCBrace();

    // Cost has no impact on negCost here
    stream << cost << " ";
    stream << variables1.size() << " " << variables2.size() << " ";
    for (string var1 : variables1)
        stream << var1 << " ";
    for (string var2 : variables2)
        stream << var2 << " ";

    if (ToulBar2::verbose >= 1)
        cout << "Stream for ssame : '" << stream.str() << "'" << endl;

    return;
}
#endif

// TB2 entry point for WCSP reading (not only wcsp format).
// Returns the global UB obtained form both the file and command line in internal Cost units
Cost WCSP::read_wcsp(const char* fileName)
{
    char* Nfile2;
    Nfile2 = strdup(fileName);
    name = string(basename(Nfile2));
    free(Nfile2);

    // Done internally by the CFN reader
    if (!ToulBar2::cfn && !ToulBar2::lp) {
        if (ToulBar2::deltaUbS.length() != 0) {
            ToulBar2::deltaUbAbsolute = string2Cost(ToulBar2::deltaUbS.c_str());
            ToulBar2::deltaUb = ToulBar2::deltaUbAbsolute;
        }

        if (ToulBar2::externalUB.size()) {
            Cost top = string2Cost(ToulBar2::externalUB.c_str());
            Double K = ToulBar2::costMultiplier;
            if (top < MAX_COST / K)
                top = top * K;
            else
                top = MAX_COST;
            ToulBar2::deltaUb = max(ToulBar2::deltaUbAbsolute, (Cost)(ToulBar2::deltaUbRelativeGap * (Double)min(top, getUb())));
            updateUb(top + ToulBar2::deltaUb);
            // as long as a true certificate as not been found we must compensate for the deltaUb in CUT
        }

        if (ToulBar2::costThresholdS.size())
            ToulBar2::costThreshold = string2Cost(ToulBar2::costThresholdS.c_str());
        if (ToulBar2::costThresholdPreS.size())
            ToulBar2::costThresholdPre = string2Cost(ToulBar2::costThresholdPreS.c_str());
        if (ToulBar2::vnsOptimumS.size())
            ToulBar2::vnsOptimum = string2Cost(ToulBar2::vnsOptimumS.c_str());
    }

    if (ToulBar2::cfn && !ToulBar2::gz && !ToulBar2::bz2 && !ToulBar2::xz) {
#ifdef BOOST
        ifstream Rfile;
        istream& stream = (ToulBar2::stdin_format.length() > 0) ? cin : Rfile;
        if (ToulBar2::stdin_format.compare("cfn") == 0) {
            CFNStreamReader fileReader(stream, this);

        } else {
            Rfile.open(fileName);
            if (!stream) {
                cerr << "Error: could not open file '" << fileName << "'." << endl;
                throw WrongFileFormat();

            } else {
                CFNStreamReader fileReader(stream, this);
            }
        }
#else
        cerr << "Error: compiling with Boost library is needed to allow to read CFN format files." << endl;
        throw WrongFileFormat();
#endif
    } else if (ToulBar2::cfn && ToulBar2::gz) {
#ifdef BOOST
        ifstream Rfile(fileName, std::ios_base::in | std::ios_base::binary);
        istream& file = Rfile;
        boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
        inbuf.push(boost::iostreams::gzip_decompressor());
        inbuf.push(file);
        std::istream stream(&inbuf);

        if (!file) {
            cerr << "Could not open cfn.gz file : " << fileName << endl;
            throw WrongFileFormat();
        } else {

            //  inbuf.push(file);
            CFNStreamReader fileReader(stream, this);
        }
#else
        cerr << "Error: compiling with Boost iostreams library is needed to allow to read gzip'd CFN format files." << endl;
        throw WrongFileFormat();
#endif
    } else if (ToulBar2::cfn && ToulBar2::bz2) {
#ifdef BOOST
#ifndef NO_BZ2
        ifstream Rfile(fileName, std::ios_base::in | std::ios_base::binary);
        istream& file = Rfile;
        boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
        inbuf.push(boost::iostreams::bzip2_decompressor());
        inbuf.push(file);
        std::istream stream(&inbuf);

        if (!file) {
            cerr << "Could not open cfn.bz2 file : " << fileName << endl;
            throw WrongFileFormat();
        } else {

            //  inbuf.push(file);
            CFNStreamReader fileReader(stream, this);
        }
#else
        cerr << "Error: Boost iostreams library requires bzip2 library to read bzip2'd CFN format files." << endl;
        throw WrongFileFormat();
#endif
#else
        cerr << "Error: compiling with Boost iostreams library is needed to allow to read bzip2'd CFN format files." << endl;
        throw WrongFileFormat();
#endif
    } else if (ToulBar2::cfn && ToulBar2::xz) {
#ifdef BOOST
#ifndef NO_LZMA
        ifstream Rfile(fileName, std::ios_base::in | std::ios_base::binary);
        istream& file = Rfile;
        boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
        inbuf.push(boost::iostreams::lzma_decompressor());
        inbuf.push(file);
        std::istream stream(&inbuf);

        if (!file) {
            cerr << "Could not open cfn.xz file : " << fileName << endl;
            throw WrongFileFormat();
        } else {

            //  inbuf.push(file);
            CFNStreamReader fileReader(stream, this);
        }
#else
        cerr << "Error: Boost iostream version 1.65 or higher requires lzma library to read xz compressed CFN format files." << endl;
        throw WrongFileFormat();
#endif
#else
        cerr << "Error: compiling with Boost iostreams library is needed to read xz compressed CFN format files." << endl;
        throw WrongFileFormat();
#endif
    } else if (ToulBar2::haplotype) {
        ToulBar2::haplotype->read(fileName, this);
    } else if (ToulBar2::pedigree) {
        if (!ToulBar2::bayesian)
            ToulBar2::pedigree->read(fileName, this);
        else
            ToulBar2::pedigree->read_bayesian(fileName, this);
    } else if (ToulBar2::uai) {
        read_uai2008(fileName);
    } else if (ToulBar2::xmlflag) {
        read_XML(fileName);
    } else if (ToulBar2::bep) {
        ToulBar2::bep->read(fileName, this);
    } else if (ToulBar2::wcnf) {
        read_wcnf(fileName);
    } else if (ToulBar2::qpbo) {
        read_qpbo(fileName);
    } else if (ToulBar2::opb) {
        read_opb(fileName);
    } else if (ToulBar2::lp) {
        read_lp(fileName);
    } else {
        read_legacy(fileName);
    }

    // common ending section for all readers

    //#ifdef BOOST
    //    if (ToulBar2::addAMOConstraints!=-1) {
    //        propagate(); // CHOOSE: initial propagation must be done beforehand (with full propagation for knapsack constraints)
    //        ToulBar2::addAMOConstraints_ = true;
    //        propagate(); // OR CHOOSE: initial propagation must be done beforehand (but only bound propagation for knapsack constraints)
    //        addAMOConstraints();
    //        ToulBar2::addAMOConstraints_ = false;
    //    }
    //#endif

    // Diverse variables structure and variables allocation and initialization
    if (ToulBar2::divNbSol > 1) {
        initDivVariables();
        // Dual variables allocation, only needed for divMethod 0 Dual or 1 Hidden
        if (ToulBar2::divMethod < 2) {
            divVarsId.resize(ToulBar2::divNbSol);
            for (unsigned int j = 0; j < ToulBar2::divNbSol - 1; j++) {
                for (Variable* x : divVariables) {
                    int xId = x->wcspIndex;
                    divVarsId[j][xId] = makeEnumeratedVariable(HIDDEN_VAR_TAG + to_string("c_sol") + to_string(j) + to_string("_") + x->getName(), 0, 2 * ToulBar2::divBound + 1);
                    EnumeratedVariable* theVar = static_cast<EnumeratedVariable*>(getVar(divVarsId[j][xId]));
                    for (unsigned int val = 0; val < theVar->getDomainInitSize(); val++) {
                        theVar->addValueName(to_string("q") + to_string(val % (ToulBar2::divBound + 1)) + to_string("_")
                            + to_string(val / (ToulBar2::divBound + 1)));
                    }
                }
            }
        }

        // Hidden variables, only needed for divMethod 1 Hidden or 2 Ternary
        if (ToulBar2::divMethod >= 1 && ToulBar2::divMethod <= 2) {
            divHVarsId.resize(ToulBar2::divNbSol); // make room for hidden state variables
            for (unsigned int j = 0; j < ToulBar2::divNbSol - 1; j++) {
                for (Variable* x : divVariables) {
                    int xId = x->wcspIndex;
                    divHVarsId[j][xId] = makeEnumeratedVariable(HIDDEN_VAR_TAG + to_string("h_sol") + to_string(j) + to_string("_") + x->getName(), 0, ToulBar2::divBound);
                    EnumeratedVariable* theVar = static_cast<EnumeratedVariable*>(getVar(divHVarsId[j][xId]));
                    for (unsigned int val = 0; val < theVar->getDomainInitSize(); val++) {
                        theVar->addValueName(to_string("q") + to_string(val));
                    }
                }
            }
        }

        // Joint DivMin MDD
        if (ToulBar2::divWidth > 0) { // add variables for relaxed constraint
            if (ToulBar2::divMethod < 2) {
                for (Variable* x : divVariables) {
                    int xId = x->wcspIndex;
                    divVarsId[ToulBar2::divNbSol - 1][xId] = makeEnumeratedVariable(HIDDEN_VAR_TAG + "c_relax_" + x->getName(), 0, ToulBar2::divWidth * ToulBar2::divWidth - 1);
                    EnumeratedVariable* theVar = static_cast<EnumeratedVariable*>(getVar(divVarsId[ToulBar2::divNbSol - 1][xId]));
                    for (unsigned int val = 0; val < theVar->getDomainInitSize(); val++) {
                        theVar->addValueName(to_string("Q") + to_string(val));
                    }
                }
            }
            if (ToulBar2::divMethod >= 1 && ToulBar2::divMethod <= 2) {
                for (Variable* x : divVariables) {
                    int xId = x->wcspIndex;
                    divHVarsId[ToulBar2::divNbSol - 1][xId] = makeEnumeratedVariable(HIDDEN_VAR_TAG + "h_relax_" + x->getName(), 0, ToulBar2::divWidth - 1);
                    EnumeratedVariable* theVar = static_cast<EnumeratedVariable*>(getVar(divHVarsId[ToulBar2::divNbSol - 1][xId]));
                    for (unsigned int val = 0; val < theVar->getDomainInitSize(); val++) {
                        theVar->addValueName(to_string("q") + to_string(val));
                    }
                }
            }
        }
    }
    return getUb();
}

// TOULBAR2 WCSP LEGACY PARSER
void WCSP::read_legacy(const char* fileName)
{
    ifstream rfile(fileName, (ToulBar2::gz || ToulBar2::bz2 || ToulBar2::xz) ? (std::ios_base::in | std::ios_base::binary) : (std::ios_base::in));
#ifdef BOOST
    boost::iostreams::filtering_streambuf<boost::iostreams::input> zfile;
    if (ToulBar2::gz) {
        zfile.push(boost::iostreams::gzip_decompressor());
    } else if (ToulBar2::bz2) {
#ifndef NO_BZ2
        zfile.push(boost::iostreams::bzip2_decompressor());
#else
        cerr << "Error: Boost requires the bzip2 library to read bz2 compressed wcsp format files." << endl;
        throw WrongFileFormat();
#endif
    } else if (ToulBar2::xz) {
#ifndef NO_LZMA
        zfile.push(boost::iostreams::lzma_decompressor());
#else
        cerr << "Error: compiling with Boost version 1.65 or higher is needed to allow to read xz compressed wcsp format files." << endl;
        throw WrongFileFormat();
#endif
    }
    zfile.push(rfile);
    istream ifile(&zfile);

    if (ToulBar2::stdin_format.length() == 0 && !rfile) {
        cerr << "Could not open wcsp file : " << fileName << endl;
        throw WrongFileFormat();
    }
    istream& file = (ToulBar2::stdin_format.length() > 0) ? cin : ifile;
#else
    if (ToulBar2::gz || ToulBar2::bz2 || ToulBar2::xz) {
        cerr << "Error: compiling with Boost iostreams library is needed to allow to read compressed wcsp format files." << endl;
        throw WrongFileFormat();
    }
    if (ToulBar2::stdin_format.length() == 0 && !rfile) {
        cerr << "Could not open wcsp file : " << fileName << endl;
        throw WrongFileFormat();
    }
    istream& file = (ToulBar2::stdin_format.length() > 0) ? cin : rfile;
#endif
    read_legacy(file);
}

void WCSP::read_legacy(istream& file)
{
    string pbname;
    unsigned int nbvar, nbval;
    int nbconstr;
    int nbvaltrue = 0;
    Cost top;
    int i, j, k, t, ic;
    string varname;
    int domsize;
    unsigned int a;
    unsigned int b;
    unsigned int c;
    Cost defval;
    Cost cost;
    Long ntuples;
    int arity;
    string funcname;
    Value funcparam1;
    Value funcparam2;
    vector<TemporaryUnaryConstraint> unaryconstrs;
    Cost inclowerbound = MIN_COST;
    int maxarity = 0;
    vector<Long> sharedSize;
    vector<vector<Cost>> sharedCosts;
    vector<vector<Tuple>> sharedTuples;
    vector<Tuple> emptyTuples;

    // ---------- PROBLEM HEADER ----------
    // read problem name and sizes
    file >> pbname;
    file >> nbvar;
    file >> nbval;
    file >> nbconstr;
    file >> top;
    if (ToulBar2::verbose >= 1)
        cout << "Read problem: " << pbname << endl;

    Double K = ToulBar2::costMultiplier;
    if (top < MAX_COST / K)
        top = top * K;
    else
        top = MAX_COST;

    ToulBar2::deltaUb = max(ToulBar2::deltaUbAbsolute, (Cost)(ToulBar2::deltaUbRelativeGap * (Double)min(top, getUb())));
    updateUb(top + ToulBar2::deltaUb);

    Tuple tup;
    vector<Tuple> tuples;
    vector<Cost> costs;

    // read variable domain sizes
    for (unsigned int i = 0; i < nbvar; i++) {
        string varname = to_string("x") + to_string(i);
        file >> domsize;
        if (domsize > nbvaltrue)
            nbvaltrue = domsize;
        if (ToulBar2::verbose >= 3)
            cout << "read " << ((i >= numberOfVariables()) ? "new" : "known") << " variable " << i << " of size " << domsize << endl;
        if (i >= numberOfVariables()) {
            DEBONLY(int theindex =)
            ((domsize >= 0) ? makeEnumeratedVariable(varname, 0, domsize - 1) : makeIntervalVariable(varname, 0, -domsize - 1));
            assert(theindex == (int)i);
        } else {
            if ((domsize >= 0) != getVar(i)->enumerated()) {
                cerr << "Variable(" << i << ") " << getVar(i)->getName() << " has a previous domain type (" << (getVar(i)->enumerated() ? ((EnumeratedVariable*)getVar(i))->getDomainInitSize() : getVar(i)->getDomainSize()) << ") different than the new one (" << domsize << ")!" << endl;
                throw WrongFileFormat();
            } else if (domsize < 0) {
                decrease(i, -domsize - 1);
            } else if (domsize >= 0 && (unsigned int)domsize != ((EnumeratedVariable*)getVar(i))->getDomainInitSize()) {
                cerr << "Variable(" << i << ") " << getVar(i)->getName() << " has a previous domain size " << (getVar(i)->enumerated() ? ((EnumeratedVariable*)getVar(i))->getDomainInitSize() : getVar(i)->getDomainSize()) << " different than the new one of " << domsize << "!" << endl;
                throw WrongFileFormat();
            }
        }
    }

    // read each constraint
    for (ic = 0; ic < nbconstr; ic++) {
        file >> arity;
        if (!file) {
            cerr << "Warning! EOF reached before reading all the cost functions (initial number of cost functions too large?)" << endl;
            break;
        }
        bool shared = (arity < 0);
        if (shared)
            arity = -arity;
        // ARITY > 3
        if (arity > NARYPROJECTIONSIZE) {
            maxarity = max(maxarity, arity);
            if (ToulBar2::verbose >= 3)
                cout << "read " << arity << "-ary cost function " << ic << " on";
            int scopeIndex[arity]; // replace arity by MAX_ARITY in case of compilation problem
            for (i = 0; i < arity; i++) {
                file >> j;
                if (ToulBar2::verbose >= 3)
                    cout << " " << j;
                scopeIndex[i] = j;
            }
            if (ToulBar2::verbose >= 3)
                cout << endl;
            file >> defval;
            if (defval == -1) {
                string gcname;
                file >> gcname;
                postGlobalFunction(scopeIndex, arity, gcname, file, &nbconstr);
            } else {
                file >> ntuples;
                int reusedconstr = -1;
                bool reused = (ntuples < 0);
                if (reused) {
                    reusedconstr = -ntuples - 1;
                    if (reusedconstr >= (int)sharedSize.size()) {
                        cerr << "Shared cost function number " << reusedconstr << " not already defined! Cannot reuse it!" << endl;
                        throw WrongFileFormat();
                    }
                    ntuples = sharedSize[reusedconstr];
                }
                if ((defval != MIN_COST) || (ntuples > 0)) {
                    Cost tmpcost = MULT(defval, K);
                    if (CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST * getUb()) && getUb() < (MAX_COST / MEDIUM_COST))
                        tmpcost *= MEDIUM_COST;
                    int naryIndex = postNaryConstraintBegin(scopeIndex, arity, tmpcost, ntuples);
                    NaryConstraint* nary = (NaryConstraint*)constrs[naryIndex];

                    tup.resize(arity);
                    tuples.clear();
                    costs.clear();
                    for (t = 0; t < ntuples; t++) {
                        if (!reused) {
                            for (i = 0; i < arity; i++) {
                                file >> tup[i]; // FIXME: why not translating from Value to tValue? (there is no difference if reading a wcsp file with all domains starting at zero)
                            }
                            file >> cost;
                            Cost tmpcost = MULT(cost, K);
                            if (CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST * getUb()) && getUb() < (MAX_COST / MEDIUM_COST))
                                tmpcost *= MEDIUM_COST;
                            if (shared) {
                                tuples.push_back(tup);
                                costs.push_back(tmpcost);
                            }
                            postNaryConstraintTupleInternal(naryIndex, tup, tmpcost);
                        } else {
                            postNaryConstraintTupleInternal(naryIndex, sharedTuples[reusedconstr][t], sharedCosts[reusedconstr][t]);
                        }
                    }
                    if (shared) {
                        assert(ntuples == (int)costs.size());
                        sharedSize.push_back(costs.size());
                        sharedCosts.push_back(costs);
                        sharedTuples.push_back(tuples);
                    }

                    if (ToulBar2::preprocessNary > 0) {
                        Cost minc = nary->getMinCost();
                        if (minc > MIN_COST) {
                            nary->addtoTuples(-minc);
                            if (ToulBar2::verbose >= 2)
                                cout << "IC0 performed for cost function " << nary << " with initial minimum cost " << minc << endl;
                            inclowerbound += minc;
                        }
                    }
                    postNaryConstraintEnd(naryIndex);
                }
            }
            // ARITY 3
        } else if (arity == 3) {
            maxarity = max(maxarity, arity);
            file >> i;
            file >> j;
            file >> k;
            if ((i == j) || (i == k) || (k == j)) {
                cerr << "Error: ternary cost function!" << endl;
                throw WrongFileFormat();
            }
            file >> defval;
            if (defval >= MIN_COST) {
                assert(vars[i]->enumerated());
                assert(vars[j]->enumerated());
                assert(vars[k]->enumerated());
                EnumeratedVariable* x = (EnumeratedVariable*)vars[i];
                EnumeratedVariable* y = (EnumeratedVariable*)vars[j];
                EnumeratedVariable* z = (EnumeratedVariable*)vars[k];
                if (ToulBar2::verbose >= 3)
                    cout << "read ternary cost function " << ic << " on " << i << "," << j << "," << k << endl;
                file >> ntuples;
                if (ntuples < 0) {
                    int reusedconstr = -ntuples - 1;
                    if (reusedconstr >= (int)sharedSize.size()) {
                        cerr << "Shared cost function number " << reusedconstr << " not already defined! Cannot reuse it!" << endl;
                        throw WrongFileFormat();
                    }
                    ntuples = sharedSize[reusedconstr];
                    assert(ntuples == (int)(x->getDomainInitSize() * y->getDomainInitSize() * z->getDomainInitSize()));
                    if ((defval != MIN_COST) || (ntuples > 0))
                        postTernaryConstraint(i, j, k, sharedCosts[reusedconstr]);
                    continue;
                }
                vector<Cost> costs;
                for (a = 0; a < x->getDomainInitSize(); a++) {
                    for (b = 0; b < y->getDomainInitSize(); b++) {
                        for (c = 0; c < z->getDomainInitSize(); c++) {
                            Cost tmpcost = MULT(defval, K);
                            if (CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST * getUb()))
                                tmpcost *= MEDIUM_COST;
                            costs.push_back(tmpcost);
                        }
                    }
                }
                for (t = 0; t < ntuples; t++) {
                    file >> a;
                    file >> b;
                    file >> c;
                    file >> cost;
                    Cost tmpcost = MULT(cost, K);
                    if (CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST * getUb()) && getUb() < (MAX_COST / MEDIUM_COST))
                        tmpcost *= MEDIUM_COST;
                    assert(a >= 0 && a < x->getDomainInitSize());
                    assert(b >= 0 && b < y->getDomainInitSize());
                    assert(c >= 0 && c < z->getDomainInitSize());
                    costs[(size_t)a * y->getDomainInitSize() * z->getDomainInitSize() + (size_t)b * z->getDomainInitSize() + c] = tmpcost;
                }
                if (shared) {
                    sharedSize.push_back(costs.size());
                    sharedCosts.push_back(costs);
                    sharedTuples.push_back(emptyTuples);
                }
                if ((defval != MIN_COST) || (ntuples > 0))
                    postTernaryConstraint(i, j, k, costs);
            } else if (defval == -1) {
                int scopeIndex[3];
                scopeIndex[0] = i;
                scopeIndex[1] = j;
                scopeIndex[2] = k;
                string gcname;
                file >> gcname;
                if (gcname == "clique") {
                    // postCliqueConstraint(scopeIndex, arity, file);
                    int skipread;
                    file >> skipread;
                    for (int a = 0; a < arity; a++) {
                        file >> skipread;
                        for (int b = skipread; b > 0; b--)
                            file >> skipread;
                    }
                    assert(vars[i]->enumerated());
                    assert(vars[j]->enumerated());
                    assert(vars[k]->enumerated());
                    EnumeratedVariable* x = (EnumeratedVariable*)vars[i];
                    EnumeratedVariable* y = (EnumeratedVariable*)vars[j];
                    EnumeratedVariable* z = (EnumeratedVariable*)vars[k];
                    vector<Cost> costs((size_t)x->getDomainInitSize() * (size_t)y->getDomainInitSize() * (size_t)z->getDomainInitSize(), MIN_COST);
                    postTernaryConstraint(i, j, k, costs); // generate a zero-cost ternary constraint instead that will absorb all its binary hard constraints
                } else {
                    postGlobalFunction(scopeIndex, 3, gcname, file, &nbconstr);
                }
            }
            // ARITY 2
        } else if (arity == 2) {
            maxarity = max(maxarity, arity);
            file >> i;
            file >> j;
            if (ToulBar2::verbose >= 3)
                cout << "read binary cost function " << ic << " on " << i << "," << j << endl;
            if (i == j) {
                cerr << "Error: binary cost function with only one variable in its scope!" << endl;
                throw WrongFileFormat();
            }
            file >> defval;
            if (defval >= MIN_COST) {
                assert(vars[i]->enumerated());
                assert(vars[j]->enumerated());
                EnumeratedVariable* x = (EnumeratedVariable*)vars[i];
                EnumeratedVariable* y = (EnumeratedVariable*)vars[j];
                file >> ntuples;
                if (ntuples < 0) {
                    int reusedconstr = -ntuples - 1;
                    if (reusedconstr >= (int)sharedSize.size()) {
                        cerr << "Shared cost function number " << reusedconstr << " not already defined! Cannot reuse it!" << endl;
                        throw WrongFileFormat();
                    }
                    ntuples = sharedSize[reusedconstr];
                    assert(ntuples == (int)(x->getDomainInitSize() * y->getDomainInitSize()));
                    if ((defval != MIN_COST) || (ntuples > 0))
                        postBinaryConstraint(i, j, sharedCosts[reusedconstr]);
                    continue;
                }
                vector<Cost> costs;
                for (a = 0; a < x->getDomainInitSize(); a++) {
                    for (b = 0; b < y->getDomainInitSize(); b++) {
                        Cost tmpcost = MULT(defval, K);
                        if (CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST * getUb()) && getUb() < (MAX_COST / MEDIUM_COST))
                            tmpcost *= MEDIUM_COST;
                        costs.push_back(tmpcost);
                    }
                }
                for (k = 0; k < ntuples; k++) {
                    file >> a;
                    file >> b;
                    file >> cost;
                    Cost tmpcost = MULT(cost, K);
                    if (CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST * getUb()) && getUb() < (MAX_COST / MEDIUM_COST))
                        tmpcost *= MEDIUM_COST;
                    assert(a >= 0 && a < x->getDomainInitSize());
                    assert(b >= 0 && b < y->getDomainInitSize());
                    costs[a * y->getDomainInitSize() + b] = tmpcost;
                }
                if (shared) {
                    sharedSize.push_back(costs.size());
                    sharedCosts.push_back(costs);
                    sharedTuples.push_back(emptyTuples);
                }
                if ((defval != MIN_COST) || (ntuples > 0))
                    postBinaryConstraint(i, j, costs);
            } else {
                file >> funcname;
                if (funcname == ">=") {
                    file >> funcparam1;
                    file >> funcparam2;
                    postSupxyc(i, j, funcparam1, funcparam2);
                } else if (funcname == ">") {
                    file >> funcparam1;
                    file >> funcparam2;
                    postSupxyc(i, j, funcparam1 + 1, funcparam2);
                } else if (funcname == "<=") {
                    file >> funcparam1;
                    file >> funcparam2;
                    postSupxyc(j, i, -funcparam1, funcparam2);
                } else if (funcname == "<") {
                    file >> funcparam1;
                    file >> funcparam2;
                    postSupxyc(j, i, -funcparam1 + 1, funcparam2);
                } else if (funcname == "=") {
                    file >> funcparam1;
                    file >> funcparam2;
                    postSupxyc(i, j, funcparam1, funcparam2);
                    postSupxyc(j, i, -funcparam1, funcparam2);
                } else if (funcname == "disj") {
                    Cost funcparam3;
                    file >> funcparam1;
                    file >> funcparam2;
                    file >> funcparam3;
                    postDisjunction(i, j, funcparam1, funcparam2, MULT(funcparam3, K));
                } else if (funcname == "sdisj") {
                    Value funcparam3;
                    Value funcparam4;
                    Cost funcparam5;
                    Cost funcparam6;
                    file >> funcparam1;
                    file >> funcparam2;
                    file >> funcparam3;
                    file >> funcparam4;
                    file >> funcparam5;
                    file >> funcparam6;
                    postSpecialDisjunction(i, j, funcparam1, funcparam2, funcparam3, funcparam4, MULT(funcparam5, K), MULT(funcparam6, K));
                } else {
                    int scopeIndex[2];
                    scopeIndex[0] = i;
                    scopeIndex[1] = j;
                    postGlobalFunction(scopeIndex, 2, funcname, file, &nbconstr);
                }
            }
            // ARITY 1
        } else if (arity == 1) {
            maxarity = max(maxarity, arity);
            file >> i;
            if (ToulBar2::verbose >= 3)
                cout << "read unary cost function " << ic << " on " << i << endl;
            if (vars[i]->enumerated()) {
                EnumeratedVariable* x = (EnumeratedVariable*)vars[i];
                file >> defval;
                if (defval == -1) {
                    int scopeIndex[1];
                    scopeIndex[0] = i;
                    string gcname;
                    file >> gcname;
                    postGlobalFunction(scopeIndex, 1, gcname, file, &nbconstr);
                } else {
                    file >> ntuples;
                    TemporaryUnaryConstraint unaryconstr;
                    unaryconstr.var = x;
                    if (ntuples < 0) {
                        int reusedconstr = -ntuples - 1;
                        if (reusedconstr >= (int)sharedSize.size()) {
                            cerr << "Shared cost function number " << reusedconstr << " not already defined! Cannot reuse it!" << endl;
                            throw WrongFileFormat();
                        }
                        ntuples = sharedSize[reusedconstr];
                        assert(ntuples == (int)x->getDomainInitSize());
                        unaryconstr.costs = sharedCosts[reusedconstr];
                        unaryconstrs.push_back(unaryconstr);
                        continue;
                    }
                    for (a = 0; a < x->getDomainInitSize(); a++) {
                        Cost tmpcost = MULT(defval, K);
                        if (CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST * getUb()) && getUb() < (MAX_COST / MEDIUM_COST))
                            tmpcost *= MEDIUM_COST;
                        unaryconstr.costs.push_back(tmpcost);
                    }
                    for (k = 0; k < ntuples; k++) {
                        file >> a;
                        file >> cost;
                        Cost tmpcost = MULT(cost, K);
                        if (CUT(tmpcost, getUb()) && (tmpcost < MEDIUM_COST * getUb()) && getUb() < (MAX_COST / MEDIUM_COST))
                            tmpcost *= MEDIUM_COST;
                        assert(a >= 0 && a < x->getDomainInitSize());
                        unaryconstr.costs[a] = tmpcost;
                    }
                    if (shared) {
                        sharedSize.push_back(x->getDomainInitSize());
                        sharedCosts.push_back(unaryconstr.costs);
                        sharedTuples.push_back(emptyTuples);
                    }
                    unaryconstrs.push_back(unaryconstr);
                }
            } else {
                file >> defval;
                if (defval == MIN_COST) {
                    cerr << "Error: unary cost function with zero penalty cost!" << endl;
                    throw WrongFileFormat();
                }
                file >> ntuples;
                Value* dom = new Value[ntuples];
                for (k = 0; k < ntuples; k++) {
                    file >> dom[k];
                    file >> cost;
                    if (cost != MIN_COST) {
                        cerr << "Error: unary cost function with non-zero cost tuple!" << endl;
                        throw WrongFileFormat();
                    }
                }
                postUnaryConstraint(i, dom, ntuples, defval);
                delete[] dom;
            }
            // ARITY 0
        } else if (arity == 0) {
            file >> defval;
            file >> ntuples;
            if (ToulBar2::verbose >= 3)
                cout << "read global lower bound contribution " << ic << " of " << defval << endl;
            if (ntuples > 1) {
                cerr << "Error: global lower bound contribution with several tuples!" << endl;
                throw WrongFileFormat();
            }
            if (ntuples == 1)
                file >> cost;
            else
                cost = defval;
            inclowerbound += MULT(cost, K);
        }
    }

    if ((getSolver() && ((Solver*)getSolver())->getWCSP() == this)) { // checks end of file only if reading the main model
        file >> funcname;
        if (file) {
            cerr << "Warning! EOF not reached after reading all the cost functions (initial number of cost functions too small?)" << endl;
        }
    }

    // merge unarycosts if they are on the same variable
    vector<int> seen(nbvar, -1);
    vector<TemporaryUnaryConstraint> newunaryconstrs;
    for (unsigned int u = 0; u < unaryconstrs.size(); u++) {
        if (seen[unaryconstrs[u].var->wcspIndex] == -1) {
            seen[unaryconstrs[u].var->wcspIndex] = newunaryconstrs.size();
            newunaryconstrs.push_back(unaryconstrs[u]);
        } else {
            for (unsigned int i = 0; i < unaryconstrs[u].var->getDomainInitSize(); i++) {
                if (newunaryconstrs[seen[unaryconstrs[u].var->wcspIndex]].costs[i] < getUb()) {
                    if (unaryconstrs[u].costs[i] < getUb())
                        newunaryconstrs[seen[unaryconstrs[u].var->wcspIndex]].costs[i] += unaryconstrs[u].costs[i];
                    else
                        newunaryconstrs[seen[unaryconstrs[u].var->wcspIndex]].costs[i] = getUb();
                }
            }
        }
    }
    unaryconstrs = newunaryconstrs;
    if (ToulBar2::sortDomains) {
        if (maxarity > 2) {
            cerr << "Error: cannot sort domains in preprocessing with non-binary cost functions." << endl;
            throw BadConfiguration();
        } else {
            ToulBar2::sortedDomains.clear();
            for (unsigned int u = 0; u < unaryconstrs.size(); u++) {
                ToulBar2::sortedDomains[unaryconstrs[u].var->wcspIndex] = unaryconstrs[u].var->sortDomain(unaryconstrs[u].costs);
            }
        }
    }

    // apply basic initial propagation AFTER complete network loading
    increaseLb(inclowerbound);

    // unary cost functions are delayed for compatibility issue (same lowerbound found) with old toolbar solver
    for (unsigned int u = 0; u < unaryconstrs.size(); u++) {
        postUnaryConstraint(unaryconstrs[u].var->wcspIndex, unaryconstrs[u].costs);
    }
    sortConstraints();

    if (ToulBar2::verbose >= 0)
        cout << "Read " << nbvar << " variables, with " << nbvaltrue << " values at most, and " << nbconstr << " cost functions, with maximum arity " << maxarity << "." << endl;
}

void WCSP::read_random(int n, int m, vector<int>& p, int seed, bool forceSubModular, string globalname)
{
    if (ToulBar2::externalUB.size()) {
        Cost top = string2Cost(ToulBar2::externalUB.c_str());
        Double K = ToulBar2::costMultiplier;
        if (top < MAX_COST / K)
            top = top * K;
        else
            top = MAX_COST;
        ToulBar2::deltaUb = max(ToulBar2::deltaUbAbsolute, (Cost)(ToulBar2::deltaUbRelativeGap * (Double)min(top, getUb())));
        updateUb(top + ToulBar2::deltaUb);
        // as long as a true certificate as not been found we must compensate for the deltaUb in CUT
    }
    naryRandom randwcsp(this, seed);
    randwcsp.Input(n, m, p, forceSubModular, globalname);

    unsigned int nbconstr = numberOfConstraints();
    sortConstraints();

    if (ToulBar2::verbose >= 0) {
        cout << "Generated random problem " << n << " variables, with " << m << " values, and " << nbconstr << " cost functions." << endl;
    }
}

void WCSP::read_uai2008(const char* fileName)
{
    // Compute the factor that enables to capture the difference in log for probability (1-10^resolution):
    if (ToulBar2::resolution == 0) {
        ToulBar2::NormFactor = 1.;
    } else {
        ToulBar2::NormFactor = (-1.0 / Log1p(-Exp10(-(TLogProb)ToulBar2::resolution)));
    }
    if (ToulBar2::NormFactor > (Pow((TProb)2., (TProb)INTEGERBITS) - 1) / (TLogProb)ToulBar2::resolution) {
        cerr << "This resolution cannot be ensured on the data type used to represent costs." << endl;
        throw BadConfiguration();
    } else if (ToulBar2::verbose >= 1) {
        cout << "NormFactor= " << ToulBar2::NormFactor << endl;
    }

    // Cost inclowerbound = MIN_COST;
    string uaitype;
    ifstream rfile(fileName, (ToulBar2::gz || ToulBar2::bz2 || ToulBar2::xz) ? (std::ios_base::in | std::ios_base::binary) : (std::ios_base::in));
#ifdef BOOST
    boost::iostreams::filtering_streambuf<boost::iostreams::input> zfile;
    if (ToulBar2::gz) {
        zfile.push(boost::iostreams::gzip_decompressor());
    } else if (ToulBar2::bz2) {
#ifndef NO_BZ2
        zfile.push(boost::iostreams::bzip2_decompressor());
#else
        cerr << "Error: Boost requires the bzip2 library to read bz2 compressed wcsp format files." << endl;
        throw WrongFileFormat();
#endif
    } else if (ToulBar2::xz) {
#ifndef NO_LZMA
        zfile.push(boost::iostreams::lzma_decompressor());
#else
        cerr << "Error: compiling with Boost version 1.65 or higher is needed to allow to read xz compressed uai/LG format files." << endl;
        throw WrongFileFormat();
#endif
    }
    zfile.push(rfile);
    istream ifile(&zfile);

    if (ToulBar2::stdin_format.length() == 0 && !rfile) {
        cerr << "Could not open uai file : " << fileName << endl;
        throw WrongFileFormat();
    }
    istream& file = (ToulBar2::stdin_format.length() > 0) ? cin : ifile;
#else
    if (ToulBar2::gz || ToulBar2::bz2 || ToulBar2::xz) {
        cerr << "Error: compiling with Boost iostreams library is needed to allow to read compressed uai format files." << endl;
        throw WrongFileFormat();
    }
    if (ToulBar2::stdin_format.length() == 0 && !rfile) {
        cerr << "Could not open uai file : " << fileName << endl;
        throw WrongFileFormat();
    }
    istream& file = (ToulBar2::stdin_format.length() > 0) ? cin : rfile;
#endif

    Cost inclowerbound = MIN_COST;
    updateUb((MAX_COST - UNIT_COST) / MEDIUM_COST / MEDIUM_COST / MEDIUM_COST / MEDIUM_COST);
    Cost upperbound = UNIT_COST;

    int nbval = 0;
    int nbvar, nbconstr;
    int i, j, k, ic;
    string varname;
    int domsize;
    Cost cost;
    Tuple s;
    int ntuples;
    int arity;
    int maxarity = 0;
    vector<TemporaryUnaryConstraint> unaryconstrs;
    vector<vector<int>> scopes;

    file >> uaitype;

    if (ToulBar2::verbose >= 3)
        cout << "Reading " << uaitype << "  file." << endl;

    bool markov = (uaitype == string("MARKOV"));
    // bool bayes = uaitype == string("BAYES");

    file >> nbvar;
    // read variable domain sizes
    for (i = 0; i < nbvar; i++) {
        string varname;
        varname = to_string("x") + to_string(i);
        file >> domsize;
        if (ToulBar2::verbose >= 1)
            cout << "read variable " << i << " of size " << domsize << endl;
        if (domsize > nbval)
            nbval = domsize;
        DEBONLY(int theindex =)
        ((domsize >= 0) ? makeEnumeratedVariable(varname, 0, domsize - 1) : makeIntervalVariable(varname, 0, -domsize - 1));
        assert(theindex == i);
    }

    file >> nbconstr;
    // read each constraint
    for (ic = 0; ic < nbconstr; ic++) {
        file >> arity;
        if (!file) {
            cerr << "Warning! EOF reached before reading all the scopes (initial number of factors too large?)" << endl;
            break;
        }
        maxarity = max(maxarity, arity);

        if (!file) {
            cerr << "Warning! EOF reached before reading all the cost functions (initial number of cost functions too large?)" << endl;
            break;
        }

        if (arity > 3) {
            vector<int> scopeIndex(arity, INT_MAX);
            if (ToulBar2::verbose >= 1)
                cout << "read nary cost function on ";

            for (i = 0; i < arity; i++) {
                file >> j;
                scopeIndex[i] = j;
                if (ToulBar2::verbose >= 1)
                    cout << j << " ";
            }
            if (ToulBar2::verbose >= 1)
                cout << endl;
            scopes.push_back(scopeIndex);
        } else if (arity == 3) {
            file >> i;
            file >> j;
            file >> k;
            if (ToulBar2::verbose >= 1)
                cout << "read ternary cost function " << ic << " on " << i << "," << j << "," << k << endl;
            scopes.push_back({ i, j, k });
            if ((i == j) || (i == k) || (k == j)) {
                cerr << "Error: ternary cost function!" << endl;
                throw WrongFileFormat();
            }
        } else if (arity == 2) {
            file >> i;
            file >> j;
            if (ToulBar2::verbose >= 1)
                cout << "read binary cost function " << ic << " on " << i << "," << j << endl;
            if (i == j) {
                cerr << "Error: binary cost function with only one variable in its scope!" << endl;
                throw WrongFileFormat();
            }
            scopes.push_back({ i, j });
        } else if (arity == 1) {
            file >> i;
            if (ToulBar2::verbose >= 1)
                cout << "read unary cost function " << ic << " on " << i << endl;
            TemporaryUnaryConstraint unaryconstr;
            unaryconstr.var = (EnumeratedVariable*)vars[i];
            unaryconstrs.push_back(unaryconstr);
            scopes.push_back({ i });
        } else if (arity == 0) {
            scopes.push_back({});
        }
        // prepare directed constraint graph used for finding a topological order
        if (!markov) {
            for (i = 0; i < arity - 1; i++) {
                getListSuccessors()->at(scopes.back()[i]).push_back(scopes.back().back());
            }
        }
    }

    ToulBar2::markov_log = 0; // for the MARKOV Case

    vector<int> ntuplesarray(scopes.size(), 0);
    vector<vector<Cost>> costs(scopes.size());
    int iunaryctr = 0;
    int ictr = 0;
    bool errorp = false;
    vector<vector<int>>::iterator it = scopes.begin();
    while (it != scopes.end()) {
        file >> ntuples;
        if (!file) {
            cerr << "Warning! EOF reached before reading all the factor tables (initial number of factors too large?)" << endl;
            break;
        }
        ntuplesarray[ictr] = ntuples;

        TProb p;
        vector<TProb> costsProb;

        TProb maxp = 0.;
        for (k = 0; k < ntuples; k++) {
            file >> p;
            if (ToulBar2::sigma > 0.) {
                Double noise = aleaGaussNoise(ToulBar2::sigma);
                if (ToulBar2::verbose >= 1)
                    cout << "add noise " << noise << " to " << p << endl;
                p = max((TProb)0., p + noise); // can create forbidden tuples
                if (!markov)
                    p = min(p, (TProb)1.); // in Bayesian networks probability cannot be greater than 1
            }
            assert(ToulBar2::uai > 1 || (p >= 0. && (markov || p <= 1.)));
            costsProb.push_back(p);
            maxp = max(maxp, p);
        }
        if (ToulBar2::uai == 1 && maxp == 0.)
            THROWCONTRADICTION;
        if (ToulBar2::uai == 2 && maxp < -1e38)
            THROWCONTRADICTION;
        // ensures conditional probabilities sum to 1 in the case of a Bayesian network
        if (!markov) {
            int domsize = getDomainInitSize((*it).back());
            for (k = 0; k < ntuples; k += domsize) {
                TProb cumul = 0.;
                int k2 = 0;
                for (; k2 < domsize - 1 && cumul + costsProb[k + k2] <= 1.; k2++) {
                    cumul += costsProb[k + k2];
                }
                assert(cumul <= 1.);
                for (; k2 < domsize ; k2++) {
                    if (!errorp && abs( cumul + costsProb[k + k2] - (TProb)1. ) > (TProb)ToulBar2::epsilon) {
                        cout << "Warning! Conditional probability table for variable " << getName((*it).back()) << " must have been normalized. Try option -epsilon with less precision to remove this warning (e.g., -epsilon=1e-4)." << endl;
                        errorp = true;
                    }
                    costsProb[k + k2] = 1. - cumul;
                    cumul = 1.;
                }
                // Do not need to update maxp because it is only used for Markov network!
            }
            assert(k == ntuples);
        }

        Cost minc = MAX_COST;
        Cost maxc = MIN_COST;
        for (k = 0; k < ntuples; k++) {
            p = costsProb[k];
            Cost cost;
            // ToulBar2::uai is 1 for .uai and 2 for .LG (log domain)
            if (markov) {
                if (ToulBar2::uai > 1) {
                    cost = LogProb2Cost((TLogProb)(p - maxp));
                } else {
                    cost = Prob2Cost(p / maxp);
                }
            } else {
                if (ToulBar2::uai > 1) {
                    cost = LogProb2Cost((TLogProb)p);
                } else {
                    cost = Prob2Cost(p);
                }
            }
            costs[ictr].push_back(cost);
            if (cost < minc)
                minc = cost;
            if (cost > maxc && cost < getUb())
                maxc = cost;
        }
        upperbound += maxc;

        if (ToulBar2::preprocessNary > 0 && minc > MIN_COST) {
            for (k = 0; k < ntuples; k++) {
                costs[ictr][k] -= minc;
            }
            if (ToulBar2::verbose >= 2)
                cout << "IC0 performed for cost function " << ictr << " with initial minimum cost " << minc << endl;
            inclowerbound += minc;
        }

        if (markov)
            ToulBar2::markov_log += ((ToulBar2::uai > 1) ? maxp : Log(maxp));

        ictr++;
        ++it;
    }

    file >> varname;
    if (file) {
        cerr << "Warning! EOF not reached after reading all the factor tables (initial number of factors too small?)" << endl;
    }

    updateUb(upperbound);

    ictr = 0;
    it = scopes.begin();
    while (it != scopes.end()) {
        ntuples = ntuplesarray[ictr];
        for (k = 0; k < ntuples; k++) {
            if (CUT(costs[ictr][k], getUb()))
                costs[ictr][k] = getUb() * MEDIUM_COST;
        }

        int arity = (*it).size();
        switch (arity) {
        case 0: {
            inclowerbound += costs[ictr][0];
            break;
        }

        case 1: {
            assert(unaryconstrs[iunaryctr].var->wcspIndex == scopes[ictr][0]);
            unaryconstrs[iunaryctr].costs.clear();
            for (unsigned int a = 0; a < unaryconstrs[iunaryctr].var->getDomainInitSize(); a++) {
                unaryconstrs[iunaryctr].costs.push_back(costs[ictr][a]);
            }
            iunaryctr++;
            if (ToulBar2::verbose >= 3)
                cout << "read unary costs." << endl;
            break;
        }

        case 2: {
            postBinaryConstraint(scopes[ictr][0], scopes[ictr][1], costs[ictr]);
            if (ToulBar2::verbose >= 3)
                cout << "read binary costs." << endl;
            break;
        }

        case 3: {
            postTernaryConstraint(scopes[ictr][0], scopes[ictr][1], scopes[ictr][2], costs[ictr]);
            if (ToulBar2::verbose >= 3)
                cout << "read ternary costs." << endl;
            break;
        }

        default: {
            int nctridx = postNaryConstraintBegin(scopes[ictr], MIN_COST, LONGLONG_MAX);
            NaryConstraint* nctr = (NaryConstraint*)getCtr(nctridx);
            assert(nctr->isNary());
            unsigned int j = 0;
            nctr->firstlex();
            while (nctr->nextlex(s, cost)) {
                //					  if (costs[j]>MIN_COST) nctr->setTuple(s, costs[j]);
                postNaryConstraintTupleInternal(nctr->wcspIndex, s, costs[ictr][j]);
                j++;
            }
            assert(j == costs[ictr].size());
            if (ToulBar2::verbose >= 3)
                cout << "read arity " << arity << " table costs." << endl;
            postNaryConstraintEnd(nctr->wcspIndex);
            break;
        }
        }
        ictr++;
        ++it;
    }
    if (ToulBar2::verbose >= 1) {
        cout << "MarkovShiftingValue= " << ToulBar2::markov_log << endl;
    }

    // apply basic initial propagation AFTER complete network loading
    increaseLb(inclowerbound);

    for (unsigned int u = 0; u < unaryconstrs.size(); u++) {
        postUnaryConstraint(unaryconstrs[u].var->wcspIndex, unaryconstrs[u].costs);
    }
    sortConstraints();
    if (ToulBar2::verbose >= 0)
        cout << "Read " << nbvar << " variables, with " << nbval << " values at most, and " << nbconstr << " cost functions, with maximum arity " << maxarity << "." << endl;

    int nevi = 0;
    ifstream fevid(ToulBar2::evidence_file.c_str());
    if (!fevid) {
        string tmpname(fileName);
        // skip compressed extensions like .xz or .gz
        if (tmpname.find(".uai") != std::string::npos) {
            tmpname.replace(tmpname.find(".uai"), tmpname.length(), ".uai");
        } else if (tmpname.find(".LG") != std::string::npos) {
            tmpname.replace(tmpname.find(".LG"), tmpname.length(), ".LG");
        }
        if (!(ToulBar2::uaieval)) {
            string strevid(string(tmpname) + string(".evid"));
            fevid.open(strevid.c_str());
            if (ToulBar2::verbose >= 0)
                cout << "No evidence file specified. Trying " << strevid << endl;
        }
        if (!fevid)
            if (ToulBar2::verbose >= 0)
                cout << "No evidence file. " << endl;
    }
    if (fevid) {
        vector<int> variables;
        vector<Value> values;
        fevid >> nevi;
        bool firstevid = true;
        if (nevi == 0)
            return;
        if (nevi == 1)
            fevid >> nevi; // UAI 2010 evidence file format assumes possible multiple evidence samples, but toulbar2 will search for the first evidence sample only!
        while (nevi) {
            if (!fevid) {
                cerr << "Error: incorrect number of evidences." << endl;
                throw WrongFileFormat();
            }
            fevid >> i;
            fevid >> j;
            if (firstevid && !fevid) { // old UAI 2008 evidence format
                variables.push_back(nevi);
                values.push_back(i);
                break;
            } else
                firstevid = false;
            variables.push_back(i);
            values.push_back(j);
            nevi--;
        }
        assignLS(variables, values);
    }
}

void WCSP::solution_UAI(Cost res)
{
    if (!ToulBar2::uai && !ToulBar2::uaieval)
        return;
    if (ToulBar2::isZ)
        return;
    if (((ToulBar2::writeSolution) ? ToulBar2::solutionFile : ToulBar2::solution_uai_file) == NULL)
        return;
    // UAI 2012 Challenge output format
    //	    ToulBar2::solution_file << "-BEGIN-" << endl;
    rewind((ToulBar2::writeSolution) ? ToulBar2::solutionFile : ToulBar2::solution_uai_file);
    fprintf((ToulBar2::writeSolution) ? ToulBar2::solutionFile : ToulBar2::solution_uai_file, "MPE\n");
    //	ToulBar2::solution_file << "1" << endl; // we assume a single evidence sample
    //    if (ToulBar2::showSolutions && !ToulBar2::uaieval) {
    //        cout << "t " << cpuTime() - ToulBar2::startCpuTime << endl;
    //        cout << "s " << -(Cost2LogProb(res) + ToulBar2::markov_log) << endl;
    //        cout << numberOfVariables();
    //        printSolution(cout);
    //    }
    fprintf((ToulBar2::writeSolution) ? ToulBar2::solutionFile : ToulBar2::solution_uai_file, "%d ", numberOfVariables() - numberOfHiddenVariables());
    printSolution((ToulBar2::writeSolution) ? ToulBar2::solutionFile : ToulBar2::solution_uai_file);
    fprintf((ToulBar2::writeSolution) ? ToulBar2::solutionFile : ToulBar2::solution_uai_file, "\n");
    //	if (opt) {
    //	  if (ToulBar2::showSolutions) cout << " LU" << endl;
    //	  ToulBar2::solution_file << " LU" << endl;
    //	} else {
    //	  if (ToulBar2::showSolutions) cout << " L" << endl;
    //	  ToulBar2::solution_file << " L" << endl;
    //	}
    fflush((ToulBar2::writeSolution) ? ToulBar2::solutionFile : ToulBar2::solution_uai_file);
}

#if defined(XMLFLAG)
#include "./xmlcsp/xmlcsp.h"
#elif defined(XMLFLAG3)
#include "./xmlcsp3/xmlcsp3.h"
#endif

void WCSP::read_XML(const char* fileName)
{
#if defined(XMLFLAG) || defined(XMLFLAG3)
    ifstream rfile(fileName, (ToulBar2::gz || ToulBar2::bz2 || ToulBar2::xz) ? (std::ios_base::in | std::ios_base::binary) : (std::ios_base::in));
#ifdef BOOST
    boost::iostreams::filtering_streambuf<boost::iostreams::input> zfile;
    if (ToulBar2::gz) {
        zfile.push(boost::iostreams::gzip_decompressor());
    } else if (ToulBar2::bz2) {
#ifndef NO_BZ2
        zfile.push(boost::iostreams::bzip2_decompressor());
#else
        cerr << "Error: Boost requires the bzip2 library to read bz2 compressed wcsp format files." << endl;
        throw WrongFileFormat();
#endif
    } else if (ToulBar2::xz) {
#ifndef NO_LZMA
        zfile.push(boost::iostreams::lzma_decompressor());
#else
        cerr << "Error: compiling with Boost version 1.65 or higher is needed to allow to read xz compressed xml format files." << endl;
        throw WrongFileFormat();
#endif
    }
    zfile.push(rfile);
    istream ifile(&zfile);

    if (ToulBar2::stdin_format.length() == 0 && !rfile) {
        cerr << "Could not open xml file : " << fileName << endl;
        throw WrongFileFormat();
    }
    istream& file = (ToulBar2::stdin_format.length() > 0) ? cin : ifile;
#else
    if (ToulBar2::gz || ToulBar2::bz2 || ToulBar2::xz) {
        cerr << "Error: compiling with Boost iostreams library is needed to allow to read compressed xml format files." << endl;
        throw WrongFileFormat();
    }
    if (ToulBar2::stdin_format.length() == 0 && !rfile) {
        cerr << "Could not open xml file : " << fileName << endl;
        throw WrongFileFormat();
    }
    istream& file = (ToulBar2::stdin_format.length() > 0) ? cin : rfile;
#endif
#endif
#if defined(XMLFLAG)
    MyCallback xmlCallBack;
    xmlCallBack.wcsp = this;
    xmlCallBack.fname = string(fileName);
    xmlCallBack.convertWCSP = true;
    try {
        XMLParser_libxml2<> parser(xmlCallBack);
        parser.setPreferredExpressionRepresentation(INFIX_C);
        parser.parse(file);
    } catch (exception& e) {
        cout.flush();
        cerr << "\n\tUnexpected exception in XML parsing\n";
        cerr << "\t" << e.what() << endl;
        throw WrongFileFormat();
    }
#elif defined(XMLFLAG3)
    MySolverCallbacks xmlCallBack; // my interface between the parser and the solver
    xmlCallBack.problem = this;
    try {
        XCSP3CoreParser parser(&xmlCallBack);
        parser.parse(file);
    } catch (exception& e) {
        cout.flush();
        cerr << "\n\tUnexpected exception in XML XCSP3 parsing\n";
        cerr << "\t" << e.what() << endl;
        throw WrongFileFormat();
    }
#else
    cerr << "\nXML format without including in Makefile flag XMLFLAG and files ./xmlcsp\n"
         << endl;
    throw WrongFileFormat();
#endif
}

void WCSP::solution_XML(bool opt)
{
#if defined(XMLFLAG) || defined(XMLFLAG3)
    if (!ToulBar2::xmlflag)
        return;

    if (opt && ToulBar2::xmlcop) {
        cout << "s OPTIMUM FOUND" << endl;
    } else {
        cout << "s SATISFIABLE" << endl;
    }

    // ofstream fsol;
    // ifstream sol;
    // sol.open(ToulBar2::writeSolution);
    // if(!sol) { cout << "cannot open solution file to translate" << endl; throw WrongFileFormat(); }
    // fsol.open("solution");
    // fsol << "SOL ";

    //    freopen(NULL, "r", ToulBar2::solutionFile);
    cout << "v ";
#ifdef XMLFLAG3
    if (ToulBar2::xmlcop) {
        cout << "<instantiation type=\"" << ((opt) ? "optimum" : "solution") << "\" cost=\"" << std::fixed << std::setprecision(0) << getSolutionValue() << std::setprecision(DECIMAL_POINT) << "\"> <list>";
    } else {
        cout << "<instantiation type=\"solution\"> <list>";
    }
    for (unsigned int i = 0; i < vars.size(); i++) {
        if (getName(i).rfind(IMPLICIT_VAR_TAG, 0) != 0 && getName(i).rfind(HIDDEN_VAR_TAG, 0) != 0) {
            cout << " " << getName(i);
        }
    }
    cout << " </list> <values>";
#endif
    vector<Value> lastsol = getSolution();
    for (unsigned int i = 0; i < vars.size(); i++) {
        int value = lastsol[i];
//        fscanf(ToulBar2::solutionFile, "%d", &value);
#ifdef XMLFLAG
        int index = ((EnumeratedVariable*)getVar(i))->toIndex(value);
        cout << Doms[varsDom[i]][index] << " ";
#else // XCSP3
        if (getName(i).rfind(IMPLICIT_VAR_TAG, 0) != 0 && getName(i).rfind(HIDDEN_VAR_TAG, 0) != 0) {
            cout << " " << value;
        }
#endif
    }
#ifdef XMLFLAG3
    cout << " </values> </instantiation>";
#endif
    cout << endl;
//    freopen(NULL, "w", ToulBar2::solutionFile);

// fsol << endl;
// fsol.close();
// sol.close();
#endif
}

void WCSP::read_wcnf(const char* fileName)
{
    ifstream rfile(fileName, (ToulBar2::gz || ToulBar2::bz2 || ToulBar2::xz) ? (std::ios_base::in | std::ios_base::binary) : (std::ios_base::in));
#ifdef BOOST
    boost::iostreams::filtering_streambuf<boost::iostreams::input> zfile;
    if (ToulBar2::gz) {
        zfile.push(boost::iostreams::gzip_decompressor());
    } else if (ToulBar2::bz2) {
#ifndef NO_BZ2
        zfile.push(boost::iostreams::bzip2_decompressor());
#else
        cerr << "Error: Boost requires the bzip2 library to read bz2 compressed wcsp format files." << endl;
        throw WrongFileFormat();
#endif
    } else if (ToulBar2::xz) {
#ifndef NO_LZMA
        zfile.push(boost::iostreams::lzma_decompressor());
#else
        cerr << "Error: compiling with Boost version 1.65 or higher is needed to allow to read xz compressed cnf/wcnf format files." << endl;
        throw WrongFileFormat();
#endif
    }
    zfile.push(rfile);
    istream ifile(&zfile);

    if (ToulBar2::stdin_format.length() == 0 && !rfile) {
        cerr << "Could not open wcnf file : " << fileName << endl;
        throw WrongFileFormat();
    }
    istream& file = (ToulBar2::stdin_format.length() > 0) ? cin : ifile;
#else
    if (ToulBar2::gz || ToulBar2::bz2 || ToulBar2::xz) {
        cerr << "Error: compiling with Boost iostreams library is needed to allow to read compressed wcnf format files." << endl;
        throw WrongFileFormat();
    }
    if (ToulBar2::stdin_format.length() == 0 && !rfile) {
        cerr << "Could not open wcnf file : " << fileName << endl;
        throw WrongFileFormat();
    }
    istream& file = (ToulBar2::stdin_format.length() > 0) ? cin : rfile;
#endif

    Double K = ToulBar2::costMultiplier;
    Cost inclowerbound = MIN_COST;
    updateUb((MAX_COST - UNIT_COST) / MEDIUM_COST / MEDIUM_COST);

    int maxarity = 0;
    vector<TemporaryUnaryConstraint> unaryconstrs;

    int nbvar, nbclauses;
    string dummy, sflag;

    file >> sflag;
    while (sflag[0] == 'c') {
        getline(file, dummy);
        file >> sflag;
    }
    if (sflag != "p") {
        cerr << "Wrong wcnf format in " << fileName << endl;
        throw WrongFileFormat();
    }

    string format, strtop;
    Cost top;
    file >> format;
    file >> nbvar;
    file >> nbclauses;
    if (format == "wcnf") {
        getline(file, strtop);
        if (strtop.size() > 0 && string2Cost((char*)strtop.c_str()) > 0) {
            if (ToulBar2::verbose >= 0)
                cout << "c (Weighted) Partial Max-SAT input format" << endl;
            top = string2Cost((char*)strtop.c_str());
            if (top < MAX_COST / K)
                top = top * K;
            else
                top = MAX_COST;
            ToulBar2::deltaUb = max(ToulBar2::deltaUbAbsolute, (Cost)(ToulBar2::deltaUbRelativeGap * (Double)min(top, getUb())));
            updateUb(top + ToulBar2::deltaUb);
        } else {
            if (ToulBar2::verbose >= 0)
                cout << "c Weighted Max-SAT input format" << endl;
        }
    } else {
        if (ToulBar2::verbose >= 0)
            cout << "c Max-SAT input format" << endl;
        Cost top = (nbclauses + 1) * K;
        ToulBar2::deltaUb = max(ToulBar2::deltaUbAbsolute, (Cost)(ToulBar2::deltaUbRelativeGap * (Double)min(top, getUb())));
        updateUb(top + ToulBar2::deltaUb);
    }

    // create Boolean variables
    for (int i = 0; i < nbvar; i++) {
        string varname;
        varname = to_string("x") + to_string(i);
        DEBONLY(int theindex =)
        makeEnumeratedVariable(varname, 0, 1);
        assert(theindex == i);
    }

    // Read each clause
    Tuple tup;
    for (int ic = 0; ic < nbclauses; ic++) {
        vector<int> scopeIndex;
        tup.clear();
        int arity = 0;
        if (ToulBar2::verbose >= 3)
            cout << "read clause on ";
        int j = 0;
        Cost cost = UNIT_COST;
        if (format == "wcnf")
            file >> cost;
        bool tautology = false;
        do {
            file >> j;
            if (j != 0 && !tautology) {
                scopeIndex.resize(arity + 1);
                scopeIndex[arity] = abs(j) - 1;
                if (arity < (int)tup.size()) {
                    tup[arity] = ((j > 0) ? 0 : 1);
                } else {
                    tup.push_back((j > 0) ? 0 : 1);
                }
                int k = 0;
                while (k < arity) {
                    if (scopeIndex[k] == scopeIndex[arity]) {
                        break;
                    }
                    k++;
                }
                if (k < arity) {
                    if (tup[k] != tup[arity]) {
                        tautology = true;
                        if (ToulBar2::verbose >= 3)
                            cout << j << " is a tautology! skipped.";
                    } else {
                        tup.pop_back();
                        scopeIndex.resize(arity);
                    }
                    continue;
                }
                arity++;
                if (ToulBar2::verbose >= 3)
                    cout << j << " ";
            }
        } while (j != 0);
        if (ToulBar2::verbose >= 3)
            cout << endl;
        if (tautology)
            continue;

        maxarity = max(maxarity, arity);

        if (arity > 3) {
            // #ifdef CLAUSE2KNAPSACK
            //             if (CUT(MULT(cost, K), getUb())) {
            //                 postKnapsackConstraint(scopeIndex,arity,file,false,false,false,tup);
            //             } else {
            // #endif
            int index = postNaryConstraintBegin(scopeIndex, MIN_COST, 1);
            postNaryConstraintTupleInternal(index, tup, MULT(cost, K));
            postNaryConstraintEnd(index);
            // #ifdef CLAUSE2KNAPSACK
            //             }
            // #endif
        } else if (arity == 3) {
            vector<Cost> costs;
            for (int a = 0; a < 2; a++) {
                for (int b = 0; b < 2; b++) {
                    for (int c = 0; c < 2; c++) {
                        costs.push_back(MIN_COST);
                    }
                }
            }
            costs[(tup[0]) * 4 + (tup[1]) * 2 + (tup[2])] = MULT(cost, K);
            postTernaryConstraint(scopeIndex[0], scopeIndex[1], scopeIndex[2], costs);
        } else if (arity == 2) {
            vector<Cost> costs;
            for (int a = 0; a < 2; a++) {
                for (int b = 0; b < 2; b++) {
                    costs.push_back(MIN_COST);
                }
            }
            costs[(tup[0]) * 2 + (tup[1])] = MULT(cost, K);
            postBinaryConstraint(scopeIndex[0], scopeIndex[1], costs);
        } else if (arity == 1) {
            EnumeratedVariable* x = (EnumeratedVariable*)vars[scopeIndex[0]];
            TemporaryUnaryConstraint unaryconstr;
            unaryconstr.var = x;
            if ((tup[0]) == 0) {
                unaryconstr.costs.push_back(MULT(cost, K));
                unaryconstr.costs.push_back(MIN_COST);
            } else {
                unaryconstr.costs.push_back(MIN_COST);
                unaryconstr.costs.push_back(MULT(cost, K));
            }
            unaryconstrs.push_back(unaryconstr);
        } else if (arity == 0) {
            inclowerbound += MULT(cost, K);
        } else {
            cerr << "Wrong clause arity " << arity << " in " << fileName << endl;
            throw WrongFileFormat();
        }
    }

    file >> dummy;
    if (file) {
        cerr << "Warning! EOF not reached after reading all the clauses (initial number of clauses too small?)" << endl;
    }

    // apply basic initial propagation AFTER complete network loading
    increaseLb(inclowerbound);

    for (unsigned int u = 0; u < unaryconstrs.size(); u++) {
        postUnaryConstraint(unaryconstrs[u].var->wcspIndex, unaryconstrs[u].costs);
    }
    sortConstraints();
    if (ToulBar2::verbose >= 0)
        cout << "c Read " << nbvar << " variables, with 2 values at most, and " << nbclauses << " clauses, with maximum arity " << maxarity << "." << endl;
}

/// \brief minimizes/maximizes \f$ X^t \times W \times X = \sum_{i=1}^N \sum_{j=1}^N W_{ij} \times X_i \times X_j \f$
/// where W is expressed by its M non-zero triangle matrix terms (W_ij, i<=j, it can be positive or negative float numbers)
/// Note : Quadratic terms for \f$ i < j \f$ are multiplied by 2 (see option -qpmult to change this value) to get a symmetric N*N squared matrix
/// Note : If N is positive, then variable domain values are {0,1}
/// Note : If N is negative, then variable domain values are {1,-1} with value 1 having index 0 and value -1 having index 1 in the output solutions
/// Note : If M is positive then minimizes the quadratic objective function, else maximizes it
/// Warning : It does not allow infinite costs (no forbidden assignments)
void WCSP::read_qpbo(const char* fileName)
{
    ifstream rfile(fileName, (ToulBar2::gz || ToulBar2::bz2 || ToulBar2::xz) ? (std::ios_base::in | std::ios_base::binary) : (std::ios_base::in));
#ifdef BOOST
    boost::iostreams::filtering_streambuf<boost::iostreams::input> zfile;
    if (ToulBar2::gz) {
        zfile.push(boost::iostreams::gzip_decompressor());
    } else if (ToulBar2::bz2) {
#ifndef NO_BZ2
        zfile.push(boost::iostreams::bzip2_decompressor());
#else
        cerr << "Error: Boost requires the bzip2 library to read bz2 compressed wcsp format files." << endl;
        throw WrongFileFormat();
#endif
    } else if (ToulBar2::xz) {
#ifndef NO_LZMA
        zfile.push(boost::iostreams::lzma_decompressor());
#else
        cerr << "Error: compiling with Boost version 1.65 or higher is needed to allow to read xz compressed qpbo format files." << endl;
        throw WrongFileFormat();
#endif
    }
    zfile.push(rfile);
    istream ifile(&zfile);

    if (ToulBar2::stdin_format.length() == 0 && !rfile) {
        cerr << "Could not open qpbo file : " << fileName << endl;
        throw WrongFileFormat();
    }
    istream& file = (ToulBar2::stdin_format.length() > 0) ? cin : ifile;
#else
    if (ToulBar2::gz || ToulBar2::bz2 || ToulBar2::xz) {
        cerr << "Error: compiling with Boost iostreams library is needed to allow to read compressed qpbo format files." << endl;
        throw WrongFileFormat();
    }
    if (ToulBar2::stdin_format.length() == 0 && !rfile) {
        cerr << "Could not open qpbo file : " << fileName << endl;
        throw WrongFileFormat();
    }
    istream& file = (ToulBar2::stdin_format.length() > 0) ? cin : rfile;
#endif

    int n = 0;
    file >> n;
    bool booldom = (n >= 0); // n positive means variable domains {0,1} else {1,-1}
    if (!booldom)
        n = -n;
    int m = 0;
    file >> m;
    if (n == 0 || m == 0)
        return;
    bool minimize = (m >= 0); // m positive means minimize the quadratic function, else maximize it
    if (!minimize)
        m = -m;
    int e = 0;
    int dummy;

    vector<int> posx(m, 0);
    vector<int> posy(m, 0);
    vector<Double> cost(m, 0.);
    for (e = 0; e < m; e++) {
        file >> posx[e];

        if (!file) {
            cerr << "Warning! EOF reached before reading all the cost sparse matrix (number of nonzero costs too large?)" << endl;
            break;
        }
        if (posx[e] > n) {
            cerr << "Warning! Variable index too large!" << endl;
            break;
        }
        file >> posy[e];
        if (posy[e] > n) {
            cerr << "Warning! Variable index too large!" << endl;
            break;
        }
        file >> cost[e];
    }
    file >> dummy;
    if (file) {
        cerr << "Warning! EOF not reached after reading all the cost sparse matrix (wrong number of nonzero costs too small?)" << endl;
    }
    m = e;

    // create Boolean variables
    for (int i = 0; i < n; i++) {
        makeEnumeratedVariable(to_string("x") + to_string(i), 0, 1);
    }

    vector<Cost> unaryCosts0(n, 0);
    vector<Cost> unaryCosts1(n, 0);

    // find total cost
    Double sumcost = 0.;
    for (int e = 0; e < m; e++) {
        sumcost += 2. * std::abs(cost[e]);
    }
    Double multiplier = Exp10((Double)ToulBar2::resolution); // warning! precision of costs depends on the chosen resolution but toulbar2 outputs costs as if it was always zero resolution.
    ToulBar2::setCostMultiplier(multiplier);
    if (!minimize)
        ToulBar2::setCostMultiplier(ToulBar2::costMultiplier * -1.0);
    if (multiplier * sumcost >= (Double)MAX_COST) {
        cerr << "This resolution cannot be ensured on the data type used to represent costs! (see option -precision)" << endl;
        throw BadConfiguration();
    }
    Cost top = (Cost)roundl(multiplier * sumcost) + UNIT_COST;
    ToulBar2::deltaUb = max(ToulBar2::deltaUbAbsolute, (Cost)(ToulBar2::deltaUbRelativeGap * (Double)min(top, getUb())));
    updateUb(top + ToulBar2::deltaUb);

    // create weighted binary clauses
    for (int e = 0; e < m; e++) {
        if (posx[e] != posy[e]) {
            vector<Cost> costs(4, 0);
            if (booldom) {
                if (cost[e] > 0) {
                    if (minimize) {
                        costs[3] = (Cost)roundl(multiplier * ToulBar2::qpboQuadraticCoefMultiplier * cost[e]);
                    } else {
                        costs[0] = (Cost)roundl(multiplier * ToulBar2::qpboQuadraticCoefMultiplier * cost[e]);
                        costs[1] = costs[0];
                        costs[2] = costs[0];
                        negCost += costs[0];
                    }
                } else {
                    if (minimize) {
                        costs[0] = (Cost)roundl(multiplier * ToulBar2::qpboQuadraticCoefMultiplier * -cost[e]);
                        costs[1] = costs[0];
                        costs[2] = costs[0];
                        negCost += costs[0];
                    } else {
                        costs[3] = (Cost)roundl(multiplier * ToulBar2::qpboQuadraticCoefMultiplier * -cost[e]);
                    }
                }
            } else {
                if (cost[e] > 0) {
                    if (minimize) {
                        costs[0] = (Cost)roundl(multiplier * ToulBar2::qpboQuadraticCoefMultiplier * 2. * cost[e]);
                        costs[3] = costs[0];
                        negCost += (Cost)roundl(multiplier * ToulBar2::qpboQuadraticCoefMultiplier * cost[e]);
                    } else {
                        costs[1] = (Cost)roundl(multiplier * ToulBar2::qpboQuadraticCoefMultiplier * 2. * cost[e]);
                        costs[2] = costs[1];
                        negCost += (Cost)roundl(multiplier * ToulBar2::qpboQuadraticCoefMultiplier * cost[e]);
                    }
                } else {
                    if (minimize) {
                        costs[1] = (Cost)roundl(multiplier * ToulBar2::qpboQuadraticCoefMultiplier * -2. * cost[e]);
                        costs[2] = costs[1];
                        negCost += (Cost)roundl(multiplier * ToulBar2::qpboQuadraticCoefMultiplier * -cost[e]);
                    } else {
                        costs[0] = (Cost)roundl(multiplier * ToulBar2::qpboQuadraticCoefMultiplier * -2. * cost[e]);
                        costs[3] = costs[0];
                        negCost += (Cost)roundl(multiplier * ToulBar2::qpboQuadraticCoefMultiplier * -cost[e]);
                    }
                }
            }
            postBinaryConstraint(posx[e] - 1, posy[e] - 1, costs);
        } else {
            if (booldom) {
                if (cost[e] > 0) {
                    if (minimize) {
                        unaryCosts1[posx[e] - 1] += (Cost)roundl(multiplier * cost[e]);
                    } else {
                        unaryCosts0[posx[e] - 1] += (Cost)roundl(multiplier * cost[e]);
                        negCost += (Cost)roundl(multiplier * cost[e]);
                    }
                } else {
                    if (minimize) {
                        unaryCosts0[posx[e] - 1] += (Cost)roundl(multiplier * -cost[e]);
                        negCost += (Cost)roundl(multiplier * -cost[e]);
                    } else {
                        unaryCosts1[posx[e] - 1] += (Cost)roundl(multiplier * -cost[e]);
                    }
                }
            } else {
                if (cost[e] > 0) {
                    if (minimize) {
                        unaryCosts0[posx[e] - 1] += (Cost)roundl(multiplier * 2. * cost[e]);
                        negCost += (Cost)roundl(multiplier * cost[e]);
                    } else {
                        unaryCosts1[posx[e] - 1] += (Cost)roundl(multiplier * 2. * cost[e]);
                        negCost += (Cost)roundl(multiplier * cost[e]);
                    }
                } else {
                    if (minimize) {
                        unaryCosts1[posx[e] - 1] += (Cost)roundl(multiplier * -2. * cost[e]);
                        negCost += (Cost)roundl(multiplier * -cost[e]);
                    } else {
                        unaryCosts0[posx[e] - 1] += (Cost)roundl(multiplier * -2. * cost[e]);
                        negCost += (Cost)roundl(multiplier * -cost[e]);
                    }
                }
            }
        }
    }

    // create weighted unary clauses
    for (int i = 0; i < n; i++) {
        if (unaryCosts0[i] > 0 || unaryCosts1[i] > 0) {
            vector<Cost> costs(2, 0);
            costs[0] = unaryCosts0[i];
            costs[1] = unaryCosts1[i];
            postUnaryConstraint(i, costs);
        }
    }
    sortConstraints();
    if (ToulBar2::verbose >= 0) {
        cout << "Read " << n << " variables, with " << 2 << " values at most, and " << m << " nonzero matrix costs (quadratic coef. multiplier: " << ToulBar2::qpboQuadraticCoefMultiplier << ", shifting value: " << -negCost << ")" << endl;
    }
}

bool isInteger(string& s) { return string("0123456789+-").find(s[0]) != string::npos; }
/// \param file: input file
/// \param token: in: previous token, out: new token (read from file or from the end of the previous token)
/// \param keep: in: relative position to start reading from the previous token (if positive and greater or equal to previous token size then reads from file else if negative subtracts from the end), out: size of the new token
/// Warning : if new token is + or - then replace to +1 or -1
/// Warning : if new token is +varname or -varname then split into +1 varname or -1 varname
void readToken(istream& file, string& token, int* keep = NULL)
{
    if (keep == NULL || *keep >= (int)token.size()) {
        file >> token;
    } else if (*keep >= 0) {
        token = token.substr(*keep);
    } else {
        assert(-(*keep) <= (int)token.size());
        token = token.substr(token.size() + (*keep));
    }
    bool twotokens = false;
    if (token == "+" || token == "-") {
        token = token + "1";
    } else if (token.size() >= 2 && (token[0] == '+' || token[0] == '-') && string("0123456789").find(token[1]) == string::npos) {
        twotokens = true;
        token = to_string(token[0]) + to_string("1") + token.substr(1);
    }
    if (keep) {
        if (twotokens) {
            *keep = 2;
        } else {
            *keep = token.size();
        }
    }
    if (ToulBar2::verbose >= 8) {
        cout << "##" << token << "##" << endl;
    }
}

void WCSP::read_opb(const char* fileName)
{
    ifstream rfile(fileName, (ToulBar2::gz || ToulBar2::bz2 || ToulBar2::xz) ? (std::ios_base::in | std::ios_base::binary) : (std::ios_base::in));
#ifdef BOOST
    boost::iostreams::filtering_streambuf<boost::iostreams::input> zfile;
    if (ToulBar2::gz) {
        zfile.push(boost::iostreams::gzip_decompressor());
    } else if (ToulBar2::bz2) {
#ifndef NO_BZ2
        zfile.push(boost::iostreams::bzip2_decompressor());
#else
        cerr << "Error: Boost requires the bzip2 library to read bz2 compressed opb format files." << endl;
        throw WrongFileFormat();
#endif
    } else if (ToulBar2::xz) {
#ifndef NO_LZMA
        zfile.push(boost::iostreams::lzma_decompressor());
#else
        cerr << "Error: compiling with Boost version 1.65 or higher is needed to allow to read xz compressed opb format files." << endl;
        throw WrongFileFormat();
#endif
    }
    zfile.push(rfile);
    istream ifile(&zfile);

    if (ToulBar2::stdin_format.length() == 0 && !rfile) {
        cerr << "Could not open opb file : " << fileName << endl;
        throw WrongFileFormat();
    }
    istream& file = (ToulBar2::stdin_format.length() > 0) ? cin : ifile;
#else
    if (ToulBar2::gz || ToulBar2::bz2 || ToulBar2::xz) {
        cerr << "Error: compiling with Boost iostreams library is needed to allow to read compressed opb format files." << endl;
        throw WrongFileFormat();
    }
    if (ToulBar2::stdin_format.length() == 0 && !rfile) {
        cerr << "Could not open opb file : " << fileName << endl;
        throw WrongFileFormat();
    }
    istream& file = (ToulBar2::stdin_format.length() > 0) ? cin : rfile;
#endif

    Cost inclowerbound = MIN_COST;
    updateUb((MAX_COST - UNIT_COST) / MEDIUM_COST / MEDIUM_COST);

    int maxarity = 0;
    int nbvar = 0;
    int nblinear = 0;
    int nbsoft = 0;
    int cursoft = 0;
    vector<TemporaryUnaryConstraint> unaryconstrs;

    map<string, int> varnames;
    string dummy, token;
    std::streampos prev;

    // skip initial comments
    readToken(file, token);
    while (token[0] == '*') {
        getline(file, dummy);
        if (dummy.find("#variable=") != string::npos) {
            int n = atoi(dummy.substr(dummy.find("#variable=") + 10).c_str());
            if (n > 0) {
                for (int i = 1; i <= n; i++) {
                    string varname = "x";
                    varname += to_string(i);
                    int var = makeEnumeratedVariable(varname, 0, 1);
                    addValueName(var, "v0");
                    addValueName(var, "v1");
                    varnames[varname] = var;
                    nbvar++;
                }
            }
        }
        if (dummy.find("#soft=") != string::npos) {
            int n = atoi(dummy.substr(dummy.find("#soft=") + 6).c_str());
            if (n > 0) {
                for (int i = 1; i <= n; i++) {
                    string varname = "s";
                    varname += to_string(i);
                    int var = makeEnumeratedVariable(varname, 0, 1);
                    addValueName(var, "v0");
                    addValueName(var, "v1");
                    varnames[varname] = var;
                    nbsoft++;
                }
            }
        }
        readToken(file, token);
    }

    // read objective function
    bool opt = true;
    int opsize = 4;
    Double multiplier = Exp10((Double)ToulBar2::resolution);
    ToulBar2::setCostMultiplier(multiplier);
    if (token.substr(0, 4) == "max:") {
        ToulBar2::setCostMultiplier(ToulBar2::costMultiplier * -1.0);
    } else if (token.substr(0, 5) == "soft:") {
        opt = false;
        opsize = 0;
        multiplier = 1.;
        ToulBar2::setCostMultiplier(1.);
        if (token.back() == ':') {
            readToken(file, token);
        }
        Cost top = string2Cost(token.substr((token[0]=='s')?5:0).c_str());
        if (top > MIN_COST) {
            updateUb(top);
        }
        if (token.back() != ';') {
            readToken(file, token);
        }
        assert(token.back() == ';');
    } else if (token.substr(0, 4) != "min:") {
        opt = false;
        opsize = 0;
        updateUb(UNIT_COST * multiplier);
    }

    if (opt) {
        do {
            Cost cost; // cost can be negative or decimal
            readToken(file, token, &opsize); // read cost or varname
            if (!file)
                break;
            if (isInteger(token)) {
                cost = string2Cost((const char*)token.c_str());
                if (opsize != (int)token.size() || token.back() != ';') {
                    readToken(file, token, &opsize); // read varname
                    assert(!isInteger(token));
                }
            } else {
                cost = UNIT_COST;
            }
            if ((cost >= MIN_COST && multiplier * cost >= (Double)(MAX_COST - UNIT_COST) / MEDIUM_COST / MEDIUM_COST / MEDIUM_COST / MEDIUM_COST) || (cost < MIN_COST && multiplier * -cost >= (Double)(MAX_COST - UNIT_COST) / MEDIUM_COST / MEDIUM_COST / MEDIUM_COST / MEDIUM_COST)) {
                cerr << "This resolution cannot be ensured on the data type used to represent costs! (see option -precision)" << endl;
                throw BadConfiguration();
            }
            cost *= ToulBar2::costMultiplier;
            if (token != ";") {
                vector<int> scopeIndex;
                while (!isInteger(token)) {
                    string varname = token.substr(0, token.size() - ((token.back() == ';') ? 1 : 0));
                    if (varname.rfind("~", 0) == 0) {
                        cerr << "Sorry, negative literals in opb format not implemented in this version!" << endl;
                        throw WrongFileFormat();
                    }
                    int var = 0;
                    if (varnames.find(varname) != varnames.end()) {
                        var = varnames[varname];
                    } else {
                        var = makeEnumeratedVariable(varname, 0, 1);
                        addValueName(var, "v0");
                        addValueName(var, "v1");
                        varnames[varname] = var;
                        nbvar++;
                    }
                    if (find(scopeIndex.begin(), scopeIndex.end(), var) == scopeIndex.end()) {
                        scopeIndex.push_back(var);
                    }
                    if (token.back() == ';')
                        break;
                    readToken(file, token, &opsize);
                    if (isInteger(token) || token == ";") {
                        opsize = 0;
                        break;
                    }
                }
                if (scopeIndex.size() > 3) {
                    Cost defval = max(MIN_COST, -cost);
                    int ctr = postNaryConstraintBegin(scopeIndex, defval, 1);
                    vector<Value> tuple(scopeIndex.size(), 1);
                    if (cost < MIN_COST) {
                        postNaryConstraintTuple(ctr, tuple, MIN_COST);
                        negCost -= cost;
                    } else {
                        postNaryConstraintTuple(ctr, tuple, cost);
                    }
                    postNaryConstraintEnd(ctr);
                } else if (scopeIndex.size() == 3) {
                    vector<Cost> costs(2 * 2 * 2, max(MIN_COST, -cost));
                    if (cost < MIN_COST) {
                        costs[7] = MIN_COST;
                        negCost -= cost;
                    } else {
                        costs[7] = cost;
                    }
                    postTernaryConstraint(scopeIndex[0], scopeIndex[1], scopeIndex[2], costs);
                } else if (scopeIndex.size() == 2) {
                    vector<Cost> costs(2 * 2, max(MIN_COST, -cost));
                    if (cost < MIN_COST) {
                        costs[3] = MIN_COST;
                        negCost -= cost;
                    } else {
                        costs[3] = cost;
                    }
                    postBinaryConstraint(scopeIndex[0], scopeIndex[1], costs);
                } else if (scopeIndex.size() == 1) {
                    EnumeratedVariable* x = (EnumeratedVariable*)vars[scopeIndex[0]];
                    TemporaryUnaryConstraint unaryconstr;
                    unaryconstr.var = x;
                    if (cost < MIN_COST) {
                        unaryconstr.costs.push_back(-cost);
                        unaryconstr.costs.push_back(MIN_COST);
                        negCost -= cost;
                    } else {
                        unaryconstr.costs.push_back(MIN_COST);
                        unaryconstr.costs.push_back(cost);
                    }
                    unaryconstrs.push_back(unaryconstr);
                } else if (scopeIndex.size() == 0) {
                    inclowerbound += cost;
                } else {
                    cerr << "Sorry! Cannot read objective function with non linear terms of arity " << scopeIndex.size() << endl;
                    throw WrongFileFormat();
                }
            }
        } while (token.back() != ';');
    }

    // read linear constraints
    while (file) {
        vector<int> scopeIndex;
        vector<Cost> coefs;
        string params;
        Cost coef; // allows long long coefficients inside linear constraints
        bool first = true;
        bool soft = false;
        Cost poscoef = MIN_COST;
        Cost negcoef = MIN_COST;
        do {
            readToken(file, token, &opsize); // read coefficient or operator or comments
            // skip comments
            while (file && token[0] == '*') {
                getline(file, dummy);
                readToken(file, token, &opsize);
            }
            if (!file || token == ";")
                break;
            if (token[0] == '[') {
                cursoft++;
                if (token.size() == 1) {
                    readToken(file, token, &opsize);
                }
                Cost cost = string2Cost((char*)token.substr((token[0] == '[')? 1 : 0).c_str());
                assert(cost >= MIN_COST);
                if (token.back() != ']') {
                    readToken(file, token, &opsize);
                }
                assert(token.back() == ']'); // TODO: allows no space after bracket
                if (cost > MIN_COST && cost < getUb()) {
                    soft = true;
                    string varname = "s" + to_string(cursoft);
                    int var = 0;
                    if (varnames.find(varname) != varnames.end()) {
                        var = varnames[varname];
                    } else {
                        var = makeEnumeratedVariable(varname, 0, 1);
                        addValueName(var, "v0");
                        addValueName(var, "v1");
                        varnames[varname] = var;
                        nbsoft++;
                    }
                    if (scopeIndex.size() == 0) {
                        scopeIndex.push_back(var);
                        coefs.push_back(0); // to be updated later
                    } else {
                        cerr << "Wrong soft constraint definition with more than one cost block! " << token << endl;
                        throw WrongFileFormat();
                    }
                    EnumeratedVariable* x = (EnumeratedVariable*)vars[var];
                    TemporaryUnaryConstraint unaryconstr;
                    unaryconstr.var = x;
                    unaryconstr.costs.push_back(MIN_COST);
                    unaryconstr.costs.push_back(cost);
                    unaryconstrs.push_back(unaryconstr);
                }
            } else if (token.substr(0, 2) == "<=" || token.substr(0, 1) == "=" || token.substr(0, 2) == ">=") {
                if (first) {
                    cerr << "Wrong constraint definition with empty left-hand side! " << token << endl;
                    throw WrongFileFormat();
                }
                opsize = (token[0] == '=') ? 1 : 2;
                string op = token.substr(0, opsize);
                readToken(file, token, &opsize); // read right coef
                assert(isInteger(token));
                coef = string2Cost((char*)token.c_str());
                maxarity = max(maxarity, (int)scopeIndex.size());
                nblinear++;
                if (op == "<=" || op == "=") {
                    params = to_string(-coef);
                    for (unsigned int i = 0; i < scopeIndex.size(); i++) {
                        if (soft && i == 0) {
                            params += to_string(" ") + to_string(poscoef - coef);
                        } else {
                            params += to_string(" ") + to_string(-coefs[i]);
                        }
                    }
                    postKnapsackConstraint(scopeIndex, params); // create <= before >= in case of equality because the first constraint may be more tight than the second one.
                }
                if (op == ">=" || op == "=") {
                    params = to_string(coef);
                    for (unsigned int i = 0; i < scopeIndex.size(); i++) {
                        if (soft && i == 0) {
                            params += to_string(" ") + to_string(-negcoef + coef);
                        } else {
                            params += to_string(" ") + to_string(coefs[i]);
                        }
                    }
                    postKnapsackConstraint(scopeIndex, params);
                }
                first = false;
            } else {
                assert(token.back() != ';');
                if (isInteger(token)) {
                    coef = string2Cost((char*)token.c_str());
                    readToken(file, token, &opsize); // read varname
                } else {
                    if (!first) {
                        cerr << "Sorry, constraints with non-linear terms in opb format not implemented in this version!" << endl;
                        throw WrongFileFormat();
                    }
                    coef = 1;
                }
                if (coef < MIN_COST) {
                    negcoef += coef;
                } else {
                    poscoef += coef;
                }
                assert(token.back() != ';');
                assert(!isInteger(token));
                if (token.back() == '=') {
                    opsize = (token[token.size() - 2] == '<' || token[token.size() - 2] == '>') ? -2 : -1;
                }
                string varname = token.substr(0, token.size() + opsize);
                if (varname.rfind("~", 0) == 0) {
                    cerr << "Sorry, negative literals in opb format not implemented in this version!!" << endl;
                    throw WrongFileFormat();
                }
                int var = 0;
                if (varnames.find(varname) != varnames.end()) {
                    var = varnames[varname];
                } else {
                    var = makeEnumeratedVariable(varname, 0, 1);
                    addValueName(var, "v0");
                    addValueName(var, "v1");
                    varnames[varname] = var;
                    nbvar++;
                }
                if (find(scopeIndex.begin(), scopeIndex.end(), var) == scopeIndex.end()) {
                    scopeIndex.push_back(var);
                    coefs.push_back(coef);
                } else {
                    coefs[find(scopeIndex.begin(), scopeIndex.end(), var) - scopeIndex.begin()] += coef;
                }
                first = false;
            }
        } while (token.back() != ';');
    }
    assert(cursoft == nbsoft);

    // apply basic initial propagation AFTER complete network loading
    postNullaryConstraint(inclowerbound);

    for (unsigned int u = 0; u < unaryconstrs.size(); u++) {
        postUnaryConstraint(unaryconstrs[u].var->wcspIndex, unaryconstrs[u].costs);
    }
    sortConstraints();
    if (ToulBar2::verbose >= 0) {
        if (nbsoft > 0) {
            cout << "c Read " << nbvar << " variables, with 2 values at most, and " << nblinear << " linear constraints, including " << nbsoft << " soft constraints, with maximum arity " << maxarity << " (cost multiplier: " << ToulBar2::costMultiplier << ", shifting value: " << -negCost << ")" << endl;
        } else {
            cout << "c Read " << nbvar << " variables, with 2 values at most, and " << nblinear << " linear constraints, with maximum arity " << maxarity << " (cost multiplier: " << ToulBar2::costMultiplier << ", shifting value: " << -negCost << ")" << endl;
        }
    }
}

void WCSP::read_lp(const char* fileName)
{
    ifstream rfile(fileName, (ToulBar2::gz || ToulBar2::bz2 || ToulBar2::xz) ? (std::ios_base::in | std::ios_base::binary) : (std::ios_base::in));
#ifdef BOOST
    boost::iostreams::filtering_streambuf<boost::iostreams::input> zfile;
    if (ToulBar2::gz) {
        zfile.push(boost::iostreams::gzip_decompressor());
    } else if (ToulBar2::bz2) {
#ifndef NO_BZ2
        zfile.push(boost::iostreams::bzip2_decompressor());
#else
        cerr << "Error: Boost requires the bzip2 library to read bz2 compressed lp format files." << endl;
        throw WrongFileFormat();
#endif
    } else if (ToulBar2::xz) {
#ifndef NO_LZMA
        zfile.push(boost::iostreams::lzma_decompressor());
#else
        cerr << "Error: compiling with Boost version 1.65 or higher is needed to allow to read xz compressed lp format files." << endl;
        throw WrongFileFormat();
#endif
    }
    zfile.push(rfile);
    istream ifile(&zfile);

    if (ToulBar2::stdin_format.length() == 0 && !rfile) {
        cerr << "Could not open lp file : " << fileName << endl;
        throw WrongFileFormat();
    }
    istream& file = (ToulBar2::stdin_format.length() > 0) ? cin : ifile;
#else
    if (ToulBar2::gz || ToulBar2::bz2 || ToulBar2::xz) {
        cerr << "Error: compiling with Boost iostreams library is needed to allow to read compressed lp format files." << endl;
        throw WrongFileFormat();
    }
    if (ToulBar2::stdin_format.length() == 0 && !rfile) {
        cerr << "Could not open lp file : " << fileName << endl;
        throw WrongFileFormat();
    }
    istream& file = (ToulBar2::stdin_format.length() > 0) ? cin : rfile;
#endif

    auto pb = baryonyx::make_problem(file);
    if (!pb) {
        cerr << "Could not read lp file : " << fileName << endl;
        throw WrongFileFormat();
    }
    if (pb.status != baryonyx::file_format_error_tag::success) {
        cerr << "Wrong lp format : " << pb.status << endl;
        throw WrongFileFormat();
    }

    updateUb((MAX_COST - UNIT_COST) / MEDIUM_COST / MEDIUM_COST);

    int maxarity = 0;
    int maxdom = 0;
    int nbvar = 0;
    int nblinear = 0;
    vector<TemporaryUnaryConstraint> unaryconstrs;

    map<string, int> varnames;
    assert(pb.vars.names.size() == pb.vars.values.size());
    for (unsigned int i = 0; i < pb.vars.names.size(); i++) {
        if (!pb.vars.values[i].touched) {
            assert(pb.vars.values[i].min == std::numeric_limits<int>::min());
            pb.vars.values[i].min = 0;
        }
        // assert(pb.vars.values[i].type == baryonyx::variable_type::binary || pb.vars.values[i].type == baryonyx::variable_type::general || (pb.vars.values[i].type == baryonyx::variable_type::real && pb.vars.values[i].min == pb.vars.values[i].max));
        if (pb.vars.values[i].max == std::numeric_limits<int>::max()) {
            cerr << "Sorry, cannot represent unbounded variable domain for " << pb.vars.names[i] << " (must add a valid upper bound)" << endl;
            throw BadConfiguration();
        }
        if (pb.vars.values[i].min == std::numeric_limits<int>::min()) {
            cerr << "Sorry, cannot represent unbounded variable domain for " << pb.vars.names[i] << " (must add a valid lower bound)" << endl;
            throw BadConfiguration();
        }
        if (pb.vars.values[i].min > pb.vars.values[i].max) {
            cerr << "Infeasible problem! " << pb.vars.names[i] << " has a lower bound " << pb.vars.values[i].min << " greater than its upper bound " << pb.vars.values[i].max << endl;
            throw Contradiction();
        }
#ifdef WCSPFORMATONLY
        int var = makeEnumeratedVariable(to_string(pb.vars.names[i]), 0, pb.vars.values[i].max - pb.vars.values[i].min);
#else
        int var = makeEnumeratedVariable(to_string(pb.vars.names[i]), pb.vars.values[i].min, pb.vars.values[i].max);
#endif
        varnames[to_string(pb.vars.names[i])] = var;
        for (int v = pb.vars.values[i].min; v <= pb.vars.values[i].max; v++) {
            string valname = "v";
            valname += to_string(v);
            addValueName(var, valname);
        }
        nbvar++;
        maxdom = max(maxdom, pb.vars.values[i].max - pb.vars.values[i].min + 1);
    }

    if (ToulBar2::resolution < 0) {
        cerr << "Sorry, cannot use a negative value for precision! (see option -precision)" << endl;
        throw BadConfiguration();
    }
    Double multiplier = Exp10((Double)min(ToulBar2::resolution, baryonyx::precision));
    ToulBar2::decimalPoint = abs(min(ToulBar2::resolution, baryonyx::precision));

    // read objective function
    ToulBar2::setCostMultiplier(1.0);
    if (pb.type == baryonyx::objective_function_type::maximize) {
        ToulBar2::setCostMultiplier(ToulBar2::costMultiplier * -1.0);
    }

    Double totalShiftCost = 0.0;
    for (auto& elem : pb.objective.elements) {
        Double mult = elem.factor * ToulBar2::costMultiplier;
        if (multiplier * std::abs(mult) * max(abs(pb.vars.values[elem.variable_index].min), abs(pb.vars.values[elem.variable_index].max)) >= (Double)(MAX_COST - UNIT_COST) / MEDIUM_COST / MEDIUM_COST / MEDIUM_COST / MEDIUM_COST) {
            cerr << "This resolution cannot be ensured on the data type used to represent costs! (see option -precision)" << endl;
            throw BadConfiguration();
        }
        mult = multiplier * min(mult * pb.vars.values[elem.variable_index].min, mult * pb.vars.values[elem.variable_index].max);
        totalShiftCost += mult;
    }
    for (auto& elem : pb.objective.qelements) {
        Double mult = elem.factor * ToulBar2::costMultiplier;
        if (multiplier * std::abs(mult) * max(abs(pb.vars.values[elem.variable_index_a].min), abs(pb.vars.values[elem.variable_index_a].max)) * max(abs(pb.vars.values[elem.variable_index_b].min), abs(pb.vars.values[elem.variable_index_b].max)) >= (Double)(MAX_COST - UNIT_COST) / MEDIUM_COST / MEDIUM_COST / MEDIUM_COST / MEDIUM_COST) {
            cerr << "This resolution cannot be ensured on the data type used to represent costs! (see option -precision)" << endl;
            throw BadConfiguration();
        }
        mult = multiplier * min(min(mult * pb.vars.values[elem.variable_index_a].min * pb.vars.values[elem.variable_index_b].min, mult * pb.vars.values[elem.variable_index_a].min * pb.vars.values[elem.variable_index_b].max), min(mult * pb.vars.values[elem.variable_index_a].max * pb.vars.values[elem.variable_index_b].min, mult * pb.vars.values[elem.variable_index_a].max * pb.vars.values[elem.variable_index_b].max));
        totalShiftCost += mult;
    }
    negCost -= (Cost)roundl(totalShiftCost);

    // linear objective
    for (auto& elem : pb.objective.elements) {
        Double mult = elem.factor * ToulBar2::costMultiplier * multiplier;
        EnumeratedVariable* x = (EnumeratedVariable*)vars[elem.variable_index];
        TemporaryUnaryConstraint unaryconstr;
        unaryconstr.var = x;
        Double shiftCost = min(mult * pb.vars.values[elem.variable_index].min, mult * pb.vars.values[elem.variable_index].max);
        for (int v = pb.vars.values[elem.variable_index].min; v <= pb.vars.values[elem.variable_index].max; v++) {
            unaryconstr.costs.push_back((Cost)roundl(mult * v - shiftCost));
        }
        unaryconstrs.push_back(unaryconstr);
    }
    // quadratic objective
    for (auto& elem : pb.objective.qelements) {
        Double mult = elem.factor * ToulBar2::costMultiplier * multiplier;
        vector<Cost> costs;
        Double shiftCost = min(min(mult * pb.vars.values[elem.variable_index_a].min * pb.vars.values[elem.variable_index_b].min, mult * pb.vars.values[elem.variable_index_a].min * pb.vars.values[elem.variable_index_b].max), min(mult * pb.vars.values[elem.variable_index_a].max * pb.vars.values[elem.variable_index_b].min, mult * pb.vars.values[elem.variable_index_a].max * pb.vars.values[elem.variable_index_b].max));
        for (int v = pb.vars.values[elem.variable_index_a].min; v <= pb.vars.values[elem.variable_index_a].max; v++) {
            for (int w = pb.vars.values[elem.variable_index_b].min; w <= pb.vars.values[elem.variable_index_b].max; w++) {
                costs.push_back((Cost)roundl(mult * v * w - shiftCost));
            }
        }
        postBinaryConstraint(elem.variable_index_a, elem.variable_index_b, costs);
    }

    // read linear constraints
    for (auto& ctr : pb.less_constraints) {
        vector<int> scopeIndex;
        string params;
        maxarity = max(maxarity, (int)ctr.elements.size());
        nblinear++;
        Double multiplier = Exp10((Double)min(ToulBar2::resolution, ctr.precision));
        if (multiplier * std::abs(ctr.value) >= (Double)(std::numeric_limits<long long int>::max())) {
            cerr << "This resolution cannot be ensured on the data type used to represent linear constraint coefficients! (see option -precision)" << endl;
            throw BadConfiguration();
        }
        params = to_string(static_cast<long long int>(-baryonyx::Ceil(ctr.value * multiplier)));
        for (unsigned int i = 0; i < ctr.elements.size(); i++) {
            int x = ctr.elements[i].variable_index;
            int min = pb.vars.values[x].min;
            int max = pb.vars.values[x].max;
            if (multiplier * std::abs(ctr.elements[i].factor) * std::max(abs(min), abs(max)) * ctr.elements.size() >= (Double)(std::numeric_limits<long long int>::max())) {
                cerr << "This resolution cannot be ensured on the data type used to represent linear constraint coefficients! (see option -precision)" << endl;
                throw BadConfiguration();
            }
            long long int coef = static_cast<long long int>(baryonyx::Floor(ctr.elements[i].factor * multiplier));
            params += to_string(" ") + to_string(max - min + 1);
            scopeIndex.push_back(x);
            for (int v = min; v <= max; v++) {
#ifdef WCSPFORMATONLY
                params += to_string(" ") + to_string(v - min);
#else
                params += to_string(" ") + to_string(v);
#endif
                params += to_string(" ") + to_string(-coef * v);
            }
        }
        postKnapsackConstraint(scopeIndex, params, false, true, false);
    }
    for (auto& ctr : pb.equal_constraints) {
        vector<int> scopeIndex;
        string params;
        maxarity = max(maxarity, (int)ctr.elements.size());
        nblinear++;
        Double multiplier = Exp10((Double)min(ToulBar2::resolution, ctr.precision));
        if (multiplier * std::abs(ctr.value) >= (Double)(std::numeric_limits<long long int>::max())) {
            cerr << "This resolution cannot be ensured on the data type used to represent linear constraint coefficients! (see option -precision)" << endl;
            throw BadConfiguration();
        }
        params = to_string(static_cast<long long int>(-baryonyx::Ceil(ctr.value * multiplier)));
        for (unsigned int i = 0; i < ctr.elements.size(); i++) {
            int x = ctr.elements[i].variable_index;
            int min = pb.vars.values[x].min;
            int max = pb.vars.values[x].max;
            if (multiplier * std::abs(ctr.elements[i].factor) * std::max(abs(min), abs(max)) * ctr.elements.size() >= (Double)(std::numeric_limits<long long int>::max())) {
                cerr << "This resolution cannot be ensured on the data type used to represent linear constraint coefficients! (see option -precision)" << endl;
                throw BadConfiguration();
            }
            long long int coef = static_cast<long long int>(baryonyx::Floor(ctr.elements[i].factor * multiplier));
            params += to_string(" ") + to_string(max - min + 1);
            for (int v = min; v <= max; v++) {
#ifdef WCSPFORMATONLY
                params += to_string(" ") + to_string(v - min);
#else
                params += to_string(" ") + to_string(v);
#endif
                params += to_string(" ") + to_string(-coef * v);
            }
        }
        postKnapsackConstraint(scopeIndex, params, false, true, false);
        params = to_string(static_cast<long long int>(baryonyx::Floor(ctr.value * multiplier)));
        for (unsigned int i = 0; i < ctr.elements.size(); i++) {
            int x = ctr.elements[i].variable_index;
            int min = pb.vars.values[x].min;
            int max = pb.vars.values[x].max;
            if (multiplier * std::abs(ctr.elements[i].factor) * std::max(abs(min), abs(max)) * ctr.elements.size() >= (Double)(std::numeric_limits<long long int>::max())) {
                cerr << "This resolution cannot be ensured on the data type used to represent linear constraint coefficients! (see option -precision)" << endl;
                throw BadConfiguration();
            }
            long long int coef = static_cast<long long int>(baryonyx::Ceil(ctr.elements[i].factor * multiplier));
            params += to_string(" ") + to_string(max - min + 1);
            scopeIndex.push_back(x);
            for (int v = min; v <= max; v++) {
#ifdef WCSPFORMATONLY
                params += to_string(" ") + to_string(v - min);
#else
                params += to_string(" ") + to_string(v);
#endif
                params += to_string(" ") + to_string(coef * v);
            }
        }
        postKnapsackConstraint(scopeIndex, params, false, true, false);
    }
    for (auto& ctr : pb.greater_constraints) {
        vector<int> scopeIndex;
        string params;
        maxarity = max(maxarity, (int)ctr.elements.size());
        nblinear++;
        Double multiplier = Exp10((Double)min(ToulBar2::resolution, ctr.precision));
        if (multiplier * std::abs(ctr.value) >= (Double)(std::numeric_limits<long long int>::max())) {
            cerr << "This resolution cannot be ensured on the data type used to represent linear constraint coefficients! (see option -precision)" << endl;
            throw BadConfiguration();
        }
        params = to_string(static_cast<long long int>(baryonyx::Floor(ctr.value * multiplier)));
        for (unsigned int i = 0; i < ctr.elements.size(); i++) {
            int x = ctr.elements[i].variable_index;
            int min = pb.vars.values[x].min;
            int max = pb.vars.values[x].max;
            if (multiplier * std::abs(ctr.elements[i].factor) * std::max(abs(min), abs(max)) * ctr.elements.size() >= (Double)(std::numeric_limits<long long int>::max())) {
                cerr << "This resolution cannot be ensured on the data type used to represent linear constraint coefficients! (see option -precision)" << endl;
                throw BadConfiguration();
            }
            long long int coef = static_cast<long long int>(baryonyx::Ceil(ctr.elements[i].factor * multiplier));
            params += to_string(" ") + to_string(max - min + 1);
            scopeIndex.push_back(x);
            for (int v = min; v <= max; v++) {
#ifdef WCSPFORMATONLY
                params += to_string(" ") + to_string(v - min);
#else
                params += to_string(" ") + to_string(v);
#endif
                params += to_string(" ") + to_string(coef * v);
            }
        }
        postKnapsackConstraint(scopeIndex, params, false, true, false);
    }

    ToulBar2::resolution = min(ToulBar2::resolution, baryonyx::precision);

    Cost upperBound = getUb();
    if (ToulBar2::externalUB.length() != 0) {
        upperBound = min(upperBound, decimalToCost(ToulBar2::externalUB, 0) + negCost);
    }
    if (ToulBar2::deltaUbS.length() != 0) {
        ToulBar2::deltaUbAbsolute = max(MIN_COST, decimalToCost(ToulBar2::deltaUbS, 0));
        ToulBar2::deltaUb = max(ToulBar2::deltaUbAbsolute, (Cost)(ToulBar2::deltaUbRelativeGap * (Double)upperBound));
        if (ToulBar2::deltaUb > MIN_COST) {
            // as long as a true certificate as not been found we must compensate for the deltaUb in CUT
            upperBound += ToulBar2::deltaUb;
        }
    }
    updateUb(upperBound);

    if (ToulBar2::costThresholdS.size())
        ToulBar2::costThreshold = decimalToCost(ToulBar2::costThresholdS, 0);
    if (ToulBar2::costThresholdPreS.size())
        ToulBar2::costThresholdPre = decimalToCost(ToulBar2::costThresholdPreS, 0);
    if (ToulBar2::vnsOptimumS.size())
        ToulBar2::vnsOptimum = decimalToCost(ToulBar2::vnsOptimumS, 0) + negCost;

    // apply basic initial propagation AFTER complete network loading
    if (multiplier * std::abs(pb.objective.value) >= (Double)(MAX_COST - UNIT_COST) / MEDIUM_COST / MEDIUM_COST / MEDIUM_COST / MEDIUM_COST) {
        cerr << "This resolution cannot be ensured on the data type used to represent costs! (see option -precision)" << endl;
        throw BadConfiguration();
    }
    postNullaryConstraint((Cost)roundl(pb.objective.value * ToulBar2::costMultiplier * multiplier));

    for (unsigned int u = 0; u < unaryconstrs.size(); u++) {
        postUnaryConstraint(unaryconstrs[u].var->wcspIndex, unaryconstrs[u].costs);
    }
    sortConstraints();
    if (ToulBar2::verbose >= 0)
        cout << "c Read " << nbvar << " variables, with " << maxdom << " values at most, and " << nblinear << " linear constraints, with maximum arity " << maxarity << " (precision: " << ToulBar2::decimalPoint << ", shifting value: " << -negCost << ")" << endl;
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
