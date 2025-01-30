/** \file tb2types.hpp
 *  \brief Macros, types, and globals.
 *
 * The main types are:
 * - ::Value : domain value
 * - ::Cost : cost value (exact type depends on compilation flag)
 * - ::Long : large integer (long long int)
 * - ::TProb : probability value (exact type depends on compilation flag)
 * - ::TLogProb : log probability value (exact type depends on compilation flag)
 * - ::Double : large float (long double)
 * - ::tValue : a short Value type for tuples
 * - ::Tuple : vector of tValues to encode tuples
 *
 * \note Compilation flag for Cost is: \c INT_COST (int), \c LONGLONG_COST (long long), or \c PARETOPAIR_COST (see ::ParetoPair)
 * \warning \c PARETOPAIR_COST is fragile.
 * \note Compilation flag for TProb is: \c DOUBLE_PROB or \c LONGDOUBLE_PROB
 * \note Compilation flag for T(Log)Prob is: \c DOUBLE_PROB or \c LONGDOUBLE_PROB
 */

#ifndef TB2TYPES_HPP_
#define TB2TYPES_HPP_

// #define INT_COST
// #define LONGLONG_COST
// #define PARETOPAIR_COST

// #define DOUBLE_PROB
// #define LONGDOUBLE_PROB

#include "utils/tb2utils.hpp"
// Must be included after tb2utils.hpp
#include "utils/tb2integer.hpp"
#ifdef QUAD_PROB
#include <boost/multiprecision/float128.hpp>
//#include <quadmath.h> // only with gcc/g++
#endif

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::ifstream;
using std::istream;
using std::istringstream;
using std::list;
using std::make_pair;
using std::make_tuple;
using std::map;
using std::max;
using std::min;
using std::numeric_limits;
using std::ofstream;
using std::ostringstream;
using std::queue;
using std::set;
using std::stack;
using std::string;
using std::stringstream;
using std::tie;
using std::tuple;

/// Special character value at the beginning of a variable's name to identify implicit variables (i.e., variables which are not decision variables)
const string IMPLICIT_VAR_TAG = "#";

/// Special character value at the beginning of a variable's name to identify hidden variables like diverse extra variables corresponding to the current sequence of diverse solutions found so far
const string HIDDEN_VAR_TAG = "^";
const string HIDDEN_VAR_TAG_HVE = "^c"; // tag for a hidden variable representing a dualized nary cost function
const string HIDDEN_VAR_TAG_HVE_PRE = "^!"; // temporary hidden variable which should disappear after dedualization in preprocessing

/// Domain value (can be positive or negative integers)
#ifdef SHORT_VALUE
typedef int16_t Value;
typedef int16_t tValue;
#else
typedef int Value;
#if defined(XMLFLAG) || defined(XMLFLAG3)
typedef int tValue;
#else
typedef int16_t tValue;
#endif
#endif

/// Maximum domain value
const Value MAX_VAL = (std::numeric_limits<Value>::max() / 2);
/// Forbidden domain value
const Value WRONG_VAL = std::numeric_limits<Value>::max();
/// Minimum domain value
const Value MIN_VAL = -(std::numeric_limits<Value>::max() / 2);
/// Maximum domain size
/// \deprecated Should use WCSP::getMaxDomainSize instead.
const Value MAX_DOMAIN_SIZE = 10000;

typedef vector<tValue> Tuple;

// For very large domains with ternary cost functions, use NARYPROJECTIONSIZE=2 instead of 3
const int NARYPROJECTIONSIZE = 3; // limit on the number of unassigned variables before nary constraints are projected to smaller-arity constraint (should be between 1 and 3)
const unsigned int NARYPROJECTION3MAXDOMSIZE = 30; // limit on the maximum initial domain size for nary to ternary projection
const Long NARYPROJECTION3PRODDOMSIZE = 10000; // limit on the cartesian product of initial domain sizes for nary to ternary projection
const Long NARYDECONNECTSIZE = 4; // maximum number of initial tuples in nary constraints in order to check for its removal (if it is always satisfied by current domains)

const int MAX_BRANCH_SIZE = 1000000;
const ptrdiff_t CHOICE_POINT_LIMIT = SIZE_MAX; // warning! converted to -1
const ptrdiff_t OPEN_NODE_LIMIT = SIZE_MAX; // warning! converted to -1

#if (defined(SHORT_COST) || defined(SHORT_VALUE))
// C++ integer promotion occurs on any arithmetic operation (i.e. int16_t ope int_16_t results to int type conversion)
inline int16_t min(int16_t x, int y)
{
    if (x < y)
        return x;
    else
        return y;
}
inline int16_t min(int x, int16_t y)
{
    if (x < y)
        return x;
    else
        return y;
}
inline int16_t max(int16_t x, int y)
{
    if (x > y)
        return x;
    else
        return y;
}
inline int16_t max(int x, int16_t y)
{
    if (x > y)
        return x;
    else
        return y;
}
#endif

#ifdef SHORT_COST
const bool PARTIALORDER = false;
typedef int16_t Cost;
const Cost MIN_COST = short{ 0 };
const Cost UNIT_COST = short{ 1 };
const Cost SMALL_COST = short{ 1 };
const Cost MEDIUM_COST = short{ 3 };
const Cost LARGE_COST = short{ 100 };
const Cost MAX_COST = (std::numeric_limits<Cost>::max() / MEDIUM_COST);

#if defined __has_builtin
#if __has_builtin(__builtin_sadd_overflow)
inline bool Add(Cost a, Cost b, Cost* c)
{
    return __builtin_sadd_overflow(a, b, c);
}
#endif
#if __has_builtin(__builtin_ssub_overflow)
inline bool Sub(Cost a, Cost b, Cost* c)
{
    return __builtin_ssub_overflow(a, b, c);
}
#endif
#if __has_builtin(__builtin_smul_overflow)
inline bool Mul(Cost a, Cost b, Cost* c)
{
    return __builtin_smul_overflow(a, b, c);
}
#endif
#else
inline bool Add(Cost a, Cost b, Cost* c)
{
    *c = a + b;
    return false;
}
inline bool Sub(Cost a, Cost b, Cost* c)
{
    *c = a - b;
    return false;
}
inline bool Mul(Cost a, Cost b, Cost* c)
{
    *c = a * b;
    return false;
}
#endif

inline Cost MIN(Cost a, Cost b)
{
    return min(a, b);
}
inline Cost MAX(Cost a, Cost b) { return max(a, b); }
inline Cost MULT(Cost a, Double b)
{
    assert(b < MAX_COST);
    if (a >= MAX_COST)
        return MAX_COST;
    else if (b <= UNIT_COST)
        return (Cost)((Double)a * b);
    else if (a < (Double)MAX_COST / b)
        return (Cost)((Double)a * b);
    else {
        cerr << "Error: cost multiplication overflow!" << endl;
        throw InternalError();
    }
}
inline Cost GLB(Cost a, Cost b) { return MIN(a, b); }
inline Cost LUB(Cost a, Cost b) { return MAX(a, b); }
inline bool GLB(Cost* a, Cost b)
{
    if (b < *a) {
        *a = b;
        return true;
    } else
        return false;
}
inline bool LUB(Cost* a, Cost b)
{
    if (b > *a) {
        *a = b;
        return true;
    } else
        return false;
}
inline bool GLBTEST(Cost a, Cost b) { return (b < a); }
inline bool LUBTEST(Cost a, Cost b) { return (b > a); }
inline bool DACTEST(Cost a, Cost b) { return (a == 0 && b > 0); }
inline bool SUPPORTTEST(Cost a, Cost b) { return false; }
inline bool SUPPORTTEST(Cost a) { return false; }
inline void initCosts() {}
#endif

#ifdef INT_COST
const bool PARTIALORDER = false;
typedef int Cost;
const Cost MIN_COST = 0;
const Cost UNIT_COST = 1;
const Cost SMALL_COST = 1;
const Cost MEDIUM_COST = 3;
const Cost LARGE_COST = 100;
const Cost MAX_COST = ((std::numeric_limits<Cost>::max() / 2) / MEDIUM_COST / MEDIUM_COST);

#if defined __has_builtin
#if __has_builtin(__builtin_saddl_overflow)
inline bool Add(Cost a, Cost b, Cost* c)
{
    return __builtin_saddl_overflow(a, b, c);
}
#endif
#if __has_builtin(__builtin_ssubl_overflow)
inline bool Sub(Cost a, Cost b, Cost* c)
{
    return __builtin_ssubl_overflow(a, b, c);
}
#endif
#if __has_builtin(__builtin_smull_overflow)
inline bool Mul(Cost a, Cost b, Cost* c)
{
    return __builtin_smull_overflow(a, b, c);
}
#endif
#else
inline bool Add(Cost a, Cost b, Cost* c)
{
    *c = a + b;
    return false;
}
inline bool Sub(Cost a, Cost b, Cost* c)
{
    *c = a - b;
    return false;
}
inline bool Mul(Cost a, Cost b, Cost* c)
{
    *c = a * b;
    return false;
}
#endif

inline Cost MIN(Cost a, Cost b)
{
    return min(a, b);
}
inline Cost MAX(Cost a, Cost b) { return max(a, b); }
inline Cost MULT(Cost a, Double b)
{
    assert(b < MAX_COST);
    if (a >= MAX_COST)
        return MAX_COST;
    else if (b <= UNIT_COST)
        return (Cost)((Double)a * b);
    else if (a < (Double)MAX_COST / b)
        return (Cost)((Double)a * b);
    else {
        cerr << "Error: cost multiplication overflow!" << endl;
        throw InternalError();
    }
}
inline Cost GLB(Cost a, Cost b) { return MIN(a, b); }
inline Cost LUB(Cost a, Cost b) { return MAX(a, b); }
inline bool GLB(Cost* a, Cost b)
{
    if (b < *a) {
        *a = b;
        return true;
    } else
        return false;
}
inline bool LUB(Cost* a, Cost b)
{
    if (b > *a) {
        *a = b;
        return true;
    } else
        return false;
}
inline bool GLBTEST(Cost a, Cost b) { return (b < a); }
inline bool LUBTEST(Cost a, Cost b) { return (b > a); }
inline bool DACTEST(Cost a, Cost b) { return (a == 0 && b > 0); }
inline bool SUPPORTTEST(Cost a, Cost b) { return false; }
inline bool SUPPORTTEST(Cost a) { return false; }
inline void initCosts() {}
#endif

#ifdef LONGLONG_COST
const bool PARTIALORDER = false;
typedef Long Cost;
struct DCost {
    Cost c;

    DCost(Cost ic) { c = ic; };
    friend ostream& operator<<(ostream& os, const DCost& r)
    {
        os << r.c;
        return os;
    }
    friend istream& operator>>(istream& is, DCost& r)
    {
        is >> r.c;
        return is;
    }
};
const Cost MIN_COST = 0;
const Cost UNIT_COST = 1;
const Cost SMALL_COST = 1;
const Cost MEDIUM_COST = 3;
const Cost LARGE_COST = 100;
const Cost MAX_COST = ((LONGLONG_MAX / 2) / MEDIUM_COST / MEDIUM_COST);

#if defined __has_builtin
#if __has_builtin(__builtin_saddll_overflow)
inline bool Add(Cost a, Cost b, Cost* c)
{
    return __builtin_saddll_overflow(a, b, c);
}
#endif
#if __has_builtin(__builtin_ssubll_overflow)
inline bool Sub(Cost a, Cost b, Cost* c)
{
    return __builtin_ssubll_overflow(a, b, c);
}
#endif
#if __has_builtin(__builtin_smulll_overflow)
inline bool Mul(Cost a, Cost b, Cost* c)
{
    return __builtin_smulll_overflow(a, b, c);
}
#endif
#else
inline bool Add(Cost a, Cost b, Cost* c)
{
    *c = a + b;
    return false;
}
inline bool Sub(Cost a, Cost b, Cost* c)
{
    *c = a - b;
    return false;
}
inline bool Mul(Cost a, Cost b, Cost* c)
{
    *c = a * b;
    return false;
}
#endif

inline Cost MIN(Cost a, Cost b)
{
    return min(a, b);
}
inline Cost MAX(Cost a, Cost b) { return max(a, b); }
inline Cost MULT(Cost a, Double b)
{
    assert(b < MAX_COST);
    if (a >= MAX_COST)
        return MAX_COST;
    else if (b <= UNIT_COST)
        return (Cost)((Double)a * b);
    else if (a < (Double)MAX_COST / b)
        return (Cost)((Double)a * b);
    else {
        cerr << "Error: cost multiplication overflow!" << endl;
        throw InternalError();
    }
}
inline Cost GLB(Cost a, Cost b) { return MIN(a, b); }
inline Cost LUB(Cost a, Cost b) { return MAX(a, b); }
inline bool GLB(Cost* a, Cost b)
{
    if (b < *a) {
        *a = b;
        return true;
    } else
        return false;
}
inline bool LUB(Cost* a, Cost b)
{
    if (b > *a) {
        *a = b;
        return true;
    } else
        return false;
}
inline bool GLBTEST(Cost a, Cost b) { return (b < a); }
inline bool LUBTEST(Cost a, Cost b) { return (b > a); }
inline bool DACTEST(Cost a, Cost b) { return (a == 0 && b > 0); }
inline bool SUPPORTTEST(Cost a, Cost b) { return false; }
inline bool SUPPORTTEST(Cost a) { return false; }
inline void initCosts() {}
#endif

#ifdef PARETOPAIR_COST
const bool PARTIALORDER = true;
#include "utils/tb2paretopair.hpp"
typedef ParetoPair Cost;
const Cost MIN_COST = PARETOPAIR_MIN;
const Cost UNIT_COST = PARETOPAIR_1;
const Cost SMALL_COST = PARETOPAIR_1;
const Cost MEDIUM_COST = PARETOPAIR_3;
const Cost LARGE_COST = PARETOPAIR_100;
const Cost MAX_COST = PARETOPAIR_MAX;
#endif

#ifdef QUAD_PROB
typedef boost::multiprecision::float128 TProb;
typedef boost::multiprecision::float128 TLogProb;
inline Cost Round(TLogProb f) { return (Cost)round(f); }
#endif

#ifdef DOUBLE_PROB
typedef double TProb;
typedef double TLogProb;
inline Cost Round(TLogProb f) { return (Cost)round(f); }
#endif

#ifdef LONGDOUBLE_PROB
typedef Double TProb;
typedef Double TLogProb;
inline Cost Round(TLogProb f) { return (Cost)roundl(f); }
#endif

inline TLogProb GLogSumExp(TLogProb logc1, TLogProb logc2) // log[exp(c1) + exp(c2)]
{
    if (logc1 == -numeric_limits<TLogProb>::infinity())
        return logc2;
    else if (logc2 == -numeric_limits<TLogProb>::infinity())
        return logc1;
    else {
        if (logc1 >= logc2)
            return logc1 + (Log1p(Exp(logc2 - logc1)));
        else
            return logc2 + (Log1p(Exp(logc1 - logc2)));
    }
}

inline TLogProb GLogSubExp(TLogProb logc1, TLogProb logc2) // log[exp(c1) - exp(c2)]
{
    if (logc1 == logc2)
        return -numeric_limits<TLogProb>::infinity();
    else if (logc1 > logc2)
        return logc1 + (Log(-Expm1(logc2 - logc1)));
    else {
        cerr << "My oh my ! Try to logarithm a negative number" << endl;
        throw InternalError();
    }
}

typedef map<int, int> TSCOPE;
typedef map<int, Value> TAssign;

typedef unsigned int uint;

/*
 * General constants (limits)
 *
 */

const int STORE_SIZE = 16;
#define INTEGERBITS (8 * sizeof(Cost) - 2)

const int MAX_ELIM_BIN = 1000000000;
const int MAX_ARITY = 1000;
/// Maximum number of tuples in n-ary cost functions
const Long MAX_NB_TUPLES = 100000LL;
const int LARGE_NB_VARS = 10000;

const int DECIMAL_POINT = 3; // default number of digits after decimal point for printing floating-point values

/*
 * Parameter settings for cost function propagation
 *
 */

// Clique constraint propagates by including unary and binary cost functions inside its scope (in practice, it can be very time-consuming)
// #define PROPAGATE_CLIQUE_WITH_BINARIES

// Transforms hard clique constraint into knapsack constraint (warning! clique of binary constraints are no more useful)
#define CLIQUE2KNAPSACK

// Transforms nary cost function on Boolean variables with a single nonzero tuple into weighted clause constraint
const bool NARY2CLAUSE = true;

// Transforms hard clause constraint into knapsack constraint
// #define CLAUSE2KNAPSACK

// Transforms knapsack constraint with unit coefficients into hard clause constraint
// #define UNITKNAPSACK2CLAUSE

// Transforms hard decomposable among constraint into knapsack constraint
// #define WAMONG2KNAPSACK

// VAC propagation becomes more incremental in pass 1 keeping removed values from previous cost threshold iterations
#define INCREMENTALVAC

// VAC propagation has optimal O(ed^2) time complexity in pass 1 (but it requires to reset support values at every cost threshold iteration)
// #define AC2001

/*
 * Abstract data type that help post a global cost function
 *
 */

// A value with weight
struct WeightedObjInt {
    int val;
    Cost weight;

    WeightedObjInt(int val_, Cost weight_ = MIN_COST)
        : val(val_)
        , weight(weight_)
    {
    }
};

// A value with upper and lower limit
struct BoundedObjValue {
    Value val;
    unsigned int upper;
    unsigned int lower;

    BoundedObjValue(Value val_, unsigned int upper_, unsigned int lower_ = 0)
        : val(val_)
        , upper(upper_)
        , lower(lower_)
    {
    }
};

// A transition in DFA
struct DFATransition {
    int start;
    int end;
    Value symbol;
    Cost weight;

    DFATransition(int start_, Value symbol_, int end_, Cost weight_ = MIN_COST)
        : start(start_)
        , end(end_)
        , symbol(symbol_)
        , weight(weight_)
    {
    }
};

// A production rule in CFG
struct CFGProductionRule {
    int from;
    Cost weight;
    int order;
    int* to;
};

// A variable-value pair with weight
struct WeightedVarValPair {
    int varIndex;
    Value val;
    Cost weight;
};

/*
 * Global variables encapsulated as static members
 *
 */

typedef void (*externalevent)(int wcspId, int varIndex, Value value, void* solver);
typedef void (*externalcostevent)(int wcspId, int varIndex, Cost cost, void* solver);
typedef void (*externalsolution)(int wcspId, void* solver);
typedef void (*externalfunc)();

typedef enum {
    ELIM_NONE = 0,
    MAX_CARD = 1,
    MIN_DEGREE = 2,
    MIN_FILL = 3,
    ELIM_MST = 4,
    CUTHILL_MCKEE = 5,
    APPROX_MIN_DEGREE = 6,
    ELIM_FILE_ORDER = 7,
    ELIM_LEXICOGRAPHIC_ORDER = 8,
    ELIM_DAG_ORDER = 9,
    ELIM_MAX
} ElimOrderType;

typedef enum {
    CONSTR_ORDER_ID = 1,
    CONSTR_ORDER_DAC = 2,
    CONSTR_ORDER_TIGHTNESS = 3,
    CONSTR_ORDER_DAC_TIGHTNESS = 4,
    CONSTR_ORDER_TIGHTNESS_DAC = 5,
    CONSTR_ORDER_RANDOM = 6,
    CONSTR_ORDER_LAG = 7,
    CONSTR_ORDER_ARITY = 8,
    CONSTR_ORDER_ARITY_DAC = 9,
    CONSTR_ORDER_THEMAX
} ConstrOrdering;

typedef enum {
    BISUPPORT_HEUR_LB = 1,
    BISUPPORT_HEUR_UB = 2,
    BISUPPORT_HEUR_MINGAP = 3,
    BISUPPORT_HEUR_MAXGAP = 4,
    BISUPPORT_HEUR_THEMAX
} BiSupportHeur;

class Pedigree;
class Haplotype;
class BEP;

typedef enum {
    LC_NC = 0,
    LC_SNIC = 0,
    LC_AC = 1,
    LC_DAC = 2,
    LC_FDAC = 3,
    LC_EDAC = 4,
    LC_THEMAX
} LcLevelType;
const int MAX_EAC_ITER = 10000;

typedef enum {
    DFBB,
    VNS,
    DGVNS,
    CPDGVNS,
    RPDGVNS,
    TREEDEC
} SearchMethod;

typedef enum {
    NOBTD,
    BTD,
    RDSBTD,
    RDS,
    ADAPTBTD
} BTDMethod;

typedef enum {
    LS_INIT_RANDOM = -1,
    LS_INIT_INF = -2,
    LS_INIT_SUP = -3,
    LS_INIT_DFBB = -4,
    LS_INIT_LDS0 = 0
} VNSSolutionInitMethod;

typedef enum {
    RANDOMVAR = 0,
    CONFLICTVARBASE = 1,
    CONNECTEDCONFLICTVAR = 2,
    CLUSTERRAND = 3,
    CONFLICTVARMAXDEG = 4,
    CLUSTERSORTED = 5,
    SEPCLUSTERSORTED = 6,
    SEPCLUSTERSORTEDV2 = 7,
    MASTERCLUSTERRAND = 8,
    PCONFLICTVAR = 9
} VNSVariableHeuristic;

typedef enum {
    VNS_ADD1 = 1,
    VNS_MULT2 = 2,
    VNS_LUBY = 3,
    VNS_ADD1JUMP = 4
} VNSInc;

typedef enum {
    WCSP_FORMAT = 1,
    CFN_FORMAT,
    WCNF_FORMAT,
    OPB_FORMAT,
    BEP_FORMAT,
    CNF_FORMAT,
    LG_FORMAT,
    MAP_FORMAT,
    PRE_FORMAT,
    QPBO_FORMAT,
    UAI_FORMAT,
    WBO_FORMAT,
    XCSP2_FORMAT
} ProblemFormat;

struct ValueCost {
    Value value;
    Cost cost;
    friend bool operator<(const ValueCost& u, const ValueCost& v) { return u.cost < v.cost; }
    friend bool operator>(const ValueCost& u, const ValueCost& v) { return u.cost > v.cost; }
    friend bool operator==(const ValueCost& u, const ValueCost& v) { return u.cost == v.cost; }
};

/**
 * It contains all toulbar2 global variables encapsulated as static class members of this class.
 *
 * Each variable may correspond to some command-line option of toulbar2 executable.
 *
 */

class ToulBar2 {
protected:
    virtual ~ToulBar2() = 0; // Trick to avoid any instantiation of ToulBar2
public:
    static string version; ///< \brief toulbar2 version number
    static int verbose; ///< \brief verbosity level (-1:no output, 0: new solutions found, 1: choice points, 2: current domains, 3: basic EPTs, 4: active cost functions, 5: detailed cost functions, 6: more EPTs, 7: detailed EPTs) (command line option -v)

    static bool FullEAC; ///< \brief VAC-integrality/Full-EAC variable ordering heuristic (command line option -vacint and optionally -A)
    static bool VACthreshold; ///< \brief automatic threshold cost value selection for VAC during search (command line option -vacthr)
    static int nbTimesIsVAC; ///< \internal do not use
    static int nbTimesIsVACitThresholdMoreThanOne; ///< \internal do not use
    static bool RASPS; ///< \internal do not use
    static int useRASPS; ///< \brief VAC-based upper bound probing heuristic (0: no rasps, 1: rasps using DFS, >1: using LDS with bounded discrepancy + 1) (command line option -raspslds or -rasps)
    static bool RASPSreset; ///< \brief reset weighted degree variable ordering heuristic after doing upper bound probing (command line option -raspsini)
    static int RASPSangle; ///< \brief automatic threshold cost value selection for probing heuristic (command line option -raspsdeg)
    static Long RASPSnbBacktracks; ///< \brief number of backtracks of VAC-based upper bound probing heuristic (command line option -rasps)
    static int RASPSnbStrictACVariables; ///< \internal do not use
    static Cost RASPSlastitThreshold; ///< \internal do not use
    static bool RASPSsaveitThresholds; ///< \internal do not use
    static vector<pair<Cost, Double>> RASPSitThresholds; ///< \internal do not use

    static int debug; ///< \brief debug mode(0: no debug, 1: current search depth and statics on nogoods for BTD, 2: idem plus some information on heuristics, 3: idem plus save problem at each node if verbose >= 1) (command line option -Z)
    static string externalUB; ///< \brief initial upper bound in CFN format
    static int showSolutions; ///< \brief shows each solution found (0: nothing, 1: value indexes, 2: value names, 3: variable&value names) (command line option -s)
    static bool showHidden; ///< \brief shows hidden variables for each solution found (command line option -s with a negative value)
    static int writeSolution; ///< \brief writes each solution found (0: nothing, 1: value indexes, 2: value names, 3: variable&value names) (command line option -w)
    static FILE* solutionFile; ///< \internal do not use
    static long solutionFileRewindPos; ///< \internal do not use
    static Long allSolutions; ///< \brief finds at most a given number of solutions with a cost strictly lower than the initial upper bound and stops (or counts the number of zero-cost satisfiable solutions in conjunction with BTD) (command line option -a)
    static int dumpWCSP; ///< \brief saves the problem in wcsp (0: do not save, 1: original or 2: after preprocessing) or cfn (3: original or 4: after preprocessing) format (command line option -z)
    static bool dumpOriginalAfterPreprocessing; ///< \brief saves the problem with initial domains after preprocessing (used in conjunction with dumpWCSP)
    static bool approximateCountingBTD; ///< \brief approximate zero-cost satisfiable solution counting using BTD (command line options -D and -a and -B=1)
    static bool binaryBranching; ///< \brief tree search using binary branching instead of n-ary branching for enumerated domains (command line option -b)
    static int dichotomicBranching; ///< \brief tree search using dichotomic branching if current domain size is strictly greater than ToulBar2::dichotomicBranchingSize (0: no dichotomic branching, 1: splitting in the middle of domain range, 2: splitting in the middle of sorted unary costs) (command line option -d)
    static unsigned int dichotomicBranchingSize; ///< \brief dichotomic branching threshold (related to command line option -d)
    static bool sortDomains; ///< \brief sorts domains in preprocessing based on increasing unary costs (command line option -sortd) \warning Works only for binary WCSPs.
    static map<int, ValueCost*> sortedDomains; ///< \internal do not use
    static bool solutionBasedPhaseSaving; ///< \brief solution-based phase saving value heuristic (command line option -solr)
    static Double bisupport; ///< \brief value heuristic in bi-objective optimization when the second objective is encapsulated by a bounding constraint (command line option -bisupport)
    static int elimDegree; ///< \brief boosting search with variable elimination of small degree (0: no variable elimination, 1: linked to at most one binary cost function, 2: linked to at most two binary cost functions, 3: linked to at most one ternary cost function and two scope-included cost functions) (command line option -e)
    static int elimDegree_preprocessing; ///< \brief  in preprocessing, generic variable elimination of degree less than or equal to a given value (0: no variable elimination) (command line option -p)
    static int elimDegree_; ///< \internal do not use
    static int elimDegree_preprocessing_; ///< \internal do not use
    static int elimSpaceMaxMB; ///< \brief maximum space size for generic variable elimination (in MegaByte) (related to command line option -p)
    static int minsumDiffusion; ///< \brief in preprocessing, applies Min Sum Diffusion algorithm a given number of iterations (command line option -M)
    static int preprocessTernaryRPC; ///< \brief in preprocessing, simulates restricted path consistency by adding ternary cost functions on most-promising triangles of binary cost functions (maximum space size in MegaByte) (command line option -t)
    static int hve; ///< \brief hidden variable encoding into a binary WCSP
    static int pwc; ///< \brief pairwise consistency by dual encoding into a binary WCSP
    static bool pwcMinimalDualGraph; ///< \brief minimizes dual intersection graph by removing redundant edges
    static int preprocessFunctional; ///< \brief in preprocessing, applies variable elimination of 0: no variable, 1: functional, or 2: bijective variables, or 3:functional before and after PWC (command line option -f)
    static bool costfuncSeparate; ///< \brief in preprocessing, applies pairwise decomposition of non-binary cost functions (command line option -dec)
    static int preprocessNary; ///< \brief in preprocessing, projects n-ary cost functions on all their scope-included binary cost functions if n is lower than a given value  (0: no projection) (command line option -n)
    static bool QueueComplexity; ///< \brief ensures optimal worst-case time complexity of DAC and EAC (command line option -o)
    static bool Static_variable_ordering; ///< \brief tree search using a static variable ordering heuristic (same order as DAC) (command line option -svo)
    static bool lastConflict; ///< \brief tree search using binary branching with last conflict backjumping variable ordering heuristic (command line options -c and -b)
    static int weightedDegree; ///< \brief weighted degree variable ordering heuristic if the number of cost functions is less than a given value (command line option -q)
    static int weightedTightness; ///< \brief in preprocessing, initializes weighted degrees associated to cost functions by their 1: average or 2: median costs (command line options -m and -q)
    static int constrOrdering; ///< \brief in preprocessing, sorts constraints based on 0: do not sort, 1: lexicographic ordering, 2: decreasing DAC ordering, 3: decreasing constraint tightness, 4: DAC then tightness, 5: tightness then DAC, 6: random order, or the opposite order if using a negative value (command line option -sortc)
    static bool MSTDAC; ///< \brief maximum spanning tree DAC ordering (command line option -mst)
    static int DEE; ///< \brief soft neighborhood substitutability, a.k.a., dead-end elimination (0: no elimination, 1: restricted form during search, 2: full in preprocessing and restricted during search, 3: full always, 4: full in preprocessing) (command line option -dee)
    static int DEE_; ///< \internal do not use
    static int nbDecisionVars; ///< \brief tree search by branching only on the first variables having a lexicographic order position below a given value, assuming the remaining variables are completely assigned by this first group of variables (0: branch on all variables) (command line option -var)
    static int lds; ///< \brief iterative limited discrepancy search (0: no LDS), use a negative value to stop the search after the given absolute number of discrepancies has been explored (command line option -l)
    static bool limited; ///< \internal do not use
    static Long restart; ///< \brief randomly breaks ties in variable ordering heuristics and Luby restarts until a given number of search nodes (command line option -L)
    static Long backtrackLimit; ///< \brief limit on the number of backtracks (command line option -bt)
    static externalevent setvalue; ///< \internal do not use
    static externalevent setmin; ///< \internal do not use
    static externalevent setmax; ///< \internal do not use
    static externalevent removevalue; ///< \internal do not use
    static externalcostevent setminobj; ///< \internal do not use
    static externalsolution newsolution; ///< \internal do not use
    static Pedigree* pedigree; ///< \internal do not use
    static Haplotype* haplotype; ///< \internal do not use
    static string map_file; ///< \internal do not use
    static bool cfn; ///< \internal do not use
    static bool gz; ///< \internal do not use
    static bool bz2; ///< \internal do not use
    static bool xz; ///< \internal do not use
    static bool bayesian; ///< \internal do not use
    static int uai; ///< \internal do not use
    static int resolution; ///< \brief defines the number of digits that should be representable in UAI/OPB/QPBO/WBO formats (command line option -precision)
    static bool resolution_Update; ///< \internal flag true when default cfn precision is modified
    static TProb errorg; ///< \internal do not use
    static TLogProb NormFactor; ///< \internal do not use
    static int foundersprob_class; ///< \internal do not use
    static vector<TProb> allelefreqdistrib; ///< \internal do not use
    static bool consecutiveAllele; ///< \internal do not use
    static bool generation; ///< \internal do not use
    static int pedigreeCorrectionMode; ///< \internal do not use
    static int pedigreePenalty; ///< \internal do not use
    static int vac; ///< \brief enforces VAC at each search node having a search depth less than the absolute value of a given value (0: no VAC, 1: VAC in preprocessing, >1: VAC during search up to a given search depth), if given a negative value then VAC is not performed inside depth-first search of hybrid best-first search method (command line option -A and possibly -hbfs)
    static int vac_prev; // allows to disconnect soft local consistency temporally
    static string costThresholdS; ///< \brief threshold cost value for VAC in CFN format (command line option -T)
    static string costThresholdPreS; ///< \brief in preprocessing, threshold cost value for VAC in CFN format (command line option -P)
    static Cost costThreshold; ///< \brief threshold cost value for VAC (command line option -T)
    static Cost costThresholdPre; ///< \brief in preprocessing, threshold cost value for VAC (command line option -P)
    static Double trwsAccuracy; ///< \brief in preprocessing , enforces TRW-S until a given accuracy is reached (command line option -trws)
    static bool trwsOrder; ///< \brief replaces DAC order by Kolmogorov's TRW-S order (command line option --trws-order)
    static unsigned int trwsNIter; ///< \brief enforces at most n iterations of TRW-S (command line option --trws-n-iters)
    static unsigned int trwsNIterNoChange; ///< \brief stops TRW-S when n iterations did not change the lower bound (command line option --trws-n-iters-no-change)
    static unsigned int trwsNIterComputeUb; ///< \brief computes an upper bound every n steps in TRW-S (command line option --trws-n-iters-compute-ub)
    static Double costMultiplier; ///< \brief multiplies all costs internally by this number when loading a problem in WCSP format (command line option -C)
    static Cost costMultiplier_; ///< \internal do not use (should be set by setCostMultiplier)
    static unsigned int decimalPoint; ///< \internal do not use
    static string deltaUbS; ///< \brief stops search if the absolute optimality gap reduces below a given value in CFN format (command line option -agap)
    static Cost deltaUb; ///< \internal do not use
    static Cost deltaUbAbsolute; ///< \brief stops search if the absolute optimality gap reduces below a given value (command line option -agap)
    static Double deltaUbRelativeGap; ///< \brief stops search if the relative optimality gap reduces below a given value (command line option -rgap)
    static int singletonConsistency; ///< \brief in preprocessing, performs singleton soft local consistency (command line option -S)
    static bool vacValueHeuristic; ///< \brief VAC-based value ordering heuristic (command line options -V and -A)
    static BEP* bep; ///< \internal do not use
    static LcLevelType LcLevel; ///< \brief soft local consistency level (0: NC, 1: AC, 2: DAC, 3: FDAC, 4: EDAC) (command line option -k)
    static LcLevelType LcLevel_prev; // allows to disconnect soft local consistency temporally
    static int maxEACIter; ///< \brief maximum number of iterations in EDAC before switching to FDAC
    static bool wcnf; ///< \internal do not use
    static bool qpbo; ///< \internal do not use
    static Double qpboQuadraticCoefMultiplier; ///< \brief defines coefficient multiplier for quadratic terms in QPBO format (command line option -qpmult)
    static bool opb; ///< \internal do not use
    static bool lp; ///< \internal do not use

    static int addAMOConstraints; ///< \brief automatically detects and adds at-most-one constraints to existing knapsack constraints
    static bool addAMOConstraints_; ///< \brief automatically detects and adds at-most-one constraints to existing knapsack constraints
    static int knapsackDP; ///< \brief solves exactly knapsack constraints using dynamic programming (at every search node or less often)
    static bool VAClin; ///< \brief solves exactly knapsack constraints using dynamic programming (at every search node or less often)

    static unsigned int divNbSol; ///< \brief upper bound on the number of diverse solutions (0: no diverse solution) (keep it small as it controls model size)
    static unsigned int divBound; ///< \brief minimum Hamming distance between diverse solutions (command line options -div and -a)
    static unsigned int divWidth; ///< \brief adds a global MDD constraint with a given maximum relaxed width for finding diverse solutions (command line option -mdd)
    static unsigned int divMethod; ///< \brief diversity encoding method (0: Dual, 1: Hidden, 2: Ternary, 3: Knapsack) (command line option -divm)
    static unsigned int divRelax; ///< \brief MDD relaxation heuristic (0: random, 1: high diversity, 2: small diversity, 3: high unary costs) (command line option -mddh)

    static char* varOrder; ///< \brief variable elimination order for DAC, BTD, and VNS methods (0: lexicographic ordering, -1: maximum cardinality search ordering, -2: minimum degree ordering, -3: minimum fill-in ordering, -4: maximum spanning tree ordering, -5: reverse Cuthill-Mckee ordering, -6: approximate minimum degree ordering, -7: same as 0, 8: lexicographic ordering using variable names, string: variable ordering filename) (command line option -O)
    static int btdMode; ///< \brief tree search exploiting tree/path decomposition (0: no tree decomposition, 1: BTD with tree decomposition, 2: RDS-BTD with tree decomposition, 3: RDS-BTD with path decomposition) (command line option -B)
    static int btdSubTree; ///< \brief in RDS-BTD, cluster index for solving only this particular rooted cluster subtree (command line option -I)
    static int btdRootCluster; ///< \brief chooses the root cluster index (command line option -R)
    static int rootHeuristic; ///< \brief root cluster heuristic (0: maximum size, 1: maximum ratio of size by height-size, 2: minimum ratio of size by height-size, 3: minimum height) (command line option -root)
    static bool reduceHeight; ///< \brief minimize cluster tree height when searching for the root cluster (command line option -minheight)

    static bool maxsateval; ///< \internal do not use
    static bool xmlflag; ///< \internal do not use
    static bool xmlcop; ///< \internal do not use
    static TLogProb markov_log; ///< \internal do not use
    static string evidence_file; ///< \internal do not use
    static FILE* solution_uai_file; ///< \internal do not use
    static string solution_uai_filename; ///< \internal do not use
    static string problemsaved_filename; ///< \internal do not use
    static bool isZ; ///< \brief computes logarithm of probability of evidence (a.k.a. log-partition function) in UAI format (command line option -logz)
    static TLogProb logZ; ///< \internal do not use (lower bound on log-partition function)
    static TLogProb logU; ///< \internal do not use (upper bound on rejected potentials)
    static TLogProb logepsilon; ///< \brief approximation factor for computing the log-partition function (command line option -epsilon)
    static Double epsilon; ///< \brief floating-point epsilon
    static bool uaieval; ///< \internal do not use
    static string stdin_format; ///< \brief file format used when reading a problem from a Unix pipe ("cfn", "wcsp", "uai", "LG", "cnf", "wcnf", "qpbo", "opb", "wbo", "lp") (command line option --stdin)

    static double startCpuTime; ///< \internal do not use
    static double startRealTime; ///< \internal do not use
    static double startRealTimeAfterPreProcessing; ///< \internal do not use (for parallel execution of HBFS and VNS)

    static int splitClusterMaxSize; ///< \brief splits large clusters into a chain of smaller embedded clusters with a number of proper variables less than a given value (command line option -j)
    static double boostingBTD; ///< \brief in BTD, merges recursively leaf clusters with their fathers if separator size smaller than ToulBar2::elimDegree, else in VNS, merges clusters if the ratio of number of separator variables by number of cluster variables is above a given threshold (command line option -E and possibly -e)
    static int maxSeparatorSize; ///< \brief merges recursively clusters with their fathers if separator size greater than a given threshold (command line option -r)
    static int minProperVarSize; ///< \brief merges recursively clusters with their fathers if the number of proper variables is less than a given threshold (command line option -X)

    static bool heuristicFreedom; ///< \brief merges clusters automatically to give more freedom to variable ordering heuristics in BTD methods (command line option -F)
    static int heuristicFreedomLimit; ///< \brief stops merging a cluster subtree during BTD search if we tried repeatedly to solve this cluster for the same separator assignment more than a given number of times (-1: no merging) (command line option -F)

    static bool Berge_Dec; ///< \internal do not use
    static bool learning; ///< \internal do not use
    static externalfunc timeOut; ///< \internal do not use
    static std::atomic<bool> interrupted; ///< \internal do not use
    static int seed; ///< \brief initial random seed value, or use current time if a negative value is given (command line option -seed)
    static Double sigma; ///< \brief initial random noise standard deviation to be added to energy values when reading UAI format files (command line option -sigma)

    static string incop_cmd; ///< \brief in preprocessing, executes INCOP local search method to produce a better initial upper bound (default parameter string value "0 1 3 idwa 100000 cv v 0 200 1 0 0", see INCOP user manual http://imagine.enpc.fr/~neveub/incop/incop1.1/usermanual.ps)  (command line option -i)
    static string pils_cmd; ///< \brief in preprocessing, executes PILS local search method to produce a better initial upper bound (default parameter string value "3 0 0.333 150 150 1500 0.1 0.5 0.1 0.1", see PILS article https://doi.org/10.1002/prot.26174)  (command line option -pils)
    static string lrBCD_cmd; ///< \brief in preprocessing, executes LR-BCD to produce a better initial upper bound (default parameter string value "5 -2 3") (command line option -lrBCD)

    static SearchMethod searchMethod; ///< \brief chooses between tree search and variable neighborhood search methods (0: tree search, 1: sequential unified VNS, 2: sequential unified decomposition guided VNS, 3: synchronous parallel UDGVNS, 4: asynchronous parallel UDGVNS, 5: tree decomposition heuristic) (command line option -vns)

    static string clusterFile; ///< \brief cluster tree decomposition filename in COV or DEC format (with or without running intersection property)
    static ofstream vnsOutput; ///< \internal do not use

    static VNSSolutionInitMethod vnsInitSol; ///< \brief initial solution for VNS-like methods (-1: random, -2: minimum domain values, -3: maximum domain values, -4: first solution found by DFS, >=0: or by LDS with at most n discrepancies (command line option -vnsini)
    static int vnsLDSmin; ///< \brief minimum discrepancy value for VNS-like methods (command line option -ldsmin)
    static int vnsLDSmax; ///< \brief maximum discrepancy value for VNS-like methods (command line option -ldsmax)
    static VNSInc vnsLDSinc; ///< \brief discrepancy increment strategy for VNS-like methods (1: Increment by 1, 2: Multiply by 2, 3: Luby operator) (command line option -ldsinc)
    static int vnsKmin; ///< \brief minimum neighborhood size for VNS-like methods (command line option -kmin)
    static int vnsKmax; ///< \brief maximum neighborhood size for VNS-like methods (command line option -kmax)
    static VNSInc vnsKinc; ///< \brief neighborhood size increment strategy for VNS-like methods (1: Increment by 1, 2: Multiply by 2, 3: Luby operator, 4: Increment by 1 until maximum cluster size then considers all variables) (command line option -kinc)

    static int vnsLDScur; ///< \internal do not use (current LDS discrepancy value, used only for debugging display)
    static int vnsKcur; ///< \internal do not use (current neighborhood size, used only for debugging display)
    static VNSVariableHeuristic vnsNeighborVarHeur; ///< \brief neighborhood heuristic method (0: random variables, 1: variables in conflict, 2: connected variables in conflict, 3: random cluster, 4: variables in conflict with maximum degree, 5: sorted cluster, 6: sorted cluster separator, 7: similar to 6, 8: randomized root cluster, 9: variables in partial conflict)
    static bool vnsNeighborChange; ///< \internal do not use // (in RADGVNS, true if change neighborhood cluster only when not improved)
    static bool vnsNeighborSizeSync; ///< \internal do not use // (in RADGVNS, true if neighborhood size is synchronized)
    static bool vnsParallelLimit; ///< \internal do not use // (in RSDGVNS and RADGVNS, true if number of parallel slaves limited by number of clusters)
    static bool vnsParallelSync; ///< \internal do not use // (true if RSGDVNS else RADGVNS)
    static string vnsOptimumS; ///< \brief stops VNS if a solution is found with a given cost (or better) in CFN format (command line option -best)
    static Cost vnsOptimum; ///< \brief stops VNS if a solution is found with a given cost (or better) (command line option -best)

    static bool parallel; ///< \brief parallel mode for tree search and VNS (see mpirun toulbar2)

    static Long hbfs; ///< \brief performs hybrid best-first search with a given limit in the number of backtracks for depth-first search before visiting another open node (0: always DFS, 1: HBFS) (related to command line option -hbfs)
    static Long hbfsGlobalLimit; ///< \brief restarts BTD-HBFS from the root cluster after a given number of backtracks (command line option -hbfs)
    static Long hbfsAlpha; ///< \brief minimum recomputation node redundancy percentage threshold value (command line option -hbfsmin)
    static Long hbfsBeta; ///< \brief  maximum recomputation node redundancy percentage threshold value (command line option -hbfsmax)
    static ptrdiff_t hbfsCPLimit; ///< \brief maximum number of stored choice points before switching to normal DFS (warning! should always be cast to std::size_t when used, so that -1 is equivalent to +infinity)
    static ptrdiff_t hbfsOpenNodeLimit; ///< \brief maximum number of stored open nodes before switching to normal DFS (command line option -open) (warning! should always be cast to std::size_t when used, so that -1 is equivalent to +infinity)
    static Long sortBFS; ///< \brief number of visited open nodes before sorting the remaining open nodes (command line option -sopen)
#ifdef OPENMPI
    static bool burst; ///< \brief in parallel HBFS, workers send open nodes as soon as possible to the master (command line option -burst)
#endif
    static Long eps; ///< \brief performs HBFS until a given number of open nodes are collected and exits (command line option -eps)
    static string epsFilename; ///< \brief a given filename to print remaining valid (lower bound less than current upper bound) open nodes as partial assignments before exits (command line option -eps)

    static bool verifyOpt; ///< \brief compiled in debug, checks if a given (optimal) solution is never pruned by propagation when the current upper bound is greater than the cost of this solution (see Solver::read_solution, related to command line option -opt)
    static Cost verifiedOptimum; ///< \brief compiled in debug, a given (optimal) solution cost (see Solver::read_solution, related to command line option -opt)

    static int bilevel; ///< \brief bilevel optimization using modified BTD with (at least) four clusters (P0 -> P1, P0 -> P2, P0 -> NegP2) corresponding to the restricted leader problem (P0 and P1), the follower problem (P2), and the negative follower problem (NegP2) (command line option -bilevel)
    static vector<unsigned int> decimalPointBLP;
    static vector<Double> costMultiplierBLP;
    static vector<Cost> negCostBLP;
    static vector<Cost> initialLbBLP;
    static vector<Cost> initialUbBLP;

    static Double getCostMultiplier() {
        assert(costMultiplier_ == max(UNIT_COST, (Cost)floorl(std::abs(costMultiplier))));
        return costMultiplier;
    }
    static Cost getCostMultiplierInt() {
        assert(costMultiplier_ == max(UNIT_COST, (Cost)floorl(std::abs(costMultiplier))));
        return costMultiplier_;
    }
    static void setCostMultiplier(Double mult) {
        costMultiplier = mult;
        costMultiplier_ = max(UNIT_COST, (Cost)floorl(std::abs(costMultiplier)));
        assert(costMultiplier_ >= 1);
    }
};

#if defined(INT_COST) || defined(SHORT_COST) || defined(LONGLONG_COST)
inline Cost rounding(Cost lb)
{
    if (ToulBar2::costMultiplier_ > 1 && ((lb % ToulBar2::costMultiplier_) != MIN_COST)) {
        return lb + ToulBar2::costMultiplier_;
    } else {
        return lb;
    }
}
inline bool CUT(Cost lb, Cost ub)
{
    return (rounding(lb) + ToulBar2::deltaUb) >= ub;
}
inline bool CSP(Cost lb, Cost ub)
{
    return CUT(lb + UNIT_COST, ub);
}
#endif

/*
 * Backtrack exception
 *
 */

#ifdef ILOGLUE
extern IloSolver IlogSolver;
#define THROWCONTRADICTION ({if (ToulBar2::verbose >= 2) cout << "... contradiction!" << endl; if (ToulBar2::weightedDegree) conflict(); IlogSolver.fail(0); })
#else
class Contradiction {
public:
    Contradiction()
    {
        if (ToulBar2::verbose >= 2)
            cout << what() << endl;
    }
    const char* what() const { return "... contradiction!"; }
};
#define THROWCONTRADICTION            \
    {                                 \
        if (ToulBar2::weightedDegree) \
            conflict();               \
        throw Contradiction();        \
    }
#endif

/*
 * Internal classes and basic data structures used everywhere
 *
 */

class Store;
class Domain;
class Variable;
class IntervalVariable;
class EnumeratedVariable;
class Constraint;
class BinaryConstraint;
class TernaryConstraint;
class NaryConstraint;
class WCSP;
class Solver;
class Cluster;
class Separator;
class TreeDecomposition;
class VACExtension;

class ConstraintLink {
public:
    Constraint* constr;
    int scopeIndex;
};

class WCSPLink {
public:
    WCSP* const wcsp;
    int wcspIndex;
    WCSPLink(WCSP* w, int index)
        : wcsp(w)
        , wcspIndex(index)
    {
    }
};

/// \brief allows one to sort pointers to WCSPLink objects (Constraints or Variables) by their wcspIndex rather than pointer values
template <class T>
bool compareWCSPIndex(const T* lhs, const T* rhs)
{
    assert(lhs);
    assert(rhs);
    int left = lhs->wcspIndex;
    int right = rhs->wcspIndex;
    if (left < 0)
        left = MAX_ELIM_BIN - left; // makes elimTernConstraints after elimBinConstraints after original constraints
    if (right < 0)
        right = MAX_ELIM_BIN - right;
    return left < right;
}

template <class T>
struct Compare {
    typedef bool (*compare_t)(const T*, const T*);
};

template <class T>
class Set : public set<T*, typename Compare<T>::compare_t> {
public:
    Set()
        : set<T*, typename Compare<T>::compare_t>(compareWCSPIndex<T>)
    {
    }
};

typedef Set<Constraint> ConstraintSet;
typedef Set<Variable> VariableSet;

// For incremental diverse solution search - relaxed constraint
typedef vector<vector<vector<vector<Cost>>>> Mdd; // mdd[layer][source][target][value] = label weight if source---val--->target exists, getUb otherwise

#endif /*TB2TYPES_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
