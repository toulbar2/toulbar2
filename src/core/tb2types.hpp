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
 * - ::String : string with extended character size (wide string) to encode tuples
 *
 * \note Compilation flag for Cost is: \c INT_COST (int), \c LONGLONG_COST (long long), or \c PARETOPAIR_COST (see ::ParetoPair)
 * \warning \c PARETOPAIR_COST is fragile.
 * \note Compilation flag for TProb is: \c DOUBLE_PROB or \c LONGDOUBLE_PROB
 * \note Compilation flag for T(Log)Prob is: \c DOUBLE_PROB or \c LONGDOUBLE_PROB
 * \note Compilation flag for String is: \c WIDE_STRING or nothing (usual C++ string)
 */

#ifndef TB2TYPES_HPP_
#define TB2TYPES_HPP_

//#define INT_COST
//#define LONGLONG_COST
//#define PARETOPAIR_COST

//#define DOUBLE_PROB
//#define LONGDOUBLE_PROB

// uncomment if using large enumerated domains with BTD or in nary cost functions
//#define WIDE_STRING

#include "utils/tb2utils.hpp"
//Must be included after tb2utils.hpp
#include "utils/tb2integer.hpp"
#ifdef QUAD_PROB
#include <quadmath.h> // only with gcc/g++
#endif

/// Domain value (can be positive or negative integers)
typedef int Value;
/// Maximum domain value
const Value MAX_VAL = (INT_MAX / 2);
/// Forbidden domain value
const Value WRONG_VAL = INT_MAX;
/// Minimum domain value
const Value MIN_VAL = -(INT_MAX / 2);
/// Maximum domain size
/// \deprecated Should use WCSP::getMaxDomainSize instead.
const Value MAX_DOMAIN_SIZE = 2000;

// For very large domains with ternary cost functions, use NARYPROJECTIONSIZE=2 instead of 3
const int NARYPROJECTIONSIZE = 3; // limit on the number of unassigned variables before nary constraints are projected to smaller-arity constraint (should be between 1 and 3)
const Long NARYDECONNECTSIZE = 4; // maximum number of initial tuples in nary constraints in order to check for its removal (if it is always satisfied by current domains)

const int MAX_BRANCH_SIZE = 1000000;
const ptrdiff_t CHOICE_POINT_LIMIT = SIZE_MAX - MAX_BRANCH_SIZE;
const ptrdiff_t OPEN_NODE_LIMIT = SIZE_MAX;

#ifdef INT_COST
const bool PARTIALORDER = false;
typedef int Cost;
const Cost MIN_COST = 0;
const Cost UNIT_COST = 1;
const Cost SMALL_COST = 1;
const Cost MEDIUM_COST = 3;
const Cost LARGE_COST = 100;
const Cost MAX_COST = ((INT_MAX / 2) / MEDIUM_COST / MEDIUM_COST);
//inline bool Add(Cost a, Cost b, Cost* c) { return __builtin_sadd_overflow(a, b, c); }
//inline bool Sub(Cost a, Cost b, Cost* c) { return __builtin_ssub_overflow(a, b, c); }
//inline bool Mul(Cost a, Cost b, Cost* c) { return __builtin_smul_overflow(a, b, c); }
inline Cost MIN(Cost a, Cost b) { return min(a, b); }
inline Cost MAX(Cost a, Cost b) { return max(a, b); }
inline Cost MULT(Cost a, double b)
{
    assert(b < MAX_COST);
    if (a >= MAX_COST)
        return MAX_COST;
    else if (b <= UNIT_COST)
        return a * b;
    else if (a < MAX_COST / b)
        return a * b;
    else {
        cerr << "Error: cost multiplication overflow!" << endl;
        exit(1);
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
//inline bool Add(Cost a, Cost b, Cost* c) { return __builtin_saddll_overflow(a, b, c); }
//inline bool Sub(Cost a, Cost b, Cost* c) { return __builtin_ssubll_overflow(a, b, c); }
//inline bool Mul(Cost a, Cost b, Cost* c) { return __builtin_smulll_overflow(a, b, c); }
inline Cost MIN(Cost a, Cost b) { return min(a, b); }
inline Cost MAX(Cost a, Cost b) { return max(a, b); }
inline Cost MULT(Cost a, double b)
{
    assert(b < MAX_COST);
    if (a >= MAX_COST)
        return MAX_COST;
    else if (b <= UNIT_COST)
        return a * b;
    else if (a < MAX_COST / b)
        return a * b;
    else {
        cerr << "Error: cost multiplication overflow!" << endl;
        exit(1);
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
typedef __float128 TProb;
typedef __float128 TLogProb;
inline Cost Round(TLogProb f) { return (Cost)roundq(f); }
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
        cerr << "My oh my ! Try to Logarithm a negative number" << endl;
        exit(0);
    }
}

const int STORE_SIZE = 16;
#define INTEGERBITS (8 * sizeof(Cost) - 2)

const int MAX_ELIM_BIN = 1000000000;
const int MAX_ARITY = 1000;
/// Maximum number of tuples in n-ary cost functions
const int MAX_NB_TUPLES = 1000000;
const int LARGE_NB_VARS = 10000;

const int DECIMAL_POINT = 3; // default number of digits after decimal point for printing floating-point values

typedef map<int, int> TSCOPE;
typedef map<int, Value> TAssign;

#ifdef NARYCHAR
#define CHAR_FIRST '0'
#else
#define CHAR_FIRST 1
#endif

typedef unsigned int uint;
/*
 * Abstract data type that help post a global cost function
 *
 */

// A value with weight
template <class Object>
struct WeightedObj {
    Object val;
    Cost weight;

    WeightedObj(const Object& val_, Cost weight_ = MIN_COST)
        : val(val_)
        , weight(weight_)
    {
    }
};

// A value with upper and lower limit
template <class Object>
struct BoundedObj {
    Object val;
    unsigned int upper;
    unsigned int lower;

    BoundedObj(const Object& val_, unsigned int upper_, unsigned int lower_ = 0)
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
    unsigned int weight;

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
    unsigned int weight;
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
    ELIM_MAX
} ElimOrderType;

class Seq;
class Pedigree;
class Haplotype;
class BEP;
class Cpd;
class Tb2ScpBranch;
class TrieNum;
#ifdef OPENMPI
class BaseJobs;
class Jobs;
#endif
class SequenceHandler;

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

struct ValueCost {
    Value value;
    Cost cost;
    friend bool operator<(const ValueCost& u, const ValueCost& v) { return u.cost < v.cost; }
    friend bool operator>(const ValueCost& u, const ValueCost& v) { return u.cost > v.cost; }
    friend bool operator==(const ValueCost& u, const ValueCost& v) { return u.cost == v.cost; }
};



///contains all global variables (mainly solver's command-line options)
namespace ToulBar2 {
    extern string version;
    extern int verbose;
    extern int debug;
    extern string externalUB;
    extern int showSolutions;
    extern char* writeSolution;
    extern FILE* solutionFile;
    extern Long allSolutions;
    extern int dumpWCSP;
    extern bool approximateCountingBTD;
    extern bool binaryBranching;
    extern int dichotomicBranching;
    extern unsigned int dichotomicBranchingSize;
    extern bool sortDomains;
    extern map<int, ValueCost*> sortedDomains;
    extern bool solutionBasedPhaseSaving;
    extern int elimDegree;
    extern int elimDegree_preprocessing;
    extern int elimDegree_;
    extern int elimDegree_preprocessing_;
    extern int elimSpaceMaxMB;
    extern int minsumDiffusion;
    extern int prodsumDiffusion;
    extern int preprocessTernaryRPC;
    extern int preprocessFunctional;
    extern bool costfuncSeparate;
    extern int preprocessNary;
    extern bool QueueComplexity;
    extern bool Static_variable_ordering; // flag for static variable ordering during search (dynamic ordering is default value)
    extern bool lastConflict;
    extern int weightedDegree;
    extern int weightedTightness;
    extern bool MSTDAC;
    extern int DEE;
    extern int DEE_;
    extern int nbDecisionVars;
    extern int lds;
    extern bool limited;
    extern Long restart;
    extern externalevent setvalue;
    extern externalevent setmin;
    extern externalevent setmax;
    extern externalevent removevalue;
    extern externalcostevent setminobj;
    extern externalsolution newsolution;
    extern Pedigree* pedigree;
    extern Haplotype* haplotype;
    extern Cpd* cpd;
    extern Tb2ScpBranch* scpbranch;
#ifdef OPENMPI
    extern BaseJobs* jobs;
#endif
    extern SequenceHandler* sequence_handler;
    extern string map_file;
    extern bool cfn;
    extern bool gz;
    extern bool xz;
    extern bool bayesian;
    extern int uai;
    extern int resolution;
    extern TProb errorg;
    extern TLogProb NormFactor;
    extern int foundersprob_class;
    extern vector<TProb> allelefreqdistrib;
    extern bool consecutiveAllele;
    extern bool generation;
    extern int pedigreeCorrectionMode;
    extern int pedigreePenalty;
    extern int vac;
    extern string costThresholdS;
    extern string costThresholdPreS;
    extern Cost costThreshold;
    extern Cost costThresholdPre;
    extern double trwsAccuracy;
    extern bool trwsOrder;
    extern unsigned int trwsNIter;
    extern unsigned int trwsNIterNoChange;
    extern unsigned int trwsNIterComputeUb;
    extern double costMultiplier;
    extern unsigned int decimalPoint;
    extern string deltaUbS;
    extern Cost deltaUb;
    extern bool singletonConsistency;
    extern bool vacValueHeuristic;
    extern BEP* bep;
    extern LcLevelType LcLevel;
    extern bool wcnf;
    extern bool qpbo;
    extern double qpboQuadraticCoefMultiplier;

    extern char* varOrder;
    extern int btdMode;
    extern int btdSubTree;
    extern int btdRootCluster;

    extern bool maxsateval;
    extern bool xmlflag;
    extern TLogProb markov_log;
    extern string evidence_file;
    extern FILE* solution_uai_file;
    extern string solution_uai_filename;
    extern string problemsaved_filename;
    extern bool isZ;
    extern float isZCelTemp;
    extern int isZUB;
    extern Cost enumUB;
    extern bool bestconf;
    extern bool stop; // STOP TB2 (TMP way)
    extern bool isGumbel;
    extern Seq* seq;
    extern TLogProb logZ;
    extern TLogProb Entropy;
    extern TLogProb Enthalpy;
    extern TLogProb GlobalLogUbZ; // Upper bound on Z.
    extern TLogProb GlobalLogLbZ; // Upper bound on Z.
    extern TLogProb logU; // upper bound on rejected potentials
    extern TLogProb logepsilon; // epsilon for Z* pruning
    extern TProb sigma; // sigma set for HBFS-Counting
    extern TrieNum* trieZ; // Trie over preprocessing Optimum Energies
    extern string Trie_File;
    extern bool isTrie_File;

    extern TLogProb ubE;
    extern TLogProb Normalizing_Constant;

    extern bool uaieval;
    extern string stdin_format; // stdin format declaration

    extern double startCpuTime;

    extern int splitClusterMaxSize;
    extern double boostingBTD;
    extern int maxSeparatorSize;
    extern int minProperVarSize;
    extern int smallSeparatorSize;

    extern int Berge_Dec; // flag for berge acyclic decomposition
    extern int nbvar; // initial number of variable (read in the file)
    extern bool learning; // if true, perform pseudoboolean learning
    extern externalfunc timeOut;
    extern bool interrupted;
    extern int seed;

    extern string incop_cmd;

    extern SearchMethod searchMethod;

    extern string clusterFile; // cluster tree decomposition file (without running intersection property)
    extern ofstream vnsOutput; // output file for VNS

    extern VNSSolutionInitMethod vnsInitSol; // initial solution strategy (search with max discrepancy limit if positive value)
    extern int vnsLDSmin; // discrepancy initial value
    extern int vnsLDSmax; // discrepancy maximum value
    extern VNSInc vnsLDSinc; // discrepancy increment strategy inside VNS
    extern int vnsKmin; // neighborhood initial size
    extern int vnsKmax; // neighborhood maximum size
    extern VNSInc vnsKinc; // neighborhood size increment strategy inside VNS

    extern int vnsLDScur; // current discrepancy (used only for debugging display)
    extern int vnsKcur; // current neighborhood size (used only for debugging display)
    extern VNSVariableHeuristic vnsNeighborVarHeur; // variable heuristic to build a neighborhood (used to differentiate VNS/DGVNS)
    extern bool vnsNeighborChange; // true if change neighborhood cluster only when not improved (only in RADGVNS)
    extern bool vnsNeighborSizeSync; // true if neighborhood size is synchronized (only in RADGVNS)
    extern bool vnsParallelLimit; // true if number of parallel slaves limited by number of clusters (only in RSDGVNS and RADGVNS)
    extern bool vnsParallelSync; // true if RSGDVNS else RADGVNS
    extern string vnsOptimumS;
    extern Cost vnsOptimum; // stops VNS if solution found with this cost (or better)
    extern bool vnsParallel; // true if in master/slaves paradigm

    extern Long hbfs; // hybrid best-first search mode (used as a limit on the number of backtracks before visiting another open search node)
    extern Long hbfsGlobalLimit; // limit on the number of nodes before stopping the search on the current cluster subtree problem
    extern Long hbfsAlpha; // inverse of minimum node redundancy goal limit
    extern Long hbfsBeta; // inverse of maximum node redundancy goal limit
    extern ptrdiff_t hbfsCPLimit; // limit on the number of choice points stored inside open node list
    extern ptrdiff_t hbfsOpenNodeLimit; // limit on the number of open nodes

    extern bool verifyOpt; // if true, for debugging purposes, checks the given optimal solution (problem.sol) is not pruned during search
    extern Cost verifiedOptimum; // for debugging purposes, cost of the given optimal solution
};

#ifdef INT_COST
inline Cost rounding(Cost lb)
{
    return (((lb % max(UNIT_COST, (Cost)floor(ToulBar2::costMultiplier))) != MIN_COST) ? (lb + (Cost)floor(ToulBar2::costMultiplier)) : lb);
}
inline bool CUT(Cost lb, Cost ub) { return rounding(lb) >= ub; }
inline bool CSP(Cost lb, Cost ub) { return CUT(lb + UNIT_COST, ub); }
#endif

#ifdef LONGLONG_COST
inline Cost rounding(Cost lb)
{
    return (((lb % max(UNIT_COST, (Cost)floor(fabs(ToulBar2::costMultiplier)))) != MIN_COST) ? (lb + (Cost)floor(fabs(ToulBar2::costMultiplier))) : lb);
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
            cout << "... contradiction!" << endl;
    }
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

/// < \brief allows one to sort pointers to WCSPLink objects (Constraints or Variables) by their wcspIndex rather than pointer values
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

#endif /*TB2TYPES_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
