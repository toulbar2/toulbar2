/** \file tb2types.hpp
 *  \brief Macros, types, and globals.
 * 
 */

#ifndef TB2TYPES_HPP_
#define TB2TYPES_HPP_

//#define INT_COST
#define LONGLONG_COST
//#define PARETOPAIR_COST

//#define DOUBLE_PROB
#define LONGDOUBLE_PROB

#include "tb2utils.hpp"

typedef int Value;

const Value MAX_VAL = (INT_MAX / 2);
const Value MIN_VAL = -(INT_MAX / 2);
const Value MAX_DOMAIN_SIZE  = 2000;

#ifdef INT_COST
const bool PARTIALORDER = false;
typedef int Cost;
const Cost MIN_COST = 0;
const Cost UNIT_COST = 1;
const Cost SMALL_COST = 1;
const Cost MEDIUM_COST = 3;
const Cost LARGE_COST = 100;
const Cost MAX_COST = ((INT_MAX / 2) / MEDIUM_COST);
inline Cost MIN(Cost a, Cost b) {return min(a,b);}
inline Cost MAX(Cost a, Cost b) {return max(a,b);}
inline Cost GLB(Cost a, Cost b) {return MIN(a,b);}
inline Cost LUB(Cost a, Cost b) {return MAX(a,b);}
inline bool GLB(Cost *a, Cost b) {if (b < *a) {*a = b; return true;} else return false;}
inline bool LUB(Cost *a, Cost b) {if (b > *a) {*a = b; return true;} else return false;}
inline bool GLBTEST(Cost a, Cost b) {return (b < a);}
inline bool LUBTEST(Cost a, Cost b) {return (b > a);}
inline bool DACTEST(Cost a, Cost b) {return (a==0 && b>0);}
inline bool SUPPORTTEST(Cost a, Cost b) {return false;}
inline bool SUPPORTTEST(Cost a) {return false;}
inline bool CUT(Cost lb, Cost ub) {return lb >= ub;}
inline bool CSP(Cost lb, Cost ub) {return (ub - lb) <= 1;}
inline void initCosts(Cost ub) {}
#endif

#ifdef LONGLONG_COST
const bool PARTIALORDER = false;
typedef Long Cost;
const Cost MIN_COST = 0;
const Cost UNIT_COST = 1;
const Cost SMALL_COST = 1;
const Cost MEDIUM_COST = 3;
const Cost LARGE_COST = 100;
const Cost MAX_COST = ((LONGLONG_MAX / 2) / MEDIUM_COST);
inline Cost MIN(Cost a, Cost b) {return min(a,b);}
inline Cost MAX(Cost a, Cost b) {return max(a,b);}
inline Cost GLB(Cost a, Cost b) {return MIN(a,b);}
inline Cost LUB(Cost a, Cost b) {return MAX(a,b);}
inline bool GLB(Cost *a, Cost b) {if (b < *a) {*a = b; return true;} else return false;}
inline bool LUB(Cost *a, Cost b) {if (b > *a) {*a = b; return true;} else return false;}
inline bool GLBTEST(Cost a, Cost b) {return (b < a);}
inline bool LUBTEST(Cost a, Cost b) {return (b > a);}
inline bool DACTEST(Cost a, Cost b) {return (a==0 && b>0);}
inline bool SUPPORTTEST(Cost a, Cost b) {return false;}
inline bool SUPPORTTEST(Cost a) {return false;}
inline bool CUT(Cost lb, Cost ub) {return lb >= ub;}
inline bool CSP(Cost lb, Cost ub) {return (ub - lb) <= 1;}
inline void initCosts(Cost ub) {}
#endif

#ifdef PARETOPAIR_COST
const bool PARTIALORDER = true;
#include "tb2paretopair.hpp"
typedef ParetoPair Cost;
const Cost MIN_COST = PARETOPAIR_MIN;
const Cost UNIT_COST = PARETOPAIR_1;
const Cost SMALL_COST = PARETOPAIR_1;
const Cost MEDIUM_COST = PARETOPAIR_3;
const Cost LARGE_COST = PARETOPAIR_100;
const Cost MAX_COST = PARETOPAIR_MAX;
#endif

#ifdef DOUBLE_PROB
typedef double TProb;
#endif

#ifdef LONGDOUBLE_PROB
typedef Double TProb;
#endif

const int STORE_SIZE = 16;
#define INTEGERBITS (8*sizeof(Cost)-2)

typedef map<int,int> TSCOPE;

#ifdef NARYCHAR
#define CHAR_FIRST 'A'
#else
#define CHAR_FIRST 1
#endif


/*
 * Global variables encapsulated as static members
 * 
 */


typedef void (*externalevent)(int wcspId, int varIndex, Value value);
typedef void (*externalcostevent)(int wcspId, int varIndex, Cost cost);

typedef enum {ELIM_NONE = 0, MAX_CARD = 1, MIN_FILL = 2, MIN_DEGREE = 3, ELIM_MAX } ElimOrderType;

class Pedigree;
class BEP;

typedef enum {LC_NC = 0, LC_AC = 1, LC_DAC = 2, LC_FDAC = 3, LC_EDAC = 4, LC_MAX } LcLevelType;

class ToulBar2
{
protected:
    virtual ~ToulBar2() = 0;	// Trick to avoid any instantiation of ToulBar2
public:
    static double version;
    static int verbose;
    static bool showSolutions;
    static bool writeSolution;
    static bool allSolutions;
    static bool binaryBranching;
    static bool dichotomicBranching;
    static unsigned int dichotomicBranchingSize;
    static int  elimDegree; 
    static int  elimDegree_preprocessing;
    static int  elimDegree_; 
    static int  elimDegree_preprocessing_;
    static int minsumDiffusion;
    static bool preprocessTernary;
    static bool preprocessTernaryHeuristic;
    static bool QueueComplexity;
    static bool lastConflict;
    static bool weightedDegree;
    static bool lds;
    static bool limited;
    static externalevent setvalue;
    static externalevent setmin;
    static externalevent setmax;
    static externalevent removevalue;
    static externalcostevent setminobj;
    static Pedigree *pedigree;
    static bool bayesian;
    static int resolution;
    static TProb errorg;
    static TProb NormFactor;
    static int foundersprob_class; 
    static vector<TProb> allelefreqdistrib;
    static bool consecutiveAllele;
    static bool generation;
    static int pedigreeCorrectionMode;
    static int vac;
    static Cost costThreshold;
    static Cost costMultiplier;
    static Cost relaxThreshold;
    static ElimOrderType elimOrderType;
    static bool singletonConsistency;
    static BEP *bep;
    static LcLevelType LcLevel;
};

/*
 * Backtrack exception
 * 
 */

#ifdef ILOGLUE
extern IloSolver IlogSolver;
#define THROWCONTRADICTION ({if (ToulBar2::verbose >= 2) cout << "... contradiction!" << endl; IlogSolver.fail(0);})
#else
class Contradiction
{
public:
    Contradiction() {if (ToulBar2::verbose >= 2) cout << "... contradiction!" << endl;}
};
#define THROWCONTRADICTION (/* conflict(), */ throw Contradiction())
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
class WCSP;
class Solver;

struct ValueCost
{
	Value value;
	Cost cost;
};

struct ConstraintLink 
{
    Constraint *constr;
    int scopeIndex;
};

class WCSPLink 
{
public:
    WCSP * const wcsp;
    int wcspIndex;
    WCSPLink(WCSP *w, int index) : wcsp(w), wcspIndex(index) {}
};

#endif /*TB2TYPES_HPP_*/
