/** \file tb2types.hpp
 *  \brief Macros, types, and globals.
 * 
 */

#ifndef TB2TYPES_HPP_
#define TB2TYPES_HPP_

//#define INT_COST
#define LONGLONG_COST

//#define DOUBLE_PROB
#define LONGDOUBLE_PROB

#include "tb2utils.hpp"

typedef int Value;

const Value MAX_VAL = (INT_MAX / 2);
const Value MIN_VAL = -(INT_MAX / 2);
const Value MAX_DOMAIN_SIZE  = 1000;

#ifdef INT_COST
typedef int Cost;
const Cost MIN_COST = 0;
const Cost MAX_COST = (INT_MAX / 2);
#endif

#ifdef LONGLONG_COST
typedef Long Cost;
const Cost MIN_COST = 0;
const Cost MAX_COST = (LONGLONG_MAX / 2);
#endif

#ifdef DOUBLE_PROB
typedef double TProb;
#endif

#ifdef LONGDOUBLE_PROB
typedef Double TProb;
#endif

const int STORE_SIZE = 16;
#define INTEGERBITS (8*sizeof(Cost)-2)


#include <map>
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

class ToulBar2
{
protected:
    virtual ~ToulBar2() = 0;	// Trick to avoid any instantiation of ToulBar2
public:
    static double version;
    static int verbose;
    static bool showSolutions;
    static bool writeSolution;
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
    static bool FDAComplexity;
    static bool FDAC;
    static bool lastConflict;
    static bool lastWConflict;
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
    static Cost relaxThreshold;
    static Cost costConstant;
    static bool makeScaleCosts;
    static ElimOrderType elimOrderType;
    static bool singletonConsistency;
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
