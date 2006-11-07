/** \file tb2types.hpp
 *  \brief Macros, types, and globals.
 * 
 */

#ifndef TB2TYPES_HPP_
#define TB2TYPES_HPP_

#ifdef ILOGLUE
#include <ilsolver/ilosolverint.h>
#else
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#endif

#include <vector>

using namespace std;

typedef long long Long;

#include "tb2utils.hpp"
#include "tb2rational.hpp"

typedef int Value;

const Value MAX_VAL = (INT_MAX / 2);
const Value MIN_VAL = -(INT_MAX / 2);

const Value MAX_DOMAIN_SIZE  = 1000;

typedef int Cost;
const Cost MAX_COST = (INT_MAX / 2);

//typedef Long Cost;
//const Cost MAX_COST = (LONG_LONG_MAX / 2);

//typedef Rational Cost;
//const Cost MAX_COST = RATIONAL_MAX;

const int STORE_SIZE = 16;

/*
 * Global variables encapsulated as static members
 * 
 */

typedef void (*externalevent)(int wcspId, int varIndex, Value value);
typedef void (*externalcostevent)(int wcspId, int varIndex, Cost cost);

class Pedigree;

class ToulBar2
{
protected:
    virtual ~ToulBar2() = 0;	// Trick to avoid any instantiation of ToulBar2
public:
    static double version;
    static int verbose;
    static bool showSolutions;
    static bool binaryBranching;
    static bool elimVarWithSmallDegree;
    static bool elimVarWithSmallDegree_; // flag activated after elimination data structures initialization
    static bool only_preprocessing;
    static bool preprocessTernary;
    static bool preprocessTernaryHeuristic;
    static bool FDAComplexity;
    static bool lastConflict;
    static externalevent setvalue;
    static externalevent setmin;
    static externalevent setmax;
    static externalevent removevalue;
    static externalcostevent setminobj;
    static Pedigree *pedigree;
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
    const int wcspIndex;
    WCSPLink(WCSP *w, int index) : wcsp(w), wcspIndex(index) {}
};

#endif /*TB2TYPES_HPP_*/
