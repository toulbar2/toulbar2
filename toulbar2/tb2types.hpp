/** \file tb2types.hpp
 *  \brief Macros, types, and globals.
 * 
 */

#ifndef TB2TYPES_HPP_
#define TB2TYPES_HPP_

#include <assert.h>

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

#include "tb2utils.hpp"

#define Value int
#define MAX_VAL (INT_MAX / 2)
#define MIN_VAL -(INT_MAX / 2)

#define Cost int
#define MAX_COST (INT_MAX / 2)

#define MAX_DOMAIN_SIZE 1000

#define STORE_SIZE 16

/*
 * Global variables encapsulated as static members
 * 
 */

typedef void (*externalevent)(int wcspId, int varIndex, Value value);

class ToulBar2
{
protected:
    virtual void dummy() = 0;	// Trick to avoid any instantiation of ToulBar2
public:
    static int verbose;
    static bool showSolutions;
    static bool binaryBranching;
    static externalevent setvalue;
    static externalevent setmin;
    static externalevent setmax;
    static externalevent removevalue;
};

/*
 * Backtrack exception
 * 
 */

class Contradiction
{
public:
    Contradiction() {if (ToulBar2::verbose >= 2) cout << "... contradiction!" << endl;}
};

/*
 * Internal classes and basic data structures used everywhere
 * 
 */

class Store;
class Domain;
class Variable;
class Constraint;
class WCSP;
class Solver;

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
