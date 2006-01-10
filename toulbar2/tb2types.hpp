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
 
class ToulBar2
{
protected:
    virtual void dummy() = 0;	// Trick to avoid any instantiation of ToulBar2
public:
    static int verbose;
    static bool showSolutions;
    static bool binaryBranching;
};

/*
 * Classes and basic data structures used everywhere
 * 
 */

class Variable;
class CostVariable;
class Constraint;
class WCSP;
class Solver;

struct ConstraintLink 
{
    Constraint *constr;
    int scopeIndex;
};

struct WCSPLink 
{
    WCSP *wcsp;
    int wcspIndex;
};

#endif /*TB2TYPES_HPP_*/
