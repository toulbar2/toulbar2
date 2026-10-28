/** 
 *  \file tb2config.hpp
 *  \brief Types definitions, automatically generated with cmake.
 */

#pragma once

// probabilities
/* #undef DOUBLE_PROB */
#define LONGDOUBLE_PROB
/* #undef QUAD_PROB */

// costs
/* #undef SHORT_COST */
/* #undef INT_COST */
#define LONGLONG_COST

// domain Values
/* #undef SHORT_VALUE */

// BOOST
#define BOOST

// BZ2 (only relevant when BOOST defined)
/* #undef NO_BZ2 */

// LZMA (only relevant when BOOST defined)
/* #undef NO_LZMA */

// XML XCSP Parser
/* #undef XMLFLAG */
#define XMLFLAG3

// MPI
/* #undef OPENMPI */

// CPLEX
/* #undef ILOGCPLEX */

// LR-BCD (Eigen)
/* #undef LR_BCD_BUILD */

// TOULBAR2_MULTITHREAD
/* #undef TOULBAR2_MULTITHREAD */

#ifdef TOULBAR2_MULTITHREAD
    #define TB2_THREAD_LOCAL thread_local
#else
    #define TB2_THREAD_LOCAL
#endif
