/** \file tb2system.hpp
 *  \brief System dependent functions.
 *
 */

#ifndef TB2SYSTEM_HPP_
#define TB2SYSTEM_HPP_

#if __cplusplus > 199711L
#define FINAL final
#else
#define FINAL
#endif

#ifdef QUAD_PROB
#include <quadmath.h>
#endif

extern const char* PrintFormatProb;

double cpuTime(); ///< \brief return CPU time in seconds with high resolution (microseconds) if available
void timeOut(int sig);
void timer(int t); ///< \brief set a timer (in seconds)
void timerStop(); ///< \brief stop a timer

#ifdef WIDE_STRING
typedef wchar_t Char;
typedef wstring String;
#define Cout wcout
#include <cwchar>
#define Strcpy wcscpy
//    #define Strncpy wcsncpy
//    #define Strcat wcscat
//    #define Strncat wcsncat
#define Strcmp wcscmp
//    #define Strncmp wcsncmp
#define Strlen wcslen
#else
typedef char Char;
typedef string String;
#define Cout cout
#define Strcpy strcpy
#define Strncpy strncpy
#define Strcat strcat
#define Strncat strncat
#define Strcmp strcmp
#define Strncmp strncmp
#define Strlen strlen
#endif

typedef long long Long;

#ifndef LONGLONG_MAX
#ifdef LINUX
#ifdef LONG_LONG_MAX
const Long LONGLONG_MAX = LONG_LONG_MAX;
#else
const Long LONGLONG_MAX = LLONG_MAX;
#endif
#endif
#ifdef WINDOWS
const Long LONGLONG_MAX = 0x7FFFFFFFFFFFFFFFLL;
#endif
#endif

typedef long double Double;

#ifdef LINUX
inline void mysrand(long seed)
{
    return srand48(seed);
}
inline int myrand() { return lrand48(); }
inline Long myrandl() { return (Long)((Long)lrand48() /**LONGLONG_MAX*/); }
#endif
#ifdef WINDOWS
inline void mysrand(long seed)
{
    return srand(seed);
}
inline int myrand() { return rand(); }
inline Long myrandl() { return (Long)((Long)rand() /**LONGLONG_MAX*/); }
//replace srand48(seed) by srand(seed)
//replace drand48() by (double(rand()) / RAND_MAX)
#endif

#ifdef DOUBLE_PROB
#ifdef LINUX
inline double Pow(double x, double y)
{
    return pow(x, y);
}
inline double Exp10(double x) { return exp10(x); }
inline double Exp(double x) { return exp(x); }
inline double Log10(double x) { return log10(x); }
inline double Expm1(double x) { return expm1(x); }
inline double Log(double x) { return log(x); }
inline double Log1p(double x) { return log1p(x); }
#endif
#ifdef WINDOWS
inline double Pow(double x, double y)
{
    return pow(x, y);
}
inline double Exp10(double x) { return pow(10., x); }
inline double Exp(double x) { return exp(x); }
inline double Log10(double x) { return log(x) / log(10.); }
inline double Log(double x) { return log(x); }
inline double Log1p(double x) { return log(1. + x); }
inline double Expm1(double x) { return exp(x) - 1.; }
#endif
#endif

#ifdef LONGDOUBLE_PROB
#ifdef LINUX
inline Double Pow(Double x, Double y)
{
    return powl(x, y);
}
inline Double Exp10(Double x) { return powl(10.l, (Double)x); }
inline Double Exp(Double x) { return expl(x); }
inline Double Log10(Double x) { return log10l(x); }
inline Double Expm1(Double x) { return expm1l(x); }
inline Double Log(Double x) { return logl(x); }
inline Double Log1p(Double x) { return log1pl(x); }
#endif
#ifdef WINDOWS
inline Double Pow(Double x, Double y)
{
    return pow(x, y);
}
inline Double Exp10(Double x) { return pow(10., x); }
inline Double Log10(Double x) { return log(x) / log(10.); }
inline Double Log(Double x) { return log(x); }
inline Double Log1p(Double x) { return log(1. + x); }
inline Double Expm1(Double x) { return exp(x) - 1.; }
#endif
#endif

#ifdef QUAD_PROB
inline std::ostream& operator<<(std::ostream& os, const __float128& f)
{
    char* y = new char[1000];
    quadmath_snprintf(y, 1000, "%.30Qg", f);
    os << y;
    delete[] y;
    return os;
}

inline std::istream& operator>>(std::istream& is, __float128& f)
{
    char* y = new char[1000];
    is >> y;
    f = strtoflt128(y, NULL);
    delete[] y;
    return is;
}

//inline __float128 abs( __float128 x ){return fabsq( x );}
//inline __float128 sqrt( __float128 x ){return sqrtq( x );}
inline __float128 Pow(__float128 x, __float128 y) { return powq(x, y); }
inline __float128 Exp(__float128 x) { return expq(x); }
inline __float128 Exp10(__float128 x) { return powq(10, x); } // Assumes 10 is representable.
inline __float128 Expm1(__float128 x) { return expm1q(x); } // Assumes 10 is representable.
inline __float128 Log(__float128 x) { return logq(x); }
inline __float128 Log10(__float128 x) { return log10q(x); }
inline __float128 Log1p(__float128 x) { return log1pq(x); }
#endif

#ifdef INT_COST
inline double to_double(const int cost)
{
    return (double)cost;
}
inline int ceil(const int e) { return e; }
inline int floor(const int e) { return e; }
inline int randomCost(int min, int max) { return min + (myrand() % (max - min + 1)); }
inline int string2Cost(const char* ptr) { return atoi(ptr); }

//cost= 0 log2= -1
//cost= 1 log2= 0
//cost= 2 log2= 1
//cost= 3 log2= 1
//cost= 4 log2= 2
//cost= 5 log2= 2
//cost= 6 log2= 2
//cost= 7 log2= 2
//cost= 8 log2= 3
//cost= 9 log2= 3
//cost= 10 log2= 3
//cost= 11 log2= 3
//cost= 12 log2= 3
//cost= 13 log2= 3
//cost= 14 log2= 3
//cost= 15 log2= 3
//cost= 16 log2= 4

// This only works for a 32bits machine
// and compiler g++ version < 4.0

/*
inline int cost2log2(int v)
{
 float x;

 if (v==0) return -1;
 x=(float) v;
 return (*(int*)&x >> 23) - 127;
}
*/

inline int cost2log2(int x)
{
    if (x == 0)
        return -1;
    int l2 = 0;
    x >>= 1;
    for (; x != 0; x >>= 1) {
        ++l2;
    }
    return (l2);
}
inline int cost2log2glb(int x) { return cost2log2(x); }
inline int cost2log2gub(int x) { return cost2log2(x); }
#endif

#ifdef LONGLONG_COST
inline double to_double(const Long cost)
{
    return (double)cost;
}
inline Long ceil(const Long e) { return e; }
inline Long floor(const Long e) { return e; }
inline Long randomCost(Long min, Long max) { return min + (myrandl() % (max - min + 1)); }

#ifdef LINUX
inline Long string2Cost(const char* ptr)
{
    return atoll(ptr);
}
#endif
#ifdef WINDOWS
inline Long string2Cost(const char* ptr)
{
    return atol(ptr);
}
#endif

inline int cost2log2(Long x)
{
    if (x == 0)
        return -1;
    int l2 = 0;
    x >>= 1;
    for (; x != 0; x >>= 1) {
        ++l2;
    }
    return (l2);
}

inline int cost2log2glb(Long x) { return cost2log2(x); }
inline int cost2log2gub(Long x) { return cost2log2(x); }
#endif

// function mkdir
#ifdef LINUX
#include <sys/stat.h>
#include <signal.h>
#endif

#ifdef WIN32
#include <direct.h> // for WINDOWS?
#endif

#ifndef SIZE_MAX
#define SIZE_MAX ((size_t)-1)
#endif

#endif /* TB2SYSTEM_HPP_ */

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
