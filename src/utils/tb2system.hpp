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

typedef long long Long;

#ifndef LONGLONG_MAX
#ifdef LONG_LONG_MAX
const Long LONGLONG_MAX = LONG_LONG_MAX;
#else
const Long LONGLONG_MAX = LLONG_MAX;
#endif
#endif

typedef long double Double;

#ifdef __WIN32__
#include <random>
inline void mysrand(long) {};
inline double mydrand() {
    static std::ranlux48 source(std::random_device{}());
    return std::uniform_real_distribution<double>(0,1)(source);
}
inline int myrand() {
    return INT_MAX * mydrand();
}
inline Long myrandl() {
    return LONG_MAX * mydrand();
}
inline Long myrandln() {
    return LONG_MAX * 2. * (mydrand() - 0.5);
}
#else
inline void mysrand(int seed)
{
    return srand48(seed);
}
inline int myrand() { return (int)lrand48(); }
inline Long myrandl() { return (Long)((Long)lrand48() /**LONGLONG_MAX*/); }
inline Long myrandln() { return (Long)((Long)mrand48() /**LONGLONG_MAX*/); }
inline double mydrand() { return drand48(); }
#endif

#ifdef DOUBLE_PROB
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

#ifdef LONGDOUBLE_PROB
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
    if (x <= 0)
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

inline Long string2Cost(const char* ptr)
{
    return atoll(ptr);
}

inline int cost2log2(Long x)
{
    if (x <= 0)
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

//luby(0)= N/A
//luby(1)= 1
//luby(2)= 1
//luby(3)= 2
//luby(4)= 1
//luby(5)= 1
//luby(6)= 2
//luby(7)= 4
//luby(8)= 1
//luby(9)= 1
//luby(10)= 2
//luby(11)= 1
//luby(12)= 1
//luby(13)= 2
//luby(14)= 4
//luby(15)= 8
//luby(16)= 1
inline Long luby(Long r)
{
    int j = cost2log2(r + 1);
    if (r + 1 == (1L << j))
        return (1L << (j - 1));
    else
        return luby(r - (1L << j) + 1);
}

// function mkdir
#include <sys/stat.h>

#ifndef __WIN32__
#include <signal.h>
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
