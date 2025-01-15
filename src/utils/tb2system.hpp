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
#include <boost/multiprecision/float128.hpp>
//#include <quadmath.h>
#endif

/*
 * Internal error exceptions
 *
 */

class InternalError : public std::exception {
public:
    InternalError()
    {
    }
    virtual const char* what() const throw() { return "... internal error found, sorry!"; }
};

class BadConfiguration : public InternalError {
public:
    const char* what() const throw() FINAL { return "... bad solver configuration!"; }
};

class WrongFileFormat : public InternalError {
public:
    const char* what() const throw() FINAL { return "... wrong problem file format!"; }
};

extern const char* PrintFormatProb;

/// \brief return real time in seconds since epoch
inline double realTime()
{
    return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count() / 1000.;
}
double cpuTime(); ///< \brief return CPU time in seconds with high resolution (microseconds) if available
void timeOut(int sig); ///< \brief set ToulBar2::interrupted to throw a TimeOut exception later
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

/* --------------------------------------------------------------------
// Random numbers
// -------------------------------------------------------------------- */

#include <random>
extern std::mt19937 myrandom_generator;

inline void mysrand(int seed_)
{
    myrandom_generator.seed(seed_);
}
inline int myrand()
{
    static std::uniform_int_distribution<int> myrandom_uidistribution(0, INT_MAX - 1);
    return myrandom_uidistribution(myrandom_generator);
}
inline Long myrandl()
{
    static std::uniform_int_distribution<Long> myrandom_uldistribution(0, LONG_MAX - 1);
    return myrandom_uldistribution(myrandom_generator);
}
inline Long myrandln()
{
    static std::uniform_int_distribution<Long> myrandom_umdistribution(-LONG_MAX, LONG_MAX - 1);
    return myrandom_umdistribution(myrandom_generator);
}
inline double mydrand()
{
    static std::uniform_real_distribution<double> myrandom_uddistribution(0.0, 1.0);
    return myrandom_uddistribution(myrandom_generator);
}
inline double mydrandl()
{
    static std::uniform_real_distribution<Double> myrandom_uddistribution(0.0, 1.0);
    return myrandom_uddistribution(myrandom_generator);
}
inline double myurand()
{
    static std::uniform_real_distribution<double> myrandom_uddistribution(-1.0, 1.0);
    return myrandom_uddistribution(myrandom_generator);
}
inline double myurandl()
{
    static std::uniform_real_distribution<Double> myrandom_uddistribution(-1.0, 1.0);
    return myrandom_uddistribution(myrandom_generator);
}
inline double mynrand()
{
    static std::normal_distribution<double> myrandom_nddistribution(0.0, 1.0);
    return myrandom_nddistribution(myrandom_generator);
}
inline double mynrandl()
{
    static std::normal_distribution<Double> myrandom_nddistribution(0.0, 1.0);
    return myrandom_nddistribution(myrandom_generator);
}

#ifdef DOUBLE_PROB
inline double Pow(double x, double y)
{
    return pow(x, y);
}
inline double Exp10(double x) { return pow(10., (double)x); }
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
//inline std::ostream& operator<<(std::ostream& os, const __float128& f)
//{
//    char y[1024];
//    quadmath_snprintf(y, 1000, "%.30Qg", f);
//    os << y;
//    return os;
//}

//inline std::istream& operator>>(std::istream& is, __float128& f)
//{
//    char y[1024];
//    is >> y;
//    f = strtoflt128(y, NULL);
//    return is;
//}

inline boost::multiprecision::float128 Pow(boost::multiprecision::float128 x, boost::multiprecision::float128 y) { return pow(x, y); }
inline boost::multiprecision::float128 Exp(boost::multiprecision::float128 x) { return exp(x); }
inline boost::multiprecision::float128 Exp10(boost::multiprecision::float128 x) { return pow(10, x); } // Assumes 10 is representable.
inline boost::multiprecision::float128 Expm1(boost::multiprecision::float128 x) { return expm1(x); }
inline boost::multiprecision::float128 Log(boost::multiprecision::float128 x) { return log(x); }
inline boost::multiprecision::float128 Log10(boost::multiprecision::float128 x) { return log10(x); }
inline boost::multiprecision::float128 Log1p(boost::multiprecision::float128 x) { return log1p(x); }
#endif

#if defined(INT_COST) || defined(SHORT_COST)
inline Double to_double(const int cost)
{
    return (Double)cost;
}
inline int ceil(const int e) { return e; }
inline int floor(const int e) { return e; }
inline int randomCost(int min, int max) { return min + (myrand() % (max - min + 1)); }
inline int string2Cost(const char* ptr) { return atoi(ptr); }

// cost= 0 log2= -1
// cost= 1 log2= 0
// cost= 2 log2= 1
// cost= 3 log2= 1
// cost= 4 log2= 2
// cost= 5 log2= 2
// cost= 6 log2= 2
// cost= 7 log2= 2
// cost= 8 log2= 3
// cost= 9 log2= 3
// cost= 10 log2= 3
// cost= 11 log2= 3
// cost= 12 log2= 3
// cost= 13 log2= 3
// cost= 14 log2= 3
// cost= 15 log2= 3
// cost= 16 log2= 4

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
inline Double to_double(const Long cost)
{
    return (Double)cost;
}
inline Long ceil(const Long e) { return e; }
inline Long floor(const Long e) { return e; }
inline Long randomCost(Long min, Long max) { return min + (myrandl() % (max - min + 1)); }

inline Long string2Cost(const char* ptr)
{
    try {
        std::string s(ptr);
        if (s.size() == 0) {
            return 0;
        }
        Double d = stold(s);
        if ((d >= 0 && d > LONGLONG_MAX) || (d < 0 && d < -LONGLONG_MAX)) {
            throw std::out_of_range("long long overflow!");
        }
    } catch (std::exception& e) {
        std::cerr << "Overflow exception: cannot convert this number " << ptr << " into a cost!" << std::endl;
        std::cerr << "\t" << e.what() << std::endl;
        throw WrongFileFormat();
    }
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

// luby(0)= N/A
// luby(1)= 1
// luby(2)= 1
// luby(3)= 2
// luby(4)= 1
// luby(5)= 1
// luby(6)= 2
// luby(7)= 4
// luby(8)= 1
// luby(9)= 1
// luby(10)= 2
// luby(11)= 1
// luby(12)= 1
// luby(13)= 2
// luby(14)= 4
// luby(15)= 8
// luby(16)= 1
inline Long luby(Long r)
{
    int j = cost2log2(r + 1);
    if (r + 1 == (1L << j))
        return (1L << (j - 1));
    else
        return luby(r - (1L << j) + 1);
}

// my aleaGauss generator
#define M_PIl 3.141592653589793238462643383279502884L /* pi */
inline double aleaGaussNoise(double s)
{
    double U1 = mydrand();
    double U2 = mydrand();
    // double U2 = myrandom_uddistribution(0.0, 1.0);
    double P = s * sqrt(-2 * log(U1)) * cos(2 * M_PIl * U2);
    return P;
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
