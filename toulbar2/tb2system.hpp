/** \file tb2system.hpp
 *  \brief System dependent functions.
 * 
 */

#ifndef TB2SYSTEM_HPP_
#define TB2SYSTEM_HPP_

typedef long long Long;
#ifndef LONGLONG_MAX
const Long LONGLONG_MAX = LONG_LONG_MAX;
#endif

typedef long double Double;

inline void mysrand(long seed) {return srand48(seed);}
inline int myrand() { return lrand48();  }
inline Long myrandl() {return (Long) ((Long)lrand48()/**LONGLONG_MAX*/);}

#ifdef DOUBLE_PROB
inline double Pow(double x, double y) {return pow(x,y);}
inline double Log10(double x) {return log10(x);}
inline double Log(double x) {return log(x);}
inline double Log1p(double x) {return log1p(x);}
#endif

#ifdef LONGDOUBLE_PROB
inline Double Pow(Double x, Double y) {return powl(x,y);}
inline Double Log10(Double x) {return log10l(x);}
inline Double Log(Double x) {return logl(x);}
inline Double Log1p(Double x) {return log1pl(x);}
#endif

#ifdef INT_COST
inline double to_double(const int cost) {return (double) cost;}
inline int ceil(const int e) {return e;}
inline int floor(const int e) {return e;}
inline int randomCost(int min, int max) { return min + (myrand() % (max - min + 1)); }
inline int string2Cost(char *ptr) {return atoi(ptr);}

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
        if (x==0) return -1;
        register int l2 = 0;
        x>>=1;
        for (; x != 0; x >>=1)
        {
                ++ l2;
        }
        return (l2);
}
inline int cost2log2glb(int x) {return cost2log2(x);}
inline int cost2log2gub(int x) {return cost2log2(x);}
#endif

#ifdef LONGLONG_COST
inline double to_double(const Long cost) {return (double) cost;}
inline Long ceil(const Long e) {return e;}
inline Long floor(const Long e) {return e;}
inline Long randomCost(Long min, Long max) { return  min + (myrandl() % (max - min + 1)); }
inline Long string2Cost(char *ptr) {return atoll(ptr);}

inline int cost2log2(Long x)
{
        if (x==0) return -1;
        register int l2 = 0;
        x>>=1;
        for (; x != 0; x >>=1)
        {
                ++ l2;
        }
        return (l2);
}
inline int cost2log2glb(Long x) {return cost2log2(x);}
inline int cost2log2gub(Long x) {return cost2log2(x);}
#endif

#endif /* TB2SYSTEM_HPP_ */
