/*
 * **************** System dependent functions **********************
 */

#include "tb2system.hpp"


#ifdef INT_COST
double Pow(double x, double y) {return pow(x,y);}
double Log10(double x) {return log10(x);}
double Log(double x) {return log(x);}
double Log1p(double x) {return log1p(x);}
#endif

#ifdef LONGLONG_COST
long double Pow(long double x, long double y) {return powl(x,y);}
long double Log10(long double x) {return log10l(x);}
long double Log(long double x) {return logl(x);}
long double Log1p(long double x) {return log1pl(x);}
#endif
