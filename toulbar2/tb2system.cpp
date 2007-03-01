/*
 * **************** System dependent functions **********************
 */

#include "tb2system.hpp"


#ifdef DOUBLE_PROB
double Pow(double x, double y) {return pow(x,y);}
double Log10(double x) {return log10(x);}
double Log(double x) {return log(x);}
double Log1p(double x) {return log1p(x);}
Cost String2Cost(char *ptr) {return atoi(ptr);}
#endif

#ifdef LONGDOUBLE_PROB
Double Pow(Double x, Double y) {return powl(x,y);}
Double Log10(Double x) {return log10l(x);}
Double Log(Double x) {return logl(x);}
Double Log1p(Double x) {return log1pl(x);}
Cost String2Cost(char *ptr) {return atoll(ptr);}
#endif
