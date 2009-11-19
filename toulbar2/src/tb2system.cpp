/*
 * ****** System dependent functions.
 */
 
#include "tb2types.hpp"
#include "tb2system.hpp"

/* --------------------------------------------------------------------
// Timer management functions
// -------------------------------------------------------------------- */
#ifdef LINUX 
#include <unistd.h> 
#include <sys/time.h>
#include <sys/times.h>

double cpuTime()
{
  static struct tms buf;

  times(&buf);
  double res = ((double) (buf.tms_utime+buf.tms_stime+buf.tms_cutime+buf.tms_cstime)) / ((double) sysconf(_SC_CLK_TCK));
  return (res>0)?res:0;
}

#else
double cpuTime()
{
  return (double) (clock() / CLOCKS_PER_SEC);
}
#endif
