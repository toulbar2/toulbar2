/*
 * **************** System dependent functions **********************
 */

#include "tb2system.hpp"

int cost2log2(int v)
{ 
  float x;

  if (v==0) return -1;
  x=(float) v; 
  return (*(int*)&x >> 23) - 127;
}
