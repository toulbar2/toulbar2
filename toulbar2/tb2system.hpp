/** \file tb2system.hpp
 *  \brief System dependent functions.
 * 
 */

#ifndef TB2SYSTEM_HPP_
#define TB2SYSTEM_HPP_

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


#endif /* TB2SYSTEM_HPP_ */
