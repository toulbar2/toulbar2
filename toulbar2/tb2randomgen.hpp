#ifndef TB2RANDOMGEN_H_
#define TB2RANDOMGEN_H_

#include "tb2wcsp.hpp"
#include <vector>
using namespace std;


class naryRandom {
public:

  WCSP& wcsp;

  naryRandom(WCSP* wcspin) : wcsp(*wcspin) { srand48(3); }
  ~naryRandom() {}

  long myrand() {  return lrand48();   }
     
  int n,m;
     
  bool connected();
  void generateTernCtr( int i, int j, int k, long p, Cost costMin = 1, Cost costMax = 1 );
  void generateBinCtr( int i, int j, long p, Cost costMin = 1, Cost costMax = 1 );
  void Input( int in_n, int in_m, vector<int>& p, int seed = 0 );  

  
  bool ini( vector<int>& index, int arity );
  long toIndex( vector<int>& index );
  int inc( vector<int>& index, int i );
  bool inc( vector<int>& index );

};






#endif /*TB2RANDOMGEN_H_*/
