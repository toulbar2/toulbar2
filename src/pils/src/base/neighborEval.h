/*
 neighborEval.h

 Author: 
  Sebastien Verel, 
  Univ. du Littoral Côte d'Opale, France.
 
  David Simoncini
  Univ. Toulouse 1 Capitole, France.

*/

#ifndef _neighborEval_h
#define _neighborEval_h

#include <iostream>
#include <vector>
#include <base/solution.h>
#include <base/costFunction.h>

class NeighborEval {
public:
  NeighborEval(CostFunction & _eval) : eval(_eval) {
  }

  virtual Cost operator()(Solution & _solution, std::pair<int, int> & _neighbor) = 0;

protected:
  // Evaluation function
  CostFunction & eval;

};

#endif
