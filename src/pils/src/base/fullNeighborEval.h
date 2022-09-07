/*
 fullNeighborEval.h

 Author: 
  Sebastien Verel, 
  Univ. du Littoral Côte d'Opale, France.
 
  David Simoncini
  Univ. Toulouse 1 Capitole, France.

*/

#ifndef _fullNeighborEval_h
#define _fullNeighborEval_h

#include <iostream>
#include <vector>
#include <base/solution.h>
#include <base/costFunction.h>
#include <base/neighborEval.h>

class FullNeighborEval : public NeighborEval {
public:
  FullNeighborEval(CostFunction & _eval) : NeighborEval(_eval) {
  }

  virtual Cost operator()(Solution & _solution, std::pair<int, int> & _neighbor) {
    Cost tmp = _solution.fitness();
    unsigned v = _solution[ _neighbor.first ];
    _solution[ _neighbor.first ] = (v + _neighbor.second) % eval.n_values[ _neighbor.first ];
    eval(_solution);
  	Cost neighFit = _solution.fitness();
	  // back
     _solution[ _neighbor.first ] = v;
	  _solution.fitness(tmp);

  	return neighFit;
  }

};

#endif
