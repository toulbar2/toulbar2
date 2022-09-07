/*
 incrNeighborEval.h

 Author: 
  Sebastien Verel, 
  Univ. du Littoral Côte d'Opale, France.
 
  David Simoncini
  Univ. Toulouse 1 Capitole, France.

*/

#ifndef _incrNeighborEval_h
#define _incrNeighborEval_h

#include <iostream>
#include <vector>
#include <base/solution.h>
#include <base/costFunction.h>
#include <base/neighborEval.h>

class IncrNeighborEval : public NeighborEval {
public:
  using NeighborEval::eval;

  IncrNeighborEval(CostFunction & _eval) : NeighborEval(_eval) {
    unsigned j ;
    interactions.resize(eval.links.size());
    for(unsigned i = 0; i < eval.links.size(); i++)
      for(unsigned k = 0; k < eval.links[i].size(); k++) {
        j = eval.links[i][k];
        interactions[j].push_back(i);
      }
  }

  virtual Cost operator()(Solution & _solution, std::pair<int, int> & _neighbor) {
    Cost delta ;

    unsigned new_value = (_solution[_neighbor.first] + _neighbor.second) % eval.n_values[_neighbor.first];
    // linear part
    delta = - eval.energy[ _neighbor.first ][ _solution[_neighbor.first] ] + eval.energy[ _neighbor.first ][ new_value ]; 
    // quadratic part
    unsigned i, j;
    for(unsigned k = 0; k < eval.links[ _neighbor.first ].size(); k++) {
      j = eval.links[ _neighbor.first ][k];
      delta += - eval.energy2[ _neighbor.first ][ j ][ _solution[_neighbor.first] ][ _solution[j] ] 
               + eval.energy2[ _neighbor.first ][ j ][ new_value ][ _solution[j] ];
    }
    for(unsigned k = 0; k < interactions[ _neighbor.first ].size(); k++) {
      i = interactions[ _neighbor.first ][k];
      delta += - eval.energy2[ i ][ _neighbor.first ][ _solution[i] ][ _solution[_neighbor.first] ] 
               + eval.energy2[ i ][ _neighbor.first ][ _solution[i] ][ new_value ];
    }
    return _solution.fitness() + delta;
  }

  // interaction between variables with j < i
  std::vector< std::vector<unsigned> > interactions;

};

#endif
