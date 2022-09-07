#ifndef __staticPerturb_h
#define __staticPerturb_h

#include <random>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <algo/perturbation.h>

#include <algo/localSearch.h>
#include <base/solution.h>
#include <base/costFunction.h>
#include <base/neighborEval.h>


/*
  Perturbation that modifies _strength variables at random
    it uses incremental evaluation neighborEval at each modification 
*/
class StaticPerturb : public Perturbation {
public:
  using LocalSearch::nEval;
  using LocalSearch::eval;

  StaticPerturb(std::mt19937 & _rng, 
                     CostFunction & _eval,
                     NeighborEval & _neighborEval, 
                     unsigned _strength) : Perturbation(_eval), rng(_rng), strength(_strength), neighborEval(_neighborEval) {
    ivar.resize(eval.size());
    for(unsigned k = 0; k < ivar.size(); k++)
        if (eval.n_values[k] > 1)
          ivar[k] = k;
  }
  
  virtual void operator()(Solution & _solution) {
    std::pair<int, int> neighbor;
    Cost neighFit;
    unsigned r;

    for(unsigned k = 0; k < strength; k++) {
      // rnd neighbor
      r = runif(rng) % (ivar.size() - k);
      neighbor.first = ivar[r];
      neighbor.second = (eval.n_values[ neighbor.first ]>1)?1 + runif(rng) % (eval.n_values[ neighbor.first ]-1):1;
      // fitness of the neighbor
      neighFit = neighborEval(_solution, neighbor);
      nEval++;
      // move
      _solution.fitness(neighFit);
      _solution[ neighbor.first ] = (_solution[ neighbor.first ] + neighbor.second) % eval.n_values[ neighbor.first ];
      // "remove" the variable
      std::swap(ivar[r], ivar[ivar.size() - 1 - k]);
    }

  }
  
  virtual void init(Solution & _solution) {
    id = std::to_string(strength) + " ";
  }

protected:
  std::mt19937 & rng;
  std::uniform_int_distribution<int> runif;
  
  // neighborhood
  std::vector< unsigned > ivar;

  // perturbation strength
  unsigned strength;

  NeighborEval & neighborEval;
};

#endif
