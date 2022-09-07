#ifndef __rndPerturb_h
#define __rndPerturb_h

#include <random>
#include <iostream>
#include <string>
#include <vector>
#include <math.h> 

#include <algo/perturbation.h>

#include <algo/localSearch.h>
#include <base/solution.h>
#include <base/costFunction.h>
#include <base/neighborEval.h>


/*
  Perturbation that modifies strength variables at random
    strength is random: between min and max values at each iteration
*/
class RandomPerturb : public Perturbation {
public:
  using LocalSearch::nEval;
  using LocalSearch::eval;

  RandomPerturb(std::mt19937 & _rng, 
                CostFunction & _eval,
                NeighborEval & _neighborEval, 
                unsigned _minStrength, unsigned _maxStrength) : Perturbation(_eval), rng(_rng),
                                                   minStrength(_minStrength), maxStrength(_maxStrength), neighborEval(_neighborEval) {
    ivar.resize(eval.size());
    for(unsigned k = 0; k < ivar.size(); k++)
        if (eval.n_values[k] > 1)
          ivar[k] = k;

    runifStrength.param( std::uniform_int_distribution<int>::param_type(minStrength, maxStrength) );
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
    strength = runifStrength(rng);
    id = std::to_string(strength) + " ";
  }

  virtual void update(Solution & _previous, Solution & _solution) {
    strength = runifStrength(rng);

    id = std::to_string(strength) + " ";
  }

protected:
  std::mt19937 & rng;
  std::uniform_int_distribution<int> runif;
  std::uniform_int_distribution<int> runifStrength;
  
  // neighborhood
  std::vector< unsigned > ivar;

  // perturbation strength
  unsigned strength;
  unsigned minStrength, maxStrength;

  NeighborEval & neighborEval;
};

#endif
