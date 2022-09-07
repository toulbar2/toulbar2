#ifndef __adaptivePerturb_h
#define __adaptivePerturb_h

#include <random>
#include <iostream>
#include <string>
#include <vector>
//#include <algorithm>
#include <math.h> 

#include <algo/perturbation.h>

#include <algo/localSearch.h>
#include <base/solution.h>
#include <base/costFunction.h>
#include <base/neighborEval.h>


/*
  Perturbation that modifies strength variables at random
    strength is adaptive: increase when remains in the same local optima, decrease when new local optima
*/
class AdaptivePerturb : public Perturbation {
public:
  using LocalSearch::nEval;
  using LocalSearch::eval;

  AdaptivePerturb(std::mt19937 & _rng, 
                CostFunction & _eval,
                NeighborEval & _neighborEval, 
                double _minStrength, double _maxStrength,
                double _incrFactor, double _decrFactor) : Perturbation(_eval), rng(_rng),
                                                   minStrength(_minStrength), maxStrength(_maxStrength),
                                                   incrFactor(_incrFactor), decrFactor(_decrFactor), neighborEval(_neighborEval) {
    ivar.resize(eval.size());    
    for(unsigned k = 0; k < ivar.size(); k++)
        if (eval.n_values[k] > 1)
          ivar[k] = k;
    if (minStrength <= 0)
      minStrength = 2;
    if (maxStrength <= 0)
      maxStrength = eval.size();
  }
  
  virtual void operator()(Solution & _solution) {
    std::pair<int, int> neighbor;
    Cost neighFit;
    unsigned r;
    unsigned n_strength;

    // rounded the strength
    double remains = strength - floor(strength);
    n_strength = int (floor(strength));
    if (remains > 0) {
      if (runifd(rng) > remains)
        n_strength++;
    }

    for(unsigned k = 0; k < n_strength; k++) {
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
    strength = minStrength;
    id = std::to_string(strength) + " ";
  }

  virtual void update(Solution & _previous, Solution & _solution) {
    if (_previous.fitness() == _solution.fitness()) 
      strength *= incrFactor;
    else
      strength *= decrFactor;

    if (strength < minStrength)
      strength = minStrength;
    else
      if (maxStrength < strength)
        strength = maxStrength;

    id = std::to_string(strength) + " ";
  }

protected:
  std::mt19937 & rng;
  std::uniform_int_distribution<int> runif;
  std::uniform_real_distribution<double> runifd;
  
  // neighborhood
  std::vector< unsigned > ivar;

  // perturbation strength
  double strength;
  double minStrength, maxStrength;
  double incrFactor, decrFactor;

  NeighborEval & neighborEval;
};

#endif
