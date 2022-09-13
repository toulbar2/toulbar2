/*
 init.h

 Author: 
  Sebastien Verel, 
  Univ. du Littoral Côte d'Opale, France.
 
  David Simoncini
  Univ. Toulouse 1 Capitole, France.

*/
#ifndef _init_h
#define _init_h

#include <random>
#include <base/solution.h>
#include <base/costFunction.h>

class Init {
public:
  Init(std::mt19937 & _rng, CostFunction & _eval) : rng(_rng), protein_size(_eval.size()) {

      runif.resize(_eval.size());

      for(unsigned i = 0; i < protein_size; i++)
        runif[i].param(std::uniform_int_distribution<int>::param_type(0, _eval.variable_size(i) - 1));
  }
  
  void operator()(Solution & _solution) {
    if (_solution.size() != protein_size) {
      _solution.resize(protein_size);
    }

    for(size_t i = 0; i < _solution.size(); i++)
      _solution[i] = runif[i](rng);
    
    _solution.invalidate();
  }

protected:
  std::mt19937 & rng;

  std::vector< std::uniform_int_distribution<int> > runif;
  
  unsigned int protein_size;
};

#endif
