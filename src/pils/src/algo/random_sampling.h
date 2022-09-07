/*
 random_sampling.h

 Author: 
  Sebastien Verel, 
  Univ. du Littoral Côte d'Opale, France.
 
  David Simoncini
  Univ. Toulouse 1 Capitole, France.

*/
#ifndef __randomSampling__h
#define __randomSampling__h

#include <vector>

#include <base/solution.h>
#include <init/init.h>


class RandomSampling {
public:
  RandomSampling(CostFunction & _eval, Init & _init, unsigned _n) : eval(_eval), init(_init), n(_n) {
  }
  
  // run the sampling
  void operator()() {
    Solution x(eval.size());

    for(unsigned i = 0; i < n; i++) {
      init(x);
      eval(x);

      //x.print(); std::cout << std::endl;

      sample.push_back( x.fitness() );
    }
  }

  // fitness values of the sample
  std::vector<Cost> sample;

protected:
  // Automata
  CostFunction & eval;

  // initiatization
  Init & init;

  // size of the sample
  unsigned n;
};

#endif
