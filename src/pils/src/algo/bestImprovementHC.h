#ifndef __bestImprovementHC_h
#define __bestImprovementHC_h

#include <random>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <algo/localSearch.h>
#include <base/solution.h>
#include <base/costFunction.h>
#include <base/neighborEval.h>

class BestImprovementHC : public LocalSearch {
public:
  using LocalSearch::nEval;
  using LocalSearch::eval;
  using LocalSearch::id;

  BestImprovementHC(std::mt19937 & _rng, 
                     CostFunction & _eval,
                     std::vector< std::pair<int, int> > & _neighborhood, NeighborEval & _neighborEval, 
                     unsigned long long _nEvalMax, 
                     std::ostream & _out) : LocalSearch(_eval), rng(_rng), neighborhood(_neighborhood), neighborEval(_neighborEval), nEvalMax(_nEvalMax), out(_out) {
    
  }
  
  virtual void operator()(Solution & _solution) {
    unsigned long long nEvalLocal = 0;
    unsigned long long nAccept = 0;
    int r, n;
    Cost neighFit, bestFit;
    bool accept = true;
    unsigned nBest;
    std::vector<unsigned> iBest(neighborhood.size());

    //out << id << nEval << " 0 " << _solution << std::endl;

    while (nEvalLocal < nEvalMax && accept) {
      out << id << nEvalLocal << " " << _solution << std::endl;

      // find best neighbor
      n = 0;
      bestFit =  neighborEval(_solution, neighborhood[n]);        
      nBest = 1;
      iBest[0] = 0;
      n++;

      while (n < neighborhood.size() && nEvalLocal < nEvalMax) {
        // test the neighbor
        neighFit =  neighborEval(_solution, neighborhood[n]);

        // update best solution
        if (neighFit < bestFit) {
          bestFit = neighFit;
          nBest = 1;
          iBest[0] = n;
        } else 
          if (neighFit == bestFit) { // when ties
            nBest++;
            iBest[nBest - 1] = n;
          }

        n++;
      }

      if (bestFit <= _solution.fitness()) {
          accept = true;

          if (nBest > 1) {
            r = iBest[runif(rng) % nBest];
          } else {
            r = iBest[0];
          }
          // move
          _solution.fitness(bestFit);
          _solution[ neighborhood[r].first ] = (_solution[ neighborhood[r].first ] + neighborhood[r].second) % eval.n_values[ neighborhood[r].first ];
//            out << id << nEval << " " << nAccept << " " << _solution << std::endl;        
      } else 
        accept = false;

      nEvalLocal++;
    }
    
    nEval += nEvalLocal;
    out << id << nEvalLocal << " " << _solution << std::endl;
  }
  
protected:
  std::mt19937 & rng;
  std::uniform_int_distribution<int> runif;

  unsigned long long nEvalMax;

  // output
  std::ostream & out;
  
  // neighborhood
  std::vector< std::pair<int, int> > neighborhood;
  NeighborEval & neighborEval;
};

#endif
