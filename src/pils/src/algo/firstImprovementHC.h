#ifndef __firstImprovementHC_h
#define __firstImprovementHC_h

#include <random>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <algo/localSearch.h>
#include <base/solution.h>
#include <base/costFunction.h>
#include <base/neighborEval.h>

class FirstImprovementHC : public LocalSearch {
public:
  using LocalSearch::nEval;
  using LocalSearch::eval;

  FirstImprovementHC(std::mt19937 & _rng, 
                     CostFunction & _eval,
                     std::vector< std::pair<int, int> > & _neighborhood, NeighborEval & _neighborEval, 
                     unsigned long long _nEvalMax, 
                     std::ostream & _out, std::string & _id) : LocalSearch(_eval), rng(_rng), neighborhood(_neighborhood), neighborEval(_neighborEval), nEvalMax(_nEvalMax), out(_out), id(_id) {
    
  }
  
  virtual void operator()(Solution & _solution) {
    unsigned long long nEvalLocal = 0;
    unsigned long long nAccept = 0;
    int r, n;
    Cost neighFit;
    bool accept = true;

    //out << id << nEval << " 0 " << _solution << std::endl;

    while (nEvalLocal < nEvalMax && accept) {
      n = 0;
      accept = false;
      while (n < neighborhood.size() && !accept && nEval < nEvalMax) {
        // random neighbor
        r = runif(rng) % (neighborhood.size() - n);
        
        // test the neighbor
        neighFit =  neighborEval(_solution, neighborhood[r]);        
        nEvalLocal++;
        nEval++;

        // acceptance
        if (neighFit <= _solution.fitness()) {
          accept = true;
          nAccept++;
          _solution.fitness(neighFit);
          _solution[ neighborhood[r].first ] = (_solution[ neighborhood[r].first ] + neighborhood[r].second) % eval.n_values[ neighborhood[r].first ];
          /*
          if (neighFit < _solution.fitness()) {
            out << id << nEval << " " << nAccept << " " << _solution << std::endl;
          }
          */
        } else {
          std::swap(neighborhood[r], neighborhood[neighborhood.size() - 1 - n]);
        }
        
        n++;
      }
      
    }
    
  }
  
protected:
  std::mt19937 & rng;
  std::uniform_int_distribution<int> runif;

  unsigned long long nEvalMax;

  // output
  std::ostream & out;
  
  // external id of the hc
  std::string & id;

  // neighborhood
  std::vector< std::pair<int, int> > neighborhood;
  NeighborEval & neighborEval;
};

#endif
