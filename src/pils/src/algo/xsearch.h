/*
  #ifndef __ils_h
  #define __ils_h
*/
#include <random>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>


#include <algo/localSearch.h>
#include <base/solution.h>
#include <base/costFunction.h>
#include <base/neighborEval.h>
#include <base/incrNeighborEval.h>
#include <algo/perturbation.h>
#include <base/Xover.h>
#include <init/init.h>


using namespace std; 

class Xsearch  : public LocalSearch {
public : 
  using LocalSearch::nEval;
  using LocalSearch::eval;
  using LocalSearch::id;

  Xsearch(CostFunction & _eval,
          IncrNeighborEval & _neighborEval,
          LocalSearch & _ls1,
          Perturbation & _pertu1,
          unsigned long long _nEvalMax,
          unsigned _flatMax) :
              LocalSearch(_eval), nEvalMax(_nEvalMax), flatMax(_flatMax), lastFit(MIN_COST),
              ls1(_ls1), pertu1(_pertu1), neighborEval(_neighborEval)
                                 {

  }

  virtual void operator()(Solution & sol1){
    Xover xo(eval);
    pertu1.init(sol1);
    ls1.nEval = 0;
    ls1(sol1);
#ifndef NDEBUG
    checkfit(sol1);
#endif
    nEval += ls1.nEval;
    unsigned flat = 0;
    while(nEval < nEvalMax && flat <= flatMax && sol1.fitness() > 0) {
      if (ToulBar2::interrupted) {
        throw TimeOut();
      }

      if ( sol1.fitness() == lastFit ){
        flat++;
      } else {
        flat=0;
      }
      lastFit = sol1.fitness();
      sol2 = sol1;
      pertu1.nEval =0;
      pertu1(sol2);
      nEval += pertu1.nEval;

      ls1.nEval = 0;
      ls1(sol2);
#ifndef NDEBUG
      checkfit(sol2);
#endif
      nEval += ls1.nEval;

      if (ToulBar2::verbose >= 1) {
          cout << eval.getLb() + sol1.fitness() << " " << eval.getLb() + sol2.fitness() << " " << nEval << " ";
      }
      xo(sol1, sol2, sol1);
#ifndef NDEBUG
      checkfit(sol1);
#endif

      ls1.nEval = 0;
      ls1(sol1);
#ifndef NDEBUG
      checkfit(sol1);
#endif
      nEval += ls1.nEval;
    }
  }


    	

protected :

  unsigned long long nEvalMax;

  unsigned flatMax;

  Cost lastFit;

  Solution sol2;

  LocalSearch & ls1;

  Perturbation & pertu1;

  IncrNeighborEval & neighborEval;

  multimap<Cost, Solution> fitnesss;

  multimap<Cost, Solution>::iterator it1;

};
