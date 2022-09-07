#ifndef __ils_h
#define __ils_h

#include <random>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <algo/localSearch.h>
#include <base/solution.h>
#include <base/costFunction.h>
#include <base/neighborEval.h>
#include <base/incrNeighborEval.h>
#include <algo/perturbation.h>

using namespace std;

class ILS : public LocalSearch {
public:
    using LocalSearch::nEval;
    using LocalSearch::eval;
    using LocalSearch::id;

    ILS(CostFunction & _eval,
        IncrNeighborEval & _neighborEval,
        LocalSearch & _ls,
        Perturbation & _perturbation,
        unsigned long long _nEvalMax, 
        unsigned _flatMax,
        std::ostream & _out) : LocalSearch(_eval), nEvalMax(_nEvalMax), flatMax(_flatMax),
        ls(_ls), perturbation(_perturbation), neighborEval(_neighborEval), out(_out) {

    }

    virtual void operator()(Solution & _solution) {
        unsigned i = 0;
        id = std::to_string(i) + " " ;
        perturbation.init(_solution);
        // first ls
        out << id << perturbation.id << nEval << " " << _solution << std::endl;
        ls.id = std::to_string(i) + " ";
        ls.nEval = 0;
        ls(_solution);
        nEval += ls.nEval;

        //checkfit(_solution);

        out << id << perturbation.id << nEval << " " << _solution << std::endl;
        out << id << perturbation.id << nEval << " " << _solution << std::endl;

        Solution previous;

        i++;
        unsigned flat = 0;

        while (nEval < nEvalMax && flat < flatMax) {

            if ( _solution.fitness() == lastFit ){
                flat++;
            } else {
                flat=0;
            }
            lastFit = _solution.fitness();

            //std::cout << i << std::endl;
            id = std::to_string(i) + " " ;

            // copy the solution before perturbation
            previous = _solution;

            // perturbation
            perturbation.nEval = 0;
            perturbation(_solution);
            nEval += perturbation.nEval;
            //checkfit(_solution);
            out << id << perturbation.id << nEval << " " << _solution << std::endl;


            // local search
            ls.nEval = 0;
            ls.id = std::to_string(i) + " ";
            ls(_solution);
            nEval += ls.nEval;

            //            checkfit(_solution);
            out << id << perturbation.id << nEval << " " << _solution << std::endl;

            // update perturbation
            perturbation.update(previous, _solution);
            // selection
            if (previous.fitness() < _solution.fitness()) {
                // then previous solution is better, so restart from that solution
                _solution = previous;
            }

            //checkfit(_solution);
            out << id << perturbation.id << nEval << " " << _solution << std::endl;

            i++;
        }

    }

protected:
  unsigned long long nEvalMax;

  unsigned flatMax;

  Cost lastFit = 0;
  
  // local search  (here HC)
  LocalSearch & ls;

  // perturbation
  Perturbation & perturbation;

  // neighborhood
  IncrNeighborEval & neighborEval;

  // output
  std::ostream & out;
};

#endif
