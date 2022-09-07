#ifndef __perturbation_h
#define __perturbation_h

#include <iostream>
#include <string>

#include <algo/localSearch.h>
#include <base/solution.h>
#include <base/costFunction.h>

/*
  Perturbation that modifies solution at random
*/
class Perturbation : public LocalSearch {
public:
    Perturbation(CostFunction & _eval) : LocalSearch(_eval) {
    }
  
    virtual void init(Solution & _solution) {

    }
  
    virtual void update(Solution & _previous, Solution & _solution) {
    
    }

};

#endif
