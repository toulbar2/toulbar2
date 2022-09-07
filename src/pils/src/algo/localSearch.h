#ifndef _localSearch_h
#define _localSearch_h

#include <string.h>

#include <base/solution.h>
#include <base/costFunction.h>

class LocalSearch {
public:
  LocalSearch(CostFunction & _eval) : eval(_eval) {
    
  }
  virtual ~LocalSearch() {}
  
  virtual void operator()(Solution & _solution) = 0;

  // number of evaluation
  unsigned long long nEval;

  // id of the search (for output)
  std::string id;
  
protected:
  CostFunction & eval;
  
};

#endif
