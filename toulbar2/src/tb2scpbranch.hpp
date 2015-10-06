#ifndef TB2SCPBRANCH_HPP
#define TB2SCPBRANCH_HPP

#include "tb2types.hpp"
#include "tb2btlist.hpp"

class Tb2ScpBranch {

public:

  // typedef struct SortItem
  // {
  //   double h;
  //   size_t i;
  // }SortItem;

  Tb2ScpBranch(){}
  //  Tb2ScpBranch(const char *filename);
  ~Tb2ScpBranch();
  //  double AA2Criterium(char c);
  //  tuple<size_t,size_t> getBounds(size_t var_index);
  tuple<size_t, size_t> getBounds(int varIndex, Value value);
  size_t moveAAFirst(ValueCost * sorted, size_t domsize, size_t left, size_t right);
  // void keep(size_t begin, size_t end, size_t var_index);
  // void remove(size_t begin, size_t end, size_t var_index);
  // void remove(size_t index, size_t var_index);

  //private:
  // vector< vector <SortItem> > sort_criterium;
  // vector< double > ub_reference;
  // vector< Value > initial_sol;
};

#endif
