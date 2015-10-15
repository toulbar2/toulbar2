#ifndef TB2GUMDISTRIB_HPP
#define TB2GUMDISTRIB__HPP

#include "toulbar2lib.hpp"
#include "tb2wcsp.hpp"
#include <iostream>
#include <vector>

class GumDistrib
{
  public:
    GumDistrib();
    ~GumDistrib();
    void Gumbel_noise(vector< vector<Cost> > costs);


};

#endif
