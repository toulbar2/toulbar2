#ifndef TB2SEQ_HPP
#define TB2SEQ_HPP

#include "tb2wcsp.hpp"
#include <iostream>
#include <vector>

class Seq {
public:
    Seq();
    ~Seq();
    Seq(std::string fname);
    std::string get_sequence() { return sequence; }
    std::vector<std::vector<bool>> get_mask() { return mask; }
    void generate_mask(std::vector<std::vector<char>> rots);
    void mask_variable(vector<vector<Cost>> m_cost);
    void show_mask();

private:
    std::string sequence;
    std::vector<std::vector<bool>> mask;
};

#endif
