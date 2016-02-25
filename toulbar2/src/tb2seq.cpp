#include <iostream>
#include <vector>
#include "tb2seq.hpp"

using namespace std;

Seq::Seq(string fname)
{
    ifstream seqdata(fname);

    if (seqdata) {
        seqdata >> sequence;
    } else {
        cout << "Could not open file: " << fname << endl;
        exit(1);
    }
}

void Seq::generate_mask(vector < vector < char > > rots)
{

    vector < bool > varmask;
    if (rots.size() == sequence.size()) {
        for (size_t i = 0; i < sequence.size(); i++) {
            for (size_t j = 0; j < rots[i].size(); j++) {
                if (rots[i][j] == sequence[i])
                    varmask.push_back(true);
                else
                    varmask.push_back(false);
            }
            mask.push_back(varmask);
            varmask.clear();
        }
    } else {
        cout << "Sequence size does not match the number of variables" << endl;
        exit(1);
    }
}

void Seq::mask_variable(vector< vector<Cost> > m_cost)
{

    for (size_t i = 0; i < m_cost.size(); i++) {
        for (size_t j = 0; j < m_cost[i].size(); j++) {
            if (!mask[i][j]) cout << "Switching : " << m_cost[i][j] << " to TOP" << endl;
        }
    }
}
