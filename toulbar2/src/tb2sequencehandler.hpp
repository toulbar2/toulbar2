#ifndef TB2SEQUENCEHANDLER_HPP
#define TB2SEQUENCEHANDLER_HPP

#include <vector>
#include <string>
#include "tb2types.hpp"

using namespace std;

class SequenceHandler {
public:
    SequenceHandler(string filename, unsigned nb_backbones);
    ~SequenceHandler();
    string get_sequence(unsigned sequence_index) { return sequences_[sequence_index]; };
    bool update_best_cost(unsigned sequence_index);
    void remove_sequence(unsigned sequence_index);
    bool distribute_bkb_and_sequence(int& bkb, unsigned& sequence_index);
    bool distribute_sequence(int backbone, unsigned& sequence_index);
    void update_sequences(int backbone, unsigned sequence_index, Cost cost);
    char getAAtype(int sequence, int varIndex) { return sequences_[sequence][varIndex]; }
    void report();

private:
    vector<string> sequences_; // all possible sequences
    vector<vector<unsigned>> bkbsequences; // sequence indexes left for each backbone
    vector<vector<Cost>> sequencecosts; // best cost found for each sequence on each backbone
    vector<unsigned> slavebkb; // backbone affected to each slave (! only accessible by master !)
    Cost bestcost; // best cost on all backbones
    unsigned bestsequence; // best sequence on all backbones
};

#endif
