#include "tb2sequencehandler.hpp"
#include <fstream>
#include <climits>
#include <time.h>
#include <stdlib.h>

using namespace std;

SequenceHandler::SequenceHandler(string filename, unsigned nb_backbones)
{
    ifstream is(filename);
    string current_seq;
    while (is) {
        is >> current_seq;
        if (is) {
            sequences_.push_back(current_seq);
        }
    }
    sequencecosts.resize(sequences_.size());
    bkbsequences.resize(nb_backbones);
    for (unsigned i = 0; i < nb_backbones; i++) {
        for (unsigned j = 0; j < sequences_.size(); j++)
            bkbsequences[i].push_back(j);
    }
    srand(time(NULL));
    bestcost = LLONG_MIN;
    bestsequence = UINT_MAX;
}

SequenceHandler::~SequenceHandler()
{
}

// checks if we have a new bestcost. Ususally we do when we call this method.
bool SequenceHandler::update_best_cost(unsigned sequence_index)
{
    Cost mincost = LLONG_MAX;
    for (unsigned i = 0; i < sequencecosts[sequence_index].size(); i++)
        if (sequencecosts[sequence_index][i] < mincost)
            mincost = sequencecosts[sequence_index][i];
    if (mincost > bestcost) {
        bestcost = mincost;
        bestsequence = sequence_index;
        return true;
    }
    return false;
}

// Removes a sequence from all possible backbones
void SequenceHandler::remove_sequence(unsigned sequence_index)
{
    for (unsigned i = 0; i < bkbsequences.size(); i++) {
        int to_remove = -1;
        for (unsigned j = 0; j < bkbsequences[i].size(); j++)
            if (bkbsequences[i][j] == sequence_index) {
                to_remove = (int)j;
                break;
            }
        if (to_remove != -1) {
            bkbsequences[i][(unsigned)to_remove] = bkbsequences[i][bkbsequences[i].size() - 1];
            bkbsequences[i].pop_back();
        }
    }
}

bool SequenceHandler::distribute_bkb_and_sequence(int& bkb, unsigned& sequence_index)
{
    for (unsigned backbone = 0; backbone < bkbsequences.size(); backbone++) {
        if (bkbsequences[backbone].size() == 0)
            continue;
        unsigned bkb_index = (unsigned)rand() % bkbsequences[backbone].size();
        sequence_index = bkbsequences[backbone][bkb_index];
        bkb = backbone + 1;
        bkbsequences[backbone][bkb_index] = bkbsequences[backbone][bkbsequences[backbone].size() - 1];
        bkbsequences[backbone].pop_back();
        return true;
    }
    return false;
}

// Randomly pick a sequence for a given backbone
bool SequenceHandler::distribute_sequence(int backbone, unsigned& sequence_index)
{
    if (bkbsequences[backbone].size() == 0)
        return false;
    unsigned bkb_index = (unsigned)rand() % bkbsequences[backbone].size();
    sequence_index = bkbsequences[backbone][bkb_index];
    bkbsequences[backbone][bkb_index] = bkbsequences[backbone][bkbsequences[backbone].size() - 1];
    bkbsequences[backbone].pop_back();
    return true;
}

// store a new cost for a sequence and check if anything can be updated
void SequenceHandler::update_sequences(int backbone, unsigned sequence_index, Cost cost)
{
    if (cost < bestcost)
        remove_sequence(sequence_index);
    else {
        sequencecosts[sequence_index].push_back(cost);
        if (sequencecosts[sequence_index].size() == bkbsequences.size()) // if this sequence has a cost for every backbone
            update_best_cost(sequence_index); // we can update bestcost to this sequence's min cost
    }
}

void SequenceHandler::report()
{
    cout << "Best candidate sequence for negative design : " << sequences_[bestsequence] << endl;
    cout << "Cost : " << bestcost << endl;
}
