#ifndef TB2SEQUENCEHANDLER_HPP
#define TB2SEQUENCEHANDLER_HPP

#include "core/tb2types.hpp"

using namespace std;

class SequenceHandler {
public:
  SequenceHandler(const string &filename, unsigned nb_backbones);
  ~SequenceHandler();
  void read_enum(const string &filename);
  string get_sequence(unsigned sequence_index) { return sequences_[sequence_index]; };
  bool update_best_cost(unsigned sequence_index);
  void remove_sequence(unsigned sequence_index);
  bool distribute_bkb_and_sequence(int& bkb, unsigned& sequence_index);
  bool distribute_sequence(int backbone, unsigned& sequence_index);
  void update_sequences(int backbone, unsigned sequence_index, Double cost);
  char getAAtype(int sequence, int varIndex) { return sequences_[sequence][varIndex]; }
  void report();

private:
  vector<string> sequences_; // all possible sequences
  vector<Double> costs_; // all positive energies
  vector<vector<unsigned>> bkbsequences; // sequence indexes left for each backbone
  vector<vector<Double>> sequencecosts; // best cost found for each sequence on each backbone
  vector<unsigned> slavebkb; // backbone affected to each slave (! only accessible by master !)
  Double bestcost; // best cost on all backbones
  unsigned bestsequence; // best sequence on all backbones
};

#endif
