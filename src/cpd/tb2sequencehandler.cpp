#include "tb2sequencehandler.hpp"
#include <fstream>
#include <climits>
#include <float.h>
#include <time.h>
#include <stdlib.h>

using namespace std;

SequenceHandler::SequenceHandler(const string &filename, unsigned nb_backbones)
{
  read_enum(filename);
  sequencecosts.resize(sequences_.size());
  bkbsequences.resize(nb_backbones);
  for (unsigned i = 0; i < nb_backbones; i++) {
    for (unsigned j = 0; j < sequences_.size(); j++)
      bkbsequences[i].push_back(j);
  }
  srand(time(NULL));
  bestcost = -LDBL_MAX;
  bestsequence = UINT_MAX;
}

SequenceHandler::~SequenceHandler()
{
}

// Reads pompd enum files. Assumes that enum file format does not change.
void SequenceHandler::read_enum(const string &filename)
{
  ifstream is(filename);
  string line;
  string word;
  while(getline(is,line))
    {
      if (line.find("sequence")!=string::npos)
        {
          istringstream iss(line);
          vector<string> words;
          while(iss >> word)
            {
              words.push_back(word);
            }
          sequences_.push_back(words[2]);
          costs_.push_back(stold(words[6]));
        }
    }
}

// checks if we have a new bestcost. Ususally we do when we call this method.
bool SequenceHandler::update_best_cost(unsigned sequence_index)
{
  Double mincost = LDBL_MAX;
  for (unsigned i = 0; i < sequencecosts[sequence_index].size(); i++)
    if (sequencecosts[sequence_index][i] < mincost)
      mincost = sequencecosts[sequence_index][i];
  if (ToulBar2::diffneg) {
    if ((mincost - costs_[sequence_index]) > (bestcost - costs_[sequence_index])) {
      bestcost = mincost;
      bestsequence = sequence_index;
      return true;
    }
  }
  else
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
void SequenceHandler::update_sequences(int backbone, unsigned sequence_index, Double cost)
{
  // if (ToulBar2::diffneg)
  //     cost = cost - costs_[sequence_index];
  sequencecosts[sequence_index].push_back(cost);
  if (cost < bestcost)
    remove_sequence(sequence_index);
  else {
    if (sequencecosts[sequence_index].size() == bkbsequences.size()) // if this sequence has a cost for every backbone
      update_best_cost(sequence_index); // we can update bestcost to this sequence's min cost
  }
}

void SequenceHandler::report()
{
  cout << "Writing report:" << endl;
  for(unsigned i=0;i<sequences_.size(); i++)
    {
      Double mincost = LDBL_MAX;
      for (unsigned j = 0; j < sequencecosts[i].size(); j++)
        {
          if (sequencecosts[i][j] < mincost)
            mincost = sequencecosts[i][j];
        }
      cout << sequences_[i] << " " << costs_[i]-mincost << endl;
    }
  cout << "Best candidate sequence for negative design : " << sequences_[bestsequence] << endl;
  if (ToulBar2::diffneg)
    {
      cout << "Energy : " << bestcost  << endl;
      cout << "Energy difference with positive states : " << bestcost - costs_[bestsequence] << endl;
    }
  else
    cout << "Energy : " << bestcost << endl;
}
