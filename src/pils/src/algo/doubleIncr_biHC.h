#ifndef __doubleIncr_biHC_h
#define __doubleIncr_biHC_h

#include <random>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <algo/localSearch.h>
#include <base/solution.h>
#include <base/costFunction.h>
#include <base/neighborEval.h>
#include <base/incrNeighborEval.h>
#include <base/bst.h>

class DoubleIncr_biHC : public LocalSearch {
public:
  using LocalSearch::nEval;
  using LocalSearch::eval;
  using LocalSearch::id;

  DoubleIncr_biHC(std::mt19937 & _rng, 
                     CostFunction & _eval,
                     IncrNeighborEval & _neighborEval, 
                     unsigned long long _nEvalMax,
                     unsigned _flatMaxLocal,
                     std::ostream & _out) : LocalSearch(_eval), rng(_rng), nEvalMax(_nEvalMax), flatMaxLocal(_flatMaxLocal), neighborEval(_neighborEval), out(_out) {
    delta.resize(eval.size());
    neighborhoodSize = 0;
    for(unsigned i = 0; i < delta.size(); i++) {
      delta[i].resize(eval.n_values[i] - 1);
      neighborhoodSize += eval.n_values[i] - 1;
    }
  }
  
  void printDelta() {
    for(unsigned i = 0; i < delta.size(); i++) {
        for(unsigned k = 0; k < delta[i].size(); k++) {
            out << i << " " << (k+1) << " " << delta[i][k] << std::endl;
        }
    }
  }

  virtual void operator()(Solution & _solution) {
    unsigned long long nEvalLocal = 0;
    btree.makeEmpty();
    bool accept = true;

    unsigned nBest;
    Cost bestDelta;
    std::pair<int, int> bestNeigh;
    std::vector< Node* > iBest(neighborhoodSize);
    unsigned r, old_value = 0;
    unsigned flatLocal = 0;

    while (nEvalLocal < nEvalMax && accept && flatLocal <= flatMaxLocal) {

      if (nEvalLocal == 0) {
        // first iteration: 
        //   compute the initial deltas
        std::pair<int, int> neighbor(0, 1);
        for(size_t i = 0; i < delta.size(); i++) {
          for(size_t k = 0; k < delta[i].size(); k++) {
            neighbor.first  = i;
            neighbor.second = k + 1;
            delta[i][k] = neighborEval(_solution, neighbor) - _solution.fitness();
            // insert in the search tree
            btree.insert(delta[i][k], neighbor);
          }
        }
      } else {
        // only update delta according to the accepted move: iBest[r]
        unsigned i_move = iBest[r]->neighbor.first;
        std::pair<int, int> neighbor(i_move, 1);

        // update indice i_move with incremental evaluation (except back_move)
        unsigned back_move = eval.n_values[i_move] - iBest[r]->neighbor.second - 1;  // minus one because indice from 0
        unsigned back_val = iBest[r]->neighbor.second - 1;
        neighbor.second = back_move + 1;
          // remove old value
        btree.remove(delta[i_move][ back_move ], neighbor);
          // compute new value
        delta[i_move][ back_move ] = - delta[i_move][ back_val ];
          // save new value
        btree.insert(delta[i_move][ back_move ], neighbor);

        unsigned k;
        for(k = 0; k < back_move; k++) {
          neighbor.second = k + 1;          
          // remove old value
          btree.remove(delta[i_move][k], neighbor);
          // compute new value
          delta[i_move][k] = neighborEval(_solution, neighbor) - _solution.fitness();
          // save new value
          btree.insert(delta[i_move][k], neighbor);
        }
        for(k = back_move + 1; k < delta[i_move].size(); k++) {
          neighbor.second = k + 1;          
          // remove old value
          btree.remove(delta[i_move][k], neighbor);
          // compute new value
          delta[i_move][k] = neighborEval(_solution, neighbor) - _solution.fitness();
          // save new value
          btree.insert(delta[i_move][k], neighbor);
        }
        // update indice i linked to i_move with i_move < i
        unsigned move_value, i;
        for(unsigned kk = 0; kk < eval.links[i_move].size(); kk++) {
          i = eval.links[i_move][ kk ];
          neighbor.first = i;
          for(k = 0; k < delta[ i ].size(); k++) {
            move_value = (_solution[i] + k + 1) % eval.n_values[i];
            neighbor.second = k + 1;
            // remove old value
            btree.remove(delta[i][k], neighbor);
            // compute new value
            delta[i][k] += eval.energy2[i_move][i][ old_value ][ _solution[i] ] - eval.energy2[i_move][i][ _solution[i_move] ][ _solution[i] ]
                         - eval.energy2[i_move][i][ old_value ][ move_value ]   + eval.energy2[i_move][i][ _solution[i_move] ][ move_value ];
            // save new value
            btree.insert(delta[i][k], neighbor);
          }
        }

        // update indice i linked to i_move with i < i_move
        for(size_t kk = 0; kk < neighborEval.interactions[i_move].size(); kk++) {
          i = neighborEval.interactions[i_move][ kk ];
          neighbor.first  = i;
          for(k = 0; k < delta[i].size(); k++) {
            move_value = (_solution[i] + k + 1) % eval.n_values[i];
            neighbor.second = k + 1;
            // remove old value
            btree.remove(delta[i][k], neighbor);
            // compute new value
            delta[i][k] += eval.energy2[i][i_move][ _solution[i] ][ old_value ]       - eval.energy2[i][i_move][ _solution[i] ][ _solution[i_move] ]
                         + eval.energy2[i][i_move][ move_value ][ _solution[i_move] ] - eval.energy2[i][i_move][ move_value ][ old_value ];
            // save new value
            btree.insert(delta[i][k], neighbor);
          }
        }
      }

      //   find best neighbor
      std::pair<int, int> bestNeigh;
      btree.minimum(bestDelta, bestNeigh);
      // find all best minima
      btree.findall(bestDelta, nBest, iBest);

      // Selection: select one of the best (if there are ties)
      if (bestDelta <= 0) { //TODO: global switch between greedy hill climbing and steepest descent
        accept = true;

        if (bestDelta == 0) {
            flatLocal++;
        } else {
            flatLocal = 0;
        }

        if (nBest > 1) {
          r = runif(rng) % nBest;
        } else {
          r = 0;
        }

        // save value before move (for double incr eval)
        old_value = _solution[ iBest[r]->neighbor.first ];

        // move
        _solution.fitness(_solution.fitness() + bestDelta);
        _solution[ iBest[r]->neighbor.first ] = (_solution[ iBest[r]->neighbor.first ] + iBest[r]->neighbor.second) % eval.n_values[ iBest[r]->neighbor.first ];
      } else 
        accept = false;

      nEvalLocal++;
    }

    nEval += nEvalLocal;

  }
  
protected:
  std::mt19937 & rng;
  std::uniform_int_distribution<int> runif;

  unsigned long long nEvalMax;
  unsigned flatMaxLocal;

  // incremental delta
  unsigned neighborhoodSize;
  std::vector< std::vector<Cost> > delta;

  // neighborhood
  IncrNeighborEval & neighborEval;

  // to compute and update the minimum with BST
  BST btree;

  // output
  std::ostream & out;
};

#endif
