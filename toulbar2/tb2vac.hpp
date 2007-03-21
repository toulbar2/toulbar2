/** \file tb2vac.hpp
 *  \brief Enforce VAC in a WCSP.
 */
 
#ifndef TB2VAC_HPP_
#define TB2VAC_HPP_

#include <stack>
#include "tb2btqueue.hpp"
#include "tb2vacutils.hpp"

/**
 * The class that enforces VAC
 */
class VACExtension {

private:

  Queue VAC;                                // non backtrackable list
  BTQueue VAC2;                             // backtrackable list
  
  /**
   * a local copy of the pointer to the network;
   */
  WCSP *wcsp;
  /**
   * whether VAC alternative should be enforced
   */
  bool alternative;
  /**
   * whether VAC decomposition should be enforced
   */
  bool decomposition;
  /**
   * the number of times the three passes have be run to have a VAC instance
   */
  int nbIterations;
  /**
   * the maximum value of unary and binary k (useful for passes 2 and 3)
   */
  Cost maxK;
  /**
   * the index of the inconsistent variable
   *   or -1 if the network is VAC
   *   or -2 if the network is not NC
   */
  int inconsistentVariable;
  /**
   * if VAC alternative is chosen, when a cost is above the threshold, it is considered as 0
   */
  Cost costThreshold;
  /**
   * the last time the threshold was updated
   */
  Long lastCostThresholdUpdate;
    
  /**
   * a queue that stores values removed by hard AC (useful for passes 1 and 2)
   */
  stack< pair<int, int> > *queueP;
  /**
   * a queue that stores values necessarily removed to increase c0 (useful for passes 2 and 3)
   */
  stack< pair<int, int> > *queueR;
  /**
   * an estimation of the lower bound of the current instance
   */
  Cost s;
  /**
   * odds about the enforcement of VAC
   */
  VACOddsRecorder *vacOddsRecorder;

  /**
   * set all the useful variables
   */
  void reset ();
  /**
   * enforce pass 1 (enforce hard AC)
   */
  void enforcePass1 ();
  /**
   * check that pass 1 actually enforced hard AC (for debugging purpose)
   */
  bool checkPass1 () const;
  /**
   * @pre should have enforced pass 1 before
   * check whether c0 could increase by enforcing VAC
   * @return true, if c0 could increase by enforcing VAC
   */
  bool isVAC () const;
  /**
   * enforce pass 2 (find minimal set to enforce hard AC)
   */
  void enforcePass2 ();
  /**
   * enforce VAC or VAC decomposition pass 3 
   */
  void enforcePass3 ();
  /**
   * actually enforce VAC pass 3 (project costs to increase c0)
   */
  void enforcePass3VAC ();
  /**
   * actually enforce VAC decomposition pass 3 (substract cost and decrease top)
   */
  void enforcePass3VACDecomposition ();
  /**
   * get an estimation of the lower bound of the current instance
   * @return the value @a s;
   */
  Cost getS () const;
  /**
   * get odds about the enforcement of VAC
   * @return a class that stores odds about the enforcement of VAC
   */
  VACOddsRecorder *getVACOddsRecorder ();

public:
  /**
   * constructor
   * @param w a pointer to the network
   * @param a whether VAC alternative should be enforced
   * @param d whether VAC decomposition should be enforced
   */
  VACExtension (WCSP *w, bool a = false, bool d = false);
  /**
   * destructor
   */
  ~VACExtension ();
  /**
   * enforce VAC (possibly several times)
   */
  void propagate ();
  /**
   * check whether VAC was actually enforced
   * @return true, if VAC was actually enforced
   */
  bool verify ();
  /**
   * check whether VAC has been enforced on all variables
   * @return empty if VAC has been enforced on all variables
   */
  bool empty ();
  /**
   * clear the propagation queue of VAC (actually, push all elements in it)
   */
  void clear ();
  /**
   * whether VAC did not change the network
   * @return true, if @a nbIterations is null
   */
  bool remainedIdle ();
  /**
   * update the queue VAC2 (possibly pops elements)
   */
  void updateQueue ();
  /**
   * @pre supposes that VAC alternative is used
   * get the threshold above which a cost is considered as 0
   * @return the threshold
   */
  Cost getCostThreshold () const;
  /**
   * check is a cost is considered as 0
   * (in VAC Alternative, a cost not greater the @a costThreshold is considered as 0)
   * @param c a cost
   * @return true if the cost is considered as 0
   */
  bool isNull(Cost c) const;
  /**
   * @pre supposes that VAC alternative is used
   * get the last the threshold was updated
   * @return the timestamp
   */
  Long getLastCostThresholdUpdate () const;
  /**
   * @pre supposes that VAC alternative is used
   * update the threshold above which a cost is considered as 0
   * @param t the new threshold
   * @param d the current timestamp
   */
  void updateThreshold (Cost t, Long d);
  /**
   * queue a variable to the VAC queue
   * @param link the link the variable
   */
  void queueVAC(DLink<VariableWithTimeStamp> *link);
  /**
   * queue a variable to the VAC2 queue
   * @param link the link the variable
   */
  void queueVAC2(DLink<Variable *> *link);
  /**
   * unqueue a variable to the VAC queue
   * @param link the link the variable
   */
  void dequeueVAC2(DLink<Variable *> *link);
};


#endif /*TB2VAC_HPP_*/
