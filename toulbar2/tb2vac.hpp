/** \file tb2vac.hpp
 *  \brief Enforce VAC in a WCSP.
 */
 
#ifndef TB2VAC_HPP_
#define TB2VAC_HPP_

#include <stack>
#include <set>
#include <list>
#include "tb2btqueue.hpp"
#include "tb2vacutils.hpp"


class tVACStat;

typedef map<Cost,int> tScale;


/**
 * The class that enforces VAC
 */
class VACExtension {

private:


  WCSP *wcsp;
  Queue VAC;                                 // non backtrackable list; the queue AC2001 used inside VAC
  Queue SeekSupport;                         // non backtrackable list; collect all variables with a value removed due to binary constraints during pass 1
  BTQueue VAC2;                              // backtrackable list; updated during AC and EDAC
  Long nbIterations;						 // updated at each pass, we will use it as timeStamp
  int inconsistentVariable;				     // is used also to check after enforcePass1() if the network is VAC


  Cost itThreshold;					         // threshold iterative descent 
  int breakCycles; 
  tScale scaleCost;
  list<Cost> scaleVAC;

  
  Cost minlambda;
    
  stack< pair<int, Value> > *queueP;  	   // values removed by hard AC (useful for passes 1 and 2)
  stack< pair<int, Value> > *queueR;  	   // values necessarily removed to increase c0 (useful for passes 2 and 3)
  
  void enforcePass1 ();		   			  				    // enforce hard AC
  bool enforcePass1( VACVariable *xj, VACConstraint* cij);  // revise and k updates
  bool checkPass1 () const;
  void enforcePass2 ();        			   // find minimal set to enforce hard AC
  bool enforcePass3 ();					   // project costs to increase c0
  void enforcePass3VACDecomposition ();    // enforce VAC decomposition pass 3 (substract cost and decrease top)

  void reset();  						 

  map<int,tVACStat*> heapAccess;
  vector<tVACStat*>   heap;	  	
  Cost sumlb;
  Long nlb;
  Long sumvars;
  int sumk;
  int theMaxK;
  
  EnumeratedVariable* nearIncVar;
  Cost 				  atThreshold;
  
public:

  VACExtension (WCSP *w);
  ~VACExtension ();

  bool firstTime() { return nbIterations == 0; } 
  
  bool isVAC () const;         						  //  should have enforced pass 1 before, if c0 could increase by enforcing VAC
  bool propagate ();
  bool isNull(Cost c) const;

  void clear();						                  // empty VAC queue 
  void queueVAC(DLink<VariableWithTimeStamp> *link);  
  void queueSeekSupport(DLink<VariableWithTimeStamp> *link);  
  void queueVAC2(DLink<Variable *> *link);
  void dequeueVAC2(DLink<Variable *> *link);
  
  void init();
  void iniThreshold();
  Cost getThreshold();
  void nextScaleCost();
  void histogram( Cost c );
  void histogram();

  set<int> singletonI;	  
  set<int> singleton;	  
  
  int varAssign; 
  Value valAssign;
  void assign(int varIndex, Value newValue) { varAssign = varIndex; valAssign = newValue; }

  void afterPreprocessing();
  void iniSingleton();
  void updateSingleton();
  void removeSingleton();

  int  getHeuristic();
  int  getVarACDom( int i );
  Cost getVarCostStat( int i );
  Long getVarTimesStat( int i );
  void updateStat(Cost lambda);
  void printStat(bool ini = false);
  void printTightMatrix();

  void minsumDiffusion();
};


#endif /*TB2VAC_HPP_*/
