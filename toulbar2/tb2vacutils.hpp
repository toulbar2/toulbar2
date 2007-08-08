/** \file tb2vacutils.hpp
 *  \brief Set of useful classes to enforce VAC
 */
 
#ifndef TB2VACUTILS_HPP_
#define TB2VACUTILS_HPP_

#include "tb2binconstr.hpp"


/**
 * A class that stores useful information on values of variables to enforce VAC
 */
class VACValue {

private:
  /**
   * the real value
   */
  Value value;
  /**
   * whether the value is removed
   */
  bool removed;
  /**
   * whether the value is marked (useful for pass 2)
   */
  bool mark;
  /**
   * when the value is removed, the index of the variable with which the value was inconstent (useful for passes 1 and 2)
   */
  int killer;
  /**
   * the temporary variable k that applies on unary constraints (useful for passes 2 and 3)
   */
  Cost k;

public:
  /**
   * constructor
   * @param v the value of this value
   */
  VACValue ();
  /**
   * destructor
   */
  ~VACValue ();

  /**
   * set all members of the class -- @a value
   * @param v the value of this element
   */
  void reset (const Value v);
  /**
   * whether the value is marked
   */
  bool isMarked () const;
  /**
   * whether the value is removed
   */
  bool isRemoved () const;
  /**
   * set @a killer
   * @return the index of the variable with which the value was inconstent, or -1 if the value is not removed
   */
  int getKiller () const;
  /**
   * get the value of the temporary variable k
   * @return the value of k
   */
  Cost getK () const;
  /**
   * accessor to @a value
   * @return the value of this instance
   */
  Value getValue () const;
  /**
   * mark the value
   */
  void setMark ();
  /**
   * remove the value
   */
  void remove ();
  /**
   * set @a killer
   * @param i the index of the variable such that this value is inconstent with it
   */
  void setKiller (const int i);
  /**
   * add some cost to the temporary variable k
   * @param c the cost
   */

  bool isNull (const Cost c);
  void setThreshold (Cost c); 
  Cost myThreshold;


  void addToK (const Cost c);
  /**
   * print a trace of this instance
   * @param os an output stream
   * @param e this instance
   * @return a trace of this instance
   */ 
  friend ostream& operator<< (ostream& os, const VACValue &e) {
    os << e.value;
    if (e.killer != -1) {
      os << "  killer = " << e.killer;
    }
    if (e.mark) {
      os << "  X";
    }
    return os;
  }
};


/**
 * A class that stores useful information on variables to enforce VAC
 */
class VACVariable : public EnumeratedVariable {

public:

private:
  /**
   * the values of the domain
   */
  VACValue **vacValues;
  /**
   * a link the VAC queue
   */
  DLink<VariableWithTimeStamp> linkVACQueue;
  /**
   * a link a secondary VAC queue
   */
  DLink<Variable *> linkVAC2Queue;

  /**
   * allocate memory
   */
  void init ();

public:
  /**
   * constructor
   * @param wcsp a link to global variables
   * @param n a name of this variable
   * @param iinf the minimum value of the domain
   * @param isup the maximum value of the domain
   */
  VACVariable (WCSP *wcsp, string n, Value iinf, Value isup);
  /**
   * constructor
   * @param wcsp a link to global variables
   * @param n a name of this variable
   * @param d the values of the domain
   * @param dsize the number of values in the domain
   */
  VACVariable (WCSP *wcsp, string n, Value *d, int dsize);
  /**
   * destructor
   */
  ~VACVariable ();

  /**
   * set all members of the class
   */
  void clear ();

  /**
   * set all members of the class
   */
  void reset ();
  /**
   * whether the variable has no possible value
   * @return true, if the variable has no value
   */
  bool isEmpty () const;
  /**
   * get the @a i th value of the domain
   * @param i the index of the value
   */
  VACValue *getValue (const unsigned int i) const;
  /**
   * whether a value is removed
   * @param i the index of this value
   * @return true if the @a i th value is removed
   */
  bool isRemoved (const unsigned int i) const;
  /**
   * remove a value
   * @param i the index of this value
   */
  void removeValue (const unsigned int i);
  /**
   * get the unary cost of a value, using VAC representation of the values
   * @param i the index of this value
   */
  Cost getIniCost (const unsigned int i); 
  Cost getVACCost (const unsigned int i);
  /**
   * set the cost of a value
   * @param i the index of this value
   * @param c the cost decrease
   */
  void setCost (const unsigned int i, const Cost c);
  /**
   * decrease a unary cost
   * @param i the index of the value
   * @param c the cost
   */
  void decreaseCost (const unsigned int i, const Cost c);
  /**
   * increase a unary cost
   * @param i the index of the value
   * @param c the cost
   */
  void increaseCost (const unsigned int i, const Cost c);
  /**
   * project some cost on some value of this variable
   * @param i the value
   * @param c the projected cost
   */
  void VACproject (const unsigned int i, const Cost c);
  /**
   * extend some cost from some value of this variable
   * @param i the value
   * @param c the extended cost
   */
  void VACextend (const unsigned int i, const Cost c);
  /**
   * enqueue in the VAC queue
   */

  /**
   * initial number of unremoved values
   */
  int nbValues;

  /**
   * the current number of unremoved values
   */
  int nbVACValues;

  Cost myThreshold;
  void setThreshold (Cost c); 
  Cost getThreshold (); 
  bool isNull (Cost c); 
 /**
   * specialized threshold per variable
   */ 
   
   
  void queueVAC();
  /**
   * enqueue in the secondary VAC queue
   */
  void queueVAC2();
  /**
   * dequeue in the secondary VAC queue
   */
  void dequeueVAC2();
  /**
   * overload the function of Variable
   * project some cost to c0
   * @param cost the projected cost
   */
  virtual void extendAll(Cost cost);
  /**
   * overload the function of EnumeratedVariable
   * notify an assignment of this variable
   * @param newValue the value of this variable
   */
  virtual void assign(Value newValue);
  /**
   * overload the function of EnumeratedVariable
   * remove a value
   * @param value the value to be removed
   */
  virtual void remove(Value value);
  /**
   * overload the function of EnumeratedVariable
   * remove a value
   * @param value the value to be removed
   */
  virtual void removeFast(Value value);
  /**
   * overload the function of EnumeratedVariable
   * increase the unary cost of a value
   * @param value the value
   * @param cost the cost increase
   */
  virtual void project (Value value, Cost cost);
  /**
   * overload the function of EnumeratedVariable
   * increase the unary cost of a value
   * @param value the value
   * @param cost the cost increase
   */
  virtual void extend (Value value, Cost cost);
  /**
   * overload the function of EnumeratedVariable
   * move the lower bound of a domain
   * @param newInf the lower bound
   */
  virtual void increase (Value newInf);
  /**
   * overload the function of EnumeratedVariable
   * move the upper bound of a domain
   * @param newSup the new upper bound
   */
  virtual void decrease (Value newSup);

  /**
   * empty @a deltaCost variable, so that they are equal to 0
   */
  void repair ();
  /**
   * set the top value of all values, considering the @a s cost computed to far
   * param s the @a s computed (estimation of the cost of the optimal solution)
   */
  void setTop (const Cost s);
  /**
   * print a trace of this instance
   * @param os an output stream
   * @param v this instance
   * @return a trace of this instance
   */
  friend ostream& operator<< (ostream& os, VACVariable &v) {
    for (unsigned int i = 0; i < v.getDomainSize(); i++) {
      if (!v.isRemoved(i)) {
        os << *(v.vacValues[i]) << "  -  ";
      }
    }
    return os;
  }
};


/**
 * A class that stores information about a binary constraint
 */
class VACConstraint : public BinaryConstraint {

private:
  /**
   * the support of the values
   */
  unsigned int *support[2];
  /**
   * the temporary variable k that applies on binary constraints, used to increase the cost of the value (useful for passes 2 and 3)
   */
  Cost *k[2];

  /**
   * allocate memory
   * @param sizeX the size of the first variable
   * @param sizeY the size of the second variable
   */
  void init (int sizeX, int sizeY);


public:
  /**
   * constructor
   * @param wcsp a link to global variables
   * @param xx the first variable
   * @param yy the second variable
   * @param tab the costs of this constraint
   * @param storeCost the delta costs
   */
  VACConstraint (WCSP *wcsp, EnumeratedVariable *xx, EnumeratedVariable *yy, vector<Cost> &tab, StoreStack<Cost, Cost> *storeCost);

  /**
   * constructor
   * @param wcsp a link to global variables
   * @param storeCost the delta costs
   */
  VACConstraint (WCSP *wcsp, StoreStack<Cost, Cost> *storeCost);
  /**
   * destructor
   */
  ~VACConstraint ();

  /**
   * set all members of the class -- especially reset @a k and @a support
   */
  void reset ();
  /**
   * get the first variable
   * @return the first variable
   */
  VACVariable *getX() const;
  /**
   * get the second variable
   * @return the second variable
   */
  VACVariable *getY() const;
  /**
   * get the @a i th variable
   * @param i the index of the variable
   * @return the @a i th variable
   */
  VACVariable *getVariable(const int i) const;

  /**
   * @warn this procedure has no side effect!
   * project some cost from this binary constraint to a unary cost
   * @param v the variable on which the cost is projected
   * @param i the value of the variable on which the cost is projected
   * @param c the projected cost
   */
  void VACproject (const int v, const unsigned int i, const Cost c);
  /**
   * @warn this procedure has no side effect!
   * extend some cost from a unary cost to this binary constraint
   * @param v the variable from which the cost is extended
   * @param i the value of which the variable from which the cost is extended
   * @param c the extended cost
   */
  void VACextend (const int v, const unsigned int i, const Cost c);
  /**
   * get the binary cost between values @a v and @a w, using VAC representation of the values
   * @param v the value of the first variable
   * @param w the value of the second variable
   * @param i is 0, if the variables are correctly ordered, 1 if one should evaluate c(@a w, @a v) instead
   * @return the binary cost
   */
    
  Cost getIniCost (const unsigned int v, const unsigned int w, const int i = 0);
  Cost getVACCost (const unsigned int v, const unsigned int w, const int i = 0);
  /**
   * set the cost of a tuple
   * @param i the index of the x variable
   * @param j the index of the y variable
   * @param i is 0, if the variables are correctly ordered, 1 if one should evaluate c(@a w, @a v) instead
   * @param c the cost decrease
   */
  void setCost (const unsigned int v, const unsigned int w, const int i, const Cost c);
  /**
   * decrease the cost of a tuple
   * @param i the index of the x variable
   * @param j the index of the y variable
   * @param i is 0, if the variables are correctly ordered, 1 if one should evaluate c(@a w, @a v) instead
   * @param c the cost decrease
   */
  void decreaseCost (const unsigned int v, const unsigned int w, const int i, const Cost c);
  /**
   * empty @a deltaCosts variable, so that they are equal to 0
   */
  void repair ();
  /**
   * set the top value of all values, considering the @a s cost computed to far
   * param s the @a s computed (estimation of the cost of the optimal solution)
   */
  void setTop (const Cost s);
  /**
   * get the temporary variable k
   * @param i the variable that should receive a cost of k to make its unary cost increase
   * @param v the value of the variable that should receive a cost of k to make its unary cost increase
   * @return temporary variable k
   */
  Cost getK (const int i, const unsigned int v) const;
  /**
   * check whether a value still has a support with respect to this constraint
   * @param i the variable that should be checked
   * @param v the value of the variable that should be checked
   * @return true, if the value @a v should be removed
   */
  bool revise (const int i, const unsigned int v);
  /**
   * set the temporary variable k
   * @param i the variable that should receive a cost of k to make its unary cost increase
   * @param v the value of the variable that should receive a cost of k to make its unary cost increase
   * @param c the cost increment
   */
  void setK (const int i, const unsigned int v, const Cost c);
  
  /**
   * print a trace of this instance
   * @param os an output stream
   * @param c this instance
   * @return a trace of this instance
   */
  friend ostream& operator<< (ostream& os, VACConstraint &c) {
    VACConstraint *cp = &c;
    BinaryConstraint *cpSuper = dynamic_cast<BinaryConstraint*>(cp);
    os << (*cpSuper);
    //return ((BinaryConstraint) c);
    return os;
  }
};


/**
 * A class that stores odds about the number of variables used to enforce VAC
 */
class VACOddsRecorder {

private:
  /**
   * number of variables in the instance
   */
  int nbVariables;
  /**
   * whether a variable has been used to enforced VAC
   */
  bool *used;
  /**
   * total number of variables used to enforce VAC
   */
  int nbTotalUsed;
  /**
   * minimum number of variables used to enforce VAC
   */
  int nbMinUsed;
  /**
   * maximum number of variables used to enforce VAC
   */
  int nbMaxUsed;
  /**
   * number of times VAC has been enforced
   */
  int nbCalls;
  /**
   * number of times the VAC c0 increase is rational
   */
  int nbRationalCalls;
  /**
   * number of times the VAC c0 increase is integer
   */
  int nbIntegerCalls;

public:
  /**
   * constructor
   * @param nv the number of variables in the instance
   */
  VACOddsRecorder (int nv);
  /**
   * destructor
   */
  ~VACOddsRecorder ();
  /**
   * set all values
   */
  void reset ();
  /**
   * set that a variable has been used to enforce VAC
   * @param i the index of the variable
   */
  void addVariable (const int i);
  /**
   * increase the number of times the VAC c0 increase is integer
   */
  void addNbRationalCalls ();
  /**
   * increase the number of times the VAC c0 increase is rational
   */
  void addNbIntegerCalls ();
  /**
   * compute the odds of the instance
   */
  void computeOdds ();
  /**
   * get the average number of variable used to enforce VAC
   * @return the average number of variable used
   */
  float getMeanUsed () const;
  /**
   * get the minimum number of variable used to enforce VAC
   * @return the minimum number of variable used
   */
  int getMinUsed () const;
  /**
   * get the minimum number of variable used to enforce VAC
   * @return the minimum number of variable used
   */
  int getMaxUsed () const;
  /**
   * get the number of times the VAC c0 increase is rational
   * @return the value @a nbRationalCalls
   */
  int getNbRationalCalls () const;
  /**
   * get the number of times the VAC c0 increase is integer
   * @return the value @a nbIntegerCalls
   */
  int getNbIntegerCalls () const;
};

#endif /*TB2VACUTILS_HPP_*/
