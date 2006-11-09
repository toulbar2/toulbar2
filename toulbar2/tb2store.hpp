/** \file tb2store.hpp
 *  \brief Generic storable data management.
 * 
 */

#ifndef TB2STORE_HPP_
#define TB2STORE_HPP_

#include "tb2types.hpp"

template <class T> class BTList;
template <class T> class DLink;
class Constraint;
class Variable;
class ConstraintLink;

/*
 * Storable stack
 * 
 */

template <class T, class V>
class StoreStack
{
    string name;
    T **pointers;
    V *content;
    int index;
    int indexMax;
    int base;

    // make it private because we don't want copy nor assignment
    StoreStack(const StoreStack &s);
    StoreStack& operator=(const StoreStack &s);
    
public:
    StoreStack(string s, int powbckmemory) : name(s) {        
        if (pow(2.,powbckmemory) >= INT_MAX) {
            cerr << "command-line initial memory size parameter " << powbckmemory << " power of two too large!" << endl;
            exit(EXIT_FAILURE);
        }
        indexMax = (int) pow(2.,powbckmemory);
        pointers = new T *[indexMax];
        content = new V[indexMax];
        index = 0;
        base = 0;
        if (ToulBar2::verbose) {
            cout << indexMax * (sizeof(V) + sizeof(T *))
                 << " Bytes allocated for " << name << " stack." << endl;
        }
    }
    
    ~StoreStack() {delete[] pointers; delete[] content;}

    void realloc() {
      T **newpointers = new T *[indexMax * 2];
      V *newcontent = new V[indexMax * 2];
      if (!newpointers || !newcontent) {
        cerr << name << " stack out of memory!" << endl;
        exit(EXIT_FAILURE);
      }
      for (int i = 0; i<indexMax; i++) {
        newpointers[i] = pointers[i];
        newcontent[i] = content[i];
      }
      delete[] pointers; delete[] content;
      pointers = newpointers;
      content = newcontent;
      indexMax *= 2;
      if (ToulBar2::verbose) {
        cout << indexMax * (sizeof(V) + sizeof(T *))
             << " Bytes allocated for " << name << " stack." << endl;
      }
    }
      
    void store(T *x, V y) {
        if (index > 0) {
            index++;
	        if (index >= indexMax) realloc();
            content[index] = y;
            pointers[index] = x;
        }
    }
    
    void store(T *x) {
        if (index > 0) {
            index++;
    	    if (index >= indexMax) realloc();
            content[index] = *x;
            pointers[index] = x;
        }
    }
    
    void store() {
	   index++;
	   if (index >= indexMax) realloc();
       pointers[index] = (T *) base;
       base = index;
    }

    void restore(int **adr, int *val, int x) {
        *adr[x] = val[x];
    }

    void restore(Long **adr, Long *val, int x) {
        *adr[x] = val[x];
    }

    void restore(Rational **adr, Rational *val, int x) {
        *adr[x] = val[x];
    }

    template <class Q> void restore(BTList<Q> **l, DLink<Q> **elt, int &x);

    void restore() {
        if (index > 0) {        // do nothing if already at depth = 0
            int x,y;
            
            x = index + 1;
            y = base;
            while (--x != y) {
                restore(pointers,content,x);
            }
            
            index = y - 1;
            base = (long) pointers[y];
        }
    }
};

/*
 * Storable basic types
 */
 
template <class T>
class StoreBasic
{
    T v;
    StoreStack<T,T> *store;

public:
    StoreBasic(T vv, StoreStack<T,T> *s) : v(vv), store(s) { }

    operator T() const {return v;}    // allows conversion from StoreBasic to T
     
    StoreBasic(const StoreBasic &elt) : v(elt.v), store(elt.store) {}

    StoreBasic &operator=(const StoreBasic &elt) {      // Warning! assignment has to be backtrackable
        if (&elt != this) {
            store->store(&v);
            v = elt.v;
        }
        return *this;
    }
    
    StoreBasic &operator=(const T vv) {
        store->store(&v);
        v = vv;
        return *this;
    }
    StoreBasic &operator+=(const T vv) {
        store->store(&v);
        v += vv;
        return *this;
    }
    StoreBasic &operator-=(const T vv) {
        store->store(&v);
        v -= vv;
        return *this;
    }
};

typedef StoreBasic<int> StoreInt;
typedef StoreBasic<Value> StoreValue;
typedef StoreBasic<Cost> StoreCost;


/*
 * Container for all storable stacks
 */
 
class Store
{
    int depth;
public:
    StoreStack<Value, Value> storeValue;
    StoreStack<BTList<Value>, DLink<Value> *> storeDomain;
    StoreStack<Cost, Cost> storeCost;
    StoreStack<BTList<ConstraintLink>, DLink<ConstraintLink> *> storeConstraint;
    StoreStack<BTList<Variable *>, DLink<Variable *> *> storeVariable;
    
    Store(int pow) : depth(0),
        storeValue("Value",pow), storeDomain("Domain",pow), storeCost("Cost",pow), 
        storeConstraint("Constraint",pow), storeVariable("Variable",pow) {}

    int getDepth() const {return depth;}
    
    void store() {
        depth++;
        storeCost.store();
        storeValue.store();
        storeDomain.store();
        storeConstraint.store();
        storeVariable.store();
    }
    
    void restore() {
        depth--;
        storeCost.restore();
        storeValue.restore();
        storeDomain.restore();
        storeConstraint.restore();
        storeVariable.restore();
    }
    
    void restore(int newDepth) {
        assert(depth >= newDepth);
        while (depth > newDepth) restore();
    }
};

#endif /*TB2STORE_HPP_*/
