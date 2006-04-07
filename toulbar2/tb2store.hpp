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
 * Warning! limitation on the second parameter of StoreStack template
 * 
 */

template <class T, class V>        // Warning! Conversions between V and int must exist (store base in content)
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
        pointers = new T *[(1 << powbckmemory)];
        content = new V[(1 << powbckmemory)];
        index = 0;
        base = 0;
        indexMax = (1 << (powbckmemory - 1));
        if (ToulBar2::verbose) {
            cout << (1 << powbckmemory) * sizeof(V) + (1 << powbckmemory) * sizeof(T *)
                 << " Bytes allocated for " << name << " stack." << endl;
        }
    }
    
    ~StoreStack() {delete[] pointers; delete[] content;}
    
    void store(T *x, V y) {
        if (index > 0) {
            index++;
            content[index] = y;
            pointers[index] = x;
        }
    }
    
    void store(T *x) {
        if (index > 0) {
            index++;
            content[index] = *x;
            pointers[index] = x;
        }
    }
    
    void store() {
        if (index >= indexMax) {
            cerr << name << " stack out of memory!" << endl;
            exit(EXIT_FAILURE);
        }
        content[++index] = (V) base;
        base = index;
    }

    void restore(int **adr, int *val, int x) {
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
            base = (int) content[y];
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
    StoreBasic(T vv, StoreStack<T,T> *s) : v(vv), store(s) {}

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
};

#endif /*TB2STORE_HPP_*/
