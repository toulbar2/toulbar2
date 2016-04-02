/** \file tb2store.hpp
 *  \brief Generic storable data management.
 *
 *	\defgroup backtrack Backtrack management
 *  Used by backtrack search methods.
 *  Allows to copy / restore the current state using Store::store and Store::restore methods.
 *  All storable data modifications are trailed into specific stacks.
 *
 *  Trailing stacks are associated to each storable type:
 *  - Store::storeValue for storable domain values ::StoreValue (value supports, etc)
 *  - Store::storeCost for storable costs ::StoreCost (inside cost functions, etc)
 *  - Store::storeDomain for enumerated domains (to manage holes inside domains)
 *  - Store::storeConstraint for backtrackable lists of constraints
 *  - Store::storeVariable for backtrackable lists of variables
 *  - Store::storeSeparator for backtrackable lists of separators (see tree decomposition methods)
 *  - Store::storeBigInteger for very large integers ::StoreBigInteger used in solution counting methods
 *
 *  Memory for each stack is dynamically allocated by part of \f$2^x\f$ with \e x initialized to ::STORE_SIZE and increased when needed.
 *  \note storable data are not trailed at depth 0.
 *  \warning ::StoreInt uses Store::storeValue stack (it assumes Value is encoded as int).
 */

#ifndef TB2STORE_HPP_
#define TB2STORE_HPP_

#include "tb2types.hpp"
#include <boost/version.hpp>
#if (BOOST_VERSION >= 105600)
#include <boost/type_index.hpp>
#else
#include <typeinfo>
#endif

template<class T> class BTList;
template<class T> class DLink;
class Constraint;
class Variable;
class Separator;
class ConstraintLink;

/*
 * Storable stack
 *
 */
template<class T, class V>
class StoreStack
{
    T **pointers;
    V *content;
    ptrdiff_t index;
    ptrdiff_t indexMax;
    ptrdiff_t base;

    // make it private because we don't want copy nor assignment
    StoreStack(const StoreStack &s);
    StoreStack &operator=(const StoreStack &s);

public:
    StoreStack(int powbckmemory = STORE_SIZE)  {
        if (pow(2., powbckmemory) >= SIZE_MAX) {
            cerr << "command-line initial memory size parameter " << powbckmemory << " power of two too large!" << endl;
            exit(EXIT_FAILURE);
        }
        indexMax = (ptrdiff_t) pow(2., powbckmemory);
        pointers = new T *[indexMax];
        content = new V[indexMax];
        index = 0;
        base = 0;
        if (ToulBar2::verbose >= 0) {
            cout << "c " << indexMax * (sizeof(V) + sizeof(T *)) << " Bytes allocated for " <<
#if (BOOST_VERSION >= 105600)
                boost::typeindex::type_id<T>().pretty_name()
#else
                typeid(T).name()
                #endif
                 << " stack." << endl;
        }
    }

    ~StoreStack() {
        delete[] pointers;
        delete[] content;
    }

    void realloc() {
        T **newpointers = new T *[indexMax * 2];
        V *newcontent = new V[indexMax * 2];
        if (!newpointers || !newcontent) {
            cerr << typeid(T).name() << " stack out of memory!" << endl;
            exit(EXIT_FAILURE);
        }
        std::copy(pointers,pointers+indexMax,newpointers);
        std::copy(content,content+indexMax,newcontent);

        delete[] pointers;
        delete[] content;
        pointers = newpointers;
        content = newcontent;
        indexMax *= 2;
        if (ToulBar2::verbose > 0) {
            cout << indexMax * (sizeof(V) + sizeof(T *)) << " Bytes allocated for " << typeid(T).name() << " stack." << endl;
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
        pointers[index] = (T *)(intptr_t) base;
        base = index;
    }

//	void restore(int **adr, int *val, ptrdiff_t x) {
//		*adr[x] = val[x];
//	}

    void restore(Value **adr, Value *val, ptrdiff_t x) {
        *adr[x] = val[x];
    }

#ifndef INT_COST
    void restore(Cost **adr, Cost *val, ptrdiff_t x) {
        *adr[x] = val[x];
    }
#endif

    void restore(BigInteger **adr, BigInteger *val, ptrdiff_t x) {
        *adr[x] = val[x];
    }
    template<class Q> void restore(BTList<Q> **l, DLink<Q> **elt, ptrdiff_t &x);

    void restore() {
        if (index > 0) { // do nothing if already at depth = 0
            ptrdiff_t x, y;

            x = index + 1;
            y = base;
            while (--x != y) {
                restore(pointers, content, x);
            }

            index = y - 1;
            base = (ptrdiff_t) pointers[y];
        }
    }
};



/*
 * Storable basic types
 */
template<class T>
class StoreBasic
{
    T v;

public:
    StoreBasic(T vv) : v(vv) {
    }

    operator T() const {
        return v;
    } ///< allows conversion from StoreBasic to T

    StoreBasic(const StoreBasic &elt) : v(elt.v) {}

    static void store() { mystore.store(); };
    static void restore() { mystore.restore(); };

    StoreBasic &operator=(const StoreBasic &elt) { ///< \note assignment has to be backtrackable
        if (&elt != this) {
            mystore.store(&v);
            v = elt.v;
        }
        return *this;
    }

    StoreBasic &operator=(const T vv) {
        mystore.store(&v);
        v = vv;
        return *this;
    }
    StoreBasic &operator+=(const T vv) {
        mystore.store(&v);
        v += vv;
        return *this;
    }
    StoreBasic &operator-=(const T vv) {
        mystore.store(&v);
        v -= vv;
        return *this;
    }

    static StoreStack<T, T> mystore;
};

template<class T>
StoreStack<T, T> StoreBasic<T>::mystore(STORE_SIZE);


//typedef StoreBasic<int> StoreInt;
typedef StoreBasic<Value> StoreValue;
typedef StoreValue StoreInt;
typedef StoreBasic<Cost> StoreCost;
typedef StoreBasic<BigInteger> StoreBigInteger;

/*
 * Container for all storable stacks
 */
class Store
{
    int depth;
    
public:
    StoreStack<BTList<Value> , DLink<Value> *> storeDomain;
    StoreStack<BTList<ConstraintLink> , DLink<ConstraintLink> *> storeConstraint;
    StoreStack<BTList<Variable *> , DLink<Variable *> *> storeVariable;
    StoreStack<BTList<Separator *> , DLink<Separator *> *> storeSeparator;

    Store(int pow) : depth(0) {
    }

    /// \return the current (backtrack / tree search) depth
    int getDepth() const {
        return depth;
    }

    /// makes a copy of the current state
    void store() {
        depth++;
        StoreCost::store();
        StoreValue::store();
        StoreBigInteger::store();
//		storeInt.store();
        storeDomain.store();
        storeConstraint.store();
        storeVariable.store();
        storeSeparator.store();
    }

    /// restores the current state to the last copy
    void restore() {
        depth--;
        StoreCost::restore();
        StoreValue::restore();
        StoreBigInteger::restore();
//		storeInt.restore();
        storeDomain.restore();
        storeConstraint.restore();
        storeVariable.restore();
        storeSeparator.restore();
    }

    /// restore the current state to the copy made at depth \c newDepth
    void restore(int newDepth) {
        assert(depth >= newDepth);
        while (depth > newDepth) restore();
    }
};

//#define storeInt storeValue
#define storeIndexList storeDomain

#endif /*TB2STORE_HPP_*/
