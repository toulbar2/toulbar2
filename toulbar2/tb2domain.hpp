/** \file tb2domain.hpp
 *  \brief Storable enumerated domain.
 * 
 */
 
#ifndef TB2DOMAIN_HPP_
#define TB2DOMAIN_HPP_

#include "tb2btlist.hpp"

extern int cmpValue(const void *v1, const void *v2);

class Domain : public BTList<Value>
{
    const int initSize;
    const Value distanceToZero;
    DLink<Value> *all;
    
    // make it private because we don't want copy nor assignment
    Domain(const Domain &s);
    Domain& operator=(const Domain &s);

    void init(Value inf, Value sup) {
        assert( sup - inf + 1 >= 1 );
        assert( sup - inf + 1 <= MAX_DOMAIN_SIZE );
        all = new DLink<Value>[sup-inf+1];
        for (int idx=0; idx<sup-inf+1; idx++) {
            all[idx].content = idx + inf;
            push_back(&all[idx], false);
        }
    }
        
public:
    typedef BTList<Value>::iterator iterator;

    // used by Variable represented by an interval only
    Domain() : BTList<Value>(NULL), initSize(0), distanceToZero(0) {}
    
    Domain(Value inf, Value sup, StoreStack<BTList<Value>, DLink<Value> *> *s)
            : BTList<Value>(s), initSize(sup - inf + 1), distanceToZero(inf) {
        init(inf, sup);
    }
    
    Domain(Value *d, int dsize, StoreStack<BTList<Value>, DLink<Value> *> *s) 
            : BTList<Value>(s), initSize(max(d,dsize)-min(d,dsize)+1), distanceToZero(min(d,dsize)) {
        assert( dsize >= 1 );
        assert( dsize <= MAX_DOMAIN_SIZE );
        qsort(d, dsize, sizeof(Value), cmpValue);
        init(d[0], d[dsize-1]);
        int i = 0;
        for (iterator iter = begin(); iter != end(); ++iter) {
            if (*iter < d[i]) BTList<Value>::erase(&all[toIndex(*iter)], false);
            else i++;
        }
    }
    
    ~Domain() {if (initSize >= 1) delete[] all;}
        
    int getInitSize() const {return initSize;}
    int toIndex(Value v) const {return v - distanceToZero;}
    Value toValue(int idx) const {return idx + distanceToZero;}
    
    bool canbe(Value v) const {return !all[toIndex(v)].removed;}
    bool cannotbe(Value v) const {return all[toIndex(v)].removed;}

    void erase(Value v) {BTList<Value>::erase(&all[toIndex(v)], true);}

    int increase(Value v) {
        iterator newInf = lower_bound(v);
        assert(canbe(*newInf));
        for (iterator iter = begin(); iter != newInf; ++iter) {
            erase(*iter);
        }
        return *newInf;
    }
    int decrease(Value v) {
        iterator newSup = upper_bound(v);
        assert(canbe(*newSup));
        for (iterator iter = rbegin(); iter != newSup; --iter) {
            erase(*iter);
        }
        return *newSup;
    }
    
    //Finds the first available element whose value is greater or equal to v
    iterator lower_bound(Value v) {
        iterator iter(&all[toIndex(v)]);
        if (cannotbe(v)) {
            ++iter;
        }
        return iter;
    }
    
    //Finds the first available element whose value is lower or equal to v
    iterator upper_bound(Value v) {
        iterator iter(&all[toIndex(v)]);
        if (cannotbe(v)) {
            --iter;
        }
        return iter;
    }
    
    friend ostream& operator<<(ostream& os, Domain &l) {
        os << "{";
        for (Domain::iterator iter = l.begin(); iter != l.end(); ++iter) {
            os << " " << *iter;
        }
        os << " }";
        return os;
    }    
};

#endif /*TB2DOMAIN_HPP_*/
