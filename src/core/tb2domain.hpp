/** \file tb2domain.hpp
 *  \brief Storable enumerated domain.
 *
 */

#ifndef TB2DOMAIN_HPP_
#define TB2DOMAIN_HPP_

#include "utils/tb2btlist.hpp"

extern int cmpValue(const void* v1, const void* v2);

class Domain : public BTList<Value> {
    bool contiguous;
    unsigned int initSize;
    Value distanceToZero;
    Value notFoundValue;
    DLink<Value>* all;
    map<Value, unsigned int> mapping;

    void init(Value inf, Value sup);
    void init(vector<Value>& dom);

    // make it private because we don't want copy nor assignment
    Domain(const Domain& s);
    Domain& operator=(const Domain& s);

public:
    typedef BTList<Value>::iterator iterator;

    Domain(Value inf, Value sup);

    Domain(vector<Value>& dom);

    ~Domain()
    {
        if (initSize >= 1)
            delete[] all;
    }

    void shrink(Value inf, Value sup);
    unsigned int getInitSize() const { return initSize; }
    unsigned int get(Value v, unsigned int notFound) const
    {
        auto it = mapping.find(v);
        return ((it == mapping.end()) ? notFound : it->second);
    }
    Value getinv(unsigned int idx, Value notFound) const { return ((idx < initSize) ? all[idx].content : notFound); }
    unsigned int toIndex(Value v) const { return ((contiguous) ? (v - distanceToZero) : get(v, initSize)); }
    Value toValue(int idx) const { return ((contiguous) ? (idx + distanceToZero) : getinv(idx, notFoundValue)); }
    unsigned int toCurrentIndex(Value v)
    {
        assert(canbe(v));
        unsigned int pos = 0;
        for (iterator iter = begin(); iter != end(); ++iter) {
            if (*iter == v)
                return pos;
            pos++;
        }
        cerr << "Bad (removed) value given as argument of toCurrentIndex function!" << endl;
        throw InternalError();
    }

    bool canbe(Value v) const { return !all[toIndex(v)].removed; }
    bool cannotbe(Value v) const { return all[toIndex(v)].removed; }

    void erase(Value v) { BTList<Value>::erase(&all[toIndex(v)], true); }

    Value increase(Value v)
    {
        iterator newInf = lower_bound(v);
        assert(canbe(*newInf));
        for (iterator iter = begin(); iter != newInf; ++iter) {
            erase(*iter);
        }
        return *newInf;
    }
    Value decrease(Value v)
    {
        iterator newSup = upper_bound(v);
        assert(canbe(*newSup));
        for (iterator iter = rbegin(); iter != newSup; --iter) {
            erase(*iter);
        }
        return *newSup;
    }

    // Finds the first available element whose value is greater or equal to v
    iterator lower_bound(Value v)
    {
        if (!contiguous) {
            auto it = mapping.lower_bound(v);
            if (it != mapping.end()) {
                v = it->first;
            }
        }
        assert(toIndex(v) >= 0 && toIndex(v) < initSize);
        iterator iter(&all[toIndex(v)]);
        if (cannotbe(v)) {
            ++iter;
        }
        return iter;
    }

    // Finds the first available element whose value is lower or equal to v
    iterator upper_bound(Value v)
    {
        if (!contiguous) {
            auto it = mapping.lower_bound(v);
            if (it != mapping.end()) {
                if (it->first != v) {
                    --it;
                    assert(it != mapping.end());
                    assert(it->first < v);
                    v = it->first;
                }
            }
        }
        assert(toIndex(v) >= 0 && toIndex(v) < initSize);
        iterator iter(&all[toIndex(v)]);
        if (cannotbe(v)) {
            --iter;
        }
        return iter;
    }

    friend ostream& operator<<(ostream& os, Domain& l);
};

#endif /*TB2DOMAIN_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
