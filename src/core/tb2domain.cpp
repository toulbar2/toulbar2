/*
 * **************** Storable enumerated domain **********************
 */

#include "tb2domain.hpp"

/*
 * Constructors and misc.
 *
 */

Domain::Domain(Value inf, Value sup)
    : BTList<Value>(&Store::storeDomain)
    , contiguous(true)
    , initSize(sup - inf + 1)
    , distanceToZero(inf)
    , notFoundValue(sup + 1)
{
    init(inf, sup);
}

Domain::Domain(vector<Value>& dom)
    : BTList<Value>(&Store::storeDomain)
    , contiguous(true)
    , initSize(*max_element(dom.begin(), dom.end()) - *min_element(dom.begin(), dom.end()) + 1)
    , distanceToZero(*min_element(dom.begin(), dom.end()))
    , notFoundValue(initSize + distanceToZero)
{
    assert(dom.size() >= 1);
    assert(dom.size() <= std::numeric_limits<tValue>::max());
    qsort(&dom[0], dom.size(), sizeof(Value), cmpValue);
    if (dom[dom.size() - 1] - dom[0] + 1 <= 2 * (int)dom.size()) {
        init(dom[0], dom[dom.size() - 1]);
        int i = 0;
        for (iterator iter = begin(); iter != end(); ++iter) {
            if (*iter < dom[i])
                BTList<Value>::erase(&all[toIndex(*iter)], false);
            else
                i++;
        }
    } else {
        contiguous = false;
        initSize = dom.size();
        init(dom);
        BTList<Value>::erase(&all[dom.size()], false);
    }
}

void Domain::init(Value inf, Value sup)
{
    assert(contiguous);
    assert(sup - inf + 1 >= 1);
    assert(sup - inf + 1 <= std::numeric_limits<tValue>::max());
#if defined(WCSPFORMATONLY) && !defined(NUMBERJACK)
    assert(distanceToZero == 0);
#endif
    all = new DLink<Value>[sup - inf + 1];
    for (int idx = 0; idx < sup - inf + 1; idx++) {
        all[idx].content = idx + inf;
        push_back(&all[idx], false);
    }
}

void Domain::init(vector<Value>& dom)
{
#if defined(WCSPFORMATONLY)
    assert(false);
#endif
    assert(!contiguous);
    assert(dom.size() <= std::numeric_limits<tValue>::max());
    all = new DLink<Value>[dom.size() + 1];
    for (unsigned int idx = 0; idx < dom.size(); idx++) {
        all[idx].content = dom[idx];
        push_back(&all[idx], false);
        mapping[dom[idx]] = idx;
    }
    all[dom.size()].content = notFoundValue;
    push_back(&all[dom.size()], false);
}

void Domain::shrink(Value inf, Value sup)
{
    assert(contiguous);
    assert(sup - inf + 1 >= 1);
    assert(sup - inf + 1 <= (Value)initSize);
#if defined(WCSPFORMATONLY)
    assert(inf == 0);
#endif
    assert(inf >= distanceToZero);
    for (int idx = 0; idx < sup - inf + 1; idx++) {
        all[idx] = all[idx + inf - distanceToZero];
    }
    initSize = sup - inf + 1;
    distanceToZero = inf;
}

int cmpValue(const void* v1, const void* v2)
{
    if (*((int*)v1) < *((int*)v2))
        return -1;
    else if (*((int*)v1) > *((int*)v2))
        return 1;
    else
        return 0;
}

ostream& operator<<(ostream& os, Domain& l)
{
    os << "{";
    for (Domain::iterator iter = l.begin(); iter != l.end(); ++iter) {
        os << " " << *iter;
    }
    os << " }";
    return os;
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
