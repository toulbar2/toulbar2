/*
 * ****** Variable with domain represented by an interval or an enumerated domain *******
 */
 
#include "tb2variable.hpp"
#include "tb2wcsp.hpp"


/*
 * Constructors and misc.
 * 
 */

Variable::Variable(string n, Value iinf, Value isup, Store *s) : enumerated(false), name(n), 
        inf(iinf, &s->storeValue), sup(isup, &s->storeValue), 
        value((iinf==isup)?iinf:(isup+1), &s->storeValue)
{
    if (s->getDepth() > 0) {
        cerr << "You cannot create a variable during the search!" << endl;
        exit(EXIT_FAILURE);
    }
}

Variable::Variable(string n, Value iinf, Value isup, Store *s, bool enumerate) : enumerated(true), name(n),
        inf(iinf, &s->storeValue), sup(isup, &s->storeValue),
        value((iinf==isup)?iinf:(isup+1), &s->storeValue),
        domain(iinf, isup, &s->storeDomain)
{
    assert(enumerate == true);
    if (s->getDepth() > 0) {
        cerr << "You cannot create a variable during the search!" << endl;
        exit(EXIT_FAILURE);
    }
}

Variable::Variable(string n, Value *d, int dsize, Store *s, bool enumerate) : enumerated(true), name(n),
        inf(min(d,dsize), &s->storeValue), sup(max(d, dsize), &s->storeValue),
        value((inf==sup)?inf:(sup+1), &s->storeValue),
        domain(d, dsize, &s->storeDomain)
{
    assert(enumerate == true);
    if (s->getDepth() > 0) {
        cerr << "You cannot create a variable during the search!" << endl;
        exit(EXIT_FAILURE);
    }
}

int Variable::getDegree()
{
    int sum = 0;
    for (unsigned int i=0; i<wcsps.size(); i++) {
    CostVariable *x = wcsps[i].wcsp->getVar(wcsps[i].wcspIndex);
    if (x != NULL) sum += x->getDegree();
    }
    return sum;
}


/*
 * Propagation methods
 * 
 */

void Variable::increase(Value newInf)
{
    if (ToulBar2::verbose >= 2) cout << "increase " << getName() << " " << inf << " -> " << newInf << endl;
    if (newInf > inf) {
        if (newInf > sup) throw Contradiction();
        else {
            if (enumerated) inf = domain.increase(newInf);
            else inf = newInf;
            if (inf == sup) assign(inf);
            else {
                for (unsigned int i=0; i<wcsps.size(); i++) {
                    wcsps[i].wcsp->increase(wcsps[i].wcspIndex);
                }
            }
        }
      }
}

void Variable::decrease(Value newSup)
{
    if (ToulBar2::verbose >= 2) cout << "decrease " << getName() << " " << sup << " -> " << newSup << endl;
    if (newSup < sup) {
        if (newSup < inf) throw Contradiction();
        else {
            if (enumerated) sup = domain.decrease(newSup);
            else sup = newSup;
            if (inf == sup) assign(sup);
            else {
                for (unsigned int i=0; i<wcsps.size(); i++) {
                    wcsps[i].wcsp->decrease(wcsps[i].wcspIndex);
                }
            }
        }
      }
}

void Variable::assign(Value newValue)
{
    if (ToulBar2::verbose >= 2) cout << "assign " << *this << " -> " << newValue << endl;
    if (unassigned() || value != newValue) {
        if (cannotbe(newValue)) throw Contradiction();
        value = newValue;
        inf = newValue;
        sup = newValue;
        for (unsigned int i=0; i<wcsps.size(); i++) {
            wcsps[i].wcsp->assign(wcsps[i].wcspIndex);
        }
    }
}

void Variable::remove(Value val)
{
    if (ToulBar2::verbose >= 2) cout << "remove " << *this << " <> " << val << endl;
    if (val == inf) increase(inf + 1);
    else if (val == sup) decrease(sup - 1);
    else if (enumerated && domain.canbe(val)) {
        domain.erase(val);
        for (unsigned int i=0; i<wcsps.size(); i++) {
            wcsps[i].wcsp->remove(wcsps[i].wcspIndex, val);
        }
    }
}
