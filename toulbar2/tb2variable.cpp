/*
 * ****** Variable with domain represented by an interval or an enumerated domain *******
 */
 
#include "tb2solver.hpp"


/*
 * Constructors and misc.
 * 
 */

Variable::Variable(string n, Value iinf, Value isup, Solver *s, VariableType t) : enumerated(false), type(t), name(n), 
        inf(iinf, &s->getStore()->storeValue), sup(isup, &s->getStore()->storeValue), 
        value((iinf==isup)?iinf:(isup+1), &s->getStore()->storeValue), solver(s)
{
    if (s->getStore()->getDepth() > 0) {
        cerr << "You cannot create a variable during the search!" << endl;
        exit(EXIT_FAILURE);
    }
    linkUnassignedVars.content = this;
    addVarToSolver();
}

Variable::Variable(string n, Value iinf, Value isup, Solver *s, VariableType t, bool enumerate) : 
        enumerated(true), type(t), name(n),
        inf(iinf, &s->getStore()->storeValue), sup(isup, &s->getStore()->storeValue),
        value((iinf==isup)?iinf:(isup+1), &s->getStore()->storeValue),
        domain(iinf, isup, &s->getStore()->storeDomain), solver(s)
{
    assert(enumerate == true);
    if (s->getStore()->getDepth() > 0) {
        cerr << "You cannot create a variable during the search!" << endl;
        exit(EXIT_FAILURE);
    }
    linkUnassignedVars.content = this;
    addVarToSolver();
}

Variable::Variable(string n, Value *d, int dsize, Solver *s, VariableType t, bool enumerate) : 
        enumerated(true), type(t), name(n),
        inf(min(d,dsize), &s->getStore()->storeValue), sup(max(d, dsize), &s->getStore()->storeValue),
        value((inf==sup)?inf:(sup+1), &s->getStore()->storeValue),
        domain(d, dsize, &s->getStore()->storeDomain), solver(s)
{
    assert(enumerate == true);
    if (s->getStore()->getDepth() > 0) {
        cerr << "You cannot create a variable during the search!" << endl;
        exit(EXIT_FAILURE);
    }
    linkUnassignedVars.content = this;
    addVarToSolver();
}

void Variable::addVarToSolver()
{
    if (type >= AUX_VAR) solver->getVars()->push_back(this);
    if (type >= DECISION_VAR && unassigned()) solver->getUnassignedVars()->push_back(&linkUnassignedVars, true);
} 

int Variable::getDegree()
{
    int sum = 0;
    for (unsigned int i=0; i<wcsps.size(); i++) {
        sum += wcsps[i].wcsp->getDegree(wcsps[i].wcspIndex);
    }
    return sum;
}


/*
 * Propagation methods
 * 
 */

void Variable::increase(WCSP *wcsp, Value newInf)
{
    if (ToulBar2::verbose >= 2) cout << "increase " << getName() << " " << inf << " -> " << newInf << endl;
    if (newInf > inf) {
        if (newInf > sup) throw Contradiction();
        else {
            if (enumerated) newInf = domain.increase(newInf);
            if (newInf == sup) assign(newInf);
            else {
            	inf = newInf;
                for (unsigned int i=0; i<wcsps.size(); i++) {
                    wcsps[i].wcsp->increase((wcsps[i].wcsp == wcsp), wcsps[i].wcspIndex);
                }
            }
        }
      }
}

void Variable::decrease(WCSP *wcsp, Value newSup)
{
    if (ToulBar2::verbose >= 2) cout << "decrease " << getName() << " " << sup << " -> " << newSup << endl;
    if (newSup < sup) {
        if (newSup < inf) throw Contradiction();
        else {
            if (enumerated) newSup = domain.decrease(newSup);
            if (inf == newSup) assign(newSup);
            else {
            	sup = newSup;
                for (unsigned int i=0; i<wcsps.size(); i++) {
                    wcsps[i].wcsp->decrease((wcsps[i].wcsp == wcsp), wcsps[i].wcspIndex);
                }
            }
        }
      }
}

void Variable::assign(Value newValue)
{
    assert(value > sup);
    if (ToulBar2::verbose >= 2) cout << "assign " << *this << " -> " << newValue << endl;
    if (unassigned() || value != newValue) {
        if (cannotbe(newValue)) throw Contradiction();
        Value lastInfBeforeAssign = inf;    // previous inf and sup values just before being assigned
        Value lastSupBeforeAssign = sup;
        value = newValue;
        inf = newValue;
        sup = newValue;
        if (type == DECISION_VAR) solver->getUnassignedVars()->erase(&linkUnassignedVars, true);
        for (unsigned int i=0; i<wcsps.size(); i++) {
            wcsps[i].wcsp->assign(wcsps[i].wcspIndex,lastInfBeforeAssign,lastSupBeforeAssign);
        }
    }
}

void Variable::remove(WCSP *wcsp, Value val)
{
    if (ToulBar2::verbose >= 2) cout << "remove " << *this << " <> " << val << endl;
    if (val == inf) increase(wcsp, inf + 1);
    else if (val == sup) decrease(wcsp, sup - 1);
    else if (enumerated && domain.canbe(val)) {
        domain.erase(val);
        for (unsigned int i=0; i<wcsps.size(); i++) {
            wcsps[i].wcsp->remove((wcsps[i].wcsp == wcsp), wcsps[i].wcspIndex, val);
        }
    }
}
