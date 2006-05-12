/*
 * ****** Variable with domain represented by an interval *******
 */
 
#include "tb2intervar.hpp"
#include "tb2wcsp.hpp"


/*
 * Constructors and misc.
 * 
 */


IntervalVariable::IntervalVariable(WCSP *w, string n, Value iinf, Value isup) : 
        Variable(w, n, iinf, isup),
        infCost(0, &w->getStore()->storeCost), supCost(0, &w->getStore()->storeCost)
{
}

void IntervalVariable::print(ostream& os)
{
    os << " [" << inf << "," << sup << "]";
    os << "/" << getDegree();
    if (unassigned()) {
        os << " < " << getInfCost() << "," << getSupCost() << " >";
    }
}


/*
 * Propagation methods
 * 
 */

void IntervalVariable::projectInfCost(Cost cost)
{
    assert(cost >= 0);
    Cost oldcost = getInfCost();
    infCost += cost;
    Cost newcost = oldcost + cost;
    if (getInf() == maxCostValue || newcost > maxCost) queueNC();
    if (newcost + wcsp->getLb() >= wcsp->getUb()) increaseFast(getInf() + 1);     // Avoid any unary cost overflow
    if (getSup() == getInf() + 1 && getInfCost() > 0 && getSupCost() > 0) {
        Cost minCost = min(getInfCost(),getSupCost());
        extendAll(minCost);
        if (ToulBar2::verbose >= 2) cout << "lower bound increased " << wcsp->getLb() << " -> " << wcsp->getLb()+minCost << endl;
        wcsp->increaseLb(wcsp->getLb() + minCost);
    }
}

void IntervalVariable::projectSupCost(Cost cost)
{
    assert(cost >= 0);
    Cost oldcost = getSupCost();
    supCost += cost;
    Cost newcost = oldcost + cost;
    if (getSup() == maxCostValue || newcost > maxCost) queueNC();
    if (newcost + wcsp->getLb() >= wcsp->getUb()) decreaseFast(getSup() - 1);     // Avoid any unary cost overflow
    if (getSup() == getInf() + 1 && getInfCost() > 0 && getSupCost() > 0) {
        Cost minCost = min(getInfCost(),getSupCost());
        extendAll(minCost);
        if (ToulBar2::verbose >= 2) cout << "lower bound increased " << wcsp->getLb() << " -> " << wcsp->getLb()+minCost << endl;
        wcsp->increaseLb(wcsp->getLb() + minCost);
    }
}

void IntervalVariable::propagateNC()
{
    if (ToulBar2::verbose >= 3) cout << "propagateNC for " << getName() << endl;
    if (getInfCost() + wcsp->getLb() >= wcsp->getUb()) increaseFast(getInf() + 1);
    if (getSupCost() + wcsp->getLb() >= wcsp->getUb()) decreaseFast(getSup() - 1);
    if (getInfCost() > getSupCost()) {
        setMaxUnaryCost(getInf(), getInfCost());
    } else {
        setMaxUnaryCost(getSup(), getSupCost());
    }
}

bool IntervalVariable::verifyNC()
{
    bool supported = false;
    if (getInfCost() + wcsp->getLb() >= wcsp->getUb()) {
        cout << *this << " has inf cost not NC!" << endl;
        return false;
    }
    if (getSupCost() + wcsp->getLb() >= wcsp->getUb()) {
        cout << *this << " has sup cost not NC!" << endl;
        return false;
    }
    supported = (getDomainSize() > 2 || getInfCost() == 0 || getSupCost() == 0);
    if (!supported) cout << *this << " not NC*!" << endl;
    return supported;
}

void IntervalVariable::increaseFast(Value newInf)
{
    if (ToulBar2::verbose >= 2) cout << "increase " << getName() << " " << inf << " -> " << newInf << endl;
    if (newInf > inf) {
        if (newInf > sup) THROWCONTRADICTION;
        else {
            if (newInf == sup) assign(newInf);
            else {
                inf = newInf;
                infCost = deltaCost;
                queueInc();
                if (ToulBar2::setmin) (*ToulBar2::setmin)(wcsp->getIndex(), wcspIndex, newInf);
            }
        }
      }
}

void IntervalVariable::increase(Value newInf)
{
    if (ToulBar2::verbose >= 2) cout << "increase " << getName() << " " << inf << " -> " << newInf << endl;
    if (newInf > inf) {
        if (newInf > sup) THROWCONTRADICTION;
        else {
            if (newInf == sup) assign(newInf);
            else {
                inf = newInf;
                infCost = deltaCost;
                if (newInf > maxCostValue) queueNC();           // single diff with increaseFast
                queueInc();
                if (ToulBar2::setmin) (*ToulBar2::setmin)(wcsp->getIndex(), wcspIndex, newInf);
            }
        }
      }
}

void IntervalVariable::decreaseFast(Value newSup)
{
    if (ToulBar2::verbose >= 2) cout << "decrease " << getName() << " " << sup << " -> " << newSup << endl;
    if (newSup < sup) {
        if (newSup < inf) THROWCONTRADICTION;
        else {
            if (inf == newSup) assign(newSup);
            else {
                sup = newSup;
                supCost = deltaCost;
                queueDec();
                if (ToulBar2::setmax) (*ToulBar2::setmax)(wcsp->getIndex(), wcspIndex, newSup);
            }
        }
      }
}

void IntervalVariable::decrease(Value newSup)
{
    if (ToulBar2::verbose >= 2) cout << "decrease " << getName() << " " << sup << " -> " << newSup << endl;
    if (newSup < sup) {
        if (newSup < inf) THROWCONTRADICTION;
        else {
            if (inf == newSup) assign(newSup);
            else {
                sup = newSup;
                supCost = deltaCost;
                if (newSup < maxCostValue) queueNC();           // single diff with decreaseFast
                queueDec();
                if (ToulBar2::setmax) (*ToulBar2::setmax)(wcsp->getIndex(), wcspIndex, newSup);
            }
        }
      }
}

void IntervalVariable::assign(Value newValue)
{
    if (ToulBar2::verbose >= 2) cout << "assign " << *this << " -> " << newValue << endl;
    if (unassigned() || getValue() != newValue) {
        if (cannotbe(newValue)) THROWCONTRADICTION;
        changeNCBucket(-1);
        maxCostValue = newValue;
        maxCost = 0;
        Cost cost = ((newValue==inf)?getInfCost():((newValue==sup)?getSupCost():0));
        inf = newValue;
        sup = newValue;
        if (cost > 0) {
            deltaCost += cost;
            if (ToulBar2::verbose >= 2) cout << "lower bound increased " << wcsp->getLb() << " -> " << wcsp->getLb()+cost << endl;
            wcsp->increaseLb(wcsp->getLb() + cost);
        }
        if (ToulBar2::setvalue) (*ToulBar2::setvalue)(wcsp->getIndex(), wcspIndex, newValue);
        for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
            (*iter).constr->assign((*iter).scopeIndex);
        }
    }
}
