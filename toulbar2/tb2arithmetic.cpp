/*
 * ****** Binary arithmetic soft constraints ******
 */

#include "tb2arithmetic.hpp"


/*
 * Constructors and misc.
 * 
 */

Supxyc::Supxyc(CostVariable *xx, CostVariable *yy, Value c, StoreStack<Cost,Cost> *storeCost, StoreStack<Value,Value> *storeValue) :
    AbstractBinaryConstraint(xx,yy), cst(c), deltaCost(0, storeCost),
    deltaValueXinf(xx->getSup()+1, storeValue), deltaValueYsup(yy->getSup()+1, storeValue),
    deltaCostXinf(0, storeCost), deltaCostYsup(0, storeCost) {}

void Supxyc::print(ostream& os)
{
    os << this << " " << x->getName() << " >= " << y->getName() << " + " << cst << endl;
}

/*
 * Propagation methods
 * 
 */

void Supxyc::propagate()
{
    if (ToulBar2::verbose >= 3) {
        print(cout);
        cout << " delta=" << deltaCost << " dxinf=" << deltaValueXinf << " dxcost=" << deltaCostXinf << " dysup=" << deltaValueYsup << " dycost=" << deltaCostYsup << endl;
    }

    // deconnect the constraint if always satisfied
    if (getX()->getInf() >= getY()->getSup() + cst) {
        deconnect();
    } else {
        Cost cost;
        
        // propagate hard constraint
        Value newInf = getY()->getInf() + cst - deltaCost - (wcsp->getUb() - wcsp->getLb());
        if (getX()->getInf() < newInf) getX()->increase(newInf);
        
        Value newSup = getX()->getSup() - cst + deltaCost + (wcsp->getUb() - wcsp->getLb());
        if (getY()->getSup() > newSup) getY()->decrease(newSup);
    
        // IC0 propagatation (increase global lower bound)
        cost = getY()->getInf() + cst - getX()->getSup() - deltaCost;
        if (getY()->getInf() == deltaValueYsup) cost -= deltaCostYsup;
        if (getX()->getSup() == deltaValueXinf) cost -= deltaCostXinf;
        if (cost > 0) {
            deltaCost += cost;
            if (ToulBar2::verbose >= 2) cout << "lower bound increased " << wcsp->getLb() << " -> " << wcsp->getLb()+cost << endl;
            wcsp->increaseLb(wcsp->getLb() + cost);
        }
        
        // BAC* propagation (increase unary costs of domain bounds)
        Value xinf = getX()->getInf();
        Value yinf = getY()->getInf();
        cost = yinf + cst - xinf - deltaCost;
        if (yinf == deltaValueYsup) cost -= deltaCostYsup;
        if (xinf == deltaValueXinf) cost -= deltaCostXinf;
        else deltaCostXinf = 0;
        if (cost > 0) {
            deltaValueXinf = xinf;
            deltaCostXinf += cost;
            getX()->projectInfCost(cost);
        }
        
        Value xsup = getX()->getSup();
        Value ysup = getY()->getSup();
        cost = ysup + cst - xsup - deltaCost;
        if (xsup == deltaValueXinf) cost -= deltaCostXinf;
        if (ysup == deltaValueYsup) cost -= deltaCostYsup;
        else deltaCostYsup = 0;
        if (cost > 0) {
            deltaValueYsup = ysup;
            deltaCostYsup += cost;
            getY()->projectSupCost(cost);
        }
    }
}
    
bool Supxyc::verify()
{
    Cost cmin,cxinf,cysup;
    
    cmin = getY()->getInf() + cst - getX()->getSup() - deltaCost
            - ((getY()->getInf() == deltaValueYsup)?deltaCostYsup:0)
            - ((getX()->getSup() == deltaValueXinf)?deltaCostXinf:0);
    if (cmin > 0) cout << "cmin=" << cmin << endl;
    cxinf = getY()->getInf() + cst - getX()->getInf() - deltaCost
            - ((getY()->getInf() == deltaValueYsup)?deltaCostYsup:0)
            - ((getX()->getInf() == deltaValueXinf)?deltaCostXinf:0);
    if (cxinf > 0) cout << "cxinf=" << cxinf << endl;
    cysup = getY()->getSup() + cst - getX()->getSup() - deltaCost
            - ((getY()->getSup() == deltaValueYsup)?deltaCostYsup:0)
            - ((getX()->getSup() == deltaValueXinf)?deltaCostXinf:0);
    if (cysup > 0) cout << "cysup=" << cysup << endl;
    bool icbac = (cmin <= 0) && (cysup <= 0) && (cxinf <= 0);
    if (!icbac) {
        print(cout);
        cout << " delta=" << deltaCost << " dxinf=" << deltaValueXinf << " dxcost=" << deltaCostXinf << " dysup=" << deltaValueYsup << " dycost=" << deltaCostYsup << endl;
    }
    return icbac;
}
