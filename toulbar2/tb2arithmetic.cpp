/*
 * ****** Binary arithmetic soft constraints ******
 */

#include "tb2arithmetic.hpp"
#include "tb2wcsp.hpp"


/*
 * Constructors and misc.
 * 
 */

Supxyc::Supxyc(WCSP *wcsp, Variable *xx, Variable *yy, Value c, 
               StoreStack<Cost, Cost> *storeCost, StoreStack<Value, Value> *storeValue) :
    AbstractBinaryConstraint<Variable,Variable>(wcsp, xx, yy), cst(c), deltaCost(0, storeCost),
    deltaValueXinf(xx->getSup()+1, storeValue), deltaValueYsup(yy->getSup()+1, storeValue),
    deltaCostXinf(0, storeCost), deltaCostYsup(0, storeCost)
{
    xx->queueInc();
    xx->queueDec();
    yy->queueInc();
    yy->queueDec();
}

void Supxyc::print(ostream& os)
{
    os << this << " " << x->getName() << " >= " << y->getName() << " + " << cst << endl;
}

void Supxyc::dump(ostream& os)
{
    os << "2 " << x->wcspIndex << " " << y->wcspIndex << " -1 >= " << cst << endl;
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
    if (x->getInf() >= y->getSup() + cst) {
        deconnect();
    } else {
        Cost cost;
        
        // propagate hard constraint
        Value newInf = ceil(y->getInf() + cst - deltaCost - (wcsp->getUb() - wcsp->getLb() - 1));
        if (x->getInf() < newInf) x->increase(newInf);
        
        Value newSup = floor(x->getSup() - cst + deltaCost + (wcsp->getUb() - wcsp->getLb() - 1));
        if (y->getSup() > newSup) y->decrease(newSup);
    
        // IC0 propagatation (increase global lower bound)
        cost = y->getInf() + cst - x->getSup() - deltaCost;
        if (y->getInf() == deltaValueYsup) cost -= deltaCostYsup;
        if (x->getSup() == deltaValueXinf) cost -= deltaCostXinf;
        if (cost > 0) {
            deltaCost += cost;
            if (ToulBar2::verbose >= 2) cout << "lower bound increased " << wcsp->getLb() << " -> " << wcsp->getLb()+cost << endl;
            wcsp->increaseLb(wcsp->getLb() + cost);
        }
        
        // BAC* propagation (increase unary costs of domain bounds)
        Value xinf = x->getInf();
        Value yinf = y->getInf();
        cost = yinf + cst - xinf - deltaCost;
        if (yinf == deltaValueYsup) cost -= deltaCostYsup;
        if (xinf == deltaValueXinf) cost -= deltaCostXinf;
        else deltaCostXinf = 0;
        if (cost > 0) {
            deltaValueXinf = xinf;
            deltaCostXinf += cost;
            x->projectInfCost(cost);
        }
        
        Value xsup = x->getSup();
        Value ysup = y->getSup();
        cost = ysup + cst - xsup - deltaCost;
        if (xsup == deltaValueXinf) cost -= deltaCostXinf;
        if (ysup == deltaValueYsup) cost -= deltaCostYsup;
        else deltaCostYsup = 0;
        if (cost > 0) {
            deltaValueYsup = ysup;
            deltaCostYsup += cost;
            y->projectSupCost(cost);
        }
    }
}
    
bool Supxyc::verify()
{
    Cost cmin,cxinf,cysup;
    
    cmin = y->getInf() + cst - x->getSup() - deltaCost
            - ((y->getInf() == deltaValueYsup)?(Cost) deltaCostYsup:0)
            - ((x->getSup() == deltaValueXinf)?(Cost) deltaCostXinf:0);
    if (cmin > 0) cout << "cmin=" << cmin << endl;
    cxinf = y->getInf() + cst - x->getInf() - deltaCost
            - ((y->getInf() == deltaValueYsup)?(Cost) deltaCostYsup:0)
            - ((x->getInf() == deltaValueXinf)?(Cost) deltaCostXinf:0);
    if (cxinf > 0) cout << "cxinf=" << cxinf << endl;
    cysup = y->getSup() + cst - x->getSup() - deltaCost
            - ((y->getSup() == deltaValueYsup)?(Cost) deltaCostYsup:0)
            - ((x->getSup() == deltaValueXinf)?(Cost) deltaCostXinf:0);
    if (cysup > 0) cout << "cysup=" << cysup << endl;
    bool icbac = (cmin <= 0) && (cysup <= 0) && (cxinf <= 0);
    if (!icbac) {
        print(cout);
        cout << " delta=" << deltaCost << " dxinf=" << deltaValueXinf << " dxcost=" << deltaCostXinf << " dysup=" << deltaValueYsup << " dycost=" << deltaCostYsup << endl;
    }
    return icbac;
}
