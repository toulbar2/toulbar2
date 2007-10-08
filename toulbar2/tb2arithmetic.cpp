/*
 * ****** Binary arithmetic soft constraints ******
 */

#include "tb2arithmetic.hpp"
#include "tb2wcsp.hpp"


/*
 * Constructors and misc.
 * 
 */

Unary::Unary(WCSP *wcsp, IntervalVariable *xx, Value *d, int dsize, Cost c, StoreStack<Value, Value> *storeValue) :
    AbstractUnaryConstraint<IntervalVariable>(wcsp, xx), permitted(d, d+dsize), penalty(c),
    deltaValueXinf(xx->getSup()+1, storeValue), deltaValueXsup(xx->getSup()+1, storeValue)
{
    xx->queueInc();
    xx->queueDec();
}

void Unary::print(ostream& os)
{
  os << this << " " << x->getName() << " var in {";
  for (set<Value>::iterator it=permitted.begin(); it != permitted.end(); ++it) {
	os << " " << *it;
  }
  os  << " } (" << penalty << ")" << endl;
}

void Unary::dump(ostream& os)
{
  os << "1 " << x->wcspIndex << " " << penalty << " " << permitted.size() << endl;
  for (set<Value>::iterator it=permitted.begin(); it != permitted.end(); ++it) {
	os << *it << " 0" << endl;
  }
}

Supxyc::Supxyc(WCSP *wcsp, Variable *xx, Variable *yy, Value c, Value delta,
               StoreStack<Cost, Cost> *storeCost, StoreStack<Value, Value> *storeValue) :
    AbstractBinaryConstraint<Variable,Variable>(wcsp, xx, yy), 
    cst(c), deltamax(delta), deltaCost(MIN_COST, storeCost),
    deltaValueXinf(xx->getSup()+1, storeValue), deltaValueYsup(yy->getSup()+1, storeValue),
    deltaCostXinf(MIN_COST, storeCost), deltaCostYsup(MIN_COST, storeCost)
{
    xx->queueInc();
    xx->queueDec();
    yy->queueInc();
    yy->queueDec();
}

void Supxyc::print(ostream& os)
{
  os << this << " " << x->getName() << " >= " << y->getName() << " + " << cst << " (" << deltamax << ")" << endl;
}

void Supxyc::dump(ostream& os)
{
  os << "2 " << x->wcspIndex << " " << y->wcspIndex << " -1 >= " << cst << " " << deltamax << endl;
}

Disjunction::Disjunction(WCSP *wcsp, Variable *xx, Variable *yy, Value cxx, Value cyy, Cost cost,
               StoreStack<Value, Value> *storeValue) :
    AbstractBinaryConstraint<Variable,Variable>(wcsp, xx, yy),
    cstx(cxx), csty(cyy), penalty(cost),
    deltaValueXinf(xx->getSup()+1, storeValue),deltaValueYinf(yy->getSup()+1, storeValue),
    deltaValueXsup(xx->getSup()+1, storeValue),deltaValueYsup(yy->getSup()+1, storeValue)
{
    xx->queueInc();
    xx->queueDec();
    yy->queueInc();
    yy->queueDec();
}

void Disjunction::print(ostream& os)
{
  os << this << " " << x->getName() << " >= " << y->getName() << " + " << csty << " or " << y->getName() << " >= " << x->getName() << " + " << cstx << " ("  << penalty<< ")" << endl;
}

void Disjunction::dump(ostream& os)
{
  os << "2 " << x->wcspIndex << " " << y->wcspIndex << " -1 disj " << cstx << " " << csty << " " << penalty << endl;
}

/*
 * Propagation methods
 * 
 */

void Unary::propagate()
{
    if (ToulBar2::verbose >= 3) {
	    print(cout);
        cout << " dxinf=" << deltaValueXinf << " dxsup=" << deltaValueXsup << endl;
    }

	assert(connected());
	set<Value>::iterator itinf = permitted.lower_bound(x->getInf());
	set<Value>::iterator itsup = permitted.upper_bound(x->getSup());
	--itsup;
	if (itinf == permitted.end() || itsup == permitted.end()) {
	  // IC0 propagatation (increase global lower bound)
	  deconnect();
	  wcsp->increaseLb(wcsp->getLb() + penalty);
	} else {
	  // propagate hard constraint
	  if (CUT(wcsp->getLb()+penalty, wcsp->getUb())) {
		if (x->getInf() < *itinf) x->increase(*itinf);
		if (x->getSup() > *itsup) x->decrease(*itsup);
	  }
	  // BAC* propagation (increase unary costs of domain bounds)
	  Value xinf = x->getInf();
	  if (xinf != deltaValueXinf && xinf != deltaValueXsup && permitted.find(xinf) == permitted.end()) {
		deltaValueXinf = xinf;
		x->projectInfCost(penalty);
	  }
	  Value xsup = x->getSup();
	  if (xsup != deltaValueXinf && xsup != deltaValueXsup && permitted.find(xsup) == permitted.end()) {
		deltaValueXsup = xsup;
		x->projectSupCost(penalty);
	  }
	}
}
     
bool Unary::verify()
{
  bool support = (permitted.lower_bound(x->getInf()) != permitted.end() || (--permitted.upper_bound(x->getSup())) != permitted.end() || deconnected());
  support = support && (permitted.find(x->getInf()) != permitted.end() || deltaValueXinf == x->getInf());
  support = support && (permitted.find(x->getSup()) != permitted.end() || deltaValueXsup == x->getSup());
  if (!support) {
	print(cout);
	x->print(cout);
	cout << endl;
	cout << " dxinf=" << deltaValueXinf << " dxsup=" << deltaValueXsup << endl;
  }
  return support;
}

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
        // propagate hard constraint
		Cost gap = wcsp->getUb() - wcsp->getLb() - 1 + deltaCost;
        Value newInf = ceil(y->getInf() + cst - ((gap < deltamax)?gap:deltamax));
        if (x->getInf() < newInf) x->increase(newInf);
        
        Value newSup = floor(x->getSup() - cst + ((gap < deltamax)?gap:deltamax));
        if (y->getSup() > newSup) y->decrease(newSup);
    
        // IC0 propagatation (increase global lower bound)
        Cost cost = y->getInf() + cst - x->getSup() - deltaCost;
		if (y->enumerated() && deltaCostYsup>MIN_COST && y->getInf() == deltaValueYsup) cost -= deltaCostYsup;
		if (x->enumerated() && deltaCostXinf>MIN_COST && x->getSup() == deltaValueXinf) cost -= deltaCostXinf;
        if (cost > MIN_COST) {
            deltaCost += cost;
            if (ToulBar2::verbose >= 2) cout << "lower bound increased " << wcsp->getLb() << " -> " << wcsp->getLb()+cost << endl;
            wcsp->increaseLb(wcsp->getLb() + cost);
        }
        
        // BAC* propagation (increase unary costs of domain bounds)
        Value xinf = x->getInf();
        Value yinf = y->getInf();
        cost = yinf + cst - xinf - deltaCost;
		if (y->enumerated() && deltaCostYsup>MIN_COST && yinf == deltaValueYsup) cost -= deltaCostYsup;
        if (deltaCostXinf>MIN_COST && xinf == deltaValueXinf) cost -= deltaCostXinf;
        else deltaCostXinf = 0;
        if (cost > MIN_COST) {
            deltaValueXinf = xinf;
			deltaCostXinf += cost;
            x->projectInfCost(cost);
        }
        
        Value xsup = x->getSup();
        Value ysup = y->getSup();
        cost = ysup + cst - xsup - deltaCost;
		if (x->enumerated() && deltaCostXinf>MIN_COST && xsup == deltaValueXinf) cost -= deltaCostXinf;
        if (deltaCostYsup>MIN_COST && ysup == deltaValueYsup) cost -= deltaCostYsup;
        else deltaCostYsup = 0;
        if (cost > MIN_COST) {
		    deltaValueYsup = ysup;
			deltaCostYsup += cost;
            y->projectSupCost(cost);
        }
    }
}
     
bool Supxyc::verify()
{
    Cost cmin=MIN_COST;
	Cost cxinf=MIN_COST;
	Cost cysup=MIN_COST;
    
    cmin = y->getInf() + cst - x->getSup() - deltaCost
            - ((y->getInf() == deltaValueYsup)?deltaCostYsup:MIN_COST)
            - ((x->getSup() == deltaValueXinf)?deltaCostXinf:MIN_COST);
    if (cmin > MIN_COST) cout << "cmin=" << cmin << endl;
    cxinf = y->getInf() + cst - x->getInf() - deltaCost
            - ((y->getInf() == deltaValueYsup)?deltaCostYsup:MIN_COST)
            - ((x->getInf() == deltaValueXinf)?deltaCostXinf:MIN_COST);
    if (cxinf > MIN_COST) cout << "cxinf=" << cxinf << endl;
    cysup = y->getSup() + cst - x->getSup() - deltaCost
            - ((y->getSup() == deltaValueYsup)?deltaCostYsup:MIN_COST)
            - ((x->getSup() == deltaValueXinf)?deltaCostXinf:MIN_COST);
    if (cysup > MIN_COST) cout << "cysup=" << cysup << endl;
    bool icbac = (cmin <= MIN_COST) && (cysup <= MIN_COST) && (cxinf <= MIN_COST);
    if (!icbac) {
	    print(cout);
		x->print(cout);
		cout << endl;
		y->print(cout);
		cout << endl;
        cout << " delta=" << deltaCost << " dxinf=" << deltaValueXinf << " dxcost=" << deltaCostXinf << " dysup=" << deltaValueYsup << " dycost=" << deltaCostYsup << endl;
    }
    return icbac;
}

void Disjunction::propagate()
{
  if (ToulBar2::verbose >= 3) {
	print(cout);
	cout << " dxinf=" << deltaValueXinf << " dxsup=" << deltaValueXsup << " dyinf=" << deltaValueYinf << " dysup=" << deltaValueYsup << endl;
  }
  assert(connected());
  assert(!x->enumerated() || !y->enumerated() || !((x->getInf()==deltaValueXinf || x->getSup()==deltaValueXsup) && (y->getInf()==deltaValueYinf || y->getSup()==deltaValueYsup)));

  // deconnect the constraint if always satisfied
  if (x->getInf() >= y->getSup() + csty || y->getInf() >= x->getSup() + cstx) {
	deconnect();
  } else if (x->getSup() < y->getInf() + csty && y->getSup() < x->getInf() + cstx) {
	// IC0 propagatation (increase global lower bound)
	deconnect();
	if (ToulBar2::verbose >= 2) cout << "lower bound increased " << wcsp->getLb() << " -> " << wcsp->getLb()+penalty << endl;
	if ((!x->enumerated() || (x->getInf()!=deltaValueXinf && x->getSup()!=deltaValueXsup)) &&
		(!y->enumerated() || (y->getInf()!=deltaValueYinf && y->getSup()!=deltaValueYsup))) {
	  wcsp->increaseLb(wcsp->getLb() + penalty);
	} else if (x->enumerated() && (x->getInf()==deltaValueXinf || x->getSup()==deltaValueXsup)) {
	  EnumeratedVariable *var = (EnumeratedVariable *) x;
	  var->extendAll(-penalty);
	  if (var->getInf()==deltaValueXinf) var->extend(deltaValueXinf,penalty);
	  if (var->getSup()==deltaValueXsup) var->extend(deltaValueXsup,penalty);
	  var->findSupport();
	} else if (y->enumerated() && (y->getInf()==deltaValueYinf || y->getSup()==deltaValueYsup)) {
	  EnumeratedVariable *var = (EnumeratedVariable *) y;
	  var->extendAll(-penalty);
	  if (var->getInf()==deltaValueYinf) var->extend(deltaValueYinf,penalty);
	  if (var->getSup()==deltaValueYsup) var->extend(deltaValueYsup,penalty);
	  var->findSupport();
	}
  } else {
	// propagate hard constraint
	if (CUT(wcsp->getLb()+penalty, wcsp->getUb())) {
	  if (x->getSup() < y->getInf() + csty) {
		Value newInf = x->getInf() + cstx;
		if (y->getInf() < newInf) y->increase(newInf);
		Value newSup = y->getSup() - cstx;
		if (x->getSup() > newSup) x->decrease(newSup);
	  } else if (y->getSup() < x->getInf() + cstx) {
		Value newInf = y->getInf() + csty;
		if (x->getInf() < newInf) x->increase(newInf);
		Value newSup = x->getSup() - csty;
		if (y->getSup() > newSup) y->decrease(newSup);			
	  }
	}
	
	// BAC* propagation (increase unary costs of domain bounds)
	Value xinf = x->getInf();
	if (xinf != deltaValueXinf && xinf != deltaValueXsup && xinf < y->getInf() + csty && xinf > y->getSup() - cstx) {
	  if (!y->enumerated() || (y->getInf()!=deltaValueYinf && y->getSup()!=deltaValueYsup)) {
		deltaValueXinf = xinf;
		x->projectInfCost(penalty);
	  }
	}
	Value yinf = y->getInf();
	if (yinf != deltaValueYinf && yinf != deltaValueYsup && yinf < x->getInf() + cstx && yinf > x->getSup() - csty) {
	  if (!x->enumerated() || (x->getInf()!=deltaValueXinf && x->getSup()!=deltaValueXsup)) {
		deltaValueYinf = yinf;
		y->projectInfCost(penalty);
	  }
	}
	Value xsup = x->getSup();
	if (xsup != deltaValueXsup && xsup != deltaValueXinf && xsup < y->getInf() + csty && xsup > y->getSup() - cstx) {
	  if (!y->enumerated() || (y->getInf()!=deltaValueYinf && y->getSup()!=deltaValueYsup)) {
		deltaValueXsup = xsup;
		x->projectSupCost(penalty);
	  }
	}
	Value ysup = y->getSup();
	if (ysup != deltaValueYsup && ysup != deltaValueYinf && ysup < x->getInf() + cstx && ysup > x->getSup() - csty) {
	  if (!x->enumerated() || (x->getInf()!=deltaValueXinf && x->getSup()!=deltaValueXsup)) {
		deltaValueYsup = ysup;
		y->projectSupCost(penalty);
	  }
	}
  }
}

bool Disjunction::verify()
{
  bool support = (x->getSup() >= y->getInf() + csty || y->getSup() >= x->getInf() + cstx || deconnected());
  Value xinf = x->getInf();
  Value yinf = y->getInf();
  Value xsup = x->getSup();
  Value ysup = y->getSup();
  support = support && (xinf >= yinf + csty || xinf <= ysup - cstx || deltaValueXinf == xinf || (y->enumerated() && (deltaValueYinf == yinf || deltaValueYsup == ysup)));
  support = support && (yinf >= xinf + cstx || yinf <= xsup - csty || deltaValueYinf == yinf || (x->enumerated() && (deltaValueXinf == xinf || deltaValueXsup == xsup)));
  support = support && (xsup >= yinf + csty || xsup <= ysup - cstx || deltaValueXsup == xsup || (y->enumerated() && (deltaValueYinf == yinf || deltaValueYsup == ysup)));
  support = support && (ysup >= xinf + cstx || ysup <= xsup - csty || deltaValueYsup == ysup || (x->enumerated() && (deltaValueXinf == xinf || deltaValueXsup == xsup)));
  if (!support) {
	print(cout);
	x->print(cout);
	cout << endl;
	y->print(cout);
	cout << endl;
	cout << " dxinf=" << deltaValueXinf << " dxsup=" << deltaValueXsup << " dyinf=" << deltaValueYinf << " dysup=" << deltaValueYsup << endl;
  }
  return support;
}
