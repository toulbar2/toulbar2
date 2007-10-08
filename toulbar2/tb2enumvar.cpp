
/*
 * ****** Variable with domain represented by an enumerated domain *******
 */
 
#include "tb2enumvar.hpp"
#include "tb2wcsp.hpp"
#include "tb2binconstr.hpp"
#include "tb2ternaryconstr.hpp"


/*
 * Constructors and misc.
 * 
 */


EnumeratedVariable::EnumeratedVariable(WCSP *w, string n, Value iinf, Value isup) : 
        Variable(w, n, iinf, isup), 
        domain(iinf, isup, &w->getStore()->storeDomain), deltaCost(MIN_COST, &w->getStore()->storeCost),
        support(iinf, &w->getStore()->storeValue)
{
    init();
}

EnumeratedVariable::EnumeratedVariable(WCSP *w, string n, Value *d, int dsize) : 
        Variable(w, n, min(d,dsize), max(d, dsize)), 
        domain(d, dsize, &w->getStore()->storeDomain), deltaCost(MIN_COST, &w->getStore()->storeCost),
        support(min(d,dsize), &w->getStore()->storeValue)
{
    init();
}

void EnumeratedVariable::init()
{
    if (wcsp->getStore()->getDepth() > 0) {
        cerr << "You cannot create a variable during the search!" << endl;
        exit(EXIT_FAILURE);
    }

    costs = vector<StoreCost>(getDomainInitSize(), StoreCost(MIN_COST, &wcsp->getStore()->storeCost));
    linkACQueue.content.var = this;
    linkACQueue.content.timeStamp = -1;
    linkDACQueue.content.var = this;
    linkDACQueue.content.timeStamp = -1;
    linkEAC1Queue.content.var = this;
    linkEAC1Queue.content.timeStamp = -1;
    linkEAC2Queue.content.var = this;
    linkEAC2Queue.content.timeStamp = -1;
}

void EnumeratedVariable::getDomain(Value *array)
{
    for (iterator iter = begin(); iter != end(); ++iter) {
    	*array = *iter;
    	++array;
    }
}

void EnumeratedVariable::getDomainAndCost(ValueCost *array)
{
    for (iterator iter = begin(); iter != end(); ++iter) {
    	array->value = *iter;
    	array->cost = getCost(*iter);
    	++array;
    }
}

Cost EnumeratedVariable::getBinaryCost(ConstraintLink c, Value myvalue, Value itsvalue) {
    return (c.scopeIndex == 0)?((BinaryConstraint *) c.constr)->getCost(myvalue, itsvalue):((BinaryConstraint *) c.constr)->getCost(itsvalue, myvalue);
}

Cost EnumeratedVariable::getBinaryCost(BinaryConstraint* c, Value myvalue, Value itsvalue) {
    return (c->getIndex(this) == 0)?c->getCost(myvalue, itsvalue):c->getCost(itsvalue, myvalue);
}


void EnumeratedVariable::print(ostream& os)
{
    if (unassigned()) {
        os << " " << domain;
    } else {
        os << " [" << inf << "," << sup << "]";
    }
    os << "/" << getDegree();
//    os << "/" << getWeightedDegree();
    if (unassigned()) {
        os << " <";
        for (iterator iter=begin(); iter != end(); ++iter) {
            os << " " << getCost(*iter);
        }
        os << " > s:" << support;
    }
}

/*
 * Propagation methods
 * 
 */

void EnumeratedVariable::queueAC()
{
    wcsp->queueAC(&linkACQueue);
}

void EnumeratedVariable::queueDAC()
{
    wcsp->queueDAC(&linkDACQueue);
}

void EnumeratedVariable::queueEAC1()
{
    wcsp->queueEAC1(&linkEAC1Queue);
}

void EnumeratedVariable::queueEAC2()
{
    wcsp->queueEAC2(&linkEAC2Queue);
}

void EnumeratedVariable::project(Value value, Cost cost)
{
    assert(cost >= MIN_COST);
    Cost oldcost = getCost(value);
    costs[toIndex(value)] += cost;
    Cost newcost = oldcost + cost;
    if (value == maxCostValue || LUBTEST(maxCost, newcost)) queueNC();
    if (DACTEST(oldcost, cost)) {
//        cout << "insert in DAC and EAC1 of " << getName() << endl;
        queueDAC();
        queueEAC1();
    }
    if (CUT(newcost + wcsp->getLb(), wcsp->getUb())) removeFast(value);     // Avoid any unary cost overflow
}

void EnumeratedVariable::projectInfCost(Cost cost)
{
    assert(cost >= MIN_COST);
    Value value = getInf();
    Cost oldcost = getCost(value);
    project(value, cost);
    if (support == value || SUPPORTTEST(oldcost, cost)) findSupport();
}

void EnumeratedVariable::projectSupCost(Cost cost)
{
    assert(cost >= MIN_COST);
    Value value = getSup();
    Cost oldcost = getCost(value);
    project(value, cost);
    if (support == value || SUPPORTTEST(oldcost, cost)) findSupport();
}

void EnumeratedVariable::extend(Value value, Cost cost)
{
    assert(cost >= MIN_COST);
	//    assert(CUT(costs[toIndex(value)], cost));
    costs[toIndex(value)] -= cost;
    assert( ToulBar2::verbose < 4 || ((cout << "extend " << getName() << " (" << value << ") -= " << cost << endl), true) );
    if (value == maxCostValue || PARTIALORDER) queueNC();
}

void EnumeratedVariable::extendAll(Cost cost)
{
  //    assert(cost > MIN_COST);
    deltaCost += cost;          // Warning! Possible overflow???
    queueNC();
}

void EnumeratedVariable::findSupport()
{
    if (cannotbe(support) || getCost(support) > MIN_COST) {
        Value newSupport = getInf();
        Cost minCost = getCost(newSupport);
        iterator iter = begin();
        for (++iter; minCost > MIN_COST && iter != end(); ++iter) {
            Cost cost = getCost(*iter);
            if (GLB(&minCost, cost)) {
                newSupport = *iter;
            }
        }
        if (minCost > MIN_COST) {
            extendAll(minCost);
            if (ToulBar2::verbose >= 2) cout << "lower bound increased " << wcsp->getLb() << " -> " << wcsp->getLb()+minCost << endl;
            wcsp->increaseLb(wcsp->getLb() + minCost);
        }
        assert(canbe(newSupport) && (getCost(newSupport) == MIN_COST || SUPPORTTEST(getCost(newSupport))));
        support = newSupport;
    }
}

void EnumeratedVariable::propagateNC()
{
    if (ToulBar2::verbose >= 3) cout << "propagateNC for " << getName() << endl;
    Value maxcostvalue = getSup()+1;
    Cost maxcost = MIN_COST;
    bool supportBroken = false;
    // Warning! the first value must be visited because it may be removed
    for (iterator iter = begin(); iter != end(); ++iter) {
        Cost cost = getCost(*iter);
        if (CUT(cost + wcsp->getLb(), wcsp->getUb())) {
            if (SUPPORTTEST(cost)) supportBroken = true;
            removeFast(*iter);
        } else if (LUB(&maxcost, cost) || cannotbe(maxcostvalue)) {
            maxcostvalue = *iter;
        }
    }
    assert(getCost(maxcostvalue) == maxcost || !LUBTEST(maxcost, getCost(maxcostvalue)));
    setMaxUnaryCost(maxcostvalue, maxcost);
    if (supportBroken) findSupport();
}

bool EnumeratedVariable::verifyNC()
{
    bool supported = true;
    Cost minCost = MIN_COST;
    for (iterator iter = begin(); iter != end(); ++iter) {
        Cost cost = getCost(*iter);
        if (CUT(cost + wcsp->getLb(), wcsp->getUb())) {
            cout << *this << " not NC!" << endl;
            return false;
        }
        GLB(&minCost, cost);
    }
    if (minCost > MIN_COST) {
        cout << *this << " not NC*!" << endl;
        supported = false;
    }
    if (cannotbe(support) || (getCost(support)>MIN_COST && !SUPPORTTEST(getCost(support)) )) {
        cout << *this << " has an unvalid NC support!" << endl;
        supported = false;
    }
    return supported;
}

void EnumeratedVariable::propagateAC()
{
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        (*iter).constr->remove((*iter).scopeIndex);
    }
}

void EnumeratedVariable::propagateDAC()
{
    for (ConstraintList::iterator iter=constrs.rbegin(); iter != constrs.rend(); --iter) {
        (*iter).constr->projectFromZero((*iter).scopeIndex);
    }
}

void EnumeratedVariable::fillEAC2(bool self)
{
  if (self) queueEAC2();
  for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
    (*iter).constr->fillEAC2((*iter).scopeIndex);
  }
}

bool EnumeratedVariable::isEAC(Value a)
{
    if (getCost(a)==MIN_COST) {
        for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
            if (!(*iter).constr->isEAC((*iter).scopeIndex, a)) {
#ifndef NDEBUG
                if (ToulBar2::verbose >=4) {
                    cout << getName() << "(" << a << ") is not EAC due to constraint " << *(*iter).constr << endl; 
                    if ((*iter).constr->arity()==3) {
                        TernaryConstraint *c = (TernaryConstraint *) (*iter).constr;
                        if (c->xy->connected()) cout << *c->xy;
                        if (c->xz->connected()) cout << *c->xz;
                        if (c->yz->connected()) cout << *c->yz;
                    }
                }
#endif
                return false;
            }
        }
        support = a;
#ifndef NDEBUG
        if (ToulBar2::verbose >=4) cout << getName() << "(" << a << ") is EAC!" << endl;
#endif
        return true;
    }
#ifndef NDEBUG
    if (ToulBar2::verbose >=4) cout << getName() << "(" << a << ") is not EAC due to unary cost " << getCost(a) << endl; 
#endif
    return false;
}

bool EnumeratedVariable::isEAC()
{
    assert(canbe(support));
    if (isEAC(support)) return true;
    else {
        for (iterator iter = begin(); iter != end(); ++iter) {
            if (*iter != support && isEAC(*iter)) return true;
        }
    }
    return false;
}

void EnumeratedVariable::propagateEAC()
{
    if (!isEAC()) {
#ifndef NDEBUG
        Cost beforeLb = wcsp->getLb();
#endif
        for (ConstraintList::iterator iter = constrs.begin(); iter != constrs.end(); ++iter) {
        	(*iter).constr->findFullSupportEAC((*iter).scopeIndex);
        }
        fillEAC2(false);
        if (unassigned()) {
            // findFullSupportEAC may have inserted current variable in EAC1
	        if (!linkEAC1Queue.removed) {
                assert(((BTList<VariableWithTimeStamp> *) wcsp->getQueueEAC1())->inBTList(&linkEAC1Queue));
                wcsp->getQueueEAC1()->remove(&linkEAC1Queue);
            }
#ifndef NDEBUG
            // check if lb has effectively been increased
            if (wcsp->getLb() == beforeLb)
                if (ToulBar2::verbose) cout << "EAC failed on " << getName() << endl;
#endif
        }
    }
}

void EnumeratedVariable::increaseFast(Value newInf)
{
    if (ToulBar2::verbose >= 2) cout << "increase " << getName() << " " << inf << " -> " << newInf << endl;
    if (newInf > inf) {
        if (newInf > sup) THROWCONTRADICTION;
        else {
            newInf = domain.increase(newInf);
            if (newInf == sup) assign(newInf);
            else {
                inf = newInf;
                queueInc();
                if (PARTIALORDER) queueDAC();
                if (ToulBar2::setmin) (*ToulBar2::setmin)(wcsp->getIndex(), wcspIndex, newInf);
            }
        }
      }
}

void EnumeratedVariable::increase(Value newInf)
{
    if (ToulBar2::verbose >= 2) cout << "increase " << getName() << " " << inf << " -> " << newInf << endl;
    if (newInf > inf) {
        if (newInf > sup) THROWCONTRADICTION;
        else {
            newInf = domain.increase(newInf);
            if (newInf == sup) assign(newInf);
            else {
                inf = newInf;
                if (newInf > maxCostValue || PARTIALORDER) queueNC();           // diff with increaseFast
                if (newInf > support || PARTIALORDER) findSupport();            // diff with increaseFast
                queueDAC();                                     // diff with increaseFast
                queueEAC1();                                     // diff with increaseFast
                queueInc();
                if (ToulBar2::setmin) (*ToulBar2::setmin)(wcsp->getIndex(), wcspIndex, newInf);
            }
        }
      }
}

void EnumeratedVariable::decreaseFast(Value newSup)
{
    if (ToulBar2::verbose >= 2) cout << "decrease " << getName() << " " << sup << " -> " << newSup << endl;
    if (newSup < sup) {
        if (newSup < inf) THROWCONTRADICTION;
        else {
            newSup = domain.decrease(newSup);
            if (inf == newSup) assign(newSup);
            else {
                sup = newSup;
                queueDec();
                if (PARTIALORDER) queueDAC();
                if (ToulBar2::setmax) (*ToulBar2::setmax)(wcsp->getIndex(), wcspIndex, newSup);
            }
        }
      }
}

void EnumeratedVariable::decrease(Value newSup)
{
    if (ToulBar2::verbose >= 2) cout << "decrease " << getName() << " " << sup << " -> " << newSup << endl;
    if (newSup < sup) {
        if (newSup < inf) THROWCONTRADICTION;
        else {
            newSup = domain.decrease(newSup);
            if (inf == newSup) assign(newSup);
            else {
                sup = newSup;
                if (newSup < maxCostValue || PARTIALORDER) queueNC();           // diff with decreaseFast
                if (newSup < support || PARTIALORDER) findSupport();            // diff with decreaseFast
                queueDAC();                                     // diff with decreaseFast
                queueEAC1();                                     // diff with decreaseFast
                queueDec();
                if (ToulBar2::setmax) (*ToulBar2::setmax)(wcsp->getIndex(), wcspIndex, newSup);
            }
        }
      }
}

void EnumeratedVariable::removeFast(Value value)
{
    if (ToulBar2::verbose >= 2) cout << "remove " << *this << " <> " << value << endl;
    if (value == inf) increaseFast(value + 1);
    else if (value == sup) decreaseFast(value - 1);
    else if (canbe(value)) {
        domain.erase(value);
        queueAC();
        if (PARTIALORDER) queueDAC();
        if (ToulBar2::removevalue) (*ToulBar2::removevalue)(wcsp->getIndex(), wcspIndex, value);
    }
}

void EnumeratedVariable::remove(Value value)
{
    if (ToulBar2::verbose >= 2) cout << "remove " << *this << " <> " << value << endl;
    if (value == inf) increase(value + 1);
    else if (value == sup) decrease(value - 1);
    else if (canbe(value)) {
        domain.erase(value);
        if (value == maxCostValue || PARTIALORDER) queueNC();
        if (value == support || PARTIALORDER) findSupport();
        queueDAC();
        queueEAC1();
        queueAC();
        if (ToulBar2::removevalue) (*ToulBar2::removevalue)(wcsp->getIndex(), wcspIndex, value);
    }
}

// this function is used ONLY for restoring the solution when
// variable elimination is tuned on
void EnumeratedVariable::assignWhenEliminated(Value newValue)
{
        assert(NCBucket == -1);
        inf = newValue;
        sup = newValue;
        support = newValue;
        maxCostValue = newValue;
        maxCost = MIN_COST;
}

void EnumeratedVariable::assign(Value newValue)
{
    if (ToulBar2::verbose >= 2) cout << "assign " << *this << " -> " << newValue << endl;
    if (unassigned() || getValue() != newValue) {
        if (cannotbe(newValue)) THROWCONTRADICTION;
        changeNCBucket(-1);
        inf = newValue;
        sup = newValue;
        support = newValue;
        maxCostValue = newValue;
        maxCost = MIN_COST;
   
        Cost cost = getCost(newValue);
        if (cost > MIN_COST) {
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


// eliminates the current (this) variable that participates
// in a single binary constraint ctr
void EnumeratedVariable::elimVar( BinaryConstraint* ctr )
{
       assert(getDegree() == 1);
       // deconnect first to be sure the current var is not involved in future propagation	
       ctr->deconnect();

       EnumeratedVariable *x = (EnumeratedVariable *) ctr->getVarDiffFrom(this);

       // to be done before propagation
       WCSP::elimInfo ei = {this, x, NULL, ctr, NULL, NULL};
       wcsp->elimInfos[wcsp->getElimOrder()] = ei;
       wcsp->elimination();

       bool supportBroken = false;
       for (iterator iter1 = x->begin(); iter1 != x->end(); ++iter1) {
           Cost mincost = MAX_COST;
           for (iterator iter = begin(); iter != end(); ++iter) {
               Cost curcost = getCost(*iter) + getBinaryCost(ctr, *iter, *iter1);
               if (curcost < mincost) mincost = curcost;
           }
           if (mincost > MIN_COST) {
		       if (x->getSupport() == *iter1) supportBroken = true;
               x->project(*iter1, mincost);
           }
       }
      if (supportBroken) {  x->findSupport();  }
}


// eliminates the current (this) variable that participates
// in two binary constraints (its links ara xylink and xzlink)
void EnumeratedVariable::elimVar( ConstraintLink  xylink,  ConstraintLink xzlink )
{
  	 assert(getDegree() == 2);
     xylink.constr->deconnect(); 	            
	 xzlink.constr->deconnect();

	 EnumeratedVariable *y = (EnumeratedVariable *) wcsp->getVar(xylink.constr->getSmallestVarIndexInScope(xylink.scopeIndex));
	 EnumeratedVariable *z = (EnumeratedVariable *) wcsp->getVar(xzlink.constr->getSmallestVarIndexInScope(xzlink.scopeIndex));
	    
	 BinaryConstraint* yz = y->getConstr(z);
     BinaryConstraint* yznew = wcsp->newBinaryConstr(y,z); 

	 for (iterator itery = y->begin(); itery != y->end(); ++itery) {
	 for (iterator iterz = z->begin(); iterz != z->end(); ++iterz) {
	    Cost mincost = MAX_COST;

	    for (iterator iter = begin(); iter != end(); ++iter) {
	        Cost curcost = getCost(*iter) +
						   getBinaryCost(xylink, *iter, *itery) +
						   getBinaryCost(xzlink, *iter, *iterz);
						   
	        if (curcost < mincost) mincost = curcost;
	    }
		yznew->setcost(*itery,*iterz,mincost);
	 }}
 
	 if(yz) {
		 yz->addCosts( yznew );
	 	 yz->reconnect();
	 } else { 
		 yz = yznew;
	 	 yz->reconnect();
	 }	

	 // to be done before propagation
	 WCSP::elimInfo ei = {this, y,z, (BinaryConstraint*) xylink.constr, (BinaryConstraint*) xzlink.constr, NULL};
	 wcsp->elimInfos[wcsp->getElimOrder()] = ei;
	 wcsp->elimination();
     yz->propagate(); 
}

// eliminates the current (this) variable that participates
// in the ternary constraint 'xyz'
// the function can fail to eliminate the current variable
// if it is linked to more than (in total) two variables.
// It returns true if the current variable was eliminated
bool EnumeratedVariable::elimVar( TernaryConstraint* xyz )
{
	BinaryConstraint* yz = NULL;
	if(xyz->xy->getIndex(this) < 0) yz = xyz->xy;
	else if(xyz->xz->getIndex(this) < 0) yz = xyz->xz;
	else if(xyz->yz->getIndex(this) < 0) yz = xyz->yz;	 
	assert(yz != NULL);

	int n2links = 0;
	int n3links = 0;

	ConstraintLink links[2] = {{NULL, 0},{NULL, 0}};
 	for(ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
 	   if((*iter).constr->arity() == 2) links[n2links++] =  (*iter);
 	   else if((*iter).constr->arity() == 3) n3links++;
 	   else return false;
 	}

	if(n3links > 1) return false;

	for(int i=0; i<n2links; i++) {		
		int idvar = links[i].constr->getSmallestVarIndexInScope(links[i].scopeIndex);
		if(xyz->getIndex( wcsp->getVar(idvar) ) < 0) return false; 
	}

	xyz->deconnect();
	if(n2links > 0) links[0].constr->deconnect();
	if(n2links > 1) links[1].constr->deconnect();

	EnumeratedVariable* y = (EnumeratedVariable*) yz->getVar(0);
	EnumeratedVariable* z = (EnumeratedVariable*) yz->getVar(1);
    
	bool flag_rev = false;
	if(n2links > 0) {
        if (links[0].constr->getSmallestVarIndexInScope(links[0].scopeIndex) != y->wcspIndex) {
            assert(links[0].constr->getSmallestVarIndexInScope(links[0].scopeIndex) == z->wcspIndex);
            flag_rev = true;
        }
	}
		
	for (iterator itery = y->begin(); itery != y->end(); ++itery) {
	for (iterator iterz = z->begin(); iterz != z->end(); ++iterz) {
	    Cost mincost = MAX_COST;
	    for (iterator iter = begin(); iter != end(); ++iter) {
	        Cost curcost = getCost(*iter) + xyz->getCost(this,y,z,*iter,*itery,*iterz); 

			if(!flag_rev) {
		        if(n2links > 0) { assert(links[0].constr->getIndex(y) >= 0);
		                          curcost += getBinaryCost(links[0], *iter, *itery); }
				if(n2links > 1) { assert(links[1].constr->getIndex(z) >= 0);
                    			  curcost += getBinaryCost(links[1], *iter, *iterz); }
			} else {
		        if(n2links > 0) { assert(links[0].constr->getIndex(z) >= 0);
   				                  curcost += getBinaryCost(links[0], *iter, *iterz); }
				if(n2links > 1) { assert(links[1].constr->getIndex(y) >= 0);
  								  curcost += getBinaryCost(links[1], *iter, *itery); }
			}
	        if (curcost < mincost) mincost = curcost;
	    }
		yz->addcost(*itery,*iterz,mincost);
	 }}

	if(!yz->connected()) yz->reconnect();
	
	// to be done before propagation
	WCSP::elimInfo ei = {this,y,z,(BinaryConstraint*) links[(flag_rev)?1:0].constr, (BinaryConstraint*) links[(flag_rev)?0:1].constr, xyz};
	wcsp->elimInfos[wcsp->getElimOrder()] = ei;
	wcsp->elimination();
    yz->propagate(); 
	return true;
}

void EnumeratedVariable::eliminate()
{
    if (ToulBar2::elimDegree_preprocessing_ >= 0 && 
        (getDegree() <= min(1,ToulBar2::elimDegree_preprocessing_) || 
         getTrueDegree() <= ToulBar2::elimDegree_preprocessing_)) {
	  wcsp->variableElimination(this); 
	  return;
    }
 	
    if (getDegree() > ToulBar2::elimDegree_) return;

    if(getDegree() > 0) {
		TernaryConstraint* ternCtr = existTernary();
				
		if(ternCtr) { if(!elimVar(ternCtr)) return; }
		else {
			if(getDegree() > 2) return;
	
			ConstraintLink xylink = *constrs.begin();
			ConstraintLink xzlink = {NULL,0};

			if(xylink.constr->arity() > 2) return;
			
			if(getDegree() == 2) {
				xzlink = *constrs.rbegin();
			    if(xzlink.constr->arity() > 2) return;

				elimVar(xylink,xzlink);		
			} else {
				BinaryConstraint* xy = (BinaryConstraint*) xylink.constr;	
				elimVar( xy );
			}		
		}
	}
	assert(getDegree() == 0);
	if (ToulBar2::verbose >= 2) cout << "Eliminate " << getName() << endl;
	assert(getCost(support) == MIN_COST); // it is ensured by previous calls to findSupport
	assign(support); // warning! dummy assigned value
}
