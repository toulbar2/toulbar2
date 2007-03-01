/*
 * ****** Variable with domain represented by an interval or an enumerated domain *******
 */
 
#include "tb2system.hpp"
#include "tb2variable.hpp"
#include "tb2wcsp.hpp"
#include "tb2binconstr.hpp"
#include "tb2ternaryconstr.hpp"

/*
 * Constructors and misc.
 * 
 */

Variable::Variable(WCSP *w, string n, Value iinf, Value isup) : 
        WCSPLink(w,w->numberOfVariables()), name(n),
        inf(iinf, &w->getStore()->storeValue), sup(isup, &w->getStore()->storeValue), 
        constrs(&w->getStore()->storeConstraint), deltaCost(0, &w->getStore()->storeCost),
        maxCost(0, &w->getStore()->storeCost), maxCostValue(iinf, &w->getStore()->storeValue), 
        NCBucket(-1, &w->getStore()->storeValue)
//        elimOrder(-1, &w->getStore()->storeValue)
{
    if (w->getStore()->getDepth() > 0) {
        cerr << "You cannot create a variable during the search!" << endl;
        exit(EXIT_FAILURE);
    }
    w->link(this);

    linkNCBucket.content = this;
    linkNCQueue.content.var = this;
    linkNCQueue.content.timeStamp = -1;
    linkIncDecQueue.content.var = this;
    linkIncDecQueue.content.timeStamp = -1;
    linkIncDecQueue.content.incdec = NOTHING_EVENT;
    linkEliminateQueue.content.var = this;
    linkEliminateQueue.content.timeStamp = -1;
}

DLink<ConstraintLink> *Variable::link(Constraint *c, int index)
{
    ConstraintLink e;
    e.constr = c;
    e.scopeIndex = index;
    DLink<ConstraintLink> *elt = new DLink<ConstraintLink>;
    elt->content = e;
    constrs.push_back(elt,true);
    return elt;
}

int cmpConstraint(const void *p1, const void *p2)
{
    DLink<ConstraintLink> *c1 = *((DLink<ConstraintLink> **) p1);
    DLink<ConstraintLink> *c2 = *((DLink<ConstraintLink> **) p2);
    int v1 = c1->content.constr->getSmallestVarIndexInScope(c1->content.scopeIndex);
    int v2 = c2->content.constr->getSmallestVarIndexInScope(c2->content.scopeIndex);
    if (v1 < v2) return -1;
    else if (v1 > v2) return 1;
    else return 0;
}

void Variable::sortConstraints()
{
    int size = constrs.getSize();
    DLink<ConstraintLink> *sorted[size]; // replace size by MAX_DOMAIN_SIZE in case of compilation problem
    int i=0;
    for (ConstraintList::iterator iter = constrs.begin(); iter != constrs.end(); ++iter) {
        sorted[i++] = iter.getElt();
    }
    qsort(sorted, size, sizeof(DLink<ConstraintLink> *), cmpConstraint);
    for (int i = 0; i < size; i++) {
        constrs.erase(sorted[i],true);
        constrs.push_back(sorted[i],true);
    }
}

void Variable::deconnect(DLink<ConstraintLink> *link)
{
	//cout << "deconnect de variable: " << *this << endl;	
    if (!link->removed) {
        getConstrs()->erase(link, true);
        if (ToulBar2::elimVarWithSmallDegree_ && getDegree() <= 3) queueEliminate();
    }
}

int Variable::getRealDegree()
{
    TSCOPE neighbors;
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        (*iter).constr->scopeUnion( neighbors, (*iter).constr );
    }
    return neighbors.size()-1;
}

Long Variable::getWeightedDegree()
{
    Long res = 0;
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        res += (*iter).constr->getConflictWeight();
    }
    return res;
}

void Variable::conflict()
{
    wcsp->conflict();
}

/*
 * Propagation methods
 * 
 */

void Variable::queueNC()
{
    wcsp->queueNC(&linkNCQueue);
}

void Variable::queueInc()
{
    wcsp->queueInc(&linkIncDecQueue);
}

void Variable::queueDec()
{
    wcsp->queueDec(&linkIncDecQueue);
}

void Variable::queueEliminate()
{
	wcsp->queueEliminate(&linkEliminateQueue);
}

void Variable::changeNCBucket(int newBucket)
{
    if (NCBucket != newBucket) {
        if (ToulBar2::verbose >= 3) cout << "changeNCbucket " << getName() << ": " << NCBucket << " -> " << newBucket << endl;
        wcsp->changeNCBucket(NCBucket, newBucket, &linkNCBucket);
        NCBucket = newBucket;
    }
}

void Variable::setMaxUnaryCost(Value a, Cost cost)
{
    assert(canbe(a));
    maxCostValue = a;
    assert(cost >= 0);
    if (maxCost != cost) {
        maxCost = cost;
        int newbucket = min(cost2log2(cost), wcsp->getNCBucketSize() - 1);
        changeNCBucket(newbucket);
    }
}

void Variable::extendAll(Cost cost)
{
    assert(cost > 0);
    deltaCost += cost;          // Warning! Possible overflow???
    queueNC();
}

void Variable::propagateIncDec(int incdec)
{
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        if (incdec & INCREASE_EVENT) (*iter).constr->increase((*iter).scopeIndex);
        if (incdec & DECREASE_EVENT) (*iter).constr->decrease((*iter).scopeIndex);
    }
}

// Looks for the constraint that links this variable with x
BinaryConstraint* Variable::getConstr( Variable* x )
{
	BinaryConstraint* ctr2;
	TernaryConstraint* ctr3;
	
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
    	if ((*iter).constr->arity() == 2) {
    		ctr2 = (BinaryConstraint*) (*iter).constr;
  			if(ctr2->getIndex(x) >= 0) return ctr2;
    	}
    	else if ((*iter).constr->arity() == 3) {
    		ctr3 = (TernaryConstraint*) (*iter).constr;
  			int idx = ctr3->getIndex(x);
  			if(idx >= 0) {
	  			int idt = (*iter).scopeIndex;
  				if((0 != idx) && (0 != idt)) return ctr3->yz;
  				else if((1 != idx) && (1 != idt)) return ctr3->xz;
  				else return ctr3->xy;
  			}
    	}
    }
    return NULL;
}     

// Looks for the ternary constraint that links this variable with x and y
TernaryConstraint* Variable::getConstr( Variable* x, Variable* y )
{
	TernaryConstraint* ctr;
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
    	if ((*iter).constr->arity() == 3) {
    		ctr = (TernaryConstraint*) (*iter).constr;
    		if((ctr->getIndex(x)  >= 0) && (ctr->getIndex(y)  >= 0)) return ctr;    		
    	}
    }
	return NULL;
}


// returns a ternary constraint if the current variable is linked to one
TernaryConstraint* Variable::existTernary()
{
	if(!wcsp->isternary) return NULL;

	TernaryConstraint* ctr;
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
    	if ((*iter).constr->arity() == 3) {
    		ctr = (TernaryConstraint*) (*iter).constr;
			return ctr;    		
    	}
    }
	return NULL;
}



double Variable::strongLinkedby( Variable* &strvar,  TernaryConstraint* &tctr1max, TernaryConstraint* &tctr2max  ) {
	double maxtight = -1;
	strvar = NULL; tctr1max = NULL; tctr2max = NULL; 

	TernaryConstraint *tctr1 = NULL; 
	
 	for(ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
 	   if((*iter).constr->arity() == 2) {
	 	  BinaryConstraint* bctr = (BinaryConstraint*) (*iter).constr;
	 	  double bintight = bctr->getTightness();
	 	  if(bintight > maxtight) { maxtight = bintight; strvar = wcsp->getVar(bctr->getSmallestVarIndexInScope((*iter).scopeIndex)); tctr1max = NULL; tctr2max = NULL; }
 	   }
 	   else if((*iter).constr->arity() == 3) {
	   	  double terntight;
	 	  tctr1 = (TernaryConstraint*) (*iter).constr;
	 	  terntight = tctr1->getTightness() + 
	 	  			  tctr1->xy->getTightness() +
	 	  			  tctr1->xz->getTightness() +
	 	  			  tctr1->yz->getTightness();
	 	  		
	 	  Variable *x1 = NULL, *x2 = NULL;			
	 	  switch((*iter).scopeIndex) {
	 	  	case 0: x1 = tctr1->getVar(1); x2 = tctr1->getVar(2); break;  
	 	  	case 1: x1 = tctr1->getVar(0); x2 = tctr1->getVar(2); break;  
	 	  	case 2: x1 = tctr1->getVar(0); x2 = tctr1->getVar(1); break; 
	 	  	default:; 
	 	  }

		  if(terntight > maxtight) { maxtight = terntight; strvar = x1; tctr1max = tctr1; tctr1max = NULL; }   	  
	 	  			  
		  for(ConstraintList::iterator iter2=iter; iter2 != constrs.end(); ++iter2) {
		     if((*iter2).constr->arity() == 3) {
				TernaryConstraint* tctr2 = (TernaryConstraint*) (*iter2).constr;
				Variable* commonvar = NULL;
				if(tctr2->getIndex(x1) >= 0) commonvar = x1;
				else if(tctr2->getIndex(x2) >= 0) commonvar = x2;
								
				if(commonvar) {
			 	  terntight += tctr2->getTightness() + 
			 	  			   tctr2->xy->getTightness() +
			 	  			   tctr2->xz->getTightness() +
			 	  			   tctr2->yz->getTightness();
			 	  			  
			 	  if(tctr1->xy->getIndex(commonvar) >= 0) terntight -= tctr1->xy->getTightness();
			 	  else if(tctr1->xz->getIndex(commonvar) >= 0) terntight -= tctr1->xz->getTightness();
			 	  else if(tctr1->yz->getIndex(commonvar) >= 0) terntight -= tctr1->yz->getTightness();

				  if(terntight > maxtight) { maxtight = terntight; strvar = commonvar; tctr1max = tctr1; tctr2max = tctr2; }   	  
				}
		     }
		  }		  
 	   }
 	}
 	
 	return maxtight;
}






ostream& operator<<(ostream& os, Variable &var) {
    os << var.name;
    var.print(os);
    if (ToulBar2::verbose >= 3) {
        for (ConstraintList::iterator iter=var.constrs.begin(); iter != var.constrs.end(); ++iter) {
            os << " (" << (*iter).constr << "," << (*iter).scopeIndex << ")";
        }
    }
    return os;
}
