/*
 * ****** Enforce VAC in a WCSP.
 */
 
#include "tb2vac.hpp"
#include <list>
#include <algorithm>
#include <vector>


class tVACStat {
public:

	int   var;
	Long  sumlb;
	Long  nlb;
	
	tVACStat(int varin) { var = varin; sumlb = 0; nlb = 0; }
};

bool cmp_function( tVACStat* v1, tVACStat* v2 )
{
	return v1->sumlb > v2->sumlb;
} 




VACExtension::VACExtension (WCSP *w) : wcsp(w), VAC2(&w->getStore()->storeVariable), nbIterations(0), inconsistentVariable(-1)
 {
  queueP = new stack< pair<int, int> >;
  queueR = new stack< pair<int, int> >;
  minlambda = MAX_COST;
  sumlb = 0;
  nlb = 0;
  varAssign = -1;
}

VACExtension::~VACExtension () {
  delete queueP;
  delete queueR;
}

void VACExtension::init()
{
  VACVariable* xi;
  for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
    xi = (VACVariable *) wcsp->getVar(i);
    xi->setThreshold(0);
  }
  iniThreshold();
  
  for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
  	tVACStat* vacinfo = new tVACStat(i);
  	heapAccess[i] = vacinfo;
  	heap.insert(heap.end(), vacinfo);
  } 
}

void VACExtension::iniThreshold()
{
	list<Cost>::iterator it = (wcsp->scaleVAC.begin());
	Cost c = *it; 
	if(wcsp->getUb() < c) c = wcsp->getUb();
	itThreshold = c; 
}

Cost VACExtension::getThreshold() { return itThreshold; }


void VACExtension::nextScaleCost() {
	Cost c = MAX_COST;
	bool done = false;
	list<Cost>::iterator it = wcsp->scaleVAC.begin();
	while((it != wcsp->scaleVAC.end()) && !done) {
		c = *it;
		done = c < itThreshold;
		++it;
	}
	if(!done)  c  = itThreshold / 2;
	if(c < ToulBar2::costThreshold) c = 0;
	
	itThreshold = c;	
}


void VACExtension::reset()
{
  VACVariable* x;
  clear();	
  while (!queueP->empty()) queueP->pop();
  while (!queueR->empty()) queueR->pop();
  int bucket = cost2log2(ToulBar2::costThreshold);
  if (bucket < 0) bucket = 0;
  for (; bucket < wcsp->getNCBucketSize(); bucket++) {
     VariableList* varlist = wcsp->getNCBucket(bucket);  
     for (VariableList::iterator iter = varlist->begin(); iter != varlist->end();) {
        x = (VACVariable*) *iter;
        if (x->unassigned() && x->getMaxCost() >= x->getThreshold())  x->queueVAC();
        ++iter;         
     }
  }
  for (BTQueue::iterator it = VAC2.begin(); it != VAC2.end(); ++it) {
    x = (VACVariable*) (*it);
    x->queueVAC();
  }
  
  
}





bool VACExtension::propagate() {	
  if (wcsp->getStore()->getDepth() >= ToulBar2::vac) {
    return false;
  }
 
  if(getVarTimesStat(varAssign) > 100) {
  	long double m = getVarCostStat(varAssign)/getVarTimesStat(varAssign);
  	if(m < ToulBar2::costMultiplier/10.) { 
  		inconsistentVariable = -1;  		
  		return false; 
  	}
  }
  
  Cost ub = wcsp->getUb();
  Cost lb = wcsp->getLb();
  minlambda = ub - lb;
 
  //bool isvac = false;	
  bool isvac = true;	
  bool util = true;
  
  breakCycles = 0;	
 
  while((!util || isvac) && itThreshold) {
  	  nbIterations++;

	  reset();  
	  wcsp->getStore()->store();
	  enforcePass1();
	  wcsp->getStore()->restore();

	  isvac = isVAC();
	  
	  if(!isvac) {
		    enforcePass2();
		    if (ToulBar2::verbose > 0) cout << "VAC Lb: " << lb << "    incvar: " << inconsistentVariable << "    minlambda: " << minlambda <<  "      itThreshold: " << itThreshold  << endl;		      
		    util = enforcePass3();
	  }
	  else nextScaleCost();	  
  }

  updateStat(wcsp->getLb() - lb);

  //if(isvac) assert(checkPass1());
  return util;
}
 



bool VACExtension::enforcePass1 (VACVariable *xj, VACConstraint *cij) {
  bool wipeout = false;
  VACVariable *xi;
  xi = (VACVariable *)cij->getVarDiffFrom(xj);	
  for (EnumeratedVariable::iterator it = xi->begin(); it != xi->end(); ++it) {
	  Value v = *it;
	  if(xi->getVACCost(v)) { xi->removeVAC(v); xi->queueVAC(); }
	  else if (cij->revise(xi,v)) {
         wipeout = xi->removeVAC(v);
         xi->setKiller(v,xj->wcspIndex);
         queueP->push(pair<int, int>(xi->wcspIndex, v));
         xi->queueVAC();
         if(wipeout) {
	        inconsistentVariable = xi->wcspIndex;
	        return true;
	     }
      }
  }
  return false;
}


void VACExtension::enforcePass1 () {
  VACVariable* xi;
  VACVariable* xj;
  VACConstraint *cij;
  //if (ToulBar2::verbose > 1) cout << "VAC Enforce Pass 1" << endl; 
  while (!VAC.empty()) {
    xj = (VACVariable*) VAC.pop_first();
	    			   
    //list<Constraint*> l;
    for (ConstraintList::iterator itc = xj->getConstrs()->begin(); itc != xj->getConstrs()->end(); ++itc) {
      cij = (VACConstraint *) (*itc).constr;
  	  if(cij->arity() == 2) {
	 	  xi = (VACVariable *)cij->getVarDiffFrom(xj);	
	      //if(xj->getMaxK(nbIterations) > 2) l.push_back(cij); else 
	      if(enforcePass1(xj,cij)) return;
	  }
    }
    /*for (list<Constraint*>::iterator itl = l.begin(); itl != l.end(); ++itl) {
      cij = (VACConstraint *) *itl;
	  if(enforcePass1(xj,cij)) return;
    }*/
  }
  inconsistentVariable = -1;
}

bool VACExtension::checkPass1 () const {
  VACConstraint *cij;
  VACVariable *xi, *xj;
  bool supportFound;
  
  for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
    xi = (VACVariable *) wcsp->getVar(i);
    for (ConstraintList::iterator iter = xi->getConstrs()->begin(); iter != xj->getConstrs()->end(); ++iter) {
      cij = (VACConstraint *) (*iter).constr;
      if(cij->arity() == 2) {
	      xj = (VACVariable *) cij->getVarDiffFrom(xi);
		  for (EnumeratedVariable::iterator iti = xi->begin(); iti != xi->end(); ++iti) {
			  Value v = *iti;	
	          supportFound = false;
			  for (EnumeratedVariable::iterator itj = xj->begin(); itj != xj->end(); ++itj) {
	              Value w = *itj;	
	              if ((!xi->getVACCost(v)) && (!xj->getVACCost(w)) && 
	                 (!cij->getVACCost(xi,xj, v, w))) {
		                supportFound = true;
		                break;
	                 }
	          }
	          if (!supportFound) {
	            return false;
	          }
	       }
	     }
     }
  }
  return true;
}

bool VACExtension::isVAC () const {
  return (inconsistentVariable == -1);
}

void VACExtension::enforcePass2 () {
  int i0 = inconsistentVariable;
  int i, j;
  VACVariable *xi0, *xi, *xj;
  VACConstraint *cij;
  Cost tmplambda;
  Value v;
  
  //if (ToulBar2::verbose > 1)  cout << "VAC Enforce Pass 2" << endl;
  
  assert(i0 >= 0);
  xi0 = (VACVariable *) wcsp->getVar(i0);
  
  	
  for (EnumeratedVariable::iterator iti0 = xi0->begin(); iti0 != xi0->end(); ++iti0) {
	v = *iti0;
    xi0->addToK(v,1, nbIterations);
    xi0->setMark(v, nbIterations);
	Cost cost = xi0->getVACCost(v);
	if (cost > 0 && cost < minlambda) minlambda = cost;
  }

  while (!queueP->empty()) {
    i = queueP->top().first;
    v = queueP->top().second;
    queueP->pop();
    xi = (VACVariable *) wcsp->getVar(i);
    if (xi->isMarked(v, nbIterations)) {
      j = xi->getKiller(v);
      xj = (VACVariable *) wcsp->getVar(j);
      queueR->push(pair<int, int>(i, v));
      cij = (VACConstraint *) xi->getConstr(xj);
	  assert(cij);	
      //if (ToulBar2::verbose > 6) cout << "x" << xi->wcspIndex << "," << v << "   killer: " << xj->wcspIndex << endl;
      
	  for (EnumeratedVariable::iterator itj = xj->begin(); itj != xj->end(); ++itj) {
		Value w = *itj;
		Cost costij = cij->getVACCost(xi,xj,v, w);
		if (costij > 0) {
          Cost tmpK = xi->getK(v, nbIterations);
		  if (xj->getKiller(w) == i) tmpK += xj->getK(w, nbIterations);
		  if (NOCUT(wcsp->getLb() + costij,wcsp->getUb())) {
			if( (costij/tmpK) < minlambda) minlambda = costij/tmpK;
		  }
		} else {
          Cost tmpK = xi->getK(v, nbIterations) - cij->getK(xj,w, nbIterations);
          if (tmpK > 0) {
            xj->addToK(w,tmpK,nbIterations);
            cij->setK(xj, w, xi->getK(v, nbIterations), nbIterations);
            if(!xj->getVACCost(w)) xj->setMark(w, nbIterations);
            else if (NOCUT(wcsp->getLb() + xj->getVACCost(w),wcsp->getUb())) {
              tmplambda = xj->getVACCost(w) / xj->getK(w, nbIterations);  	
			  if(tmplambda < minlambda) minlambda = tmplambda;
            }
          }  
        }
      }
    }
  }
  //if (maxK == 0) {
  //  maxK = wcsp->getUb() - wcsp->getLb();
  //}
  if (ToulBar2::verbose > 1) cout << "minLambda: " << minlambda << "\t\t (lb = " << wcsp->getLb() << ", ub = " << wcsp->getUb() << ")" << endl;
}

bool VACExtension::enforcePass3 () {
  bool util = (minlambda >= 1);
  /*if(util) {
	  Cost ub = wcsp->getUb();
	  Cost lb = wcsp->getLb();
	  util = ( (ub - lb)/(10000*ToulBar2::costMultiplier) ) < minlambda;
  }*/
  Cost lambda = minlambda;	

  //if (ToulBar2::verbose > 2) cout << "VAC Enforce Pass 3.   minlambda " << minlambda << " , var: " << inconsistentVariable << endl;
  int i, j;
  VACVariable *xi, *xj;
  VACConstraint *cij;
  int i0 = inconsistentVariable;
  VACVariable *xi0 = (VACVariable *) wcsp->getVar(i0);
  Value w;
	  
  int minvar = -1;
  int minval = -1;
  Cost minc = MAX_VAL;
  Cost maxk = 0;
    	  
  if(!util) {
	    while (!queueR->empty()) { 
		    j = queueR->top().first;
		    w = queueR->top().second; 
		    xj = (VACVariable *) wcsp->getVar(j);
			i = xj->getKiller(w);
			xi = (VACVariable *) wcsp->getVar(i);
			cij = (VACConstraint *) xi->getConstr(xj);
		    queueR->pop();
		    for (EnumeratedVariable::iterator iti = xi->begin(); iti != xi->end(); ++iti) {
				Value v = *iti;
			    if(xi->getK(v, nbIterations)) { 
			    	if(xi->getVACCost(v)) {
				    	if(minc > xi->getCost(v)) {
				    		minc = xi->getCost(v);
				    		minvar = i;
				    		minval = v;
				    	}
				    	if(maxk < xi->getK(v, nbIterations)) maxk = xi->getK(v, nbIterations);
				    } 
			    }
		    }
	    }
	    assert(minvar >= 0);
	 	xi = (VACVariable *) wcsp->getVar(minvar);	
		xi->setThreshold( minc + 1 );
	    breakCycles++;
	    //if (ToulBar2::verbose > 1) cout << "BreakCycle   (var: " << minvar << ", val= " << minval << ")   thr: " <<  xi->getThreshold() << endl;    
		if(breakCycles > 5) { inconsistentVariable = -1; itThreshold = 0; }
	    return false;
  }

  while (!queueR->empty()) {
    j = queueR->top().first;
    w = queueR->top().second;
    queueR->pop();
    xj = (VACVariable *) wcsp->getVar(j);
    i = xj->getKiller(w);
    xi = (VACVariable *) wcsp->getVar(i);
    cij = (VACConstraint *) xi->getConstr(xj);
    assert(cij);

    for (EnumeratedVariable::iterator iti = xi->begin(); iti != xi->end(); ++iti) {
	  Value v = *iti;	
      if ((NOCUT(wcsp->getLb() + cij->getCost(xi,xj,v,w), wcsp->getUb())) && 
      	  (cij->getCost(xi,xj, v, w) < lambda * xj->getK(w, nbIterations))) {
	          Cost ecost = lambda * xj->getK(w, nbIterations) - cij->getCost(xi,xj,v, w);
		      cij->VACextend(xi, v, ecost);
      }
    }
    cij->VACproject(xj, w, lambda * xj->getK(w, nbIterations));
  }
  xi0->extendAll(lambda);
  wcsp->increaseLb(wcsp->getLb() + lambda);
  return true;
}

void VACExtension::updateStat(Cost lambda)
{
  sumlb += lambda;
  nlb++;
  //tVACStat* v = heapAccess[inconsistentVariable];	

  if(varAssign >= 0) {
	  tVACStat* v = heapAccess[varAssign];	
	  v->sumlb += lambda;
	  v->nlb++;
  }
}

Cost VACExtension::getVarCostStat( int i )
{
  tVACStat* v = heap[i];
  return v->sumlb;	
}

Long VACExtension::getVarTimesStat( int i )
{
  if(i < 0) return 0; 
  tVACStat* v = heap[i];
  return v->nlb;	
}

void VACExtension::afterPreprocessing() 
{
	int discarded = 0;
	for (unsigned int i=0; i<wcsp->numberOfConstraints(); i++) {
		Constraint* c = wcsp->getCtr(i);
        if (c->connected()) c->computeTightness();
        if(c->getTightness() < ToulBar2::relaxThreshold) {
        	c->deconnect();
        	discarded++;
        }
    }
    if(discarded) cout << "WARNING num of discarded ctrs: " << discarded << endl;
}





bool VACExtension::isNull (Cost c) const {
  return (c < itThreshold);
}

void VACExtension::clear () {
  while (!VAC.empty()) VAC.pop();
}


void VACExtension::queueVAC(DLink<VariableWithTimeStamp> *link) {
  assert(ToulBar2::vac);
  VAC.push(link, wcsp->getNbNodes());
}

void VACExtension::queueVAC2(DLink<Variable *> *link) {
  assert(ToulBar2::vac);
  VAC2.push(link);
}

void VACExtension::dequeueVAC2(DLink<Variable *> *link) {
  assert(ToulBar2::vac);
  VAC2.remove(link);
}


void VACExtension::printStat(bool ini)
{
	long double mean = (long double) sumlb / (long double) nlb;
	cout << "VAC mean lb/incr: " << mean << "     total increments: " << nlb << endl; 
	if(ini) cout << "Lb after VAC: " << wcsp->getLb() << endl; 
	
	//sort(heap.begin(), heap.end(), cmp_function);
	/*cout << "Vars: ";
	vector<tVACStat*>::iterator it = heap.begin();
	while(it != heap.end()) {
		tVACStat* v = *it;
		if(v->sumlb) cout << "(" << v->var << "," << v->sumlb << ") "; 
		++it;
	}
	cout << endl;*/
}




