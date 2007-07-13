/*
 * ****** Enforce VAC in a WCSP.
 */
 
#include "tb2vac.hpp"

VACExtension::VACExtension (WCSP *w, bool a, bool d) : VAC2(&w->getStore()->storeVariable), wcsp(w), alternative(a), decomposition(d), nbIterations(0), inconsistentVariable(-1), lastCostThresholdUpdate(-1), s(0), vacOddsRecorder(NULL) {
  // temporary disable threshold trick
  queueP = new stack< pair<int, int> >;
  queueR = new stack< pair<int, int> >;
  minlambda = MAX_COST;
}

VACExtension::~VACExtension () {
  delete queueP;
  delete queueR;
  delete vacOddsRecorder;
}

void VACExtension::clean()
{
  unsigned int i;
  VACVariable *x;
  for (i = 0; i < wcsp->numberOfVariables(); i++) {
    x = (VACVariable *) wcsp->getVar(i);
	x->clear();
  }
}


bool VACExtension::propagate() {
  if (vacOddsRecorder == NULL) {
    vacOddsRecorder = new VACOddsRecorder(wcsp->numberOfVariables());
  }
  // SHOULD CHANGE THAT TOO!!!
  //if (storeData->getDepth() < 5) {
  //  return;
  //}
  
  Cost lb = wcsp->getLb();
  nbIterations = 0;
  vacOddsRecorder->reset();

  bool util = true;

  reset();
  enforcePass1();
  if(!isVAC()) {
  	nbIterations++;
    enforcePass2();
    util = enforcePass3();
  }
 
  if (ToulBar2::verbose > 0) {
    if (wcsp->getLb() > lb) {
      cout << "VAC lb :" << lb << " -> " << wcsp->getLb() << endl;
      /*if (decomposition) {
        cout << "\tS is " << getS() << endl;
      }
      cout << "\tnb integral increases: " << getVACOddsRecorder()->getNbIntegerCalls() << endl;
      cout << "\tnb rational increases: " << getVACOddsRecorder()->getNbRationalCalls();
      if (getVACOddsRecorder()->getNbRationalCalls() > 0) {
        cout <<  ", where: " << endl;
        cout << "\t\tmin nb variables: " << getVACOddsRecorder()->getMinUsed() << endl;
        cout << "\t\taverage nb variables: " << getVACOddsRecorder()->getMeanUsed() << endl;
        cout << "\t\tmax nb variables: " << getVACOddsRecorder()->getMaxUsed() << endl;
      }*/
    }
  }
  assert(empty());
  return util;
}

bool VACExtension::verify () {
  return (checkPass1());
}

bool VACExtension::empty () {
  return (true);
}

void VACExtension::clear () {
  VACVariable *xi;
  while (!VAC.empty()) {
    xi = (VACVariable *) VAC.pop();
  }
}

void VACExtension::reset () {
  unsigned int i, j;
  VACConstraint *cij;
  VACVariable *xi;
  
  // reset the queues
  while (!queueP->empty()) {
    queueP->pop();
  }
  while (!queueR->empty()) {
    queueR->pop();
  }
  for (BTQueue::iterator it = VAC2.begin(); it != VAC2.end(); ++it) {
    xi = (VACVariable*) (*it);
    xi->queueVAC();
  }

  // reset the variables and constraints
  for (i = 0; i < wcsp->numberOfVariables(); i++) {
    xi = (VACVariable *) wcsp->getVar(i);
    xi->propagateNC();
  }
  for (i = 0; i < wcsp->numberOfVariables(); i++) {
    xi = (VACVariable *) wcsp->getVar(i);
    xi->reset();
    if (xi->isEmpty()) {
      inconsistentVariable = -2;
      if (ToulBar2::verbose > 1) {
        cout << "Instance is not NC*!" << endl;
      }
      return;
    }
    for (ConstraintList::iterator iter = xi->getConstrs()->begin(); iter != xi->getConstrs()->end(); ++iter) {
      cij = (VACConstraint *) (*iter).constr;
      j = cij->getVar(1 - (*iter).scopeIndex)->wcspIndex;
      if (i > j) {
        cij->reset();
      }
    }
  }
  
  minlambda = wcsp->getUb() - wcsp->getLb();
}

void VACExtension::enforcePass1 () {
  unsigned int i, j;
  VACConstraint *cij;
  VACVariable *xi, *xj;
  if (ToulBar2::verbose > 1) {
	  cout << "------------------" << endl;
	  cout << "VAC Enforce Pass 1" << endl; 
	  cout << "------------------" << endl; 
  }
  while (!VAC.empty()) {
    xj = (VACVariable*) VAC.pop_first();
    j = xj->wcspIndex;
    for (ConstraintList::iterator iter = xj->getConstrs()->begin(); iter != xj->getConstrs()->end(); ++iter) {
      cij = (VACConstraint *) (*iter).constr;
      i = cij->getVar(1 - (*iter).scopeIndex)->wcspIndex;
      xi = (VACVariable *) wcsp->getVar(i);
      for (unsigned int v = 0; v < xi->getDomainSize(); v++) {
        if (!xi->isRemoved(v)) {
          if (((VACConstraint *) cij)->revise(1 - (*iter).scopeIndex, v)) {
            if (ToulBar2::verbose > 1) cout << "value (x" << i << ", " << xi->getValue(v)->getValue() << ") removed because of c" << i << "," << j << endl;
            xi->removeValue(v);
            xi->getValue(v)->setKiller(j);
            queueP->push(pair<int, int>(i, v));
            xi->queueVAC();
            xi->queueVAC2();
            xj->queueVAC2();
          }
        }
      }
      if (xi->isEmpty()) {
        inconsistentVariable = i;
	    if (ToulBar2::verbose > 1) cout << "Inconsistent var: " << i << endl;
        return;
      }
    }
  }
  inconsistentVariable = -1;
}

bool VACExtension::checkPass1 () const {
  unsigned int i, j;
  VACConstraint *cij;
  VACVariable *xi, *xj;
  bool supportFound;
  
  for (i = 0; i < wcsp->numberOfVariables(); i++) {
    xi = (VACVariable *) wcsp->getVar(i);
    for (ConstraintList::iterator iter = xi->getConstrs()->begin(); iter != xj->getConstrs()->end(); ++iter) {
      cij = (VACConstraint *) (*iter).constr;
      j = cij->getVar(1 - (*iter).scopeIndex)->wcspIndex;
      xj = (VACVariable *) wcsp->getVar(j);
      for (unsigned int v = 0; v < xi->getDomainSize(); v++) {
        if (!xi->isRemoved(v)) {
          supportFound = false;
          for (unsigned int w = 0; w < xj->getDomainSize(); w++) {
            if (!xj->isRemoved(w)) {
              if ((!xi->getVACCost(v)) && 
                 (!cij->getVACCost(v, w, (*iter).scopeIndex)) && 
                 (!xj->getVACCost(w))) {
                supportFound = true;
                break;
              }
            }
          }
          if (!supportFound) {
            if (ToulBar2::verbose > 1) {
              cout << "(x" << i << ", " << xi->getValue(v)->getValue() << ") has no support wrt x" << j << "\t\t(c" << i << "(" << xi->getValue(v)->getValue() << ") = " << xi->getVACCost(v) << ")" << endl;
              cout << *cij << endl;
            }
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
  unsigned int v;
  if (ToulBar2::verbose > 1) {
	  cout << "------------------" << endl; 
	  cout << "VAC Enforce Pass 2" << endl; 
	  cout << "------------------" << endl; 
  }
  assert(i0 >= 0);
  xi0 = (VACVariable *) wcsp->getVar(i0);
  for (unsigned int v = 0; v < xi0->getDomainSize(); v++) {
    xi0->getValue(v)->addToK(1);
    xi0->getValue(v)->setMark();
	Cost cost = xi0->getVACCost(v);
	if (cost > 0 && cost < minlambda) minlambda = cost;
  }

  while (!queueP->empty()) {
    i = queueP->top().first;
    v = queueP->top().second;
    queueP->pop();
    xi = (VACVariable *) wcsp->getVar(i);
    if (xi->getValue(v)->isMarked()) {
      j = xi->getValue(v)->getKiller();
      xj = (VACVariable *) wcsp->getVar(j);
      queueR->push(pair<int, int>(i, v));
      cij = (VACConstraint *) xi->getConstr(xj);

	  if (ToulBar2::verbose > 2) {
	        cout << "VAC pops x" << i << "(" << v << ") = " << xi->getVACCost(v) << "      k: " << xi->getValue(v)->getK() << "            killed by: " << j << endl;
	  }
      for (unsigned int w = 0; w < xj->getDomainSize(); w++) {
		Cost costij = cij->getVACCost(v, w, cij->getIndex(xi));
		if (costij > 0) {
          Cost tmpK = xi->getValue(v)->getK();
		  if (xj->getValue(w)->getKiller() == i) tmpK += xj->getValue(w)->getK();
		  
		  if (wcsp->getLb() + costij < wcsp->getUb()) {
			if( (costij/tmpK) < minlambda) minlambda = costij/tmpK;
		  }
		} else {
          Cost tmpK = xi->getValue(v)->getK() - cij->getK(cij->getIndex(xj), w);
          if (tmpK > 0) {
            xj->getValue(w)->addToK(tmpK);
            cij->setK(cij->getIndex(xj), w, xi->getValue(v)->getK());
            if(!xj->getVACCost(w)) xj->getValue(w)->setMark();
            else if (wcsp->getLb() + xj->getVACCost(w) < wcsp->getUb()) {
              tmplambda = xj->getVACCost(w) / xj->getValue(w)->getK();  	
			  if(tmplambda < minlambda) minlambda = tmplambda;
            }
          }  
        }
  		if (ToulBar2::verbose > 4) {
	        cout << "x" << j << "(" << w << ") = " << xj->getVACCost(w) << "      k: " << xj->getValue(w)->getK() << endl;
	        cout << "c" << i << "," << j << "(" << v << "," << w << ") = " << cij->getVACCost(v, w, cij->getIndex(xi)) << "      k: " << cij->getK(cij->getIndex(xj), w) << endl;
		}      
      }
    }
  }
  //if (maxK == 0) {
  //  maxK = wcsp->getUb() - wcsp->getLb();
  //}
  if (ToulBar2::verbose > 1) cout << "---------minLambda: " << minlambda << "\t\t (lb = " << wcsp->getLb() << ", ub = " << wcsp->getUb() << ")" << endl;
}

bool VACExtension::enforcePass3 () {
  bool done = false;
  vacOddsRecorder->addNbIntegerCalls(); 

  if (decomposition) enforcePass3VACDecomposition();
  else {
    //enforcePass3VAC();
    done = enforcePass3VACint();
  }
  return done;
}
  
void VACExtension::enforcePass3VAC () {
  Cost lambda = (Cost) minlambda;
  int i, j;
  VACVariable *xi, *xj;
  VACConstraint *cij;
  int i0 = inconsistentVariable;
  VACVariable *xi0 = (VACVariable *) wcsp->getVar(i0);
  unsigned int w;

  while (!queueR->empty()) {
    j = queueR->top().first;
    w = queueR->top().second;
    queueR->pop();
    xj = (VACVariable *) wcsp->getVar(j);
    vacOddsRecorder->addVariable(j);
    i = xj->getValue(w)->getKiller();
    xi = (VACVariable *) wcsp->getVar(i);
    vacOddsRecorder->addVariable(i);
    cij = (VACConstraint *) xi->getConstr(xj);
    for (unsigned int v = 0; v < xi->getDomainSize(); v++) {
      if ((wcsp->getLb() + cij->getVACCost(v, w, cij->getIndex(xi)) < wcsp->getUb()) && (cij->getVACCost(v, w, cij->getIndex(xi)) < lambda * xj->getValue(w)->getK())) {
        cij->VACextend(cij->getIndex(xi), v, lambda * xj->getValue(w)->getK() - cij->getVACCost(v, w, cij->getIndex(xi)));
      }
    }
    cij->VACproject(cij->getIndex(xj), w, lambda * xj->getValue(w)->getK());
  }
  if (ToulBar2::verbose > 5) {
    cout << "lower bound increased " << wcsp->getLb() << " -> " << wcsp->getLb()+lambda << endl;
  }

  xi0->extendAll(lambda);
  wcsp->increaseLb(wcsp->getLb() + lambda);
  vacOddsRecorder->computeOdds();
}



bool VACExtension::enforcePass3VACint () {
  bool util = (minlambda >= 1);
  Cost lambda = minlambda;	
  if (ToulBar2::verbose > 1) {
	  cout << "------------------" << endl; 
	  cout << "VAC Enforce Pass 3" << endl; 
	  cout << "------------------" << endl; 
  }
  int i, j;
  VACVariable *xi, *xj;
  VACConstraint *cij;
  int i0 = inconsistentVariable;
  VACVariable *xi0 = (VACVariable *) wcsp->getVar(i0);
  unsigned int w;

	  
  int minvar = -1;
  Cost minc = MAX_VAL;
  Cost maxk = 0;
  	  
  if(!util) {
    while (!queueR->empty()) { 
	    j = queueR->top().first;
	    w = queueR->top().second; 
	    xj = (VACVariable *) wcsp->getVar(j);
		i = xj->getValue(w)->getKiller();
		xi = (VACVariable *) wcsp->getVar(i);
		cij = (VACConstraint *) xi->getConstr(xj);
	    queueR->pop();

	    for (unsigned int v = 0; v < xi->getDomainSize(); v++) {
		    if(xi->getValue(v)->getK()) { 
		    	if(xi->getVACCost(v)) {
			    	if(minc > xi->getIniCost(v)) {
			    		minc = xi->getIniCost(v);
			    		minvar = i;
			    	}
			    	if(maxk < xi->getValue(v)->getK()) maxk = xi->getValue(v)->getK();
			    	
			    } 
		    }
	    }
    }
	xi = (VACVariable *) wcsp->getVar(minvar);	
	xi->setThreshold( minc + 1 );
    if (ToulBar2::verbose > 1) cout << "BreakCycle   var: " << minvar << "   thr: " <<  xi->getThreshold() << "          maxk: " << maxk << endl; 	
    return false;
  }

  while (!queueR->empty()) {
    j = queueR->top().first;
    w = queueR->top().second;
    queueR->pop();
    xj = (VACVariable *) wcsp->getVar(j);
    i = xj->getValue(w)->getKiller();
    xi = (VACVariable *) wcsp->getVar(i);
    cij = (VACConstraint *) xi->getConstr(xj);

    vacOddsRecorder->addVariable(i);
    vacOddsRecorder->addVariable(j);
    for (unsigned int v = 0; v < xi->getDomainSize(); v++) {
      if ((wcsp->getLb() + cij->getIniCost(v, w, cij->getIndex(xi)) < wcsp->getUb()) && 
      	  (cij->getIniCost(v, w, cij->getIndex(xi)) < lambda * xj->getValue(w)->getK())) {
          if (ToulBar2::verbose > 1) {
            cout << "extending " << (lambda * xj->getValue(w)->getK() - cij->getIniCost(v, w, cij->getIndex(xi))) << " from (x" << i << ", " << xi->getValue(v)->getValue() << ") to c" << i << "," << j << "   (" << xi->getIniCost(v) << " is available) and " << xi->myThreshold << endl;
          }
          Cost ecost = lambda * xj->getValue(w)->getK() - cij->getIniCost(v, w, cij->getIndex(xi));
	      cij->VACextend(cij->getIndex(xi), v, ecost);
      }
    }
    if (ToulBar2::verbose > 1) {
      cout << "projecting " << (lambda * xj->getValue(w)->getK()) << " from c" << j << "," << i << " to (x" << j << ", " << xj->getValue(w)->getValue() << ")" << endl;
    }
    cij->VACproject(cij->getIndex(xj), w, lambda * xj->getValue(w)->getK());
	if (w == (unsigned int) xj->getSupport()) {
        Value newSupport = xj->getInf();
        Cost minCost = xj->getCost(newSupport);
        EnumeratedVariable::iterator iter = xj->begin();
        for (++iter; minCost > 0 && iter != xj->end(); ++iter) {
            Cost cost = xj->getCost(*iter);
            if (cost < minCost) {
                minCost = cost;
                newSupport = *iter;
            }
        }
		xj->setSupport(newSupport);
	}
  }
  xi0->extendAll(lambda);
  wcsp->increaseLb(wcsp->getLb() + lambda);
  vacOddsRecorder->computeOdds();
  return true;
}





void VACExtension::enforcePass3VACDecomposition () {
  unsigned int i, j;
  VACConstraint *cij;
  VACVariable *xi, *xj;
  Cost lambda = (Cost) minlambda;
  Cost cost;

  for (i = 0; i < wcsp->numberOfVariables(); i++) {
    xi = (VACVariable *) wcsp->getVar(i);
    for (unsigned int v = 0; v < xi->getDomainSize(); v++) {
      if (xi->getValue(v)->getK() != 0) {
        if (xi->getVACCost(v) != 0) {
          vacOddsRecorder->addVariable(i);
          cost = ceil(lambda * xi->getValue(v)->getK());
          if (ToulBar2::verbose > 5) {
            cout << "decreasing " << cost << " from (x" << i << ", " << xi->getValue(v)->getValue() << ")     (" << xi->getVACCost(v) << " is available)" << endl;
          }
          xi->decreaseCost(v, cost);
          assert(xi->getVACCost(v) >= 0);
        }
        else {
          if (xi->getValue(v)->getKiller() != -1) {
            j = xi->getValue(v)->getKiller();
            xj = (VACVariable *) wcsp->getVar(j);
            cij = (VACConstraint *) xi->getConstr(xj);
            for (unsigned int w = 0; w < xj->getDomainSize(); w++) {
              if (cij->getVACCost(v, w, cij->getIndex(xi))) {
                vacOddsRecorder->addVariable(i);
                vacOddsRecorder->addVariable(j);
                cost = ceil(lambda * xi->getValue(v)->getK());
                if (ToulBar2::verbose > 5) {
                  cout << "decreasing " << cost << " from c" << i << "," << j << "(" << xi->getValue(v)->getValue() << ", " << xj->getValue(w)->getValue() << ") to x" << i << "   (" << cij->getVACCost(v, w, cij->getIndex(xi)) << " is available)" << endl;
                }
                cij->decreaseCost(v, w, cij->getIndex(xi), cost);
                assert(cij->getVACCost(v, w, cij->getIndex(xi)) >= 0);
              }
            }
          }
        }
      }
    }
  }
  s += ceil(lambda);
  vacOddsRecorder->computeOdds();

  // SetTop
  //   Not used because too time consuming
  /*
  for (i = 0; i < wcsp->numberOfVariables(); i++) {
    xi = (VACVariable *) wcsp->getVar(i);
    xi->setTop(s);
    for (ConstraintList::iterator iter = xi->getConstrs()->begin(); iter != xi->getConstrs()->end(); ++iter) {
      cij = (VACConstraint *) (*iter).constr;
      j = cij->getVar(1 - (*iter).scopeIndex)->wcspIndex;
      if (i < j) {
        cij->setTop(s);
      }
    }
  }
  */
}

void VACExtension::updateQueue() {
  unsigned int i, j;
  VACConstraint *cij;
  VACVariable *xi, *xj;
  bool fullSupport = false;
  int size = VAC2.getSize();
  
  for (int s = 0; s < size; s++) {
    xj = (VACVariable*) VAC2.pop_first();
    j = xj->wcspIndex;
    for (ConstraintList::iterator iter = xj->getConstrs()->begin(); iter != xj->getConstrs()->end(); ++iter) {
      cij = (VACConstraint *) (*iter).constr;
      i = cij->getVar(1 - (*iter).scopeIndex)->wcspIndex;
      xi = (VACVariable *) wcsp->getVar(i);
      for (unsigned int v = 0; v < xi->getDomainSize(); v++) {
        fullSupport = false;
        if (!xi->getVACCost(v)) {
          for (unsigned int w = 0; w < xj->getDomainSize(); w++) {
            if ((!cij->getVACCost(v, w, 1 - (*iter).scopeIndex)) &&  (!xj->getVACCost(w))) {
              fullSupport = true;
              break;
            }
          }
        }
        if (!fullSupport) {
          xj->queueVAC2();
          break;
        }
      }
      if (!fullSupport) {
        break;
      }
    }
  }
}

Cost VACExtension::getS () const {
  return s;
}

VACOddsRecorder *VACExtension::getVACOddsRecorder () {
  return vacOddsRecorder;
}

bool VACExtension::remainedIdle () {
  return (nbIterations == 0);
}

Cost VACExtension::getCostThreshold () const {
  return ToulBar2::costThreshold;
}

bool VACExtension::isNull (Cost c) const {
  if (ToulBar2::vacAlternative) {
    return (c < ToulBar2::costThreshold);
  }
  else {
    return (c == 0);
  }
}

Long VACExtension::getLastCostThresholdUpdate () const {
  return lastCostThresholdUpdate;
}

void VACExtension::updateThreshold (Cost t, Long d) {
  // temporary disable threshold trick
  /*
  VACVariable *xi;
  lastCostThresholdUpdate = d;
  if (t == costThreshold) {
    return;
  }
  costThreshold = t;
  for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
    xi = (VACVariable *) wcsp->getVar(i);
    xi->queueVAC2();
  }
  */
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
