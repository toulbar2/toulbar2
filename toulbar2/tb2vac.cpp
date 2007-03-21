/*
 * ****** Enforce VAC in a WCSP.
 */
 
#include "tb2vac.hpp"

VACExtension::VACExtension (WCSP *w, bool a, bool d) : VAC2(&w->getStore()->storeVariable), wcsp(w), alternative(a), decomposition(d), nbIterations(0), maxK(0), inconsistentVariable(-1), costThreshold(0), lastCostThresholdUpdate(-1), s(0), vacOddsRecorder(NULL) {
  queueP = new stack< pair<int, int> >;
  queueR = new stack< pair<int, int> >;
}

VACExtension::~VACExtension () {
  delete queueP;
  delete queueR;
  delete vacOddsRecorder;
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

void VACExtension::propagate() {
  if (vacOddsRecorder == NULL) {
    vacOddsRecorder = new VACOddsRecorder(wcsp->numberOfVariables());
  }

  Cost c;
  
  // SHOULD CHANGE THAT!!!
  if ((getLastCostThresholdUpdate() == -1) || (getLastCostThresholdUpdate() < wcsp->getNbNodes() + 5)) {
	c = (Long) pow(10., floor(log10(to_double(wcsp->getUb() - wcsp->getLb())))-1.);
	if (c == 0) {
	  c = 1;
	}
	updateThreshold(c, wcsp->getNbNodes());
  }

  // SHOULD CHANGE THAT TOO!!!
  //if (storeData->getDepth() < 5) {
  //  return;
  //}
  
  Cost lb = wcsp->getLb();
  nbIterations = 0;
  vacOddsRecorder->reset();

  if (ToulBar2::verbose) {
    cout << "Threshold: " << costThreshold << endl;
    cout << "VAC2";
    VAC2.print(cout);
  }

  reset();

  enforcePass1();
  while (!isVAC()) {
    nbIterations++;
    enforcePass2();
    enforcePass3();
    reset();
    enforcePass1();
  }
  updateQueue();
  assert(checkPass1());

  if (ToulBar2::verbose) {
    if (wcsp->getLb() > lb) {
      cout << "VAC: lb :" << lb << " -> " << wcsp->getLb() << endl;
      if (decomposition) {
        cout << "\tS is " << getS() << endl;
      }
      cout << "\tnb integral increases: " << getVACOddsRecorder()->getNbIntegerCalls() << endl;
      cout << "\tnb rational increases: " << getVACOddsRecorder()->getNbRationalCalls();
      if (getVACOddsRecorder()->getNbRationalCalls() > 0) {
        cout <<  ", where: " << endl;
        cout << "\t\tmin nb variables: " << getVACOddsRecorder()->getMinUsed() << endl;
        cout << "\t\taverage nb variables: " << getVACOddsRecorder()->getMeanUsed() << endl;
        cout << "\t\tmax nb variables: " << getVACOddsRecorder()->getMaxUsed() << endl;
      }
      cout << endl;
        }
  }

  assert(empty());
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
      if (ToulBar2::verbose > 5) {
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
  
  maxK = 0;
}

void VACExtension::enforcePass1 () {
  unsigned int i, j;
  VACConstraint *cij;
  VACVariable *xi, *xj;
  
  while (!VAC.empty()) {
    if (ToulBar2::verbose > 5) {
      cout << "VAC";
      VAC.print(cout);
    }
    
    xj = (VACVariable*) VAC.pop_first();
    j = xj->wcspIndex;
    for (ConstraintList::iterator iter = xj->getConstrs()->begin(); iter != xj->getConstrs()->end(); ++iter) {
      cij = (VACConstraint *) (*iter).constr;
      i = cij->getVar(1 - (*iter).scopeIndex)->wcspIndex;
      xi = (VACVariable *) wcsp->getVar(i);
      for (unsigned int v = 0; v < xi->getDomainSize(); v++) {
        if (!xi->isRemoved(v)) {
          if (((VACConstraint *) cij)->revise(1 - (*iter).scopeIndex, v)) {
            if (ToulBar2::verbose > 5) {
              cout << "value (x" << i << ", " << xi->getValue(v)->getValue() << ") removed because of c" << i << "," << j << endl;
            }
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
              if ((xi->getVACCost(v) <= costThreshold) && (cij->getVACCost(v, w, (*iter).scopeIndex) <= costThreshold) && (xj->getVACCost(w) <= costThreshold)) {
                supportFound = true;
                break;
              }
            }
          }
          if (!supportFound) {
            if (ToulBar2::verbose > 5) {
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
  Cost tmpK;
  unsigned int v;

  assert(i0 >= 0);

  xi0 = (VACVariable *) wcsp->getVar(i0);
  for (unsigned int v = 0; v < xi0->getDomainSize(); v++) {
    xi0->getValue(v)->addToK(1);
    xi0->getValue(v)->setMark();
  }
  if (ToulBar2::verbose > 5) {
    cout << "x" << i0 << " marked" << endl;
  }

  while (!queueP->empty()) {
    i = queueP->top().first;
    v = queueP->top().second;
    queueP->pop();
    xi = (VACVariable *) wcsp->getVar(i);
    if (xi->getValue(v)->isMarked()) {
      j = xi->getValue(v)->getKiller();
      if (ToulBar2::verbose > 5) {
        cout << "projection of " << xi->getValue(v)->getK() << " from x" << j << " to (x" << i << ", " << xi->getValue(v)->getValue() << ")" << endl;
      }
      xj = (VACVariable *) wcsp->getVar(j);
      queueR->push(pair<int, int>(i, v));
      cij = (VACConstraint *) xi->getConstr(xj);
      for (unsigned int w = 0; w < xj->getDomainSize(); w++) {
        if (ToulBar2::verbose > 5) {
          cout << "\t\tvalue " << xj->getValue(w)->getValue() << ": ";
        }
        if (((alternative) && (cij->getVACCost(v, w, cij->getIndex(xi)) >= 1)) ||
             ((!alternative) && (cij->getVACCost(v, w, cij->getIndex(xi)) != 0))) {
          if (ToulBar2::verbose > 5) {
            cout << "cost provided by binary constraint (" << cij->getVACCost(v, w, cij->getIndex(xi)) << " available)" << endl;
          }
          tmpK = xi->getValue(v)->getK() / cij->getVACCost(v, w, cij->getIndex(xi));
          if (xj->getValue(w)->getKiller() == i) {
            tmpK += xj->getValue(w)->getK() / cij->getVACCost(v, w, cij->getIndex(xi));
          }
          if (wcsp->getLb() + cij->getVACCost(v, w, cij->getIndex(xi)) < wcsp->getUb()) {
            maxK = (tmpK > maxK)? tmpK: maxK;
          }
        }
        else {
          if (ToulBar2::verbose > 5) {
            cout << "cost provided by unary constraint (" << xj->getVACCost(w) << " available)";
          }
          tmpK = xi->getValue(v)->getK() - cij->getK(cij->getIndex(xj), w);
          if (ToulBar2::verbose > 5) {
            cout << "     (" << cij->getK(cij->getIndex(xj), w) << " has been extended)" << endl;
          }
          if (tmpK > 0) {
            if (ToulBar2::verbose > 5) {
              cout << "\t\t\textending " << tmpK << endl;
            }
            xj->getValue(w)->addToK(tmpK);
            cij->setK(cij->getIndex(xj), w, xi->getValue(v)->getK());
            if (((alternative) && (xj->getVACCost(w) < 1)) ||
                ((!alternative) && (xj->getVACCost(w) == 0))) {
              xj->getValue(w)->setMark();
            }
            else if (wcsp->getLb() + xj->getVACCost(w) < wcsp->getUb()) {
              maxK = ((xj->getValue(w)->getK() / xj->getVACCost(w)) > maxK)? (xj->getValue(w)->getK() / xj->getVACCost(w)): maxK;
            }
          }
        }
      }
    }
  }
  if (maxK == 0) {
    maxK = wcsp->getUb() - wcsp->getLb();
  }
  if (ToulBar2::verbose > 5) {
    cout << "maxK: " << maxK << "\t\t (lb = " << wcsp->getLb() << ", ub = " << wcsp->getUb() << ")" << endl;
  }
}

void VACExtension::enforcePass3 () {
  Cost lambda = 1/maxK;
  double lambdaD = to_double(lambda);
  if (floor(lambdaD) == lambdaD) {
    vacOddsRecorder->addNbIntegerCalls();
  }
  else {
    vacOddsRecorder->addNbRationalCalls();
  }
  if (decomposition) {
    if (floor(lambdaD) == lambdaD) {
      enforcePass3VAC();
    }
    else {
      enforcePass3VACDecomposition();
    }
  }
  else {
    enforcePass3VAC();
  }
}
  
void VACExtension::enforcePass3VAC () {
  Cost lambda = 1/maxK;
  int i, j;
  VACVariable *xi, *xj;
  VACConstraint *cij;
  Cost tmpK;
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
        if (ToulBar2::verbose > 5) {
          cout << "extending " << (lambda * xj->getValue(w)->getK() - cij->getVACCost(v, w, cij->getIndex(xi))) << " from (x" << i << ", " << xi->getValue(v)->getValue() << ") to c" << i << "," << j << "   (" << xi->getVACCost(v) << " is available)" << endl;
        }
        cij->VACextend(cij->getIndex(xi), v, lambda * xj->getValue(w)->getK() - cij->getVACCost(v, w, cij->getIndex(xi)));
      }
    }
    if (ToulBar2::verbose > 5) {
      cout << "projecting " << (lambda * xj->getValue(w)->getK()) << " from c" << j << "," << i << " to (x" << j << ", " << xj->getValue(w)->getValue() << ")" << endl;
    }
    cij->VACproject(cij->getIndex(xj), w, lambda * xj->getValue(w)->getK());
  }
  xi0->extendAll(lambda);
  if (ToulBar2::verbose > 5) {
    cout << "lower bound increased " << wcsp->getLb() << " -> " << wcsp->getLb()+lambda << endl;
  }
  wcsp->increaseLb(wcsp->getLb() + lambda);
  vacOddsRecorder->computeOdds();
}

/*
void VACExtension::enforcePass3VACDecomposition () {
  Cost lambda = 1/maxK;
  int i, j;
  VACVariable *xi, *xj;
  VACConstraint *cij;
  unsigned int w;
  Cost cost;

  while (!queueR->empty()) {
    j = queueR->top().first;
    w = queueR->top().second;
    queueR->pop();
    vacOddsRecorder->addVariable(j);
    xj = (VACVariable *) wcsp->getVar(j);
    i = xj->getValue(w)->getKiller();
    xi = (VACVariable *) wcsp->getVar(i);
    cij = (VACConstraint *) xi->getConstr(xj);
    for (unsigned int v = 0; v < xi->getDomainSize(); v++) {
      if (cij->getExtended(cij->getIndex(xi), v) < lambda * xj->getValue(w)->getK()) {
        cost = ceil(lambda * xj->getValue(w)->getK() - cij->getExtended(cij->getIndex(xi), v));
        if (xi->getVACCost(v) != 0) {
          if (ToulBar2::verbose > 5) {
            cout << "decreasing " << cost << " from (x" << i << ", " << v << ") to c" << i << "," << j << "   (" << xi->getVACCost(v) << " is available)" << endl;
          }
          xi->decreaseCost(xi->getValue(v)->getValue(), cost);
          assert(xi->getVACCost(v) >= 0);
          cij->addExtended(cij->getIndex(xi), v, cost);
        }
        else if (cij->getVACCost(v, w, cij->getIndex(xi)) != 0) {
          cost = ceil(lambda * xj->getValue(w)->getK());
          if (ToulBar2::verbose > 5) {
            cout << "decreasing " << cost << " from c" << j << "," << i << "(" << w << ", " << v << ") to x" << j << "   (" << cij->getVACCost(v, w, cij->getIndex(xi)) << " is available)" << endl;
          }
          cij->decreaseCost(v, w, cij->getIndex(xi), cost);
          assert(cij->getVACCost(v, w, cij->getIndex(xi)) >= 0);
        }
      }
    }
  }
  s += ceil(lambda);
  vacOddsRecorder->computeOdds();
}
*/

void VACExtension::enforcePass3VACDecomposition () {
  unsigned int i, j;
  VACConstraint *cij;
  VACVariable *xi, *xj;
  Cost lambda = 1/maxK;
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
              if (cij->getVACCost(v, w, cij->getIndex(xi)) != 0) {
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

/*
void VACExtension::enforcePass3VAC () {
  Cost lambda = 1/maxK;
  unsigned int i, j;
  VACVariable *xi, *xj;
  VACConstraint *cij;
  Cost cost;
  unsigned int w;
  int i0 = inconsistentVariable;
  VACVariable *xi0 = (VACVariable *) wcsp->getVar(i0);

  while (!queueR->empty()) {
    j = queueR->top().first;
    w = queueR->top().second;
    queueR->pop();
    xj = (VACVariable *) wcsp->getVar(j);
    i = xj->getValue(w)->getKiller();
    xi = (VACVariable *) wcsp->getVar(i);
    cij = (VACConstraint *) xi->getConstr(xj);
    for (unsigned int v = 0; v < xi->getDomainSize(); v++) {
      if (cij->getVACCost(v, w, cij->getIndex(xi)) < lambda * xj->getValue(w)->getK()) {
        if (ToulBar2::verbose > 5) {
          cout << "extending " << (lambda * xj->getValue(w)->getK() - cij->getVACCost(v, w, cij->getIndex(xi))) << " from (x" << i << ", " << v << ") to c" << i << "," << j << "   (" << xi->getVACCost(v) << " is available)" << endl;
        }
        cij->VACextend(cij->getIndex(xi), v, lambda * xj->getValue(w)->getK() - cij->getVACCost(v, w, cij->getIndex(xi)));
      }
    }
    if (ToulBar2::verbose > 5) {
      cout << "projecting " << (lambda * xj->getValue(w)->getK()) << " from c" << j << "," << i << " to (x" << j << ", " << w << ")" << endl;
    }
    cij->VACproject(cij->getIndex(xj), w, lambda * xj->getValue(w)->getK());
  }
  xi0->extendAll(lambda);
  if (ToulBar2::verbose > 5) {
    cout << "lower bound increased " << wcsp->getLb() << " -> " << wcsp->getLb()+lambda << endl;
  }
  wcsp->increaseLb(wcsp->getLb() + lambda);
  
  for (i = 0; i < wcsp->numberOfVariables(); i++) {
    xi = (VACVariable *) wcsp->getVar(i);
    xi->repair();
    for (ConstraintList::iterator iter = xi->getConstrs()->begin(); iter != xi->getConstrs()->end(); ++iter) {
      cij = (VACConstraint *) (*iter).constr;
      j = cij->getVar(1 - (*iter).scopeIndex)->wcspIndex;
      if (i < j) {
        cij->repair();
      }
    }
  }

  for (i = 0; i < wcsp->numberOfVariables(); i++) {
    xi = (VACVariable *) wcsp->getVar(i);
    for (unsigned int v = 0; v < xi->getDomainSize(); v++) {
      cost = floor(xi->getVACCost(v));
      if (xi->getVACCost(v) != cost) {
        if (ToulBar2::verbose > 5) {
          cout << "decreasing (x" << i << ", " << v << ") to " << cost << "    (" << xi->getVACCost(v) << " is stored)" << endl;
        }
        xi->setCost(xi->getValue(v)->getValue(), cost);
        assert(xi->getVACCost(v) >= 0);
      }
      for (ConstraintList::iterator iter = xi->getConstrs()->begin(); iter != xi->getConstrs()->end(); ++iter) {
        cij = (VACConstraint *) (*iter).constr;
        j = cij->getVar(1 - (*iter).scopeIndex)->wcspIndex;
        if (j > i) {
          xj = (VACVariable *) wcsp->getVar(j);
          for (unsigned int w = 0; w < xj->getDomainSize(); w++) {
            cost = (floor(cij->getVACCost(v, w, cij->getIndex(xi))));
            if (cij->getVACCost(v, w, cij->getIndex(xi)) != cost) {
              if (ToulBar2::verbose > 5) {
                cout << "decreasing c" << i << "," << j << "(" << v << ", " << w << ") to " << cost << "   (" << cij->getVACCost(v, w, cij->getIndex(xi)) << " is stored)" << endl;
              }
              cij->setCost(v, w, cij->getIndex(xi), cost);
              assert(cij->getVACCost(v, w, cij->getIndex(xi)) >= 0);
            }
          }
        }
      }
    }
  }
}
*/

void VACExtension::updateQueue() {
  unsigned int i, j;
  VACConstraint *cij;
  VACVariable *xi, *xj;
  bool fullSupport;
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
        if (xi->getVACCost(v) <= getCostThreshold()) {
          for (unsigned int w = 0; w < xj->getDomainSize(); w++) {
            if ((cij->getVACCost(v, w, 1 - (*iter).scopeIndex) <= getCostThreshold()) && (xj->getVACCost(w) <= getCostThreshold())) {
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

Cost VACExtension::getCostThreshold () {
  return costThreshold;
}

Long VACExtension::getLastCostThresholdUpdate () {
  return lastCostThresholdUpdate;
}

void VACExtension::updateThreshold (Cost t, Long d) {
  VACVariable *xi;
  lastCostThresholdUpdate = d;
  if (t == costThreshold) {
    return;
  }
  costThreshold = 1; // t;
  for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
    xi = (VACVariable *) wcsp->getVar(i);
    xi->queueVAC2();
  }
}
