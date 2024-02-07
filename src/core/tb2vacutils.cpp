/*
 * ****** Set of useful classes to enforce VAC
 */

#include "search/tb2clusters.hpp"
#include "tb2vacutils.hpp"

VACVariable::VACVariable(WCSP* wcsp, string n, Value iinf, Value isup)
    : EnumeratedVariable(wcsp, n, iinf, isup)
    , vac(wcsp->vac)
    , myThreshold(MIN_COST)
{
    init();
}

VACVariable::VACVariable(WCSP* wcsp, string n, vector<Value>& dom)
    : EnumeratedVariable(wcsp, n, dom)
    , vac(wcsp->vac)
    , myThreshold(MIN_COST)
{
    init();
}

VACVariable::~VACVariable()
{
}

void VACVariable::init()
{
    maxk_timeStamp = 0;
    maxk = 0;
    for (unsigned int a = 0; a < getDomainInitSize(); a++) {
        mark.push_back(0);
        k_timeStamp.push_back(0);
        k.push_back(0);
        killer.push_back(0);
        PBkillers.emplace_back();
    }
    linkVACQueue.content.var = this;
    linkVACQueue.content.timeStamp = -1;
#ifdef INCREMENTALVAC
    linkVAC2Queue.content.var = this;
    linkVAC2Queue.content.timeStamp = -1;
#endif
    linkSeekSupport.content.var = this;
    linkSeekSupport.content.timeStamp = -1;
}

// void VACVariable::remove(Value value)
//{
//     if (ToulBar2::singletonConsistency)
//         vac->singleton.insert(wcsp->getMaxDomainSize() * wcspIndex + value);
//     EnumeratedVariable::remove(value);
// }
//
// void VACVariable::removeFast(Value value)
//{
//     if (ToulBar2::singletonConsistency)
//         vac->singleton.insert(wcsp->getMaxDomainSize() * wcspIndex + value);
//     EnumeratedVariable::removeFast(value);
// }
//
// void VACVariable::increase(Value newInf)
//{
//     if (ToulBar2::singletonConsistency)
//         for (int i = inf; i < newInf; i++)
//             vac->singleton.insert(wcsp->getMaxDomainSize() * wcspIndex + i);
//     EnumeratedVariable::increase(newInf);
// }
//
// void VACVariable::decrease(Value newSup)
//{
//     if (ToulBar2::singletonConsistency)
//         for (int i = sup; i > newSup; i--)
//             vac->singleton.insert(wcsp->getMaxDomainSize() * wcspIndex + i);
//     EnumeratedVariable::decrease(newSup);
// }

/******************************
 * min-sum diffusion algorithm
 */

bool VACVariable::averaging()
{
    Tuple tuple;
    Cost Top = wcsp->getUb();
    bool change = false;
    EnumeratedVariable* x;
    EnumeratedVariable* y;
    Constraint* ctr = NULL;
    ConstraintList::iterator itc = getConstrs()->begin();
    if (itc != getConstrs()->end())
        ctr = (*itc).constr;
    while (ctr) {
        if (ctr->isBinary()) {
            BinaryConstraint* bctr = (BinaryConstraint*)ctr;
            x = (EnumeratedVariable*)bctr->getVarDiffFrom((Variable*)this);
            for (iterator it = begin(); it != end(); ++it) {
                Cost cu = getCost(*it);
                Cost cmin = Top;
                for (iterator itx = x->begin(); itx != x->end(); ++itx) {
                    Cost cbin = bctr->getCost(this, x, *it, *itx);
                    if (cbin < cmin)
                        cmin = cbin;
                }
                assert(cmin < Top);
                Double mean = to_double(cmin + cu) / 2.;
                Double extc = to_double(cu) - mean;
                if (std::abs(extc) >= 1) {
                    Cost costi = (Long)extc;
                    for (iterator itx = x->begin(); itx != x->end(); ++itx) {
                        bctr->addcost(this, x, *it, *itx, costi);
                    }
                    if (mean > to_double(cu))
                        project(*it, -costi);
                    else
                        extend(*it, costi);
                    change = true;
                }
            }
        } else if (ctr->isTernary() && !ctr->isSep()) {
            TernaryConstraint* tctr = (TernaryConstraint*)ctr;
            x = (EnumeratedVariable*)tctr->getVar(0);
            if (x == this)
                x = (EnumeratedVariable*)tctr->getVar(1);
            y = (EnumeratedVariable*)tctr->getVarDiffFrom((Variable*)this, (Variable*)x);
            for (iterator it = begin(); it != end(); ++it) {
                Cost cu = getCost(*it);
                Cost cmin = Top;
                for (iterator itx = x->begin(); itx != x->end(); ++itx) {
                    for (iterator ity = y->begin(); ity != y->end(); ++ity) {
                        Cost ctern = tctr->getCost(this, x, y, *it, *itx, *ity);
                        if (ctern < cmin)
                            cmin = ctern;
                    }
                }
                assert(cmin < Top);
                Double mean = to_double(cmin + cu) / 2.;
                Double extc = to_double(cu) - mean;
                if (std::abs(extc) >= 1) {
                    Cost costi = (Long)extc;
                    for (iterator itx = x->begin(); itx != x->end(); ++itx) {
                        for (iterator ity = y->begin(); ity != y->end(); ++ity) {
                            tctr->addCost(this, x, y, *it, *itx, *ity, costi);
                        }
                    }
                    if (mean > to_double(cu))
                        project(*it, -costi);
                    else
                        extend(*it, costi);
                    change = true;
                }
            }
        } else if (ctr->isNary() && !ctr->isSep()) {
            NaryConstraint* nctr = (NaryConstraint*)ctr;
            for (iterator it = begin(); it != end(); ++it) {
                Cost cu = getCost(*it);
                Cost cmin = Top;
                int tindex = nctr->getIndex(this);
                Cost cost;
                Long nbtuples = 0;
                nctr->first();
                while (nctr->next(tuple, cost)) {
                    nbtuples++;
                    if (toValue(tuple[tindex]) == (*it) && cost < cmin)
                        cmin = cost;
                }
                if (nctr->getDefCost() < cmin && nbtuples < nctr->getDomainSizeProduct() / getDomainSize())
                    cmin = nctr->getDefCost();
                //				assert(cmin < Top);
                Double mean = to_double(cmin + cu) / 2.;
                Double extc = to_double(cu) - mean;
                if (std::abs(extc) >= 1) {
                    Cost costi = (Cost)extc;
                    nctr->addtoTuples(this, *it, costi);
                    if (mean > to_double(cu))
                        project(*it, -costi);
                    else
                        extend(*it, costi);
                    change = true;
                }
            }
        }
        ++itc;
        if (itc != getConstrs()->end())
            ctr = (*itc).constr;
        else
            ctr = NULL;
    }
    return change;
}

/************************************************************
 * VACBinaryConstraint:
 *   A class that stores information about a binary cost function
 */

VACBinaryConstraint::VACBinaryConstraint(WCSP* wcsp, EnumeratedVariable* xx, EnumeratedVariable* yy, vector<Cost>& tab)
    : BinaryConstraint(wcsp, xx, yy, tab)
    , myThreshold(MIN_COST)
{
    for (unsigned int a = 0; a < xx->getDomainInitSize(); a++) {
        kX.push_back(0);
        kX_timeStamp.push_back(0);
    }
    for (unsigned int b = 0; b < yy->getDomainInitSize(); b++) {
        kY.push_back(0);
        kY_timeStamp.push_back(0);
    }
}

VACBinaryConstraint::VACBinaryConstraint(WCSP* wcsp)
    : BinaryConstraint(wcsp)
    , myThreshold(MIN_COST)
{
}

VACBinaryConstraint::~VACBinaryConstraint()
{
}

void VACBinaryConstraint::VACproject(VACVariable* x, Value v, Cost c)
{
    assert(ToulBar2::verbose < 4 || ((cout << "project(C" << getVar(0)->getName() << "," << getVar(1)->getName() << ", (" << x->getName() << "," << v << "), " << c << ")" << endl), true));
    wcsp->revise(this);
    TreeDecomposition* td = wcsp->getTreeDec();
    if (td)
        td->addDelta(cluster, x, v, c);

    int index = x->toIndex(v);
    // TO BE REPLACED BY A LOOP ON THE DOMAIN IN ORDER TO AVOID SUBTRACTING TOP???
    if (!getIndex(x))
        deltaCostsX[index] += c;
    else
        deltaCostsY[index] += c;
    assert(getCost((EnumeratedVariable*)x, (EnumeratedVariable*)getVarDiffFrom(x), v, getVarDiffFrom(x)->getInf()) >= MIN_COST);
    assert(getCost((EnumeratedVariable*)x, (EnumeratedVariable*)getVarDiffFrom(x), v, getVarDiffFrom(x)->getSup()) >= MIN_COST);
    x->VACproject(v, c);
}

void VACBinaryConstraint::VACextend(VACVariable* x, Value v, Cost c)
{
    assert(ToulBar2::verbose < 4 || ((cout << "extend(C" << getVar(0)->getName() << "," << getVar(1)->getName() << ", (" << x->getName() << "," << v << "), " << c << ")" << endl), true));

    TreeDecomposition* td = wcsp->getTreeDec();
    if (td)
        td->addDelta(cluster, x, v, -c);

    int index = x->toIndex(v);
    // TO BE REPLACED BY A LOOP ON THE DOMAIN IN ORDER TO AVOID SUBTRACTING TOP???
    if (!getIndex(x))
        deltaCostsX[index] -= c;
    else
        deltaCostsY[index] -= c;
    x->VACextend(v, c);
}

bool VACBinaryConstraint::revise(VACVariable* var, Value v)
{
    wcsp->revise(this);
    VACVariable* xi = (VACVariable*)getVar(0);
    VACVariable* xj = (VACVariable*)getVar(1);
    Value sup = getSupport(var, v);
    if (var != xi) {
        xi = (VACVariable*)getVar(1);
        xj = (VACVariable*)getVar(0);
    }

    if (xj->canbe(sup)) {
        if (xj->getVACCost(sup) > MIN_COST) {
            xj->removeVAC(sup);
        } else {
            if (getVACCost(xi, xj, v, sup) == MIN_COST) {
                return false;
            }
        }
    }

    for (EnumeratedVariable::iterator it = xj->lower_bound(sup + 1); it != xj->end(); ++it) {
        Value w = *it;
        if (xj->getVACCost(w) > MIN_COST) {
            xj->removeVAC(w);
        } else {
            if (getVACCost(xi, xj, v, w) == MIN_COST) {
                setSupport(xi, v, w);
                return false;
            }
        }
    }

#ifndef AC2001
    for (EnumeratedVariable::iterator it = xj->upper_bound(sup - 1); it != xj->rend(); --it) {
        Value w = *it;
        if (xj->getVACCost(w) > MIN_COST) {
            xj->removeVAC(w);
        } else {
            if (getVACCost(xi, xj, v, w) == MIN_COST) {
                setSupport(xi, v, w);
                return false;
            }
        }
    }
#endif

    return true;
}

/************************************************************
 * VACTernaryConstraint:
 *   A class that stores information about a ternary cost function
 */
/*
VACTernaryConstraint::VACTernaryConstraint (WCSP *wcsp, EnumeratedVariable *xx, EnumeratedVariable *yy, EnumeratedVariable *zz, BinaryConstraint *xy, BinaryConstraint *xz, BinaryConstraint *yz, vector<Cost> &tab) :  TernaryConstraint(wcsp, xx, yy, zz, xy, xz, yz, tab, storeCost), myThreshold(MIN_COST)
{
   for (int a = 0; a < xx->getDomainInitSize(); a++) {
                kX.push_back(0);
                kX_timeStamp.push_back(0);
   }
   for (int b = 0; b < yy->getDomainInitSize(); b++) {
                kY.push_back(0);
                kY_timeStamp.push_back(0);
   }
   for (int c = 0; c < zz->getDomainInitSize(); c++) {
                kZ.push_back(0);
                kZ_timeStamp.push_back(0);
   }
}

VACTernaryConstraint::VACTernaryConstraint (WCSP *wcsp) : TernaryConstraint(wcsp) , myThreshold(MIN_COST)
{
   for (int a = 0; a < wcsp->getMaxDomainSize(); a++) {
                kX.push_back(0);
                kX_timeStamp.push_back(0);
   }
   for (int b = 0; b < wcsp->getMaxDomainSize(); b++) {
                kY.push_back(0);
                kY_timeStamp.push_back(0);
   }
   for (int c = 0; c < wcsp->getDomainInitSize(); c++) {
                kZ.push_back(0);
                kZ_timeStamp.push_back(0);
   }
}

VACTernaryConstraint::~VACTernaryConstraint ()
{
}

void VACTernaryConstraint::VACproject (VACVariable* x, Value v, Cost c) {
  assert(ToulBar2::verbose < 4 || ((cout << "project(C" << getVar(0)->getName() << "," << getVar(1)->getName() << "," << getVar(2)->getName() << ", (" << x->getName() << "," << v << "), " << c << ")" << endl), true));

  TreeDecomposition* td = wcsp->getTreeDec();
  if(td) td->addDelta(cluster,x,v,c);

  int index = x->toIndex(v);
  // TO BE REPLACED BY A LOOP ON THE DOMAIN IN ORDER TO AVOID SUBTRACTING TOP???
  if(getIndex(x)==0) deltaCostsX[index] += c;
  else if(getIndex(x)==1) deltaCostsY[index] += c;
  else deltaCostsZ[index] += c;
  x->VACproject(v, c);
}

void VACTernaryConstraint::VACextend(VACVariable* x, Value v, Cost c) {
  assert(ToulBar2::verbose < 4 || ((cout << "extend(C" << getVar(0)->getName() << "," << getVar(1)->getName() << "," << getVar(2)->getName() << ", (" << x->getName() << "," << v << "), " << c << ")" << endl), true));

  TreeDecomposition* td = wcsp->getTreeDec();
  if(td) td->addDelta(cluster,x,v,-c);

  int index = x->toIndex(v);
  // TO BE REPLACED BY A LOOP ON THE DOMAIN IN ORDER TO AVOID SUBTRACTING TOP???
  if(getIndex(x)==0) deltaCostsX[index] -= c;
  else if(getIndex(x)==1) deltaCostsY[index] -= c;
  else deltaCostsZ[index] -= c;
  x->VACextend(v, c);
}

int VACTernaryConstraint::getK (VACVariable* var, Value v, Long timeStamp) {
  if(var == (VACVariable*) getVar(0)) {
        if(kX_timeStamp[var->toIndex(v)] < timeStamp) return 0;
        else return kX[var->toIndex(v)];
  }  else if(var == (VACVariable*) getVar(1)) {
        if(kY_timeStamp[var->toIndex(v)] < timeStamp) return 0;
        else return kY[var->toIndex(v)];
  }  else {
        if(kZ_timeStamp[var->toIndex(v)] < timeStamp) return 0;
        else return kZ[var->toIndex(v)];
  }
}

void VACTernaryConstraint::setK (VACVariable* var, Value v, int c, Long timeStamp) {
  if(var == getVar(0)) {
        kX[var->toIndex(v)] = c;
        kX_timeStamp[var->toIndex(v)] = timeStamp;
  } else if(var == getVar(1)) {
    kY[var->toIndex(v)] = c;
        kY_timeStamp[var->toIndex(v)] = timeStamp;
  }	 else {
        kZ[var->toIndex(v)] = c;
        kZ_timeStamp[var->toIndex(v)] = timeStamp;
  }
}

bool VACTernaryConstraint::isNull (Cost c)
{
  VACVariable* xi = (VACVariable*) getVar(0);
  return (xi->isSimplyNull(c) || (c < myThreshold));
}

bool VACTernaryConstraint::revise (VACVariable* var, Value v) {
  bool wipeout = false;
  VACVariable* xi = (VACVariable*) getVar(0);
  VACVariable* xj = (VACVariable*) getVar(1);
  VACVariable* xl = (VACVariable*) getVar(2);
  pair<Value,Value> sup = getSupport(var,v);
  pair<Value,Value> minsup = sup;
  if(var != xi) {
          if (var != xj) {
                 xi = (VACVariable*)getVar(2); xj = (VACVariable*)getVar(0); xl = (VACVariable*)getVar(1);
          } else {
                 xi = (VACVariable*)getVar(1); xj = (VACVariable*)getVar(0);
          }
  }
  Cost cost, minCost = wcsp->getUb();

  assert(getindex(xj) < getindex(xl)); // check support is correctly oriented w.r.t. xj/first and xl/second
  if(xj->canbe(sup.first) && xl->canbe(sup.second)) {
          bool unarytest = true;
          if(xj->getVACCost(sup.first) != MIN_COST) { wipeout = xj->removeVAC(sup.first);  unarytest= false;}
          if(xl->getVACCost(sup.second) != MIN_COST) { wipeout = xl->removeVAC(sup.second); unarytest= false;}
          if (unarytest) {
                  if (getVACCost(xi,xj,xl,v,sup.first,sup.second) == MIN_COST) {
                    return false;
                  }
          }
  }

  for (EnumeratedVariable::iterator it = xj->lower_bound(sup); it != xj->end(); ++it) {
          Value w = *it;
          if(xj->getVACCost(w) != MIN_COST) { wipeout = xj->removeVAC(w); xj->queueVAC(); }
          else {
              cost = getVACCost(xi,xj,v, w);
              if (cost == MIN_COST) {
                setSupport(xi,v,w);
                return false;
              } else if (cost < minCost) {
                  minCost = cost;
                  minsup = w;
              }
          }
  }
  for (EnumeratedVariable::iterator it = xj->begin(); it != xj->lower_bound(sup); ++it) {
          Value w = *it;
          if(xj->getVACCost(w) != MIN_COST) { wipeout = xj->removeVAC(w); xj->queueVAC(); }
          else {
              cost = getVACCost(xi,xj,v, w);
              if (cost == MIN_COST) {
                setSupport(xi,v,w);
                return false;
              } else if (cost < minCost) {
                  minCost = cost;
                  minsup = w;
              }
          }
  }

  setSupport(xi,v,minsup);
  return true;
}
 */

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
