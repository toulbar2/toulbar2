/*
 * ****** Binary constraint applied on variables with enumerated domains ******
 */

#include <functional>
#include "tb2binconstr.hpp"
#include "tb2wcsp.hpp"
#include "tb2clusters.hpp"

/*
 * Constructors and misc.
 *
 */

BinaryConstraint::BinaryConstraint(WCSP *wcsp, EnumeratedVariable *xx, EnumeratedVariable *yy, vector<Cost> &tab, StoreStack<Cost, Cost> *storeCost) :
        AbstractBinaryConstraint<EnumeratedVariable,EnumeratedVariable>(wcsp, xx, yy), sizeX(xx->getDomainInitSize()), sizeY(yy->getDomainInitSize())
{
    deltaCostsX = vector<StoreCost>(sizeX,StoreCost(MIN_COST,storeCost));
    deltaCostsY = vector<StoreCost>(sizeY,StoreCost(MIN_COST,storeCost));
    assert(tab.size() == sizeX * sizeY);
    supportX = vector<Value>(sizeX,y->getInf());
    supportY = vector<Value>(sizeY,x->getInf());

	costs = vector<StoreCost>(sizeX*sizeY,StoreCost(MIN_COST,storeCost));

    for (unsigned int a = 0; a < x->getDomainInitSize(); a++)
         for (unsigned int b = 0; b < y->getDomainInitSize(); b++)
                costs[a * sizeY + b] = tab[a * sizeY + b];

    propagate();
}

BinaryConstraint::BinaryConstraint(WCSP *wcsp, StoreStack<Cost, Cost> *storeCost)
 	: AbstractBinaryConstraint<EnumeratedVariable,EnumeratedVariable>(wcsp), sizeX(0), sizeY(0)
{
//	unsigned int maxdomainsize = wcsp->getMaxDomainSize();
//    deltaCostsX = vector<StoreCost>(maxdomainsize,StoreCost(MIN_COST,storeCost));
//    deltaCostsY = vector<StoreCost>(maxdomainsize,StoreCost(MIN_COST,storeCost));
//    supportX = vector<Value>(maxdomainsize,0);
//    supportY = vector<Value>(maxdomainsize,0);
    linkX = new DLink<ConstraintLink>;
    linkY = new DLink<ConstraintLink>;

//    costs = vector<StoreCost>(maxdomainsize*maxdomainsize,StoreCost(MIN_COST,storeCost));
//    for (unsigned int a = 0; a < maxdomainsize; a++)
//         for (unsigned int b = 0; b < maxdomainsize; b++)
//                costs[a * maxdomainsize + b] = MIN_COST;
}

void BinaryConstraint::print(ostream& os)
{
    os << this << " BinaryConstraint(" << x->getName() << "," << y->getName() << ")";
	if (ToulBar2::weightedDegree) os << "/" << getConflictWeight();
    if(wcsp->getTreeDec()) os << "   cluster: " << getCluster() << endl;
    else os << endl;
    if (ToulBar2::verbose >= 5) {
        for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
            for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                os << " " << getCost(*iterX, *iterY);
            }
            os << endl;
        }
    }
}

void BinaryConstraint::dump(ostream& os, bool original)
{
  os << "2 " << ((original)?(x->wcspIndex):x->getCurrentVarId()) << " " << ((original)?(y->wcspIndex):y->getCurrentVarId()) << " " << MIN_COST << " " << x->getDomainSize() * y->getDomainSize() << endl;
  int i=0;
  for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX, i++) {
	int j=0;
	for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY, j++) {
	  os << ((original)?(*iterX):i) << " " << ((original)?(*iterY):j) << " " << ((original)?getCost(*iterX, *iterY):min(wcsp->getUb(),getCost(*iterX, *iterY))) << endl;
	}
  }
}

/*
 * Propagation methods
 *
 */
bool BinaryConstraint::project(EnumeratedVariable *x, Value value, Cost cost, vector<StoreCost> &deltaCostsX)
{
	assert(ToulBar2::verbose < 4 || ((cout << "project(C" << getVar(0)->getName() << "," << getVar(1)->getName() << ", (" << x->getName() << "," << value << "), " << cost << ")" << endl), true));
    
    // hard binary constraint costs are not changed
    if (!CUT(cost + wcsp->getLb(), wcsp->getUb())) {
	    TreeDecomposition* td = wcsp->getTreeDec();
	    if(td) td->addDelta(cluster,x,value,cost);
    	deltaCostsX[x->toIndex(value)] += cost;  // Warning! Possible overflow???
    }
    	
    Cost oldcost = x->getCost(value);
    x->project(value, cost);
#ifdef DEECOMPLETE
    getVarDiffFrom(x)->queueDEE();
#endif
    return (x->getSupport() == value || SUPPORTTEST(oldcost, cost));
}


void BinaryConstraint::extend(EnumeratedVariable *x, Value value, Cost cost, vector<StoreCost> &deltaCostsX)
{
	assert(ToulBar2::verbose < 4 || ((cout << "extend(C" << getVar(0)->getName() << "," << getVar(1)->getName() << ", (" << x->getName() << "," << value << "), " << cost << ")" << endl), true));

    TreeDecomposition* td = wcsp->getTreeDec();
    if(td) td->addDelta(cluster,x,value,-cost);

    deltaCostsX[x->toIndex(value)] -= cost;  // Warning! Possible overflow???
    x->extend(value, cost);
}

void BinaryConstraint::permute(EnumeratedVariable *xin, Value a, Value b)
{
	EnumeratedVariable *yin = y;
	if(xin != x) yin = x;

	vector<Cost> aux;
    for (EnumeratedVariable::iterator ity = yin->begin(); ity != yin->end(); ++ity)  aux.push_back( getCost(xin, yin, a, *ity) );
    for (EnumeratedVariable::iterator ity = yin->begin(); ity != yin->end(); ++ity)  setcost(xin, yin, a, *ity, getCost(xin, yin, b, *ity));

	vector<Cost>::iterator itc = aux.begin();
    for (EnumeratedVariable::iterator ity = yin->begin(); ity != yin->end(); ++ity) {
    	setcost(xin, yin, b, *ity, *itc);
    	++itc;
    }
}

bool BinaryConstraint::isFunctional(EnumeratedVariable* xin, EnumeratedVariable* yin, map<Value, Value> &functional)
{
  assert(xin != yin);
  assert(getIndex(xin) >= 0);
  assert(getIndex(yin) >= 0);
  bool isfunctional = true;
  functional.clear();
  for (EnumeratedVariable::iterator itx = xin->begin(); isfunctional && itx != xin->end(); ++itx) {
	bool first = true;
	for (EnumeratedVariable::iterator ity = yin->begin(); isfunctional && ity != yin->end(); ++ity) {
	  if (!CUT(getCost(xin,yin,*itx,*ity) + wcsp->getLb(), wcsp->getUb())) {
		if (first) {
		  functional[*itx] = *ity;
		  first = false;
		} else {
		  isfunctional = false;
		}
	  }
	}
	assert(!first); // assumes it is SAC already
  }
  if (isfunctional) return true;
  functional.clear();
  return false;
}

pair< pair<Cost,Cost>, pair<Cost,Cost> > BinaryConstraint::getMaxCost(int varIndex, Value a, Value b)
{
//    	cout << "getMaxCost(" << getVar(varIndex)->getName() << ") " << a << " <-> " << b << endl << *this << endl;
	Cost maxcosta = MIN_COST;
	Cost diffcosta = MIN_COST;
	Cost maxcostb = MIN_COST;
	Cost diffcostb = MIN_COST;
	if (varIndex == 0) {
		Cost ucosta = x->getCost(a);
		Cost ucostb = x->getCost(b);
		for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
			Cost costa = getCost(a, *iterY);
			Cost costb = getCost(b, *iterY);
			if (costa > maxcosta) maxcosta = costa;
			if (costb > maxcostb) maxcostb = costb;
			Cost ucosty = y->getCost(*iterY);
			if (!CUT(ucostb + costb + ucosty + wcsp->getLb(), wcsp->getUb())) {
				if (costa-costb > diffcosta) diffcosta = costa-costb;
			}
			if (!CUT(ucosta + costa + ucosty + wcsp->getLb(), wcsp->getUb())) {
				if (costb-costa > diffcostb) diffcostb = costb-costa;
			}
		}
	} else {
		Cost ucosta = y->getCost(a);
		Cost ucostb = y->getCost(b);
		for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
   			Cost costa = getCost(*iterX, a);
   			Cost costb = getCost(*iterX, b);
			if (costa > maxcosta) maxcosta = costa;
			if (costb > maxcostb) maxcostb = costb;
			Cost ucostx = x->getCost(*iterX);
			if (!CUT(ucostb + costb + ucostx + wcsp->getLb(), wcsp->getUb())) {
				if (costa-costb > diffcosta) diffcosta = costa-costb;
			}
			if (!CUT(ucosta + costa + ucostx + wcsp->getLb(), wcsp->getUb())) {
				if (costb-costa > diffcostb) diffcostb = costb-costa;
			}
		}
	}
	assert(maxcosta >= diffcosta);
	assert(maxcostb >= diffcostb);
	return make_pair(make_pair(maxcosta,diffcosta), make_pair(maxcostb,diffcostb));
}

void BinaryConstraint::preprocessTRWS () {
  for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
    unsigned int indexX = x->toIndex(*iterX);
		for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
      unsigned int indexY = y->toIndex(*iterY);
      costs[indexX * sizeY + indexY] -= deltaCostsX[indexX] + deltaCostsY[indexY];
    }
  }
  for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
    deltaCostsX[x->toIndex(*iterX)] = 0;
  }
  for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
    deltaCostsY[y->toIndex(*iterY)] = 0;
  }
}

Cost BinaryConstraint::getCostTRWS (Value v1, Value v2) {
  unsigned int indexX = x->toIndex(v1);
  unsigned int indexY = y->toIndex(v2);
  return costs[indexX * sizeY + indexY] - (deltaCostsX[indexX] + deltaCostsY[indexY]);
}

void BinaryConstraint::projectTRWS (EnumeratedVariable *var, Value value, Cost cost) {
  vector<StoreCost> &deltaCosts = (var == x)? deltaCostsX: deltaCostsY;
  deltaCosts[var->toIndex(value)] += cost;
  var->project(value, cost);
}

Cost BinaryConstraint::projectTRWS (EnumeratedVariable *s, EnumeratedVariable *t, Cost deltaUnary, bool last) {
  unsigned int size = sizeY;
  std::function<unsigned int (unsigned int, unsigned int)> getCostIndexX = [size](unsigned int xi, unsigned int yi) { return xi * size + yi; };
  std::function<unsigned int (unsigned int, unsigned int)> getCostIndexY = [size](unsigned int yi, unsigned int xi) { return xi * size + yi; };
  auto getCostIndex = (t == x)? getCostIndexX: getCostIndexY;
  vector<StoreCost> &tDeltaCosts = (t == x)? deltaCostsX: deltaCostsY;
  vector<StoreCost> &sDeltaCosts = (t == x)? deltaCostsY: deltaCostsX;
  // step 1: update message
  Cost delta = numeric_limits<Cost>::max();
  for (EnumeratedVariable::iterator tIter = t->begin(); tIter != t->end(); ++tIter) {
    unsigned int tIndex  = t->toIndex(*tIter);
    Cost         minCost = numeric_limits<Cost>::max();
    for (EnumeratedVariable::iterator sIter = s->begin(); sIter != s->end(); ++sIter) {
      unsigned int sIndex = s->toIndex(*sIter);
      Cost         cost   = static_cast<Cost>(trunc(s->getTRWSGamma() * (s->getCost(*sIter)-deltaUnary))) - sDeltaCosts[sIndex] + costs[getCostIndex(tIndex, sIndex)];
      minCost             = min<Cost>(cost, minCost);
      //cout << "\t\t\tcost: " << cost << " (" << (s->getTRWSGamma() * s->getCost(*sIter)) << " - " << sDeltaCosts[sIndex] << " + " << costs[getCostIndex(tIndex, sIndex)] << ")" << endl;
    }
    //cout << "\t\tmin cost: " << minCost << endl;
    if (minCost > tDeltaCosts[tIndex]) t->project(*tIter, minCost - tDeltaCosts[tIndex], true);
   else                               t->extend(*tIter, tDeltaCosts[tIndex] - minCost);
    tDeltaCosts[tIndex] = minCost;
    delta = min<Cost>(delta, minCost);
  }
  if (! last) {
    // step 2: normalize messages
    for (EnumeratedVariable::iterator tIter = t->begin(); tIter != t->end(); ++tIter) {
      tDeltaCosts[t->toIndex(*tIter)] -= delta;
    }
  }
  //cout << "\t\tprojecting " << minMessage << endl;
  return delta;
}

Cost BinaryConstraint::normalizeTRWS () {
  Cost minCost = numeric_limits<Cost>::max();
  for (EnumeratedVariable::iterator xIter = x->begin(); xIter != x->end(); ++xIter) {
    for (EnumeratedVariable::iterator yIter = y->begin(); yIter != y->end(); ++yIter) {
      minCost = min<Cost>(minCost, getCostTRWS(*xIter, *yIter));
    }
  }
  for (EnumeratedVariable::iterator xIter = x->begin(); xIter != x->end(); ++xIter) {
    unsigned int ix = x->toIndex(*xIter);
    for (EnumeratedVariable::iterator yIter = y->begin(); yIter != y->end(); ++yIter) {
      unsigned int iy = y->toIndex(*yIter);
      costs[ix * sizeY + iy] -= minCost;
    }
  }
  return minCost;
}

void BinaryConstraint::addDeltaTRWS (vector < Cost > &delta, Cost minCost, unsigned int side) {
  preprocessTRWS();
  vector<StoreCost> &deltaCosts = (side == 0)? deltaCostsX: deltaCostsY;
  EnumeratedVariable *var = dynamic_cast<EnumeratedVariable *>(getVar(side));
  for (EnumeratedVariable::iterator iter = var->begin(); iter != var->end(); ++iter) {
    unsigned int i = var->toIndex(*iter);
    if (delta[i] != 0) cout << "\t\tC" << x->wcspIndex << "," << y->wcspIndex << "[" << i << "] (side " << side << ") += " << delta[i] << "\n";
    deltaCosts[i] += delta[i];
  }
  for (EnumeratedVariable::iterator xIter = x->begin(); xIter != x->end(); ++xIter) {
    unsigned int ix = x->toIndex(*xIter);
    for (EnumeratedVariable::iterator yIter = y->begin(); yIter != y->end(); ++yIter) {
      unsigned int iy = y->toIndex(*yIter);
      costs[ix * sizeY + iy] -= minCost;
    }
  }
  for (EnumeratedVariable::iterator xIter = x->begin(); xIter != x->end(); ++xIter) {
    unsigned int ix = x->toIndex(*xIter);
    for (EnumeratedVariable::iterator yIter = y->begin(); yIter != y->end(); ++yIter) {
      unsigned int iy = y->toIndex(*yIter);
      if (costs[ix * sizeY + iy] - (deltaCostsX[ix] + deltaCostsY[iy]) < 0) cout << "\t\tProblem!  C" << x->wcspIndex << "," << y->wcspIndex << "[" << ix << ", " << iy << "] = " << costs[ix * sizeY + iy] << " - (" << deltaCostsX[ix] << " + " << deltaCostsY[iy] << ") = " << (costs[ix * sizeY + iy] - (deltaCostsX[ix] + deltaCostsY[iy])) << " < 0\n";
    }
  }
} 

/*
Cost BinaryConstraint::addDeltaTRWS (std::array < vector < Cost >, 2 > &delta) {
  for (int side = 0; side < 2; ++side) {
    vector<StoreCost> &deltaCosts = (side == 0)? deltaCostsX: deltaCostsY;
    EnumeratedVariable *var = dynamic_cast<EnumeratedVariable *>(getVar(side));
    for (EnumeratedVariable::iterator iter = var->begin(); iter != var->end(); ++iter) {
      unsigned int i = var->toIndex(*iter);
      if (delta[side][i] != 0) cout << "\t\tC" << x->wcspIndex << "," << y->wcspIndex << "[" << i << "] (side " << side << ") += " << delta[side][i] << "\n";
      deltaCosts[i] += delta[side][i];
    }
  }
  Cost minCost = numeric_limits<Cost>::max();
  for (EnumeratedVariable::iterator xIter = x->begin(); xIter != x->end(); ++xIter) {
    unsigned int ix = x->toIndex(*xIter);
    for (EnumeratedVariable::iterator yIter = y->begin(); yIter != y->end(); ++yIter) {
      unsigned int iy = y->toIndex(*yIter);
      minCost = min<Cost>(minCost, costs[ix * sizeY + iy] - (deltaCostsX[ix] + deltaCostsY[iy]));
    }
  }
  for (EnumeratedVariable::iterator xIter = x->begin(); xIter != x->end(); ++xIter) {
    unsigned int ix = x->toIndex(*xIter);
    for (EnumeratedVariable::iterator yIter = y->begin(); yIter != y->end(); ++yIter) {
      unsigned int iy = y->toIndex(*yIter);
      costs[ix * sizeY + iy] -= minCost;
    }
  }
  return minCost;
} 
*/

// Project TRW-S cost to variable
/*
void BinaryConstraint::trwsUpdateMessage(EnumeratedVariable *t, vector < Cost > &M)
{
  assert(connected());
  if (ToulBar2::verbose >= 3) cout << "\ttrwsUpdateMessage C" << t->getName() << "," << getVarDiffFrom(t)->getName() << endl;
  // step 0: be sure of the direction
  EnumeratedVariable *s = (EnumeratedVariable *) getVarDiffFrom(t);
  unsigned int size = sizeY;
  std::function<unsigned int (unsigned int, unsigned int)> getCostIndexX = [size](unsigned int xi, unsigned int yi) { return xi * size + yi; };
  std::function<unsigned int (unsigned int, unsigned int)> getCostIndexY = [size](unsigned int yi, unsigned int xi) { return xi * size + yi; };
  auto getCostIndex = (t == x)? getCostIndexX: getCostIndexY;
  vector<StoreCost> &thisDeltaCosts = (t == x)? deltaCostsX: deltaCostsY;
  vector<StoreCost> &thatDeltaCosts = (t == x)? deltaCostsY: deltaCostsX;
  // step 1: update message
  Cost minMessage = numeric_limits<Cost>::max();
  for (EnumeratedVariable::iterator thisIter = t->begin(); thisIter != t->end(); ++thisIter) {
    unsigned int thisIndex = t->toIndex(*thisIter);
    Cost         minCost   = numeric_limits<Cost>::max();
    //cout << "\t\tvalue: " << thisIndex << endl;
    for (EnumeratedVariable::iterator thatIter = s->begin(); thatIter != s->end(); ++thatIter) {
      unsigned int thatIndex = s->toIndex(*thatIter);
      Cost         cost      = static_cast<Cost>(trunc(s->getTRWSGamma() * s->getCost(*thatIter))) - thatDeltaCosts[thatIndex] + costs[getCostIndex(thisIndex, thatIndex)];
      minCost                = min<Cost>(cost, minCost);
      //cout << "\t\t\tcost: " << cost << " (" << (s->getTRWSGamma() * s->getCost(*thatIter)) << " - " << thatDeltaCosts[thatIndex] << " + " << costs[getCostIndex(thisIndex, thatIndex)] << ")" << endl;
      if (minCost <= MIN_COST) break;
    }
    //cout << "\t\tmin cost: " << minCost << endl;
    if (minCost > MIN_COST) {
    //if ((minCost > MIN_COST) && (minCost > thisDeltaCosts[thisIndex])) {
      //cout << "\t\t\tcost move " << thisDeltaCosts[thisIndex] << " -> " << minCost << endl;
      t->project(*thisIter, minCost - thisDeltaCosts[thisIndex], true);
      thisDeltaCosts[thisIndex] = minCost;
    }
    minMessage = min<Cost>(minCost, minMessage);
  }
  // step 2: normalize messages
  if (minMessage > MIN_COST) {
    for (EnumeratedVariable::iterator thisIter = t->begin(); thisIter != t->end(); ++thisIter) {
      thisDeltaCosts[t->toIndex(*thisIter)] -= minMessage;
    }
    //cout << "\t\tprojecting " << minMessage << endl;
    t->extendAll(minMessage);
    t->projectLB(minMessage);
  }
}
*/


