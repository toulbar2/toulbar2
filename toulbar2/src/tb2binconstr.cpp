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

// Project TRW-S cost to variable
void BinaryConstraint::trwsUpdateMessage(EnumeratedVariable *variable)
{
  assert(connected());
  //if (ToulBar2::verbose >= 3) cout << "\ttrwsUpdateMessage C" << variable->getName() << "," << getVarDiffFrom(variable)->getName() << endl;
  // step 0: be sure of the direction
  EnumeratedVariable *othervar = (EnumeratedVariable *) getVarDiffFrom(variable);
  unsigned int size = sizeY;
  std::function<unsigned int (unsigned int, unsigned int)> getCostIndexX = [size](unsigned int xi, unsigned int yi) { return xi * size + yi; };
  std::function<unsigned int (unsigned int, unsigned int)> getCostIndexY = [size](unsigned int yi, unsigned int xi) { return xi * size + yi; };
  auto getCostIndex = (variable == x)? getCostIndexX: getCostIndexY;
  vector<StoreCost> &thisDeltaCosts = (variable == x)? deltaCostsX: deltaCostsY;
  vector<StoreCost> &thatDeltaCosts = (variable == x)? deltaCostsY: deltaCostsX;
  // step 1: update message
  Cost minMessage = numeric_limits<Cost>::max();
  for (EnumeratedVariable::iterator iterX = variable->begin(); iterX != variable->end(); ++iterX) {
    unsigned int xIndex = variable->toIndex(*iterX);
    Cost minCost        = numeric_limits<Cost>::max();
    //cout << "\t\tvalue: " << xIndex << endl;
    for (EnumeratedVariable::iterator iterY = othervar->begin(); iterY != othervar->end(); ++iterY) {
      unsigned int yIndex = othervar->toIndex(*iterY);
      Cost         cost   = static_cast<Cost>(trunc(othervar->getTRWSGamma() * othervar->getCost(*iterY) - thatDeltaCosts[yIndex] + costs[getCostIndex(xIndex, yIndex)]));
      //cout << "\t\t\tcost: " << cost << " (" << (othervar->getTRWSGamma() * othervar->getCost(*iterY)) << " - " << thatDeltaCosts[yIndex] << " + " << costs[getCostIndex(xIndex, yIndex)] << ")" << endl;
      minCost             = min<double>(cost, minCost);
      if (minCost <= MIN_COST) break;
    }
    //cout << "\t\tmin cost: " << minCost << endl;
    if (thisDeltaCosts[xIndex] != minCost) {
      //cout << "\t\t\tcost move " << thisDeltaCosts[xIndex] << " -> " << minCost << endl;
      variable->project(*iterX, minCost - thisDeltaCosts[xIndex], true);
      thisDeltaCosts[xIndex] = minCost;
    }
    minMessage = min<double>(minCost, minMessage);
  }
  // step 2: normalize messages
  if (minMessage > MIN_COST) {
    for (EnumeratedVariable::iterator iterX = variable->begin(); iterX != variable->end(); ++iterX) {
      thisDeltaCosts[variable->toIndex(*iterX)] -= minMessage;
    }
    //cout << "\t\tprojecting " << minMessage << endl;
    projectLB(minMessage);
  }
}


