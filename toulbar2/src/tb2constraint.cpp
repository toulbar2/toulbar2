/*
 * **************** Abstract constraint **************
 */

#include "tb2constraint.hpp"
#include "tb2wcsp.hpp"
#include "tb2clusters.hpp"

/*
 * Constructor
 *
 */

Constraint::Constraint(WCSP *w) : WCSPLink(w,w->numberOfConstraints()), conflictWeight(1), fromElim1(NULL), fromElim2(NULL)
{
    w->link(this);
    tight = -1;
    isSep_ = false;
    isDuplicate_ = false;
    cluster = -1;
}

Constraint::Constraint(WCSP *w, int elimCtrIndex) : WCSPLink(w,elimCtrIndex), conflictWeight(1), fromElim1(NULL), fromElim2(NULL)
{
    tight = -1;
    isSep_ = false;
    isDuplicate_ = false;
    cluster = -1;
}

Long Constraint::getDomainSizeProduct()
{
  if (arity()==0) return 0;
  Long cartesianProduct = 1;
  for (int i=0; i<arity(); i++) {
	// trap overflow numbers
	if (cartesianProduct > LONGLONG_MAX / MAX_DOMAIN_SIZE) return LONGLONG_MAX;
	cartesianProduct *= getVar(i)->getDomainSize();
  }
  return cartesianProduct;
}

void Constraint::conflict()
{
    wcsp->conflict();
}

void Constraint::projectLB(Cost cost)
{
  if (ToulBar2::verbose >= 2) cout << "lower bound increased " << wcsp->getLb() << " -> " << wcsp->getLb()+cost;
  if (wcsp->td) {
	if (ToulBar2::verbose >= 2) cout << " in cluster C" << getCluster() << " (from " << wcsp->td->getCluster(getCluster())->getLb() << " to " << wcsp->td->getCluster(getCluster())->getLb() + cost << ")";
	wcsp->td->getCluster(getCluster())->increaseLb(cost);
  }
  if (ToulBar2::verbose >= 2) cout << endl;
  wcsp->increaseLb(cost);
}

void Constraint::sumScopeIncluded( Constraint* ctr )
{
	int ar = arity();
	EnumeratedVariable** scopethis = new EnumeratedVariable * [arity()];
	for(int i=0;i<ar;i++) scopethis[i] = (EnumeratedVariable*) getVar(i);

	Cost Top = wcsp->getUb();
	Cost c;
	String t;

	if(getDefCost() < Top) {      // enumeration case
		firstlex();
		while( nextlex(t,c) ) {
			Cost cplus = ctr->evalsubstr(t, this);
			if(c + cplus < Top) addtoTuple( t, cplus, scopethis);
			else setTuple( t, Top, scopethis);
		}
		setDefCost(Top);
	} else {
		first();
		while( next(t,c) ) {
			Cost cplus = ctr->evalsubstr(t, this);
			if(c + cplus < Top) addtoTuple( t, cplus, scopethis);
			else setTuple( t, Top, scopethis);
		}
	}

	delete [] scopethis;
}


void Constraint::assignCluster() {
	TreeDecomposition* td = wcsp->getTreeDec();
	if(!td) return;
	Cluster* lowest = td->getRoot();
	for(int i=0;i<arity();i++) if (getVar(i)->unassigned()) {
		Variable* x = getVar(i);
		Cluster* c = td->getCluster( x->getCluster() );
		if(lowest->isDescendant(c)) lowest = c;
	}
	cluster = lowest->getId();
}

Cost Constraint::getMinCost()
{
// 	    Cost minc = MAX_COST;
//         String tuple;
//         Cost cost;
//         firstlex();
//         while (nextlex(tuple,cost)) {
//             if (cost < minc) minc = cost;
//         }
//         return minc;

  Cost minc = MAX_COST;
  String tuple;
  Cost cost;
  Long nbtuples = 0;
  first();
  while (next(tuple,cost)) {
	nbtuples++;
	if (cost < minc) minc = cost;
  }
  if (getDefCost() < minc && nbtuples < getDomainSizeProduct()) minc = getDefCost();
  return minc;
}

bool Constraint::universal()
{
//   String tuple;
//   Cost cost;
//   firstlex();
//   while (nextlex(tuple,cost)) {
// 	if (cost > MIN_COST) return false;
//   }
//   return true;

  String tuple;
  Cost cost;
  Long nbtuples = 0;
  first();
  while (next(tuple,cost)) {
	nbtuples++;
	if (cost > MIN_COST) return false;
  }
  if (getDefCost() > MIN_COST && nbtuples < getDomainSizeProduct()) return false;
  return true;
}


