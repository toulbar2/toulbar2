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

Constraint::Constraint(WCSP *w) : WCSPLink(w,w->numberOfConstraints()), conflictWeight(1)
{
    w->link(this);
    tight = -1;
}

Constraint::Constraint(WCSP *w, int elimCtrIndex) : WCSPLink(w,elimCtrIndex), conflictWeight(1)
{
    tight = -1;
}

void Constraint::projectLB(Cost cost)
{
  if (ToulBar2::verbose >= 2) cout << "lower bound increased " << wcsp->getLb() << " -> " << wcsp->getLb()+cost << endl;
  if (wcsp->vac) {
	wcsp->td->getCluster(getCluster())->increaseLb(wcsp->td->getCluster(getCluster())->getLb() + cost);
  }
  wcsp->increaseLb(wcsp->getLb() + cost);  
}

void Constraint::sumScopeIncluded( Constraint* ctr ) 
{
	int ar = arity();
	EnumeratedVariable** scopethis = new EnumeratedVariable * [arity()];
	for(int i=0;i<ar;i++) scopethis[i] = (EnumeratedVariable*) getVar(i);
	
	Cost Top = wcsp->getUb();
	Cost c;
	string t;

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
