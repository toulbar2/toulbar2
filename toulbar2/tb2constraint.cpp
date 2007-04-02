/*
 * **************** Abstract constraint **************
 */

#include "tb2constraint.hpp"
#include "tb2wcsp.hpp"

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


void Constraint::sumScopeIncluded( Constraint* ctr ) 
{
	int ar = arity();
	EnumeratedVariable** scopethis = new EnumeratedVariable * [arity()];
	for(int i=0;i<ar;i++) scopethis[i] = (EnumeratedVariable*) getVar(i);
	
	Cost Top = wcsp->getUb();
	string t;
	Cost c;
			
	first();
	while( next(t,c) ) {
		Cost cplus = ctr->evalsubstr(t, this);
		if(c + cplus < Top) setTuple( t, c + cplus, scopethis);
		else setTuple( t, Top, scopethis);
	}
	
	delete [] scopethis;
}
