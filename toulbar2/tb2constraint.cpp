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
