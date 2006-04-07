/*
 * **************** Abstract constraint **************
 */

#include "tb2constraint.hpp"
#include "tb2wcsp.hpp"

/*
 * Constructor
 * 
 */

Constraint::Constraint(WCSP *w) : WCSPLink(w,w->numberOfConstraints())
{
    w->link(this);
}
