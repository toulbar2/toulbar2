/*
 * **************** Abstract constraints of predefined arities **************
 */

#include "tb2abstractconstr.hpp"

/*
 * Constructors and misc.
 * 
 */

AbstractBinaryConstraint::AbstractBinaryConstraint(CostVariable *xx, CostVariable *yy) :
    x(xx), y(yy), linkX(NULL), linkY(NULL)
{
    assert(xx != yy);
    linkX = xx->postConstraint(this,0);
    linkY = yy->postConstraint(this,1);
}

int AbstractBinaryConstraint::getSmallestVarIndexInScope(int forbiddenScopeIndex)
{
    return (forbiddenScopeIndex)?x->wcspIndex:y->wcspIndex;
}
