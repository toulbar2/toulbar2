/*
 * **************** Abstract constraints of predefined arities **************
 */

#include "tb2abstractconstr.hpp"

/*
 * Constructors and misc.
 *
 */

/// \return size of the cartesian product of all initial domains in the constraint scope.
/// \warning use deprecated MAX_DOMAIN_SIZE for performance.
Long AbstractNaryConstraint::getDomainInitSizeProduct()
{
    if (arity()==0) return 0;
    Long cartesianProduct = 1;
    for (int i=0; i<arity(); i++) {
        // trap overflow numbers
        if (cartesianProduct > LONGLONG_MAX / MAX_DOMAIN_SIZE) return LONGLONG_MAX;
        cartesianProduct *= scope[i]->getDomainInitSize();
    }
    return cartesianProduct;
}

// sorts scope variables by increasing DAC order
int cmpDAC(const void *p1, const void *p2)
{
    EnumeratedVariable *var1 = *((EnumeratedVariable **) p1);
    EnumeratedVariable *var2 = *((EnumeratedVariable **) p2);
    int v1 = var1->getDACOrder();
    int v2 = var2->getDACOrder();
    if (v1 > v2) return 1;
    else if (v1 < v2) return -1;
    else return 0;
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

