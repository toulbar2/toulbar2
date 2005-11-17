/*
 * **************** Storable enumerated domain **********************
 */

#include "tb2domain.hpp"

int cmpValue(const void *v1, const void *v2)
{
    if (*((int *) v1) < *((int *) v2)) return -1;
    else if (*((int *)v1) > *((int *) v2)) return 1;
    else return 0;
}
