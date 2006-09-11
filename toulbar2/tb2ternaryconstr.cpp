/*
 * ****** Binary constraint applied on variables with enumerated domains ******
 */

#include "tb2ternaryconstr.hpp"
#include "tb2wcsp.hpp"


/*
 * Constructors and misc.
 * 
 */

TernaryConstraint::TernaryConstraint(WCSP *wcsp, 
									 EnumeratedVariable *xx, 
									 EnumeratedVariable *yy,
									 EnumeratedVariable *zz,
									 vector<Cost> &tab, 
									 StoreStack<Cost, Cost> *storeCost) 
									 
					: AbstractTernaryConstraint<EnumeratedVariable,EnumeratedVariable,EnumeratedVariable>(wcsp, xx, yy, zz), 
					  sizeX(xx->getDomainInitSize()), sizeY(yy->getDomainInitSize()), sizeZ(zz->getDomainInitSize()), costs(tab)
{
}


void TernaryConstraint::propagate()
{
	/* if(x->assigned() && y->assigned() && z->assigned()) {
   		deconnect();
   		wcsp->increaseLb( getCost(x->getValue(), y->getValue(), z->getValue()  ) );   
   	}*/
}



void TernaryConstraint::print(ostream& os)
{
    os << this << " TernaryConstraint(" << x->getName() << "," << y->getName() << "," << z->getName() << ")" << endl;
    if (ToulBar2::verbose >= 5) {
        for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
            for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
			   for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
		         os << " " << getCost(*iterX, *iterY, *iterZ);
			   }
            }
            os << endl;
        }
    }
}




