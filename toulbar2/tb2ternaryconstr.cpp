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

double TernaryConstraint::computeTightness()
{
   int count = 0;
   double sum = 0;
   for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
      for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
	      for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
			sum += getCost(*iterX, *iterY, *iterZ);
			count++;
       }
     }
   }
   
   tight = sum / (double) count;
   return tight;
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

void TernaryConstraint::projectTernaryBinary( BinaryConstraint* yzin )
{
    bool flag = false;
    EnumeratedVariable* xin = NULL;
    if(yzin->getIndex(x) < 0)      xin = x;
    else if(yzin->getIndex(y) < 0) xin = y;
    else if(yzin->getIndex(z) < 0) xin = z;
    
    EnumeratedVariable* xx = xin;
    EnumeratedVariable* yy = (EnumeratedVariable*) yzin->getVar(0);
    EnumeratedVariable* zz = (EnumeratedVariable*) yzin->getVar(1);
     
    for (EnumeratedVariable::iterator itery = yy->begin(); itery != yy->end(); ++itery) {
    for (EnumeratedVariable::iterator iterz = zz->begin(); iterz != zz->end(); ++iterz) {
        Cost mincost = MAX_COST;
        for (EnumeratedVariable::iterator iterx = xx->begin(); iterx != xx->end(); ++iterx) {
            Cost curcost = getCost(xx, yy, zz, *iterx, *itery, *iterz);                            
            if (curcost < mincost) mincost = curcost;
        }
        assert(mincost != MAX_COST);
        if (mincost > 0) {
            flag = true;
            yzin->addcost(*itery,*iterz,mincost);   // project mincost to binary
//            if (ToulBar2::verbose >= 2) cout << "ternary projection of " << mincost << " from C" << xx->getName() << "," << yy->getName() << "," << zz->getName() << " to C" << yy->getName() << "," << zz->getName() << "(" << *itery << "," << *iterz << ")" << endl;
            
            if (connected() && mincost < wcsp->getUb()) {  // substract mincost from ternary constraint
                for (EnumeratedVariable::iterator iterx = xx->begin(); iterx != xx->end(); ++iterx) {
                    addcost(xx, yy, zz, *iterx, *itery, *iterz, -mincost);                            
                }
            }               
        }
    }}
    
    if (flag) {
        if(!yzin->connected()) yzin->reconnect();
        yzin->propagate();
//      yy->queueAC();  yy->queueDAC();
//      zz->queueAC();  zz->queueDAC();
    }
}
