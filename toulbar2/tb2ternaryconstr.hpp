/** \file tb2ternaryconstr.hpp
 *  \brief Ternary constraint applied on variables with enumerated domains.
 *
 */

#ifndef TB2TERNARYCONSTR_HPP_
#define TB2TERNARYCONSTR_HPP_


#include "tb2abstractconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2binconstr.hpp"

//class TernaryConstraint;
//typedef Cost (TernaryConstraint::*GetCostMember)(Value vx, Value vy, Value vz);

class TernaryConstraint : public AbstractTernaryConstraint<EnumeratedVariable,EnumeratedVariable,EnumeratedVariable>
{
    unsigned int sizeX;
    unsigned int sizeY;
    unsigned int sizeZ;
    vector<Cost> costs;
    
   
public:
    TernaryConstraint(WCSP *wcsp, 
					  EnumeratedVariable *xx, 
					  EnumeratedVariable *yy, 
					  EnumeratedVariable *zz, 
					  vector<Cost> &tab, 
					  StoreStack<Cost, Cost> *storeCost);

    void setBinaries( BinaryConstraint* xyin, BinaryConstraint* xzin, BinaryConstraint* yzin ) { xy = xyin; xz = xzin; yz = yzin; }

    BinaryConstraint* xy;
    BinaryConstraint* xz;
    BinaryConstraint* yz;

	~TernaryConstraint() {}
    
    Cost getCost(Value vx, Value vy, Value vz) {
        int ix = x->toIndex(vx);
        int iy = y->toIndex(vy);
        int iz = z->toIndex(vz);
        Cost res = costs[ix*sizeX*sizeY + iy*sizeY + iz];
        assert(res >= 0);
        return res;
    }

    Cost getCost(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz) {  
        int vindex[3];
        vindex[ getIndex(xx) ] = xx->toIndex(vx);
        vindex[ getIndex(yy) ] = yy->toIndex(vy);
        vindex[ getIndex(zz) ] = zz->toIndex(vz);
      
        Cost res = costs[vindex[0]*sizeX*sizeY + vindex[1]*sizeY + vindex[2]];
        assert(res >= 0);
        return res;
    }

    void addCosts( EnumeratedVariable* xin, EnumeratedVariable* yin, EnumeratedVariable* zin, vector<Cost>& costsin ) {
		assert(costsin.size() == costs.size());

        int vindex[3];
        int sizeXin = xin->getDomainInitSize();
        int sizeYin = yin->getDomainInitSize();
       
		for (EnumeratedVariable::iterator iterx = xin->begin(); iterx != xin->end(); ++iterx) {
		for (EnumeratedVariable::iterator itery = yin->begin(); itery != yin->end(); ++itery) {
		for (EnumeratedVariable::iterator iterz = zin->begin(); iterz != zin->end(); ++iterz) {

		    int vxin = xin->toIndex(*iterx);
		    int vyin = yin->toIndex(*itery);
		    int vzin = zin->toIndex(*iterz);

	        vindex[ getIndex(xin) ] = vxin;
	        vindex[ getIndex(yin) ] = vyin;
	        vindex[ getIndex(zin) ] = vzin;
	        
			costs[vindex[0]*sizeX*sizeY + vindex[1]*sizeY + vindex[2]] += costsin[vxin*sizeXin*sizeYin + vyin*sizeYin + vzin];
	    }}}
    }


    void propagate();
     
    void remove(int varIndex) {}

    void projectFromZero(int varIndex) {} 

    void increase(int index) {}
    
	void decrease(int index) {}  

	void projectTernaryBinary( EnumeratedVariable* xin, BinaryConstraint* yzin )
	{
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
			yzin->addcost(*itery,*iterz,mincost);   
		}}
		

		if(!yzin->connected()) yzin->reconnect();
		yzin->propagate();
		yy->queueAC();  yy->queueDAC();
    	zz->queueAC();  zz->queueDAC();
		
	}
    
	void assign(int varIndex) {
	    deconnect();                    
		switch(varIndex) {
			case 0: projectTernaryBinary(x, yz);  break;
			case 1: projectTernaryBinary(y, xz);  break;
			case 2: projectTernaryBinary(z, xy);  break;
			default:;
		}
    }
        
    bool verify() { return true; }
    
    void print(ostream& os);


};

#endif /*TB2TERNARYCONSTR_HPP_*/
