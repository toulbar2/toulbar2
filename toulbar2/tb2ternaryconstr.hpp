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

    void addcost( EnumeratedVariable* xin, EnumeratedVariable* yin, EnumeratedVariable* zin, int vxi, int vyi, int vzi, Cost c ) {

        int vindex[3];
	    int vx = xin->toIndex(vxi);
	    int vy = yin->toIndex(vyi);
	    int vz = zin->toIndex(vzi);

        vindex[ getIndex(xin) ] = vx;
        vindex[ getIndex(yin) ] = vy;
        vindex[ getIndex(zin) ] = vz;
	        
		costs[vindex[0]*sizeX*sizeY + vindex[1]*sizeY + vindex[2]] += c;
    }


    void propagate();
     
    void remove(int varIndex) {}

    void projectFromZero(int varIndex) {} 

    void increase(int index) {}
    
	void decrease(int index) {}
	  
	void projectTernary() {
		projectTernaryBinary(xy); //if(!xy->connected()) xy->reconnect();
		projectTernaryBinary(xz); //if(!xz->connected()) xz->reconnect();
		projectTernaryBinary(yz); //if(!yz->connected()) yz->reconnect();
		
		//xy->propagate();
		//xz->propagate();
		//yz->propagate();
	}

    void projectTernaryBinary( BinaryConstraint* yzin );

	void extendTernary()
	{
		Cost c; 
		for(EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
		for(EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
		for(EnumeratedVariable::iterator iterz = z->begin(); iterz != z->end(); ++iterz) {
			if(xy->connected()) { c = xy->getCost(*iterx, *itery); addcost(x,y,z,*iterx,*itery,*iterz, c); xy->addcost(*iterx, *itery, -c); xy->deconnect(); }
			if(xz->connected()) { c = xz->getCost(*iterx, *itery); addcost(x,y,z,*iterx,*itery,*iterz, c); xz->addcost(*iterx, *iterz, -c); xz->deconnect(); }
			if(yz->connected()) { c = yz->getCost(*iterx, *itery); addcost(x,y,z,*iterx,*itery,*iterz, c); yz->addcost(*itery, *iterz, -c); yz->deconnect(); }
		}}}
	}

	BinaryConstraint* commonBinary( TernaryConstraint* t )
	{
		if( (t->getIndex(xy->getVar(0)) >= 0) && (t->getIndex(xy->getVar(1)) >= 0) ) return xy;
		else if((t->getIndex(xz->getVar(0)) >= 0) && (t->getIndex(xz->getVar(1)) >= 0)) return xz;
		else if((t->getIndex(yz->getVar(0)) >= 0) && (t->getIndex(yz->getVar(1)) >= 0)) return yz;
		return NULL;
	}
    
	void assign(int varIndex) {
	    deconnect();                    
		switch(varIndex) {
			case 0: projectTernaryBinary(yz);  break;
			case 1: projectTernaryBinary(xz);  break;
			case 2: projectTernaryBinary(xy);  break;
			default:;
		}
    }
        
    bool verify() { return true; }

    double computeTightness();

    void print(ostream& os);
};

#endif /*TB2TERNARYCONSTR_HPP_*/
