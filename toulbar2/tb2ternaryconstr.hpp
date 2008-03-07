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
    vector<StoreCost> costs;
    vector<StoreCost> deltaCostsX;
    vector<StoreCost> deltaCostsY;
    vector<StoreCost> deltaCostsZ;
    vector< pair<Value,Value> > supportX;
    vector< pair<Value,Value> > supportY;
    vector< pair<Value,Value> > supportZ;

    // return true if unary support of x is broken
    bool project(EnumeratedVariable *x, Value value, Cost cost, vector<StoreCost> &deltaCostsX)
    {
		assert(ToulBar2::verbose < 1 || ((cout << "project(C" << getVar(0)->getName() << "," << getVar(1)->getName() << "," << getVar(2)->getName() << ", (" << x->getName() << "," << value << "), " << cost << ")" << endl), true));
        // hard binary constraint costs are not changed
        if (!CUT(cost + wcsp->getLb(), wcsp->getUb())) deltaCostsX[x->toIndex(value)] += cost;  // Warning! Possible overflow???
        Cost oldcost = x->getCost(value);
        x->project(value, cost);
        return (x->getSupport() == value || SUPPORTTEST(oldcost, cost));
    }
    
    void extend(EnumeratedVariable *x, Value value, Cost cost, vector<StoreCost> &deltaCostsX)
    {
		assert(ToulBar2::verbose < 4 || ((cout << "extend(C" << getVar(0)->getName() << "," << getVar(1)->getName() << "," << getVar(2)->getName() << ", (" << x->getName() << "," << value << "), " << cost << ")" << endl), true));
        deltaCostsX[x->toIndex(value)] -= cost;  // Warning! Possible overflow???
        x->extend(value, cost);
    }

    void project(BinaryConstraint* xy, EnumeratedVariable *x, EnumeratedVariable *y, EnumeratedVariable *z, Value valx, Value valy, Cost cost) {
		assert(ToulBar2::verbose < 4 || ((cout << "project(C" << getVar(0)->getName() << "," << getVar(1)->getName() << "," << getVar(2)->getName() << ", ((" << x->getName() << "," << valx << "),(" << y->getName() << "," << valy << ")), " << cost << ")" << endl), true));
        for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
            addcost(x,y,z,valx,valy,*iterZ,-cost);
        }
        xy->addcost(x,y,valx,valy,cost);
    }

    void extend(BinaryConstraint* xy, EnumeratedVariable *x, EnumeratedVariable *y, EnumeratedVariable *z, Value valx, Value valy, Cost cost) {
		assert(ToulBar2::verbose < 4 || ((cout << "extend(C" << getVar(0)->getName() << "," << getVar(1)->getName() << "," << getVar(2)->getName() << ", ((" << x->getName() << "," << valx << "),(" << y->getName() << "," << valy << ")), " << cost << ")" << endl), true));
        for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
            addcost(x,y,z,valx,valy,*iterZ,cost);
        }
        assert(xy->connected());
        xy->addcost(x,y,valx,valy,-cost);
    }
    
    void findSupport(EnumeratedVariable *x, EnumeratedVariable *y,  EnumeratedVariable *z,
            vector< pair<Value,Value> > &supportX, vector<StoreCost> &deltaCostsX, 
            vector< pair<Value,Value> > &supportY, vector< pair<Value,Value> > &supportZ);
    void findFullSupport(EnumeratedVariable *x, EnumeratedVariable *y,  EnumeratedVariable *z,
            vector< pair<Value,Value> > &supportX, vector<StoreCost> &deltaCostsX, 
            vector< pair<Value,Value> > &supportY, vector<StoreCost> &deltaCostsY, 
            vector< pair<Value,Value> > &supportZ, vector<StoreCost> &deltaCostsZ,
            BinaryConstraint* xy, BinaryConstraint* xz, BinaryConstraint* yz);
    bool verify(EnumeratedVariable *x, EnumeratedVariable *y, EnumeratedVariable *z);

    bool isEAC(EnumeratedVariable *x, Value a, EnumeratedVariable *y, EnumeratedVariable *z,
               vector< pair<Value,Value> > &supportX);
               
    void findSupportX() {findSupport(x,y,z,supportX,deltaCostsX,supportY,supportZ);}
    void findSupportY() {findSupport(y,x,z,supportY,deltaCostsY,supportX,supportZ);}
    void findSupportZ() {findSupport(z,x,y,supportZ,deltaCostsZ,supportX,supportY);}
    void findFullSupportX() {if (y->wcspIndex < z->wcspIndex) findFullSupport(x,y,z,supportX,deltaCostsX,supportY,deltaCostsY,supportZ,deltaCostsZ,xy,xz,yz);
        else findFullSupport(x,z,y,supportX,deltaCostsX,supportZ,deltaCostsZ,supportY,deltaCostsY,xz,xy,yz);}
    void findFullSupportY() {if (x->wcspIndex < z->wcspIndex) findFullSupport(y,x,z,supportY,deltaCostsY,supportX,deltaCostsX,supportZ,deltaCostsZ,xy,yz,xz);
        else findFullSupport(y,z,x,supportY,deltaCostsY,supportZ,deltaCostsZ,supportX,deltaCostsX,yz,xy,xz);}
    void findFullSupportZ() {if (x->wcspIndex < y->wcspIndex) findFullSupport(z,x,y,supportZ,deltaCostsZ,supportX,deltaCostsX,supportY,deltaCostsY,xz,yz,xy);
        else findFullSupport(z,y,x,supportZ,deltaCostsZ,supportY,deltaCostsY,supportX,deltaCostsX,yz,xz,xy);}
    bool verifyX() {return verify(x,y,z);}
    bool verifyY() {return verify(y,x,z);}
    bool verifyZ() {return verify(z,x,y);}
    
public:
    TernaryConstraint(WCSP *wcsp, 
					  EnumeratedVariable *xx, 
					  EnumeratedVariable *yy, 
					  EnumeratedVariable *zz, 
					  BinaryConstraint* xy,
					  BinaryConstraint* xz,
					  BinaryConstraint* yz,
					  vector<Cost> &tab, 
					  StoreStack<Cost, Cost> *storeCost);

	TernaryConstraint(WCSP *wcsp, StoreStack<Cost, Cost> *storeCost); 

    void setBinaries( BinaryConstraint* xyin, BinaryConstraint* xzin, BinaryConstraint* yzin ) { xy = xyin; xz = xzin; yz = yzin; }

    BinaryConstraint* xy;
    BinaryConstraint* xz;
    BinaryConstraint* yz;

	~TernaryConstraint() {}
    
    Cost getCost(Value vx, Value vy, Value vz) {
        int ix = x->toIndex(vx);
        int iy = y->toIndex(vy);
        int iz = z->toIndex(vz);
        Cost res = costs[ix*sizeY*sizeZ + iy*sizeZ + iz] - deltaCostsX[ix] - deltaCostsY[iy] - deltaCostsZ[iz];
        assert(res >= MIN_COST);
        return res;
    }

    Cost getCost(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz) {  
        int vindex[3];
        vindex[ getIndex(xx) ] = xx->toIndex(vx);
        vindex[ getIndex(yy) ] = yy->toIndex(vy);
        vindex[ getIndex(zz) ] = zz->toIndex(vz);
      
        Cost res = costs[vindex[0]*sizeY*sizeZ + vindex[1]*sizeZ + vindex[2]] - deltaCostsX[vindex[0]] - deltaCostsY[vindex[1]] - deltaCostsZ[vindex[2]];
        assert(res >= MIN_COST);
        return res;
    }

    Cost getCostWithBinaries(EnumeratedVariable* xx, EnumeratedVariable* yy, EnumeratedVariable* zz, Value vx, Value vy, Value vz) {
        int vindex[3];
        vindex[ getIndex(xx) ] = xx->toIndex(vx);
        vindex[ getIndex(yy) ] = yy->toIndex(vy);
        vindex[ getIndex(zz) ] = zz->toIndex(vz);
      
        Cost res = costs[vindex[0]*sizeY*sizeZ + vindex[1]*sizeZ + vindex[2]] - deltaCostsX[vindex[0]] - deltaCostsY[vindex[1]] - deltaCostsZ[vindex[2]];
        if (xy->connected()) res += xy->getCost(x,y,vindex[0],vindex[1]);
        if (xz->connected()) res += xz->getCost(x,z,vindex[0],vindex[2]);
        if (yz->connected()) res += yz->getCost(y,z,vindex[1],vindex[2]);
        assert(res >= MIN_COST);
        return res;
    }

   void addCosts( TernaryConstraint* xyz ) {
		int ix, iy, iz;
		for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
		for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
		for (EnumeratedVariable::iterator iterz = z->begin(); iterz != z->end(); ++iterz) {
        	ix = x->toIndex(*iterx); iy = y->toIndex(*itery); iz = z->toIndex(*iterz);
        	if(costs[ix*sizeY*sizeZ + iy*sizeZ + iz] < wcsp->getUb())
        		costs[ix*sizeY*sizeZ + iy*sizeZ + iz] += xyz->getCost(x,y,z, *iterx,*itery,*iterz);
	    }}}
    }


    void addCosts( EnumeratedVariable* xin, EnumeratedVariable* yin, EnumeratedVariable* zin, vector<Cost>& costsin ) {
		assert(costsin.size() == costs.size());

        int vindex[3];
        int sizeYin = yin->getDomainInitSize();
        int sizeZin = zin->getDomainInitSize();
       
		for (EnumeratedVariable::iterator iterx = xin->begin(); iterx != xin->end(); ++iterx) {
		for (EnumeratedVariable::iterator itery = yin->begin(); itery != yin->end(); ++itery) {
		for (EnumeratedVariable::iterator iterz = zin->begin(); iterz != zin->end(); ++iterz) {

		    int vxin = xin->toIndex(*iterx);
		    int vyin = yin->toIndex(*itery);
		    int vzin = zin->toIndex(*iterz);

	        vindex[ getIndex(xin) ] = vxin;
	        vindex[ getIndex(yin) ] = vyin;
	        vindex[ getIndex(zin) ] = vzin;
	        
			if(costs[vindex[0]*sizeY*sizeZ + vindex[1]*sizeZ + vindex[2]]  < wcsp->getUb())
				costs[vindex[0]*sizeY*sizeZ + vindex[1]*sizeZ + vindex[2]] += costsin[vxin*sizeYin*sizeZin + vyin*sizeZin + vzin];
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
	        
		costs[vindex[0]*sizeY*sizeZ + vindex[1]*sizeZ + vindex[2]] += c;
    }

    void setcost( EnumeratedVariable* xin, EnumeratedVariable* yin, EnumeratedVariable* zin, int vxi, int vyi, int vzi, Cost c ) {
        int vindex[3];
	    int vx = xin->toIndex(vxi);
	    int vy = yin->toIndex(vyi);
	    int vz = zin->toIndex(vzi);
        vindex[ getIndex(xin) ] = vx;
        vindex[ getIndex(yin) ] = vy;
        vindex[ getIndex(zin) ] = vz;
		costs[vindex[0]*sizeY*sizeZ + vindex[1]*sizeZ + vindex[2]] = c;
    }

    
    void propagate() {
        if (x->assigned()) {
            assign(0);
            return;
        }
        if (y->assigned()) {
            assign(1);
            return;
        }
        if (z->assigned()) {
            assign(2);
            return;
        }
        x->queueAC(); 
        y->queueAC(); 
        z->queueAC(); 
        x->queueDAC();
        y->queueDAC();
        z->queueDAC();
        x->queueEAC1();
        y->queueEAC1();
        z->queueEAC1();
    }
    
    void remove(int varIndex) {
        switch(varIndex) {
            case 0: if (ToulBar2::LcLevel==LC_AC || getDACScopeIndex()!=1) findSupportY(); if(connected()&&(ToulBar2::LcLevel==LC_AC || getDACScopeIndex()!=2)) findSupportZ(); break;
            case 1: if (ToulBar2::LcLevel==LC_AC || getDACScopeIndex()!=0) findSupportX(); if(connected()&&(ToulBar2::LcLevel==LC_AC || getDACScopeIndex()!=2)) findSupportZ(); break;
            case 2: if (ToulBar2::LcLevel==LC_AC || getDACScopeIndex()!=0) findSupportX(); if(connected()&&(ToulBar2::LcLevel==LC_AC || getDACScopeIndex()!=1)) findSupportY(); break;
        }
    }

    void projectFromZero(int varIndex) {
        switch(varIndex) {
            case 0: if (getDACScopeIndex()==1) findFullSupportY(); else if (getDACScopeIndex()==2) findFullSupportZ(); break;
            case 1: if (getDACScopeIndex()==0) findFullSupportX(); else if (getDACScopeIndex()==2) findFullSupportZ(); break;
            case 2: if (getDACScopeIndex()==0) findFullSupportX(); else if (getDACScopeIndex()==1) findFullSupportY(); break;
        }
    }

    //Trick! instead of doing remove(index) now, let AC queue do the job. 
    //So several incdec events on the same constraint can be merged into one AC event
    void increase(int index) {if (index==0) x->queueAC(); else if (index==1) y->queueAC(); else z->queueAC();}
    void decrease(int index) {if (index==0) x->queueAC(); else if (index==1) y->queueAC(); else z->queueAC();}
    
    void assign(int varIndex) {
        deconnect();                    
        switch(varIndex) {
            case 0: projectTernaryBinary(yz);  break;
            case 1: projectTernaryBinary(xz);  break;
            case 2: projectTernaryBinary(xy);  break;
        }
    }

  void fillEAC2(EnumeratedVariable *x, EnumeratedVariable *y, EnumeratedVariable *z,
                vector< pair<Value,Value> > &supportX) {
    assert(x->canbe(x->getSupport()));
    assert(getIndex(y) < getIndex(z));
    int xindex = x->toIndex(x->getSupport());
    Value ysupport = supportX[xindex].first;
    Value zsupport = supportX[xindex].second;
    if (y->cannotbe(ysupport) || z->cannotbe(zsupport) || 
        getCostWithBinaries(x,y,z,x->getSupport(),ysupport,zsupport) + y->getCost(ysupport) + z->getCost(zsupport) > MIN_COST) {
            x->queueEAC2();
    }
  }
                
  void fillEAC2(int varIndex) {
        switch(varIndex) {
            case 0: fillEAC2(y,x,z,supportY); fillEAC2(z,x,y,supportZ); break;
            case 1: fillEAC2(x,y,z,supportX); fillEAC2(z,x,y,supportZ); break;
            case 2: fillEAC2(x,y,z,supportX); fillEAC2(y,x,z,supportY); break;
        }
  }
   
  bool isEAC(int varIndex, Value a) {
//      if (varIndex==getDACScopeIndex()) return true;
    switch(varIndex) {
            case 0: return isEAC(x,a,y,z,supportX); break;
            case 1: return isEAC(y,a,x,z,supportY); break;
            case 2: return isEAC(z,a,x,y,supportZ); break;
            default: exit(EXIT_FAILURE);
    }
    return true;
  }

    void findFullSupportEAC(int varIndex) {
//          if (varIndex==getDACScopeIndex()) return;
        switch(varIndex) {
            case 0: findFullSupportX(); break;
            case 1: findFullSupportY(); break;
            case 2: findFullSupportZ(); break;
        }
    }
            
    bool verify() {
	  if (ToulBar2::LcLevel==LC_DAC) {
        switch(getDACScopeIndex()) {
            case 0: return verifyX(); break;
            case 1: return verifyY(); break;
            case 2: return verifyZ(); break;
		    default: return false;
        }
	  } else {
		return verifyX() && verifyY() && verifyZ();
	  }
	}
	  
    void projectTernaryBinary( BinaryConstraint* yzin );

	void projectTernary() {
		projectTernaryBinary(xy);
		projectTernaryBinary(xz);
		projectTernaryBinary(yz);
	}

	void extendTernary()
	{
		Cost c; 
		for(EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
		for(EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
		for(EnumeratedVariable::iterator iterz = z->begin(); iterz != z->end(); ++iterz) {
			if(xy->connected()) { c = xy->getCost(x, y, *iterx, *itery); addcost(x,y,z,*iterx,*itery,*iterz, c); }
			if(xz->connected()) { c = xz->getCost(x, z, *iterx, *iterz); addcost(x,y,z,*iterx,*itery,*iterz, c); }
			if(yz->connected()) { c = yz->getCost(y, z, *itery, *iterz); addcost(x,y,z,*iterx,*itery,*iterz, c); }
		}}}

        xy->clearCosts();
        xz->clearCosts();
        yz->clearCosts();

        xy->deconnect();
        xz->deconnect();
        yz->deconnect();         
	}

	BinaryConstraint* commonBinary( TernaryConstraint* t )
	{
		if( (t->getIndex(xy->getVar(0)) >= 0) && (t->getIndex(xy->getVar(1)) >= 0) ) return xy;
		else if((t->getIndex(xz->getVar(0)) >= 0) && (t->getIndex(xz->getVar(1)) >= 0)) return xz;
		else if((t->getIndex(yz->getVar(0)) >= 0) && (t->getIndex(yz->getVar(1)) >= 0)) return yz;
		return NULL;
	}

    double computeTightness();


    EnumeratedVariable* xvar;
    EnumeratedVariable* yvar;
    EnumeratedVariable* zvar;
    EnumeratedVariable::iterator itvx;
    EnumeratedVariable::iterator itvy;
    EnumeratedVariable::iterator itvz;

    void first() {
    	itvx = x->begin(); 
    	itvy = y->begin(); 
    	itvz = z->begin(); 
    	xvar = x;
    	yvar = y;
    	zvar = z;
    }
    
    bool next( string& t, Cost& c) 
    { 
    	char tch[4];
    	if(itvx != xvar->end()) {
    		int ix = xvar->toIndex(*itvx);
	    	tch[0] = ix + CHAR_FIRST;
	    	if(itvy != yvar->end()) {
	    		int iy = yvar->toIndex(*itvy);
		    	tch[1] = iy + CHAR_FIRST;
		    	if(itvz != zvar->end()) {
		    		int iz = zvar->toIndex(*itvz);
			    	tch[2] = iz + CHAR_FIRST;
		 	    	tch[3] = '\0';
			    	t = tch;
			    	c = getCost(xvar,yvar,zvar,*itvx, *itvy, *itvz); 
					++itvz;
			    	return true;
	    		} else {
		    		++itvy;
		    		itvz = zvar->begin();
		    		return next(t,c);
		    	}
	    	} else {
		    		++itvx;
		    		itvy = yvar->begin();
		    		return next(t,c);
		    }
    	}
    	return false; 
    }

    void firstlex() { first(); }
    bool nextlex( string& t, Cost& c) { return next(t,c); } 


	void setTuple( string& t, Cost c, EnumeratedVariable** scope_in ) {
		Value v0 = scope_in[0]->toValue(t[0]-CHAR_FIRST);
		Value v1 = scope_in[1]->toValue(t[1]-CHAR_FIRST);
		Value v2 = scope_in[2]->toValue(t[2]-CHAR_FIRST);
		setcost( scope_in[0], scope_in[1], scope_in[2], v0, v1, v2, c );		
	}

	void addtoTuple( string& t, Cost c, EnumeratedVariable** scope_in ) {
		Value v0 = scope_in[0]->toValue(t[0]-CHAR_FIRST);
		Value v1 = scope_in[1]->toValue(t[1]-CHAR_FIRST);
		Value v2 = scope_in[2]->toValue(t[2]-CHAR_FIRST);
		addcost( scope_in[0], scope_in[1], scope_in[2], v0, v1, v2, c );		
	}

    Cost getDefCost() { return wcsp->getUb(); }
	void setDefCost( Cost df ) {}


	Cost evalsubstr( string& s, Constraint* ctr )
	{
		Value vals[3];
		int count = 0;
		
		for(int i=0;i<arity();i++) {
			EnumeratedVariable* var = (EnumeratedVariable*) getVar(i); 
			int ind = ctr->getIndex(var);
			if(ind >= 0) { vals[i] = var->toValue(s[ind] - CHAR_FIRST); count++; }	
		}
		if(count == 3) return getCost(vals[0], vals[1], vals[2]);
		else return MIN_COST;
	}    
	

    void fillElimConstr( EnumeratedVariable* xin, EnumeratedVariable* yin, EnumeratedVariable* zin)
	{
		x = xin;
		y = yin;
		z = zin;
		sizeX = x->getDomainInitSize();
		sizeY = y->getDomainInitSize();
		sizeZ = z->getDomainInitSize();
		linkX->removed = true;
		linkY->removed = true;
		linkZ->removed = true;
		linkX->content.constr = this;
		linkY->content.constr = this; 
		linkZ->content.constr = this; 
		linkX->content.scopeIndex = 0;
		linkY->content.scopeIndex = 1;
		linkZ->content.scopeIndex = 2;
        if (x->wcspIndex < y->wcspIndex && x->wcspIndex < z->wcspIndex) dacvar = 0;
        else if (y->wcspIndex < x->wcspIndex && y->wcspIndex < z->wcspIndex) dacvar = 1;
        else dacvar = 2;
	}
   
    void fillElimConstrBinaries()
	{
        BinaryConstraint* xy_ = x->getConstr(y); 
		if(xy_)  xy = xy_;  else { xy = wcsp->newBinaryConstr(x,y); wcsp->elimBinOrderInc(); }

        BinaryConstraint* xz_ = x->getConstr(z); 
		if(xz_)  xz = xz_;  else { xz = wcsp->newBinaryConstr(x,z); wcsp->elimBinOrderInc(); }

        BinaryConstraint* yz_ = y->getConstr(z); 
		if(yz_)  yz = yz_;  else { yz = wcsp->newBinaryConstr(y,z); wcsp->elimBinOrderInc(); }
	}


    void print(ostream& os);
    void dump(ostream& os);
};

#endif /*TB2TERNARYCONSTR_HPP_*/
