/** \file tb2abstractconstr.hpp
 *  \brief Abstract constraints of predefined arities
 *
 */

#ifndef TB2ABSTRACTCONSTR_HPP_
#define TB2ABSTRACTCONSTR_HPP_

#include "tb2constraint.hpp"
#include "tb2variable.hpp"
#include "tb2enumvar.hpp"
#include "tb2wcsp.hpp"


template<class T1, class T2>
class AbstractBinaryConstraint : public Constraint
{
protected:
    T1 *x;
    T2 *y;
    DLink<ConstraintLink> *linkX;
    DLink<ConstraintLink> *linkY;
    int dacvar;
    
public:
    AbstractBinaryConstraint(WCSP *wcspin, T1 *xx, T2 *yy) : Constraint(wcspin), x(xx), y(yy), linkX(NULL), linkY(NULL) {
        assert(xx != yy);
        linkX = xx->link(this,0);
        linkY = yy->link(this,1);
        if (xx->wcspIndex < yy->wcspIndex) dacvar = 0; else dacvar = 1;
    }

    AbstractBinaryConstraint(WCSP *wcspin) : Constraint(wcspin,0), x(NULL), y(NULL), linkX(NULL), linkY(NULL) 
    { }

    virtual ~AbstractBinaryConstraint() {delete linkX; delete linkY;}

    bool connected() {return !linkX->removed && !linkY->removed;}
    bool deconnected() {return linkX->removed || linkY->removed;}
    void deconnect() {
        if (connected()) {
            if (ToulBar2::verbose >= 3) cout << "deconnect " << this << endl; 
            x->deconnect(linkX);
            y->deconnect(linkY);
        }
    }
    void reconnect() {
        if (deconnected()) {
            x->getConstrs()->push_back(linkX, true); 
            y->getConstrs()->push_back(linkY, true);
        }
    }

    int arity() const {return 2;}
    
    Variable *getVar(int varCtrIndex) const {return (varCtrIndex == 0)?x:y;}

    Variable *getVarDiffFrom( Variable* v ) const  {
		if(v == x) return y;
		else if(v == y) return x;
		else abort();
	}
		    
    int getIndex(Variable* var) const 
    {
    	if(var == x) return 0;
    	else if(var == y) return 1;
    	return -1;
    }

    int getSmallestVarIndexInScope(int forbiddenScopeIndex) {assert(forbiddenScopeIndex >= 0); assert(forbiddenScopeIndex < 2); return (forbiddenScopeIndex)?x->wcspIndex:y->wcspIndex;}
    int getDACScopeIndex() {return dacvar;}
};


template<class T1, class T2, class T3>
class AbstractTernaryConstraint : public Constraint
{
protected:
    T1 *x;
    T2 *y;
    T3 *z;
    DLink<ConstraintLink> *linkX;
    DLink<ConstraintLink> *linkY;
    DLink<ConstraintLink> *linkZ;
    int dacvar;
    
public:
    AbstractTernaryConstraint(WCSP *wcsp, T1 *xx, T2 *yy, T2 *zz) : Constraint(wcsp), x(xx), y(yy), z(zz), linkX(NULL), linkY(NULL), linkZ(NULL) {
        assert(xx != yy);
        assert(xx != zz);
        assert(yy != zz);
        linkX = xx->link(this,0);
        linkY = yy->link(this,1);
        linkZ = zz->link(this,2);
        if (xx->wcspIndex < yy->wcspIndex && xx->wcspIndex < zz->wcspIndex) dacvar = 0;
        else if (yy->wcspIndex < xx->wcspIndex && yy->wcspIndex < zz->wcspIndex) dacvar = 1;
        else dacvar = 2;
    }

    virtual ~AbstractTernaryConstraint() {delete linkX; delete linkY; delete linkZ;}

    bool connected() {return !linkX->removed && !linkY->removed && !linkZ->removed;}
    bool deconnected() {return linkX->removed || linkY->removed || linkZ->removed;}
    void deconnect() {
        if (connected()) {
            if (ToulBar2::verbose >= 3) cout << "deconnect " << this << endl; 
            x->deconnect(linkX);
            y->deconnect(linkY);
            z->deconnect(linkZ);
        }
    }
    void reconnect() {
        if (deconnected()) {
            x->getConstrs()->push_back(linkX, true); 
            y->getConstrs()->push_back(linkY, true);
            z->getConstrs()->push_back(linkZ, true);
        }
    }

    int arity() const {return 3;}
    
    Variable *getVar(int varCtrIndex) const 
	{
		switch(varCtrIndex) { case 0: return x; break;
						      case 1: return y; break;
						      case 2: return z; break;
						      default: abort(); }
	}

    Variable *getVarDiffFrom( Variable* v1, Variable* v2 ) const 
	{
		if      ((x == v1) && (y == v2)) return z;
		else if ((x == v2) && (y == v1)) return z;
		else if ((x == v1) && (z == v2)) return y;
		else if ((x == v2) && (z == v1)) return y;
		else if ((y == v1) && (z == v2)) return x;
		else if ((y == v2) && (z == v1)) return x;
		else abort();
	}

    int getIndex(Variable* var) const 
    {
    	if(var == x) return 0;
    	else if(var == y) return 1;
    	else if(var == z) return 2;
    	return -1;
    }


    int getSmallestVarIndexInScope(int forbiddenScopeIndex) 
	{
        assert(forbiddenScopeIndex >= 0); assert(forbiddenScopeIndex < 3);
		switch (forbiddenScopeIndex) {
            case 0: return min(y->wcspIndex,z->wcspIndex); break;
            case 1: return min(x->wcspIndex,z->wcspIndex); break;
            case 2: return min(x->wcspIndex,y->wcspIndex); break;
            default: abort();
        }
	}
    int getDACScopeIndex() {return dacvar;}
};








#include <set>

class AbstractNaryConstraint : public Constraint
{
protected:
	
	int arity_;
    int arity() const {return arity_;}
	
	typedef struct { 
		EnumeratedVariable* var;
		int pos;
	} varElem;

	struct indexCmp {
		bool operator()(const varElem* e1, const varElem* e2) const { return e1->var->wcspIndex > e2->var->wcspIndex; }
	};
	
	EnumeratedVariable** scope;
    typedef set<varElem*, indexCmp> SCOPE;
    SCOPE scope_inv;
	
	int nconnected_links;
    DLink<ConstraintLink>** links;
    
public:
    AbstractNaryConstraint(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in) : Constraint(wcsp), arity_(arity_in) 
    {
    	nconnected_links = 0;
    	scope = new EnumeratedVariable* [arity_];
    	links = new DLink<ConstraintLink>* [arity_];
    	
		varElem* ve;
		for(int i=0; i < arity_; i++) {
			ve = new varElem;
			ve->var = scope_in[i];
			ve->pos = i;
			scope_inv.insert( ve );
			scope[i] = ve->var;
			links[i] = ve->var->link(this,i);
			if(!links[i]->removed) nconnected_links++;
		}
    }

    virtual ~AbstractNaryConstraint() {}

    Variable *getVar(int varCtrIndex) const { 
    	assert(varCtrIndex < arity_); 
    	return scope[varCtrIndex];
    }

	int getIndex(Variable* var) const { 
		varElem ve;
		ve.var = (EnumeratedVariable*) var;
		SCOPE::iterator it = scope_inv.find(&ve);
		if(it == scope_inv.end()) return -1;
		else return (*it)->pos;
    }


	bool connected() {
       for(int i=0;i<arity_;i++) if(!links[i]->removed) return true;
       return false;
	}
	
    bool deconnected() {
       for(int i=0;i<arity_;i++) if(!links[i]->removed) return false;
       return true;
    }

    void deconnect() {
        if (connected()) {
            if (ToulBar2::verbose >= 3) cout << "deconnect " << this << endl; 
            for(int i=0;i<arity_;i++) {
            	scope[i]->deconnect( links[i] );
            }
        }
    }

    void reconnect() {
        if (deconnected()) {
            for(int i=0;i<arity_;i++) {
	            scope[i]->getConstrs()->push_back(links[i], true); 
	        }
        }
    }

	virtual Cost eval( string& t ) { return -1; }
	virtual void insertTuple( string t, Cost c, EnumeratedVariable** scope_in ) { }

    int getSmallestVarIndexInScope(int forbiddenScopeIndex) 
    {
		return -1;
	}

    int getDACScopeIndex() {return -1;}



	virtual void sum( AbstractNaryConstraint* nary ) {}
	virtual void project( EnumeratedVariable* x ) {}
    
};



#endif /*TB2ABSTRACTCONSTR_HPP_*/
