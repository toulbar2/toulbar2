/** \file tb2constraint.hpp
 *  \brief Abstract constraint.
 *
 */

#ifndef TB2CONSTRAINT_HPP_
#define TB2CONSTRAINT_HPP_

#include "tb2types.hpp"


using namespace std;


class Constraint : public WCSPLink
{
    Long conflictWeight;
    
    // make it private because we don't want copy nor assignment
    Constraint(const Constraint &c);
    Constraint& operator=(const Constraint &c);

public:
    Constraint(WCSP *wcsp);
    Constraint(WCSP *wcsp, int elimCtrIndex);
        
    virtual ~Constraint() {}

    // remove a constraint from the set of active constraints
    virtual bool connected() {cout << "dummy connected on (" << this << ")!" << endl;return true;}
    virtual bool deconnected() {cout << "dummy deconnected on (" << this << ")!" << endl;return false;}
    virtual void deconnect() {cout << "dummy deconnect on (" << this << ")!" << endl;}
    virtual void reconnect() {cout << "dummy reconnect on (" << this << ")!" << endl;}

    virtual int arity() const = 0;
    virtual Variable *getVar(int scopeIndex) const = 0;
    virtual int getIndex(Variable* var) const = 0;

    Long getConflictWeight() const {return conflictWeight;}
    void incConflictWeight() {conflictWeight++;}
    
	double tight;
    double getTightness() { if(tight < 0) computeTightness(); return tight; }
    virtual double  computeTightness() = 0;
    
    
    // return the smallest wcsp index in the constraint scope except for one variable having a forbidden scope index
    virtual int getSmallestVarIndexInScope(int forbiddenScopeIndex) = 0;
    virtual int getDACScopeIndex() = 0;

    virtual void propagate() = 0;
    virtual void increase(int index) {propagate();}
    virtual void decrease(int index) {propagate();}
    virtual void remove(int index) {propagate();}
    virtual void projectFromZero(int index) {}
    virtual void assign(int index) {propagate();}

    virtual void fillEAC2(int index) {}
    virtual bool isEAC(int index, Value a) {return true;}
    virtual void findFullSupportEAC(int index) {}

    virtual bool verify() {return true;};
    
    virtual void print(ostream& os) {os << this << " Unknown constraint!";}

    virtual void dump(ostream& os) {os << this << " Unknown constraint!";}


    virtual void first() {}
    virtual void first(EnumeratedVariable** scope_in) {}
    virtual bool next( string& t, Cost& c) { return false; }
	virtual void setTuple( string t, Cost c, EnumeratedVariable** scope_in ) {}
	virtual void scopeSet( TSCOPE& scope_inv ) {}	

	void scopeCommon( TSCOPE& scope_in, Constraint* ctr ) 
	{
		TSCOPE scope1,scope2;
		scopeSet( scope1 );
		ctr->scopeSet( scope2 );
		
		set_intersection( scope1.begin(), scope1.end(),
				  	   	  scope2.begin(), scope2.end(),
					  	  inserter(scope_in, scope_in.begin()) );		
	}
	
		
	void scopeUnion( TSCOPE& scope_in, Constraint* ctr ) 
	{
		TSCOPE scope1,scope2;
		scopeSet( scope1 );
		ctr->scopeSet( scope2 );
		
		set_union( scope1.begin(), scope1.end(),
		  	   	   scope2.begin(), scope2.end(),
			  	   inserter(scope_in, scope_in.begin()) );		
	}
		
	void scopeDifference( TSCOPE& scope_in, Constraint* ctr )
	{
		TSCOPE scope1,scope2;
		scopeSet( scope1 );
		ctr->scopeSet( scope2 );
		
		set_difference( scope1.begin(), scope1.end(),
			  	   	    scope2.begin(), scope2.end(),
				  	    inserter(scope_in, scope_in.begin()) );		
	}


    friend ostream& operator<<(ostream& os, Constraint &c) {
        c.print(os);
        return os;
    }
};

#endif /*TB2CONSTRAINT_HPP_*/
