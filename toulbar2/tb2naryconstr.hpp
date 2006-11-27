#ifndef TB2NARYCONSTR_HPP_
#define TB2NARYCONSTR_HPP_

#include "tb2abstractconstr.hpp"
#include "tb2ternaryconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2wcsp.hpp"

#include <map>
using namespace std;

#define CHAR_FIRST 'a'

class NaryConstraint : public AbstractNaryConstraint
{
	typedef map<string,Cost> TUPLES;
    TUPLES f;

	Cost default_cost;          // default cost returned when tuple t is not found in TUPLES (used by function eval(t)
	bool store_top; 		    // this is true when default_cost < getUb() meaning that tuples with cost greater than ub must be stored
	
	StoreInt nonassigned;       // nonassigned variables during search, must be backtrackable (storeint) !

public:
	NaryConstraint(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in, Cost defval);


	BinaryConstraint* xy;      	// xy is an empty constraint that is created at start time for 
								// projecting the nary constraint when all variables but 2 are assigned  	
	void projectNaryBinary();

	void insertTuple( string t, Cost c, EnumeratedVariable** scope_in );  
    Cost eval( string s );


	void assign(int varIndex);


    double computeTightness() { return 0; }


    void propagate() {
        for(int i=0;connected() && i<arity_;i++) {         
            if (getVar(i)->assigned()) assign(i);
        }
    };
    void increase(int index) {}
    void decrease(int index) {}
    void remove(int index) {}
    void projectFromZero(int index) {}
    void fillEAC2(int index) {}
    bool isEAC(int index, Value a) {return true;}
    void findFullSupport(int index) {}

    bool verify() {return true;}
    

	void project( EnumeratedVariable* x );
	void sum( AbstractNaryConstraint* nary );


	void fillRandom();
    void print(ostream& os);
    void dump(ostream& os);
};



#endif /*TB2NARYCONSTR_HPP_*/
