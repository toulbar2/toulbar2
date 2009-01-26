#ifndef TB2NARYCONSTR_HPP_
#define TB2NARYCONSTR_HPP_

#include "tb2abstractconstr.hpp"
#include "tb2ternaryconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2wcsp.hpp"

#include <set>


class NaryConstraint : public AbstractNaryConstraint
{
	public:
	
	BinaryConstraint* xy;      	// xy is an empty constraint that is created at start time for 
								// projecting the nary constraint when all variables but 2 are assigned  	

	TernaryConstraint* xyz;
	
	     
	Cost default_cost;          // default cost returned when tuple t is not found in TUPLES (used by function eval(t)
	bool store_top; 		    // this is true when default_cost < getUb() meaning that tuples with cost greater than ub must be stored
	StoreInt nonassigned;       // nonassigned variables during search, must be backtrackable (storeint) !

	String iterTuple;
	String evalTuple;

	vector<EnumeratedVariable::iterator> it_values;
	void firstlex();
    bool nextlex( String& t, Cost& c);
	
	NaryConstraint(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in, Cost defval);
	NaryConstraint(WCSP *wcsp);
	
	virtual void setTuple( String& tin, Cost c, EnumeratedVariable** scope_in = NULL ) = 0;

	virtual void addtoTuple( String& tin, Cost c, EnumeratedVariable** scope_in = NULL ) = 0;
	
	virtual void setTuple( int* tin, Cost c, EnumeratedVariable** scope_in ) 
	{
		Char* buf = new Char [arity_ + 1];
		for(int i=0;i<arity_;i++) buf[i] = tin[i]+CHAR_FIRST; 
		buf[arity_] = '\0';
		String str = String(buf); 
		setTuple( str, c, scope_in );
		delete [] buf;
	}

	virtual void addtoTuple( int* tin, Cost c, EnumeratedVariable** scope_in ) 
	{
		Char* buf = new Char [arity_ + 1];
		for(int i=0;i<arity_;i++) buf[i] = tin[i]+CHAR_FIRST; 
		buf[arity_] = '\0';
		String str = String(buf); 
		addtoTuple( str, c, scope_in );
		delete [] buf;
	}

	
	
    virtual Cost eval( String& s ) = 0; 
	Cost evalsubstr( String& s, Constraint* ctr );
	
	void assign(int varIndex);
	
	void projectNary();
	void projectNaryTernary(TernaryConstraint* xyz);	
	void projectNaryBinary(BinaryConstraint* xy);	
	
    void propagate() {
        for(int i=0;connected() && i<arity_;i++) {         
            if (getVar(i)->assigned()) assign(i);
        }
    };

    virtual void project( EnumeratedVariable* x, bool addUnaryCtr = true ) = 0;

    double computeTightness() { return 0; }
    bool   verify() {return true;}
    void   increase(int index) {}
    void   decrease(int index) {}
    void  remove(int index) {}
 
    void starrule(String& t, Cost minc);	    
    void projectFromZero(int index);

    virtual void print(ostream& os) {}
};


class NaryConstraintMap : public NaryConstraint
{
	typedef map<String,Cost> TUPLES;
    TUPLES* pf;


public:

	NaryConstraintMap(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in, Cost defval);
	NaryConstraintMap(WCSP *wcsp);
	virtual ~NaryConstraintMap();


	bool consistent( String& t );
    Cost eval( String& s );
	
	Cost getDefCost() { return default_cost; }
	void setDefCost( Cost df ) { default_cost = df; }
	void changeDefCost( Cost df );

  
    
    set<Constraint*>* filters;
	void resetFilters();
	void fillFilters();
	 
	void project( EnumeratedVariable* x, bool addUnaryCtr = true );
	void sum( NaryConstraintMap* nary );

	TUPLES::iterator  tuple_it;

	void first();
    bool next( String& t, Cost& c);
    
	void setTuple( String& tin, Cost c, EnumeratedVariable** scope_in = NULL );
	void addtoTuple( String& tin, Cost c, EnumeratedVariable** scope_in = NULL );
    void insertSum( String& t1, Cost c1, Constraint* ctr1, String t2, Cost c2, Constraint* ctr2, bool bFilters = false );  
	void permute( EnumeratedVariable** scope_in );
	
	void projectxy( EnumeratedVariable* x, EnumeratedVariable* y, TUPLES& fproj);
	void projectxyz( EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z, TUPLES& fproj);
	void preproject3();
	void preprojectall2();

	void fillRandom();
    void print(ostream& os);
    void dump(ostream& os);
    
};



class Trie;

class TrieNode {

  public:
    TrieNode();

	void iniLeaf(Char *w);
	void iniNonLeaf(Char ch);

    Cost c;
	
  private:
    bool leaf, endOfWord;
    Char *letters;
    Char *word;

    TrieNode **ptrs;
    friend class Trie;
};


class Trie {
  public:
    Trie() : notFound(-1) {}
    Trie(Char*, Cost c);
    void insert(Char*, Cost c);
    TrieNode* find(const Char*);
     
    void printTrie(); 
    
  private:
    TrieNode *root;
    const int notFound;
    Char prefix[80];
    int  position(TrieNode*,Char);
    void addCell(Char,TrieNode*,int);
    TrieNode* createLeaf(Char,Char*,TrieNode*);
    void printTrie(int,TrieNode*,Char*);
};





class NaryConstrie : public NaryConstraint
{
	
  public:
	Trie* f;

	
	NaryConstrie(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in, Cost defval);
	NaryConstrie(WCSP *wcsp);
	virtual ~NaryConstrie();
	
	void setTuple( String& tin, Cost c, EnumeratedVariable** scope_in = NULL );
	void addtoTuple( String& tin, Cost c, EnumeratedVariable** scope_in = NULL );

    void project( EnumeratedVariable* x, bool addUnaryCtr = true ) {};

    Cost eval( String& s );

    void print(ostream& os);

};



#endif /*TB2NARYCONSTR_HPP_*/
