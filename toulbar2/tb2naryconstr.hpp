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

	string iterTuple;
	string evalTuple;

	vector<EnumeratedVariable::iterator> it_values;
	void firstlex();
    bool nextlex( string& t, Cost& c);
	
	NaryConstraint(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in, Cost defval);
	NaryConstraint(WCSP *wcsp);
	
	virtual void setTuple( string& tin, Cost c, EnumeratedVariable** scope_in = NULL ) = 0;

	virtual void addtoTuple( string& tin, Cost c, EnumeratedVariable** scope_in = NULL ) = 0;
	
	virtual void setTuple( int* tin, Cost c, EnumeratedVariable** scope_in ) 
	{
		char* buf = new char [arity_];
		for(int i=0;i<arity_;i++) buf[i] = tin[i]+CHAR_FIRST; 
		buf[arity_] = '\0';
		string str = string(buf); 
		setTuple( str, c, scope_in );
		delete [] buf;
	}

	virtual void addtoTuple( int* tin, Cost c, EnumeratedVariable** scope_in ) 
	{
		char* buf = new char [arity_];
		for(int i=0;i<arity_;i++) buf[i] = tin[i]+CHAR_FIRST; 
		buf[arity_] = '\0';
		string str = string(buf); 
		addtoTuple( str, c, scope_in );
		delete [] buf;
	}

	
	
    virtual Cost eval( string& s ) = 0; 
	Cost evalsubstr( string& s, Constraint* ctr );
	
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
 
    void starrule(string& t, Cost minc);	    
    void projectFromZero(int index);

    virtual void print(ostream& os) {}
};


class NaryConstraintMap : public NaryConstraint
{
	typedef map<string,Cost> TUPLES;
    TUPLES* pf;


public:

	NaryConstraintMap(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in, Cost defval);
	NaryConstraintMap(WCSP *wcsp);
	virtual ~NaryConstraintMap();


	bool consistent( string& t );
    Cost eval( string& s );
	
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
    bool next( string& t, Cost& c);
    
	void setTuple( string& tin, Cost c, EnumeratedVariable** scope_in = NULL );
	void addtoTuple( string& tin, Cost c, EnumeratedVariable** scope_in = NULL );
    void insertSum( string& t1, Cost c1, Constraint* ctr1, string t2, Cost c2, Constraint* ctr2, bool bFilters = false );  
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

	void iniLeaf(char *w);
	void iniNonLeaf(char ch);

    Cost c;
	
  private:
    bool leaf, endOfWord;
    char *letters;
    char *word;

    TrieNode **ptrs;
    friend class Trie;
};


class Trie {
  public:
    Trie() : notFound(-1) {}
    Trie(char*, Cost c);
    void insert(char*, Cost c);
    TrieNode* find(const char*);
     
    void printTrie(); 
    
  private:
    TrieNode *root;
    const int notFound;
    char prefix[80];
    int  position(TrieNode*,char);
    void addCell(char,TrieNode*,int);
    TrieNode* createLeaf(char,char*,TrieNode*);
    void printTrie(int,TrieNode*,char*);
};





class NaryConstrie : public NaryConstraint
{
	
  public:
	Trie* f;

	
	NaryConstrie(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in, Cost defval);
	NaryConstrie(WCSP *wcsp);
	virtual ~NaryConstrie();
	
	void setTuple( string& tin, Cost c, EnumeratedVariable** scope_in = NULL );
	void addtoTuple( string& tin, Cost c, EnumeratedVariable** scope_in = NULL );

    void project( EnumeratedVariable* x, bool addUnaryCtr = true ) {};

    Cost eval( string& s );

    void print(ostream& os);

};



#endif /*TB2NARYCONSTR_HPP_*/
