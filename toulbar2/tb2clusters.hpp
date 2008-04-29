 
#ifndef TB2CLUSTERS_HPP_
#define TB2CLUSTERS_HPP_

#include "tb2wcsp.hpp"
#include "tb2solver.hpp"
#include "tb2enumvar.hpp"
#include "tb2naryconstr.hpp"

#include <set>
#include <list>


class Cluster;
class Nogood;

typedef set<int>	       TVars;
typedef set<Constraint*>   TCtrs;
typedef map<int,Value>     TAssign;
typedef set<Cluster*>	   TClusters;

typedef pair <Cost,bool>     TPairNG;
typedef pair <Cost,string>   TPairSol;

typedef map<string, TPairNG>  TNoGoods;
typedef map<string, TPairSol> TSols;


class Separator : public AbstractNaryConstraint
{
  private:

	Cluster*					  cluster;
	TVars   					  vars;
    vector< vector<StoreCost> >   delta;    // structure to record the costs that leave the cluster
	  										// inicialized with iniDelta()

	StoreInt nonassigned;       			// nonassigned variables during search, must be backtrackable (storeint) !
    StoreInt isUsed;
    StoreCost lbPrevious;
    StoreInt optPrevious;

    TNoGoods  					  nogoods;
    TSols  						  solutions;
    DLink<Separator *>            linkSep;

	string t;    // buffer for a separator tuple                   
	string s;    // buffer for a solution tuple                   
	
	
  public:

	Separator(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in);
	Separator(WCSP *wcsp);

    void queueSep() { wcsp->queueSeparator(&linkSep); }
    void unqueueSep() { wcsp->unqueueSeparator(&linkSep); }

	void assign(int varIndex);
    void propagate();

    bool used() { return isUsed; }

    void setup(Cluster* cluster_in);

	TVars& 			 getVars() { return vars; }    
    
    void addDelta( unsigned int posvar, Value value, Cost cost ) {
    	assert( posvar < vars.size() ); 
    	delta[posvar][value] += cost; 
    }

    void set( Cost c, bool opt );
    bool get( Cost& res, bool& opt );

	void solRec( Cost ub );
	bool solGet( TAssign& a, string& sol ); 

	void resetOpt();
        
    TVars::iterator  begin() { return vars.begin(); } 
    TVars::iterator  end()   { return vars.end(); } 
    bool  is(int i) { return vars.find(i) != vars.end(); } 


	// Obliged to include these methods; if not abstract class problem
    double computeTightness() { return 0; }
    bool   verify() {return true;}
    void   increase(int index) {}
    void   decrease(int index) {}
    void   remove(int index) {}
    
    void   print(ostream& os);
};



#define NOCLUSTER -1

class Cluster {

 private:
  	  TreeDecomposition*  td;
	  WCSP*				  wcsp;
	  list<TAssign*>      assignments;
	  TVars				  vars;
	  TCtrs			      ctrs;
	  TClusters           edges;              // adjacent clusters 

	  StoreCost           lb;	
	  Cost				  lb_opt;
	  Cost				  ub;
	  
	  StoreInt			  active;
	  Cluster*			  parent;             // parent cluster

 public:
	  Separator* 		  sep;

	  Cluster (TreeDecomposition* tdin);
	  ~Cluster();

	  int id;								  // the id corresponds to the vector index of the cluster in ClusteredWCSP
	
	  int getId() { return id; }

      // ----------------------------------------- Interface Functions
	  WCSP* 		getWCSP() { return wcsp; }

	  bool 			isVar( int i );
	  bool 			isSepVar( int i );
	  bool			isActive() { int a = active; return a == 1; }

	  Cost	        getLbRec();
	  Cost	        getLbRecNoGoods();
	  Cost	        getLbRecNoGood(bool& opt);

	  Cost		    getOpt() { return lb_opt; }
	  Cost			getUb()  { return ub; }
	  Cost			getLb()  { return lb; }
	  Cost			getLb_opt()  { return lb_opt; }
	  void			setUb(Cost c)  {ub = c; }
	  void			setLb(Cost c)  {lb = c; }
	  void			setLb_opt(Cost c)  {lb_opt = c; }
	  int			getNbVars() { return vars.size(); }
	  TVars&		getVars() { return vars; }	
	  TCtrs&		getCtrs() { return ctrs; }	
	  TClusters&	getEdges() { return edges; }
	  void 			setSep( Separator* sepin ) { sep = sepin; } 
	  void 			addVar( Variable* x );
	  void 			removeVar( Variable* x );
	  void 			addVars( TVars& vars );
	  void 			addEdge( Cluster* c );
	  void 			addEdges( TClusters& cls );
	  void 			removeEdge( Cluster* c );
	  void 			addCtrs( TCtrs& ctrsin );
	  void 			addCtr( Constraint* c );
	  void 			addAssign( TAssign* a );
	  void 		    updateUb();
	  

	  TClusters::iterator removeEdge( TClusters::iterator it );
	  


  	  void 		    getSolution( TAssign& sol );
	  
	  void 			setParent(Cluster* p);
	  Cluster*		getParent();
	  TClusters&	getDescendants();
	  Cluster*		nextSep( Variable* v ); 
	  bool			isDescendant( Cluster* c2 );
	  bool  	    isEdge( Cluster* c );

	  int			sepSize() { if(sep) return sep->arity(); else return 0; }

	  
      TClusters     descendants;  	   // set of descendants	
	  vector<bool>  quickdescendants;	
	  
   	  void deconnectSep();

	  	
	
      // ----------------------------------------------------------
	  Cost 			eval(TAssign* a);
	  void 			setWCSP();                              // sets the WCSP to the cluster problem, deconnecting the rest
	  void 			reactivate();
	  void 			deactivate();
	  void 			increaseLb( Cost newlb );

	  void setup() { if(sep) sep->setup(this); }

	  
	  void addDelta( int posvar, Value value, Cost cost ) { if(sep) sep->addDelta(posvar,value,cost); }
	  void nogoodRec( Cost c, bool opt ) { if(sep) sep->set(c,opt); }	
      Cost nogoodGet( bool& opt ) { Cost c = MIN_COST; sep->get(c,opt); return c; }	
      void resetOpt() { if(sep) sep->resetOpt(); }
	  
	  void solutionRec(Cost ub) { if(sep) sep->solRec(ub); }

	  TVars::iterator beginVars() { return vars.begin(); }
	  TVars::iterator endVars()   { return vars.end(); }
	  TVars::iterator beginSep() { return sep->begin(); }
	  TVars::iterator endSep()   { return sep->end(); }
	  TCtrs::iterator beginCtrs() { return ctrs.begin(); }
	  TCtrs::iterator endCtrs()   { return ctrs.end(); }
	  TClusters::iterator beginEdges() { return edges.begin(); }
	  TClusters::iterator endEdges()   { return edges.end(); }
	  TClusters::iterator beginDescendants() { return descendants.begin(); }
	  TClusters::iterator endDescendants()   { return descendants.end(); }

	  void print();	  
	  void printStats() { if(!sep) return; sep->print(cout); }
	  	
};



class TreeDecomposition  {
private:
	WCSP*			  wcsp;	
	vector<Cluster*>  clusters; 
	StoreInt   		  currentCluster;
	list<Cluster*> 	  roots;

public:

	TreeDecomposition(WCSP* wcsp_in);

    WCSP* 		getWCSP() { return wcsp; }
	
	Cluster*	getCluster( unsigned int i ) { assert( 0 <= i && i < clusters.size() ); return clusters[i]; }
	Cluster*   	var2Cluster( int v );	
	
	void setCurrentCluster(Cluster *c) {currentCluster = c->getId();}
    Cluster* getCurrentCluster() {return getCluster(currentCluster);}

    bool isInCurrentClusterSubTree(int idc); 
    bool isActiveAndInCurrentClusterSubTree(int idc); 
	
	void buildFromOrder();	     			    // builds the tree cluster of clusters from a given order
	void fusions();                  			// fusions all redundant clusters after build_from_order is called
	bool fusion();          		            // one fusion step

    void newSolution();

	bool reduceHeight( Cluster* c );
	void makeDescendants( Cluster* c );
	
	int      makeRooted();
	void     makeRootedRec( Cluster* c,  TClusters& visited );
	Cluster* getRoot() {return roots.front();}
	 
	int height( Cluster* r, Cluster* father );
	int height( Cluster* r );

	bool verify();

	void intersection( TVars& v1, TVars& v2, TVars& vout );
	void difference( TVars& v1, TVars& v2, TVars& vout );
	void sum( TVars& v1, TVars& v2, TVars& vout ); 
	bool included( TVars& v1, TVars& v2 ); 	
	void clusterSum( TClusters& v1, TClusters& v2, TClusters& vout );	
	
    void addDelta(int c, EnumeratedVariable *x, Value value, Cost cost);

	Cluster* getBiggerCluster( TClusters& visited );	
	
	
	TClusters c1edges;
	TClusters c2edges;


	void print( Cluster* c = NULL, int recnum = 0);
};






#endif 
