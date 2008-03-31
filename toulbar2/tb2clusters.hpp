 
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

typedef pair <Cost,bool>   TPair;
typedef map<string, TPair> TNoGoods;


class Separator 
{
  private:

	Cluster*					  cluster;
	TVars   					  vars;
    vector< vector<StoreCost> >   delta;    // structure to record the costs that leave the cluster
	  										// inicialized with iniDelta()

    TNoGoods  					  nogoods;

  public:

	Separator();

    void setup(Cluster* cluster_in);

	TVars& 			 getVars() { return vars; }    
    
    void addDelta( int posvar, Value value, Cost cost ) { delta[posvar][value] += cost; }


    void set( Cost c, bool opt );
    Cost get( bool& opt );
    
    TVars::iterator  begin() { return vars.begin(); } 
    TVars::iterator  end()   { return vars.end(); } 
    bool  is(int i) { return vars.find(i) != vars.end(); } 

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
	  Separator 		  sep;

	  StoreCost           lb;	
	  Cost				  lb_opt;
	  Cost				  ub;
	  
	  
	  Cluster*			  parent;             // parent cluster
      TClusters    	      descendants;  	 // set of descendants	

 public:
	  Cluster (TreeDecomposition* tdin);
	  ~Cluster();

	  int id;								  // the id corresponds to the vector index of the cluster in ClusteredWCSP
	
	  int getId() { return id; }

      // ----------------------------------------- Interface Functions
	  WCSP* 		getWCSP() { return wcsp; }

	  bool 			isVar( int i );
	  bool 			isSepVar( int i );

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
	  
	  void 			setParent(Cluster* p);
	  Cluster*		getParent();
	  TVars&		getSep();
	  TClusters&	getDescendants();
	  Cluster*		nextSep( Variable* v ); 
	  bool			isDescendant( Cluster* c2 );
	  
	  
      // ----------------------------------------------------------
	  Cost 			eval(TAssign* a);
	  void 			setWCSP();                              // sets the WCSP to the cluster problem, deconnecting the rest
	  void 			activate();
	  void 			deactivate();
	  void 			increaseLb( Cost newlb );

	  void setup() { sep.setup(this); }

	  
	  void addDelta( int posvar, Value value, Cost cost ) { sep.addDelta(posvar,value,cost); }
	  void nogoodRec( Cost c, bool opt ) { sep.set(c,opt); }	
      Cost nogoodGet( bool& opt ) { return sep.get(opt); }	


	  TVars::iterator beginVars() { return vars.begin(); }
	  TVars::iterator endVars()   { return vars.end(); }
	  TVars::iterator beginSep() { return sep.begin(); }
	  TVars::iterator endSep()   { return sep.end(); }
	  TCtrs::iterator beginCtrs() { return ctrs.begin(); }
	  TCtrs::iterator endCtrs()   { return ctrs.end(); }
	  TClusters::iterator beginEdges() { return edges.begin(); }
	  TClusters::iterator endEdges()   { return edges.end(); }
	  TClusters::iterator beginDescendants() { return descendants.begin(); }
	  TClusters::iterator endDescendants()   { return descendants.end(); }

	  void print();	  
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
	
	Cluster*	getCluster( int i ) { return clusters[i]; }
	Cluster*   	var2Cluster( int v );	
	
	void setCurrentCluster(Cluster *c) {currentCluster = c->getId();}
    Cluster* getCurrentCluster() {return getCluster(currentCluster);}
    bool isInCurrentClusterSubTree(int idc) {return getCurrentCluster()->isDescendant(getCluster(idc));}
	
	void buildFromOrder();	     			    // builds the tree cluster of clusters from a given order
	void fusions();                  			// fusions all redundant clusters after build_from_order is called
	bool fusion();          		            // one fusion step

	int      makeRooted( int icluster );
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

	void print( Cluster* c = NULL, int recnum = 0);
	
};



class ClusteredSolver : public Solver {	
public:

	ClusteredSolver(int storeSize, Cost initUpperBound);


private:


};




#endif 
