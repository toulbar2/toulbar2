/** \file tb2clusters.hpp
 *  \brief Cluster Tree Decomposition data-structures.
 * 
 */
 
#ifndef TB2CLUSTERS_HPP_
#define TB2CLUSTERS_HPP_

#include "tb2wcsp.hpp"
#include "tb2solver.hpp"
#include "tb2enumvar.hpp"
#include "tb2naryconstr.hpp"

#include <set>
#include <list>


class Cluster;

typedef set<int>	       TVars;
typedef set<Constraint*>   TCtrs;
typedef map<int,Value>     TAssign;
typedef set<Cluster*>	   TClusters;

typedef pair<Cost,bool>     TPairNG;
typedef pair<Cost,string>   TPairSol;

typedef map<string, TPairNG>  TNoGoods;
typedef map<string, TPairSol> TSols;


class Separator : public AbstractNaryConstraint
{
  private:
	Cluster*					  cluster;
	TVars   					  vars;
    vector< vector<StoreCost> >   delta;    // structure to record the costs that leave the cluster
	StoreInt nonassigned;       			// number of non assigned variables during search
    StoreInt isUsed;
    StoreCost lbPrevious;
    StoreInt optPrevious;

    TNoGoods  					  nogoods;
    TSols  						  solutions;
    DLink<Separator *>            linkSep; // link to insert the separator in PendingSeparator list

	string t;    // temporary buffer for a separator tuple                   
	string s;    // temporary buffer for a solution tuple                   
	
  public:
	Separator(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in);
	Separator(WCSP *wcsp);

    void setup(Cluster* cluster_in);

	TVars& getVars() { return vars; }    
    int getNbVars() { return vars.size();  }
    bool  is(int i) { return vars.find(i) != vars.end(); } 
    TVars::iterator  begin() { return vars.begin(); } 
    TVars::iterator  end()   { return vars.end(); } 

    void set( Cost c, bool opt );
    bool get( Cost& res, bool& opt );

	void solRec( Cost ub );
	bool solGet( TAssign& a, string& sol ); 

	void resetOpt(); 

    void queueSep() { wcsp->queueSeparator(&linkSep); }
    void unqueueSep() { wcsp->unqueueSeparator(&linkSep); }
    
    void addDelta( unsigned int posvar, Value value, Cost cost ) {
    	assert( posvar < vars.size() ); 
    	delta[posvar][value] += cost; 
    }
    Cost getCurrentDelta(); // separator variables may be unassigned

    bool used() { return isUsed; }

	void assign(int varIndex);
    void propagate();

	// NaryConstraint methods not used by Separator class
    double computeTightness() { return 0; }
    bool   verify() {return true;}
    void   increase(int index) {}
    void   decrease(int index) {}
    void   remove(int index) {}
    
    void   print(ostream& os);
};


class Cluster
{
 private:
  TreeDecomposition*  td;
  WCSP*				  wcsp;
  int                 id; // corresponding to the vector index of the cluster in the tree decomposition
  TVars				  vars; // contains all variables inside a cluster including separator variables
  TCtrs			      ctrs; // intermediate usage by bucket elimination (DO NOT USE)
  TClusters           edges; // adjacent clusters (includes parent cluster before makeRooted is done)

  Cluster*			  parent; // parent cluster
  TClusters           descendants; // set of cluster descendants (including itself)
  TVars               varsTree; // set of variables in cluster descendants (including itself)
  vector<bool>        quickdescendants;

  Separator* 	      sep; // associated separator with parent cluster
  StoreCost           lb; // current cluster lower bound deduced by propagation
  Cost				  lbRDS; // global cluster lower bound found by RDS
  StoreInt			  active; // unactive if a nogood including this cluster has been used by propagation

 public:
	  Cluster (TreeDecomposition* tdin);
	  ~Cluster();

	  void          setup(); 

	  int           getId() { return id; }
	  void          setId(int iid) { id=iid; }

	  WCSP* 		getWCSP() { return wcsp; }

      Separator*    getSep() { return sep; }
	  void 			setSep( Separator* sepin ) { sep = sepin; }
	  int			sepSize() { if(sep) return sep->arity(); else return 0; }
      void          deconnectSep(); // deconnect all the constraints on separator variables and assigns separator variables to their support value
      bool 			isSepVar( int i ) { if(!sep) return false; return sep->is(i); }

      bool 			isVar( int i ) { TVars::iterator it = vars.find(i); return it != vars.end(); }
	  int			getNbVars() { return vars.size(); }
	  TVars&		getVars() { return vars; }	
	  TVars&		getVarsTree() { return varsTree; }	
	  void 			addVars( TVars& vars );
	  void 			addVar( Variable* x );
	  void 			removeVar( Variable* x );

	  void 			setParent(Cluster* p) { parent = p; }
	  Cluster*		getParent() { return parent; }
	  TClusters&	getEdges() { return edges; }
	  void 			addEdges( TClusters& cls );
	  void 			addEdge( Cluster* c );
	  void 			removeEdge( Cluster* c );
	  TClusters::iterator removeEdge( TClusters::iterator it );

	  TClusters&	getDescendants() { return descendants; }
	  bool  	    isEdge( Cluster* c );
      void          accelerateDescendants(); 
	  bool			isDescendant( Cluster* c ) { return quickdescendants[c->getId()]; }

	  TCtrs&		getCtrs() { return ctrs; }	
	  void 			addCtrs( TCtrs& ctrsin );
	  void 			addCtr( Constraint* c );

	  bool			isActive() { int a = active; return a == 1; }
	  void 			deactivate();
	  void 			reactivate();

	  Cost			getLb()  { return lb; }
	  void			setLb(Cost c)  { lb = c; }
      void 			increaseLb( Cost addToLb ) { lb += addToLb; }
	  Cost		    getLbRDS() { Cost delta = (sep)?sep->getCurrentDelta():MIN_COST; return max(lbRDS - delta, MIN_COST); }
	  void			setLbRDS(Cost c)  {assert(!sep || sep->getCurrentDelta()==MIN_COST); lbRDS = c; }
	  Cost	        getLbRec();
	  Cost	        getLbRecRDS();
	  
	  void          addDelta( int posvar, Value value, Cost cost ) {assert(sep); sep->addDelta(posvar,value,cost); }

	  void          nogoodRec( Cost c, bool opt ) { if(sep) sep->set(c,opt); }	
      Cost          nogoodGet( bool& opt ) { Cost c = MIN_COST; sep->get(c,opt); return c; }	

      void          resetOptRec(Cluster *rootCluster);

	  void          solutionRec(Cost ub) { if(sep) sep->solRec(ub); }
      void 		    getSolution( TAssign& sol ); // updates sol by the recorded solution found for a separator assignment also given in sol

	  void 			setWCSP2Cluster();   // sets the WCSP to the cluster problem, deconnecting the rest

	  TVars::iterator beginVars() { return vars.begin(); }
	  TVars::iterator endVars()   { return vars.end(); }
	  TVars::iterator beginVarsTree() { return varsTree.begin(); }
	  TVars::iterator endVarsTree()   { return varsTree.end(); }
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

	  void printStatsRec() { 
	  		TClusters::iterator it = beginEdges();
			while(it != endEdges()) {
				(*it)->sep->print(cout);
				(*it)->printStatsRec();
				++it;
			}
	  }
};


class TreeDecomposition  {
private:
    WCSP*			  wcsp;	
    vector<Cluster*>  clusters; 
    list<Cluster*> 	  roots; // intermediate list used by makeRooted method, only one root at the end

    Cluster*          rootRDS; // root cluster of the current RDS iteration

    StoreInt   		  currentCluster; // used to restrict local propagation (NC) and boosting by variable elimination to the current cluster's subtree
    vector<StoreInt>  deltaModified; // accelerator avoiding unecessary checks to delta structure if it is empty (Boolean value)

public:

	TreeDecomposition(WCSP* wcsp_in);

    WCSP* 		getWCSP() { return wcsp; }
	
	int			getNbOfClusters() { return clusters.size(); }
	Cluster*	getCluster( unsigned int i ) { assert( 0 <= i && i < clusters.size() ); return clusters[i]; }
	Cluster*   	var2Cluster( int v );	
	
	void setCurrentCluster(Cluster *c) {currentCluster = c->getId();}
    Cluster* getCurrentCluster() {return getCluster(currentCluster);}

    bool isInCurrentClusterSubTree(int idc); 
    bool isActiveAndInCurrentClusterSubTree(int idc); 
	
    //main function to build a cluster tree/path decomposition:
    // - builds a bucket for each variable following a given variable elimination order
    // - builds a tree/path decomposition from the buckets
    // - associate constraints to clusters, with special treatment for ternary constraints (duplicate flag)
	void buildFromOrder();	     			
	void treeFusions();                  	// merges all redundant clusters
	bool treeFusion();          		    // one fusion step
    void pathFusions(vector<int> &order);   // builds a path decomposition of clusters from a given order

	void fusion( Cluster* ci, Cluster* cj );
	bool reduceHeight( Cluster* c, Cluster *father );
    void splitClusterRec( Cluster* c,  Cluster* father, unsigned int maxsize );
    TVars boostingVarElimRec( Cluster* c,  Cluster* father,  Cluster* grandfather, unsigned int maxsize );
    void mergeClusterRec( Cluster* c,  Cluster* father, unsigned int maxsepsize );
	void heuristicFusionRec( Cluster* c, Cluster* noc );

	void makeDescendants( Cluster* c );
    bool isDescendant( Variable* x, Variable* y ) { return getCluster(x->getCluster())->isDescendant(getCluster(y->getCluster())); }
	
    int      makeRooted(); // defines a rooted cluster tree decomposition from an undirected one
	void     makeRootedRec( Cluster* c,  TClusters& visited );
	Cluster* getBiggerCluster( TClusters& visited );	
	Cluster* getRoot() {return roots.front();}
	Cluster* getRootRDS() {return rootRDS; }
	void setRootRDS(Cluster* rdsroot) { rootRDS = rdsroot; }
	 
	int height( Cluster* r );
	int height( Cluster* r, Cluster* father );

	void intersection( TVars& v1, TVars& v2, TVars& vout );
	void difference( TVars& v1, TVars& v2, TVars& vout );
	void sum( TVars& v1, TVars& v2, TVars& vout ); 
    bool included( TVars& v1, TVars& v2 ); 	// true if v1 is included in v2
	void clusterSum( TClusters& v1, TClusters& v2, TClusters& vout );

    bool isDeltaModified(int varIndex) { return deltaModified[varIndex]; }
	Cost getLbRecRDS() { 
		Cluster* c = getCluster(currentCluster);
		Cost res = c->getLbRecRDS(); 
		return max(res,c->getLbRDS());	
	}
    void addDelta(int c, EnumeratedVariable *x, Value value, Cost cost);
    void newSolution( Cost lb );

	bool verify();
	
	void print( Cluster* c = NULL, int recnum = 0);
    void printStats( Cluster* c = NULL );
};

#endif 
