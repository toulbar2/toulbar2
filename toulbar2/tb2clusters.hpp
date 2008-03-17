 
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
typedef list<Constraint*>  TCtrs;
typedef map<int,Value>     TAssign;
typedef set<Cluster*>	   TClusters;


#define NOCLUSTER -1


class NaryNogood : public NaryConstrie
{
  private:

    TVars implication;
	void use(); 
	Cost lb;

  public:

	NaryNogood(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in);
	NaryNogood(WCSP *wcsp);
	void assign(int varIndex);
	TVars& getImplication() { return implication; }
};



class Cluster {

 private:
  	  TreeDecomposition*  td;
	  WCSP*				  wcsp;
	  list<TAssign*>      assignments;
	  TVars				  vars;
	  TCtrs			      ctrs;
	  TClusters           edges;
	  StoreCost           lb;	
	  Cost				  lb_opt;
	  Cost				  ub;
	  
	  vector< vector<StoreCost> >   delta;    // structure to record the costs that leave the cluster
	  										  // inicialized with iniDelta()
	  
	  int				  parent;
	  TVars				  sep;                // separator vars with parent cluster
      TClusters    	      ancestors;  	      // set of ancestors	

 public:
	  Cluster (TreeDecomposition* tdin);
	  ~Cluster();

	  int id;								  // the id corresponds to the vector index of the cluster in ClusteredWCSP
	
	  int getId() { return id; }

      // ----------------------------------------- Interface Functions
	  WCSP* 		getWCSP() { return wcsp; }

	  bool 			isVar( int i );
	  bool 			isSepVar( int i );


	  Cost		    getOpt() { return lb_opt; }
	  Cost			getUb()  { return ub; }
	  Cost			getLb()  { return lb; }
	  Cost			getLbRec();
	  Cost			getLb_opt()  { return lb_opt; }
	  void			setUb(Cost c)  {ub = c; }
	  void			setLb(Cost c)  {lb = c; }
	  void			setLb_opt(Cost c)  {lb_opt = c; }
	  int			getNbVars() { return vars.size(); }
	  TVars&		getVars() { return vars; }	
	  TCtrs&		getCtrs() { return ctrs; }	
	  TClusters&	getEdges() { return edges; }
	  void 			addVar( Variable* x );
	  void 			addVars( TVars& vars );
	  void 			addEdge( Cluster* c );
	  void 			addEdges( TClusters& cls );
	  void 			removeEdge( Cluster* c );
	  void 			addCtrs( TCtrs& ctrsin );
	  void 			addCtr( Constraint* c );
	  void 			addAssign( TAssign* a );
	  
	  void 			iniDelta();
	  
	  void 		    updateUb();
	  
	  
	  
	  void 			setParent(int p);
	  Cluster*		getParent();
	  TVars&		getSep();
	  TClusters&	getAncestors();
	  Cluster*		nextSep( Variable* v ); 
	  bool			isAncestor( Cluster* c2 );
      // ----------------------------------------------------------

	  Cost 			eval(TAssign* a);
	  void 			set();                              // sets the WCSP to the cluster problem, deconnecting the rest

	  void 			activate();
	  void 			deactivate();

	  void 			increaseLb( Cost newlb );


	  TVars::iterator beginVars() { return vars.begin(); }
	  TVars::iterator endVars()   { return vars.end(); }
	  TVars::iterator beginSep() { return sep.begin(); }
	  TVars::iterator endSep()   { return sep.end(); }
	  TCtrs::iterator beginCtrs() { return ctrs.begin(); }
	  TCtrs::iterator endCtrs()   { return ctrs.end(); }
	  TClusters::iterator beginEdges() { return edges.begin(); }
	  TClusters::iterator endEdges()   { return edges.end(); }

	  void print();	  
};



class TreeDecomposition  {
private:
	WCSP*			  wcsp;	
    Cluster*    	  currentCluster;
	vector<Cluster*>  clusters; 
	list<Cluster*> 	  roots;

public:

	TreeDecomposition(WCSP* wcsp_in);

    WCSP* 		getWCSP() { return wcsp; }
	
	Cluster*	getCluster( int i ) { return clusters[i]; }
	Cluster*   	var2Cluster( int v );	
	
  void setCurrentCluster(Cluster *c) {currentCluster = c;}
  Cluster* getCurrentCluster() {return currentCluster;}
  bool isInCurrentClusterSubTree(Cluster* c, Cluster* rec);
  bool isInCurrentClusterSubTree(int idc) {return isInCurrentClusterSubTree(clusters[idc], currentCluster);}
	
	void buildFromOrder();	     			    // builds the tree cluster of clusters from a given order
	void fusions();                  			// fusions all redundant clusters after build_from_order is called
	bool fusion();                   		    // one fusion step

	int  makeRooted( int icluster );
	void makeRootedRec( Cluster* c,  TClusters& visited );
  Cluster* getRoot() {return roots.front();}

	int height( Cluster* r, Cluster* father );
	int height( Cluster* r );



	void intersection( TVars& v1, TVars& v2, TVars& vout );
	void difference( TVars& v1, TVars& v2, TVars& vout );
	void sum( TVars& v1, TVars& v2, TVars& vout ); 
	bool included( TVars& v1, TVars& v2 ); 	

	void clusterSum( TClusters& v1, TClusters& v2, TClusters& vout );	
	
	void print( Cluster* c = NULL);
	
};



class ClusteredSolver : public Solver {	
public:

	ClusteredSolver(int storeSize, Cost initUpperBound);


private:


};




#endif 
