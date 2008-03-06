 
#ifndef TB2CLUSTERS_HPP_
#define TB2CLUSTERS_HPP_

#include "tb2wcsp.hpp"
#include "tb2solver.hpp"
#include "tb2enumvar.hpp"
#include "tb2naryconstr.hpp"

#include <set>
#include <list>

class Cluster;
class ClusteredWCSP;

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
	  ClusteredWCSP*	wcsp;
	  list<TAssign*>    assignments;
	  TVars				vars;
	  TCtrs			    ctrs;
	  TClusters         edges;
	  StoreCost         lb;	
	  Cost				lb_opt;
	  Cost				ub;

 public:
	  Cluster (ClusteredWCSP *w);
	  Cluster ( Cluster& c );
	  ~Cluster();

	  int id;											// the id corresponds to the vector index of the cluster in ClusteredWCSP
	
	  int getId() { return id; }

      // ----------------------------------------- Interface Functions
	  ClusteredWCSP* getWCSP() { return wcsp; }

	  Cost		    getOpt() { return lb_opt; }
	  Cost			getUb()  { return ub; }
	  Cost			getLb()  { return lb; }
	  Cost			getLb_opt()  { return lb_opt; }
	  void			setUb(Cost c)  {ub = c; }
	  void			setLb(Cost c)  {lb = c; }
	  void			setLb_opt(Cost c)  {lb_opt = c; }
	  int			getNbVars() { return vars.size(); }
	  TVars&		getVars() { return vars; }	
	  TClusters&	getEdges() { return edges; }
	  void 			addVar( Variable* x );
	  void 			addVars( TVars& vars );
	  bool 			isVar( int i );
	  void 			addEdge( Cluster* c );
	  void 			addEdges( TClusters& cls );
	  void 			removeEdge( Cluster* c );
	  void 			addCtr( Constraint* c );
	  void 			addAssign( TAssign* a );
	  void 		    updateUb();
      // ----------------------------------------------------------

	  Cost 			eval(TAssign* a);
	  void 			set();                              // sets the WCSP to the cluster problem, deconnecting the rest

	  void 			activate();
	  void 			deactivate();

	  void 			increaseLb( Cost newlb );


	  TVars::iterator beginVars() { return vars.begin(); }
	  TVars::iterator endVars()   { return vars.end(); }
	  TClusters::iterator beginEdges() { return edges.begin(); }
	  TClusters::iterator endEdges()   { return edges.end(); }

	  void print();	  
};



class ClusteredWCSP : public WCSP {
private:
	vector<Cluster*> 		clusters; 

public:

	ClusteredWCSP(Store *s, Cost upperBound);
	
	void build_from_order();	     			// builds the tree cluster of clusters from a given order
	void fusions();                  			// fusions all redundant clusters after build_from_order is called
	bool fusion();                   		    // one fusion step


	void intersection( TVars& v1, TVars& v2, TVars& vout );
	void difference( TVars& v1, TVars& v2, TVars& vout );
	void sum( TVars& v1, TVars& v2, TVars& vout ); 
	bool included( TVars& v1, TVars& v2 ); 	
};



class ClusteredSolver : public Solver {	
public:

	ClusteredSolver(int storeSize, Cost initUpperBound);


private:


};




#endif 
