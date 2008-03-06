
 
#include "tb2clusters.hpp"
#include "tb2naryconstr.hpp"

#include <list>
#include <algorithm>


/************************************************************************************************/
// NaryNogood

NaryNogood::NaryNogood(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in)
			: NaryConstrie(wcsp, scope_in, arity_in, 0)
{
	lb = 0;
}


NaryNogood::NaryNogood(WCSP *wcsp)
			: NaryConstrie(wcsp)
{
	lb = 0;
}

void NaryNogood::use() {
	int i = 0;
	char* t = new char [arity_ + 1];			
	for(i = 0; i < arity_;i++) t[i] = CHAR_FIRST + scope[i]->toIndex(scope[i]->getValue());
	t[i] = '\0';
	deconnect();
    wcsp->increaseLb(wcsp->getLb() + eval(string(t)));
	delete [] t;
}


void NaryNogood::assign(int varIndex) {
	int i = 0;
	if (connected(varIndex)) {
        deconnect(varIndex);	
	    nonassigned = nonassigned - 1;	   
	    if(nonassigned == 0) {
			TVars::iterator it = implication.begin();
			bool allunassigned = true;
			while(allunassigned  && (it != implication.end())) {
				allunassigned = allunassigned && scope[i]->unassigned();
				++it;
				i++;
			}

			if(allunassigned) use();
	    }
    }
}



/************************************************************************************************/
// Cluster

Cluster::Cluster(ClusteredWCSP *w) : wcsp(w), lb(MIN_COST, &w->getStore()->storeCost)		 
{
	lb     = MIN_COST;
	lb_opt = MIN_COST;
	ub	   = wcsp->getUb();
}

Cluster::~Cluster() {
}

Cluster::Cluster( Cluster& c ) : lb(MIN_COST, & ((WCSP*)c.getWCSP())->getStore()->storeCost)
{
    wcsp = c.wcsp;
	vars = c.vars;
	ctrs = c.ctrs;
	ub = c.getUb();
	lb = c.getLb();
	assignments = c.assignments;
}


void Cluster::addVar( Variable* x ) { vars.insert(x->wcspIndex); }

void Cluster::addVars( TVars& morevars ) { 
	set_union( vars.begin(), vars.end(),
			   morevars.begin(), morevars.end(),
			   inserter(vars, vars.begin()) );			 	  
}


void Cluster::addCtr( Constraint* c ) { ctrs.push_back(c); }

void Cluster::addEdge( Cluster* c ) { edges.insert( c ); }

void Cluster::addEdges( TClusters& cls ) 
{
	set_union( edges.begin(), edges.end(),
			   cls.begin(), cls.end(),
			   inserter(edges, edges.begin()) );			 	  
}

void Cluster::removeEdge( Cluster* c ) 
{
	TClusters::iterator it = edges.find(c); 
	if(it != edges.end()) edges.erase(it);
}	

void Cluster::addAssign( TAssign* a ) { 
	assignments.push_back(a); 
}


bool Cluster::isVar( int i ) {
	TVars::iterator it = vars.find(i);
	return it != vars.end();
}

void Cluster::updateUb( ) {
	TAssign& a =  * new TAssign;
	for(unsigned int i=0; i < wcsp->numberOfVariables(); i++) {
		EnumeratedVariable* x = (EnumeratedVariable*) wcsp->getVar(i);
		if(isVar(i)) {
			Value v;
			if(x->assigned()) v = x->getValue();
			else			  v = x->getSupport();	
			a[i] = v;
        } 
	}
	addAssign(&a);
}


void Cluster::increaseLb( Cost newlb ) {
	lb = newlb;
	if(newlb > lb_opt) {
		deactivate();
	}
}


void Cluster::activate() {
}



void Cluster::deactivate() {
}



Cost Cluster::eval( TAssign* a  ) {
	Cost ubold = wcsp->getUb();
	wcsp->getStore()->store();
	wcsp->setLb(0);
	set();
	bool valid = true;
	TAssign::iterator it = a->begin();
	while(it != a->end()) {
		EnumeratedVariable* x = (EnumeratedVariable*) wcsp->getVar( it->first );	
		Value v = it->second;
		if(x->unassigned()) {
			 try { 
				 wcsp->enforceUb();
				 x->assign(v);
				 wcsp->propagate();
			 } catch(Contradiction) { valid = false; }  	 
		}
		++it;
	}
	Cost ubcluster = wcsp->getLb();
	wcsp->getStore()->restore();
	wcsp->setUb(ubold);
	if(valid) return ubcluster;
	else return ubold;
}


void Cluster::set() {
	for(unsigned int i=0;i<wcsp->numberOfVariables();i++) {
		if(!isVar(i)) {	
			EnumeratedVariable* x = (EnumeratedVariable*) wcsp->getVar(i);	
			for (ConstraintList::iterator it=x->getConstrs()->begin(); it != x->getConstrs()->end(); ++it) {
				Constraint* ctr = (*it).constr;
				ctr->deconnect();
			}
		}
	}
	for(unsigned int i=0;i<wcsp->numberOfVariables();i++) {
		if(!isVar(i)) {	
			EnumeratedVariable* x = (EnumeratedVariable*) wcsp->getVar(i);	
			x->assign( x->getSupport() );
		}
	}
}

void Cluster::print() {
	cout << "(" << id << ",n:" << getNbVars() << ",lb:" << getLb() << ") ";
}

/*****************************************************************************************/
/* ClusteredWCSP																		 */
/*****************************************************************************************/

ClusteredWCSP::ClusteredWCSP(Store *s, Cost upperBound) : WCSP(s,upperBound) {
}


void ClusteredWCSP::fusions()
{
	while(fusion()) {}


	list<Cluster*> lclusters;
	for(unsigned int i=0; i < clusters.size(); i++) {
		if(clusters[i])	lclusters.push_back( clusters[i] );
	}

	clusters.clear();
		
	int treewidth = 0;
	
	int i = 0;
	list<Cluster*>::iterator it = lclusters.begin();
	while(it != lclusters.end()) {
		Cluster* c = *it;
		clusters.push_back( c );
		c->id = i; 
		++it;
		i++;
		
		if(c->getNbVars() > treewidth) treewidth = c->getNbVars();
	}
	
	cout << "Tree decomposition width: " << treewidth << endl; 
}

bool ClusteredWCSP::fusion()
{
	list<Cluster*> lclusters;

	for(unsigned int i = 0; i < numberOfVariables(); i++) {
		Cluster* c = clusters[i];
		if(c)	lclusters.push_back( c );
	}
	
	list<Cluster*>::iterator it = lclusters.begin();
	while(it != lclusters.end()) {
		Cluster* c = *it;
	
		TClusters::iterator itj =  c->beginEdges();
		while(itj != c->endEdges()) {
			Cluster* cj = *itj;
			if(included(c->getVars(), cj->getVars())) {
				c->addVars(cj->getVars());
				c->removeEdge(cj);
				c->addEdges(cj->getEdges());
				clusters[ cj->getId() ] = NULL;
				return true;
			}
			++itj;	
		}
		++it;
	}

	return false;
}



void ClusteredWCSP::build_from_order()
{
	vector<int> order;

 	ifstream file(ToulBar2::varOrder);
    if (!file) {
        cerr << "Could not open file " << ToulBar2::varOrder << endl;
        exit(EXIT_FAILURE);
    }	
    

    while(file) {
    	int ix;
    	file >> ix;
    	if(file) order.push_back(ix);
    } 
	 
	clusters.clear();
	for(unsigned int i=0;i<numberOfVariables();i++) {
		Cluster* c = new Cluster( this );
		c->id = i;
		c->addVar( getVar( order[i] ) );
		clusters.push_back( c );
	}
	
	
	set<Constraint*> usedctrs;
	
	for(unsigned int i=0;i<numberOfVariables();i++) {
		Variable* x = getVar( order[i] );
		Cluster* c  = clusters[i];
    
	    ConstraintList* xctrs = x->getConstrs();		
	    for (ConstraintList::iterator it=xctrs->begin(); it != xctrs->end(); ++it) {
            Constraint* ctr = (*it).constr;
			bool used = usedctrs.find( ctr ) != usedctrs.end();
            if(!used) {
            	usedctrs.insert( ctr );
            	c->addCtr(ctr);
            	for(int k=0; k < ctr->arity(); k++) c->addVar( ctr->getVar(k) );
            }
	    }
	    
		for(unsigned int j=i+1;j<numberOfVariables();j++) 
		{
			if(c->isVar(order[j])) {
				Cluster* cj  = clusters[j];
				TVars::iterator it = c->beginVars();
				while(it != c->endVars()) { 
					int ivar = *it;
					if(ivar != order[i]) cj->addVar( getVar(ivar) ); 
					++it; 
				}
				c->addEdge( cj );
				break;
			}
		} 
		
		   
	}	

	fusions();
}








void ClusteredWCSP::intersection( TVars& v1, TVars& v2, TVars& vout ) {	
	set_intersection( v1.begin(), v1.end(),
			  	   	  v2.begin(), v2.end(),
				  	  inserter(vout, vout.begin()) );			 	  
}

void ClusteredWCSP::difference( TVars& v1, TVars& v2, TVars& vout ) {	
	set_difference( v1.begin(), v1.end(),
			  	   	v2.begin(), v2.end(),
				  	inserter(vout, vout.begin()) );			 	  
}

void ClusteredWCSP::sum( TVars& v1, TVars& v2, TVars& vout ) {	
	set_union( v1.begin(), v1.end(),
			   v2.begin(), v2.end(),
			   inserter(vout, vout.begin()) );			 	  
}

bool ClusteredWCSP::included( TVars& v1, TVars& v2 ) {	
	TVars vout;
	set_union( v1.begin(), v1.end(),
			   v2.begin(), v2.end(),
			   inserter(vout, vout.begin()) );			 	  
	return vout.size() == v1.size() || vout.size() == v2.size(); 
}





/*****************************************************************************************/
/* Clustered SOLVER																		 */
/*****************************************************************************************/



ClusteredSolver::ClusteredSolver(int storeSize, Cost initUpperBound) : Solver(storeSize, initUpperBound) {
}
