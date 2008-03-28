
 
#include "tb2clusters.hpp"
#include "tb2naryconstr.hpp"

#include <list>
#include <algorithm>


/************************************************************************************************/
// NaryNogood

Separator::Separator() 
{
	cluster = NULL;
}

void Separator::setup(Cluster* cluster_in) {
	cluster = cluster_in;
	delta.clear();
	TVars::iterator it = vars.begin();
	while(it != vars.end()) {
		EnumeratedVariable* var = (EnumeratedVariable*) cluster->getWCSP()->getVar(*it);
	    delta.push_back( vector<StoreCost>(var->getDomainInitSize(), StoreCost(MIN_COST, &cluster->getWCSP()->getStore()->storeCost)) );
		++it;
	}
}

void Separator::set( Cost c, bool opt ) { 
	int arity = vars.size();
	char* t = new char [arity+1];			
	int i = 0;
	Cost deltares = 0;
	TVars::iterator it = vars.begin();
	while(it != vars.end()) {
	    assert(cluster->getWCSP()->assigned(*it));
	    Value val = cluster->getWCSP()->getValue(*it);
		t[i] = val;
	    deltares += delta[i][val];   // delta structrue
		++it;
		i++;
	}
	t[i] = '\0';
	assert(!opt || c + deltares >= 0);
	if (ToulBar2::verbose >= 1) cout << "Learn nogood " << c << " + " << deltares << " on cluster " << cluster->getId() << endl;
	nogoods[string(t)] = TPair(max(0, c + deltares), opt); 
	delete [] t;
}    

Cost Separator::get( bool& opt ) {
	int arity = vars.size();
	char* t = new char [arity+1];			
	int i = 0;
	Cost res = 0;
	TVars::iterator it = vars.begin();
	while(it != vars.end()) {
	  assert(cluster->getWCSP()->assigned(*it));
	  Value val = cluster->getWCSP()->getValue(*it);
	  t[i] = (char) val;	  // build the tuple
	  res -= delta[i][val];   // delta structrue
	  ++it;
	  i++;
	}
	t[i] = '\0';
	TPair p = nogoods[string(t)];
	if (ToulBar2::verbose >= 1) cout << "Use nogood " << res << "," << p.first << "," << p.second << " on cluster " << cluster->getId() << endl;
	res += p.first;
    opt = p.second;
	delete [] t;
	return max(0,res);
}








/************************************************************************************************/
// Cluster

Cluster::Cluster(TreeDecomposition *tdin) : td(tdin), wcsp(tdin->getWCSP()), lb(MIN_COST, &wcsp->getStore()->storeCost)	 
{
	lb     = MIN_COST;
	lb_opt = MIN_COST;
	ub	   = wcsp->getUb();
}

Cluster::~Cluster() {
}



void Cluster::addVar( Variable* x ) { vars.insert(x->wcspIndex); }
void Cluster::removeVar( Variable* x ) { vars.erase(x->wcspIndex); }

void Cluster::addVars( TVars& morevars ) { 
	set_union( vars.begin(), vars.end(),
			   morevars.begin(), morevars.end(),
			   inserter(vars, vars.begin()) );			 	  
}


void Cluster::addCtr( Constraint* c ) { ctrs.insert(c); }

void Cluster::addEdge( Cluster* c ) { edges.insert(c); }

void Cluster::addEdges( TClusters& cls ) 
{
	set_union( edges.begin(), edges.end(),
			   cls.begin(), cls.end(),
			   inserter(edges, edges.begin()) );			 	  
}

void Cluster::addCtrs( TCtrs& ctrsin ) {
	set_union( ctrs.begin(), ctrs.end(),
			   ctrsin.begin(), ctrsin.end(),
			   inserter(ctrs, ctrs.begin()) );			 	  
}

void Cluster::removeEdge( Cluster* c )  {
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

bool Cluster::isSepVar( int i ) {
	return sep.is(i);
}


void 	    Cluster::setParent(Cluster* p) { parent = p; }
Cluster*    Cluster::getParent() { return parent; }
TVars&	    Cluster::getSep() { return sep.getVars(); }
TClusters&	Cluster::getDescendants() { return descendants; }

bool Cluster::isDescendant( Cluster* c ) {
	return descendants.find( c ) != descendants.end(); 
}

Cluster* Cluster::nextSep( Variable* v ) { 
	Cluster* c = getParent();
	if(!c) return NULL;
	if(c->isVar( v->wcspIndex ) ) return c;
	else return NULL;
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


Cost Cluster::getLbRec() {
  Cost res = lb; 
  for (TClusters::iterator iter = beginEdges(); iter!= endEdges(); ++iter) {
	res += (*iter)->getLbRec();
  } 
  return res; 
}

Cost Cluster::getLbRecNoGoods() {
  Cost res = lb;
  for (TClusters::iterator iter = beginEdges(); iter!= endEdges(); ++iter) {
	Cost propalb = (*iter)->getLbRec();
    bool dummy_opt;
	Cost nogood = (*iter)->nogoodGet(dummy_opt);
	res += max(propalb, nogood);
  } 
  return res;	
}

Cost Cluster::getLbRecNoGood(bool& opt) {
  Cost propalb = getLbRec();
  Cost nogood = nogoodGet(opt);
  return max(propalb, nogood);
}

  
void Cluster::activate() {
}



void Cluster::deactivate() {
}



Cost Cluster::eval( TAssign* a  ) {
	Cost ubold = wcsp->getUb();
	wcsp->getStore()->store();
	wcsp->setLb(0);
	setWCSP();
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


void Cluster::setWCSP() {
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
	//cout << "(" << id << ",n:" << getNbVars() << ",lb:" << getLb() << ") ";
	
	cout << "cluster " << getId();

	cout << " vars {";
	TVars::iterator itp = beginVars();
	while(itp != endVars()) {
		if (!isSepVar(*itp)) {
		  cout << *itp << ",";
		  assert(wcsp->getVar(*itp)->getCluster() == getId());
		}
		++itp;
	} 
	cout << "\b}";

	cout << " U sep {";
	TVars::iterator its = beginSep();
	while(its != endSep()) {
		cout << *its;
		++its;
		if(its != endSep()) cout << ",";
	}
	cout << "}";


	if (!edges.empty()) {
	  cout << " sons {";
	  TClusters::iterator itc = beginEdges();
	  while(itc != endEdges()) {
		cout << (*itc)->getId();
		++itc;
		if(itc != endEdges()) cout << ",";
	  }
	  cout << "}";
	}

//  	cout << " ctrs {";
//  	TCtrs::iterator itctr = beginCtrs();
//  	while(itctr != endCtrs()) {
//  	  Constraint* ctr = *itctr;
//  	  cout << "( "; 
//  	  for(int i=0;i<ctr->arity();i++) cout << ctr->getVar(i)->wcspIndex << " ";
//  	  cout << ")"; 
//  	  ++itctr;
//  	}
//  	cout << "}";

//    	cout << " descendants {";
//    	TClusters::iterator itd = beginDescendants();
//    	while(itd != endDescendants()) {
//  	  cout << (*itd)->getId();
//  	  ++itd;
//  	  if(itd != endDescendants()) cout << ",";
//    	}
//    	cout << "}";
	cout << endl;
}

/*****************************************************************************************/
/* ClusteredWCSP																		 */
/*****************************************************************************************/

TreeDecomposition::TreeDecomposition(WCSP* wcsp_in) : wcsp(wcsp_in), 
  currentCluster(-1, &wcsp_in->getStore()->storeValue) {
}

void TreeDecomposition::fusions()
{
	while(fusion());

 	int treewidth = 0;	
	set<Cluster*> sclu;
	for(unsigned int i=0; i < clusters.size(); i++) {
		if(clusters[i])	{
			Cluster* c = clusters[i];
			sclu.insert( c );
	   	    if(c->getNbVars() > treewidth) treewidth = c->getNbVars();
  
		}
	}
	
	int i = 0;
	clusters.clear();
	set<Cluster*>::iterator it = sclu.begin();
	while(it != sclu.end()) {
		Cluster* c = *it;
		c->id = i++;
		clusters.push_back(*it);
		++it;
	}

	cout << "Tree decomposition width: " << treewidth - 1 << endl;
}

bool TreeDecomposition::fusion( )
{
	bool done = false;
	for(unsigned int i=0; i < clusters.size(); i++) {
		if(!clusters[i]) continue;
		Cluster* c = clusters[i];
		TClusters::iterator it =  c->beginEdges();
		while(it != c->endEdges()) {
			Cluster* cj = *it;
			if((c->getId() < cj->getId()) && included(c->getVars(), cj->getVars())) {
				c->addVars(cj->getVars());
				c->addCtrs(cj->getCtrs());
				c->addEdges(cj->getEdges());
				TClusters::iterator itk =  cj->beginEdges();
				while(itk != cj->endEdges()) {
					Cluster* ck = *itk;
					ck->removeEdge(cj);
					ck->addEdge(c);
					++itk;
				}
				c->removeEdge(c);
				clusters[ cj->getId() ] = NULL;
				delete cj;
				done = true;
			}
			++it;	
		}
	}
	return done;
}



void TreeDecomposition::buildFromOrder()
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
	 
	if(clusters.size() > 0) {
		for(unsigned int i=0;i<clusters.size();i++) {
			Cluster* c = clusters[i];
			if(c) delete c;
		}
	} 
	clusters.clear();

	for(unsigned int i=0;i<wcsp->numberOfVariables();i++) {
		Cluster* c = new Cluster( this );
		c->id = i;
		c->addVar( wcsp->getVar( order[i] ) );
		clusters.push_back( c );
	}
	
	
	set<Constraint*> usedctrs;
	
	for(unsigned int i=0;i<wcsp->numberOfVariables();i++) {
		Variable* x = wcsp->getVar( order[i] );
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
	    
		for(unsigned int j=i+1;j<wcsp->numberOfVariables();j++) 
		{
			if(c->isVar(order[j])) {
				Cluster* cj  = clusters[j];
				TVars::iterator it = c->beginVars();
				while(it != c->endVars()) { 
					cj->addVar( wcsp->getVar(*it) ); 
					++it; 
				}
				cj->removeVar(x);
				c->addEdge( cj );
				cj->addEdge( c );
				break;
			}
		} 
		
		   
	}	
	fusions();
	int h = makeRooted(0);
	cout << "tree height: " << h << endl;
	print();

	assert(verify());
}





void TreeDecomposition::makeRootedRec( Cluster* c,  TClusters& visited )
{
	TClusters::iterator itj =  c->beginEdges();
	c->getDescendants().insert(c);
	while(itj != c->endEdges()) {
		Cluster* cj = *itj;
		cj->removeEdge(c);
		cj->setParent(c);
		visited.insert(cj);
		intersection(c->getVars(), cj->getVars(), cj->getSep());
		makeRootedRec( cj, visited );
		clusterSum(c->getDescendants(), cj->getDescendants(), c->getDescendants());
		++itj;	
	}	
}

int TreeDecomposition::makeRooted( int icluster )
{
	Cluster* root = clusters[icluster];
	roots.clear();
	roots.push_back(root);

	for(unsigned int i = 0; i < clusters.size(); i++) {
		Cluster* c = clusters[i];
		c->getDescendants().clear();
		c->getSep().clear();
	}

	TClusters visited;
	root->setParent( NULL );
	visited.insert(root);
	makeRootedRec(root, visited);


	for(unsigned int i = 0; i < clusters.size(); i++) {
		Cluster* c = clusters[i];

		TCtrs::iterator itctr = c->beginCtrs();
		while(itctr != c->endCtrs()) {
			Constraint* ctr = *itctr;
			ctr->setCluster(c->getId());
			++itctr;
		}

		int posx = 0;
		TVars::iterator itv = c->beginVars();
		while(itv != c->endVars()) {
			Variable* var = wcsp->getVar(*itv);
			if(!c->isSepVar( var->wcspIndex )) var->setCluster(c->getId());
			else var->addCluster(c->getId(), posx++);  // we add the cluster and also the position of the variable for the delta stucture
			++itv;
		}
		
		c->setup();
	}


	if(visited.size() < clusters.size()) {
		// it was a forest and not a tree
	}

	return height(root);
}



int TreeDecomposition::height( Cluster* r, Cluster* father )
{
	int maxh = 0;
	TClusters::iterator it = r->beginEdges();
	while(it != r->endEdges()) {
		Cluster* adjr = *it;
		if(adjr != father) {
			int h = height(adjr,r);
			if(h > maxh) maxh = h;
		}
		++it;
	}
	return maxh + 1;
}


int TreeDecomposition::height( Cluster* r )
{
	int maxh = 0;
	TClusters::iterator it = r->beginEdges();
	while(it != r->endEdges()) {
		int h = height(*it,r);
		if(h > maxh) maxh = h;
		++it;
	}
	return maxh + 1;		
}


bool TreeDecomposition::verify()
{
	for(unsigned int i=0;i<wcsp->numberOfVariables();i++) {
		Variable* x = wcsp->getVar( i );
		Cluster* ci  = clusters[x->getCluster()];
    
    	if( ! (ci->isVar(x->wcspIndex) && !ci->isSepVar(x->wcspIndex)) ) return false;
    
//  	    ConstraintList* xctrs = x->getConstrs();		
//  	    for (ConstraintList::iterator it=xctrs->begin(); it != xctrs->end(); ++it) {
//              Constraint* ctr = (*it).constr;
//  			Cluster* cj  = clusters[ctr->getCluster()];
//              int arity = ctr->arity();
//              for(i=0;i<arity;i++) {
//          		Variable* x = ctr->getVar(i);
        		    	
//              }
//  	    }
	}
	return true;
}



void TreeDecomposition::intersection( TVars& v1, TVars& v2, TVars& vout ) {	
	set_intersection( v1.begin(), v1.end(),
			  	   	  v2.begin(), v2.end(),
				  	  inserter(vout, vout.begin()) );			 	  
}

void TreeDecomposition::difference( TVars& v1, TVars& v2, TVars& vout ) {	
	set_difference( v1.begin(), v1.end(),
			  	   	v2.begin(), v2.end(),
				  	inserter(vout, vout.begin()) );			 	  
}

void TreeDecomposition::sum( TVars& v1, TVars& v2, TVars& vout ) {	
	set_union( v1.begin(), v1.end(),
			   v2.begin(), v2.end(),
			   inserter(vout, vout.begin()) );			 	  
}

bool TreeDecomposition::included( TVars& v1, TVars& v2 ) {	
	TVars vout;
	set_union( v1.begin(), v1.end(),
			   v2.begin(), v2.end(),
			   inserter(vout, vout.begin()) );			 	  
	return vout.size() == v1.size() || vout.size() == v2.size(); 
}

void TreeDecomposition::clusterSum( TClusters& v1, TClusters& v2, TClusters& vout ) {	
	set_union( v1.begin(), v1.end(),
			   v2.begin(), v2.end(),
			   inserter(vout, vout.begin()) );			 	  
}

void TreeDecomposition::addDelta(int cyid, EnumeratedVariable *x, Value value, Cost cost)
{
  Cluster* cy = getCluster( cyid );
  Cluster* cx = getCluster( x->getCluster() );
  if(! cy->isDescendant( cx ) ) {
	int ckid,posx;
	x->beginCluster();        	
	while( x->nextCluster(ckid,posx) ) {
	  Cluster* ck = getCluster( ckid );
	  if(ck->isDescendant(cy)) {
		if (ToulBar2::verbose >= 1) cout << "add delta " << cost << " to " << x->wcspIndex << "," << value << " on cluster " << ck->getId() << endl;
		ck->addDelta(posx, value, cost);
	  }
	}
  }
}

void TreeDecomposition::print( Cluster* c, int recnum )
{
	if(!c) {
//  		for(unsigned int i=0;i<wcsp->numberOfVariables();i++) {
//  			Variable* x = wcsp->getVar(i);
//  			x->beginCluster();
//  			int c,posx;
//  			cout << x->wcspIndex << " appears in sep {";
//  			while(x->nextCluster(c,posx)) {
//  				cout << c << " ";
//  			}
//  			cout << "}" << endl;
//  		} 
		if(roots.empty()) return;
		c = * roots.begin();
	}

	for (int i=0; i<recnum; i++) cout << "  ";
    c->print();	
	

	TClusters::iterator ita = c->beginEdges();
	while(ita != c->endEdges()) {
		print( *ita, recnum+1 );
		++ita;
	}

}

/*****************************************************************************************/
/* Clustered SOLVER																		 */
/*****************************************************************************************/



ClusteredSolver::ClusteredSolver(int storeSize, Cost initUpperBound) : Solver(storeSize, initUpperBound) {
}
