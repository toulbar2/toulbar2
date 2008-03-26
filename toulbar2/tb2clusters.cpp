
 
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


void Cluster::addCtr( Constraint* c ) { ctrs.push_back(c); }

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
	TVars::iterator it = sep.find(i);
	return it != sep.end();
}


void 	    Cluster::setParent(Cluster* p) { parent = p; }
Cluster*    Cluster::getParent() { return parent; }
TVars&	    Cluster::getSep() { return sep; }
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

void Cluster::iniDelta() {
	delta.clear();
	TVars::iterator it = beginSep();
	while(it != endSep()) {
		EnumeratedVariable* var = (EnumeratedVariable*) wcsp->getVar(*it);
	    delta.push_back( vector<StoreCost>(var->getDomainInitSize(), StoreCost(MIN_COST, &wcsp->getStore()->storeCost)) );
		++it;
	}
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
	//cout << "(" << id << ",n:" << getNbVars() << ",lb:" << getLb() << ") ";
	
	cout << "cluster " << getId();
	cout << "    vars {";
	
	TVars::iterator it = beginVars();
	while(it != endVars()) {
		cout << *it;
		++it;
		if(it != endVars()) cout << ",";
	} 
	cout << "}   sep {";

	TVars::iterator its = beginSep();
	while(its != endSep()) {
		cout << *its;
		++its;
		if(its != endSep()) cout << ",";
	}
	
	cout << "}";

	cout << "    proper {";
	TVars::iterator itp = beginVars();
	while(itp != endVars()) {
		if (!isSepVar(*itp)) cout << *itp;
		++itp;
		if(itp != endVars()) cout << ",";
	} 
	cout << "}";

	cout << "    edges {";
	TClusters::iterator itc = beginEdges();
	while(itc != endEdges()) {
		cout << (*itc)->getId();
		++itc;
		if(itc != endEdges()) cout << ",";
	}
	cout << "}";

	cout << "    descendants {";
	TClusters::iterator itd = beginDescendants();
	while(itd != endDescendants()) {
		cout << (*itd)->getId();
		++itd;
		if(itd != endDescendants()) cout << ",";
	}
	cout << "}" << endl;
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
		c->iniDelta();

		int ids = 0;
		TVars::iterator itv = c->beginVars();
		while(itv != c->endVars()) {
			Variable* var = wcsp->getVar(*itv);
			if(!c->isSepVar( var->wcspIndex )) var->setCluster(c->getId());
			else { var->addCluster(c->getId(), ids);
				   ids++;
			}
			++itv;
		}
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


void TreeDecomposition::print( Cluster* c )
{
	if(!c) {
		if(roots.empty()) return;
		c = * roots.begin();
	}

    c->print();	
	

	TClusters::iterator ita = c->beginEdges();
	while(ita != c->endEdges()) {
		print( *ita );
		++ita;
	}

}

/*****************************************************************************************/
/* Clustered SOLVER																		 */
/*****************************************************************************************/



ClusteredSolver::ClusteredSolver(int storeSize, Cost initUpperBound) : Solver(storeSize, initUpperBound) {
}
