
 
#include "tb2clusters.hpp"
#include "tb2naryconstr.hpp"

#include <list>
#include <algorithm>


/************************************************************************************************/
// NaryNogood

Separator::Separator(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in)
  : AbstractNaryConstraint(wcsp, scope_in, arity_in), nonassigned(arity_in, &wcsp->getStore()->storeValue), isUsed(false, &wcsp->getStore()->storeValue), lbPrevious(MIN_COST, &wcsp->getStore()->storeCost), optPrevious(false, &wcsp->getStore()->storeValue)
{
	cluster = NULL;
	char* tbuf = new char [arity_in+1];
	tbuf[arity_in] = '\0';
    for(int i=0;i<arity_in;i++) {
    	tbuf[i] = CHAR_FIRST;
    	int domsize = scope_in[i]->getDomainInitSize();
    	vars.insert(  scope_in[i]->wcspIndex );
        if(domsize + CHAR_FIRST > 125) {
		  cerr << "Nary constraints overflow. Try undefine NARYCHAR in makefile." << endl; 
		  exit(EXIT_FAILURE);
		}
    } 	
	t = string(tbuf);
	delete [] tbuf; 
	
    linkSep.content = this;
    
    // initial "delayed" propagation
	for(int i=0;connected() && i<arity_;i++) {         
	  if (getVar(i)->assigned()) assign(i);
	}              
}


Separator::Separator(WCSP *wcsp)
			: AbstractNaryConstraint(wcsp), nonassigned(0, &wcsp->getStore()->storeValue), isUsed(false, &wcsp->getStore()->storeValue), lbPrevious(MIN_COST, &wcsp->getStore()->storeCost), optPrevious(false, &wcsp->getStore()->storeValue)
{
}


void Separator::assign(int varIndex)
{
  if (connected(varIndex)) {
	deconnect(varIndex);	
	nonassigned = nonassigned - 1;
	assert(nonassigned >= 0);
	if(nonassigned == 0) {
	  assert(cluster->isActive());
	  queueSep();
	}
  }
}


void Separator::propagate()
{


  if(nonassigned == 0 && wcsp->getTreeDec()->isInCurrentClusterSubTree(cluster->getParent()->getId())) {
	Cost res = MIN_COST; 
	bool opt = false;
	get(res,opt);
	if (cluster->isActive()) {
	  Cost lbpropa = cluster->getLbRec();
	  Cost lb = res - lbpropa;
	  if(opt || lb>MIN_COST) {
		if (ToulBar2::verbose >= 1) cout << "nogood C" << cluster->getId() << " used in advance (lbpropa=" << lbpropa << " ,lb=" << lb << ")" << endl; 
		assert(lb >= MIN_COST);
		if (opt) unqueueSep();
		// atomic operations:
		isUsed = true;
		cluster->deactivate();
		cluster->getParent()->increaseLb(cluster->getParent()->getLb()+lbpropa);
		if (lb>MIN_COST) projectLB(lb); // project into global lb and into parent cluster
		lbPrevious = res;
		optPrevious = opt;
		// end of atomic operations.
	  }
	} else if (isUsed && cluster->getParent()->isActive()) {
	  if (res > lbPrevious || (opt==true && optPrevious==false)) {
		if (ToulBar2::verbose >= 1) cout << "nogood C" << cluster->getId() << " used in advance (lbPrevious=" << lbPrevious << " ,optPrevious=" << optPrevious << " ,res=" << res << " ,opt=" << opt << ")" << endl; 
		if (opt) unqueueSep();
		// atomic operations:
		if (res > lbPrevious) projectLB(res - lbPrevious); // project into global lb and into parent cluster
		lbPrevious = res;
		optPrevious = opt;
		// end of atomic operations.		
	  }
	}
  }
}



void Separator::setup(Cluster* cluster_in) {
	cluster = cluster_in;
	AbstractNaryConstraint::cluster = cluster_in->getParent()->getId();
	delta.clear();
	TVars::iterator it = vars.begin();
	while(it != vars.end()) {
		EnumeratedVariable* var = (EnumeratedVariable*) cluster->getWCSP()->getVar(*it);
	    delta.push_back( vector<StoreCost>(var->getDomainInitSize(), StoreCost(MIN_COST, &cluster->getWCSP()->getStore()->storeCost)) );
		++it;
	}
}

void Separator::set( Cost c, bool opt ) { 
	int i = 0;
	Cost deltares = MIN_COST;
	if (ToulBar2::verbose >= 1) cout << "( ";
	TVars::iterator it = vars.begin();
	while(it != vars.end()) {
	    assert(cluster->getWCSP()->assigned(*it));
	    Value val = cluster->getWCSP()->getValue(*it);
		if (ToulBar2::verbose >= 1) cout << "(" << *it << "," << val << ") ";
		t[i] = val + CHAR_FIRST;
	    deltares += delta[i][val];   // delta structrue
		++it;
		i++;
	}
	assert(!opt || c + deltares >= MIN_COST);
	if (ToulBar2::verbose >= 1) cout << ") Learn nogood " << c << " + delta=" << deltares << "(opt=" << opt << ")" << " on cluster " << cluster->getId() << endl;
	assert(nogoods.find(string(t)) == nogoods.end() || nogoods[string(t)].second <= max(MIN_COST, c + deltares));
	nogoods[t] = TPair(max(MIN_COST, c + deltares), opt); 
}    

bool Separator::get( Cost& res, bool& opt ) {
	int i = 0;
	res = MIN_COST;
	opt = false;
	
	if (ToulBar2::verbose >= 1) cout << "( ";
	TVars::iterator it = vars.begin();
	while(it != vars.end()) {
	  assert(cluster->getWCSP()->assigned(*it));
	  Value val = cluster->getWCSP()->getValue(*it);
	  if (ToulBar2::verbose >= 1) cout << "(" << *it << "," << val << ") ";
	  t[i] = val + CHAR_FIRST;	 // build the tuple
	  res -= delta[i][val];      // delta structure
	  ++it;
	  i++;
	}
	TNoGoods::iterator itng = nogoods.find(t);
	if(itng != nogoods.end()) {
		TPair p = itng->second;
		if (ToulBar2::verbose >= 1) cout << ") Use nogood " << p.first << ", delta=" << res << " (opt=" << p.second << ") on cluster " << cluster->getId() << " (active=" << cluster->isActive() << ")" << endl;
		assert(!p.second || res + p.first >= MIN_COST);
		res += p.first;
	    opt = p.second;
		res = max(MIN_COST,res);
		return true;
	} else {
	    res = MIN_COST;
		if (ToulBar2::verbose >= 1) cout << ") NOT FOUND for cluser " <<  cluster->getId() << endl;
		return false;
	}
}








/************************************************************************************************/
// Cluster

Cluster::Cluster(TreeDecomposition *tdin) : td(tdin), wcsp(tdin->getWCSP()), lb(MIN_COST, &wcsp->getStore()->storeCost), active( true, &wcsp->getStore()->storeValue )	 
{
	lb     = MIN_COST;
	lb_opt = MIN_COST;
	ub	   = wcsp->getUb();
	sep    = NULL;
	parent = NULL;
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
	if(!sep) return false;
	return sep->is(i);
}


void 	    Cluster::setParent(Cluster* p) { parent = p; }
Cluster*    Cluster::getParent() { return parent; }
TClusters&	Cluster::getDescendants() { return descendants; }

bool Cluster::isDescendant( Cluster* c ) {
	return quickdescendants[c->getId()];
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
}


Cost Cluster::getLbRec() {
  assert(isActive());
  Cost res = lb; 
  for (TClusters::iterator iter = beginEdges(); iter!= endEdges(); ++iter) {
	if((*iter)->isActive()) res += (*iter)->getLbRec();
  } 
  return res; 
}

Cost Cluster::getLbRecNoGoods() {
  assert(isActive());
  Cost res = lb;
  for (TClusters::iterator iter = beginEdges(); iter!= endEdges(); ++iter) {
	if(!(*iter)->isActive()) continue;
	Cost propalb = (*iter)->getLbRec();
    bool dummy_opt;
	Cost nogood = MIN_COST;
	nogood = (*iter)->nogoodGet(dummy_opt);
	res += max(propalb, nogood);
  } 
  return res;	
}

Cost Cluster::getLbRecNoGood(bool& opt) {
  assert(isActive());
  Cost propalb = getLbRec();
  Cost nogood = MIN_COST;
  nogood = nogoodGet(opt);
  return max(propalb, nogood);
}

  
void Cluster::reactivate() {
  assert(!isActive());
  active = true;
  for (TClusters::iterator iter = beginEdges(); iter!= endEdges(); ++iter) {
	assert(!(*iter)->isActive());
	if (!(*iter)->sep->used()) (*iter)->reactivate();  	
  }
}


void Cluster::deactivate() {
  if (isActive()) {
	if (ToulBar2::verbose >= 1) cout << "deactive cluster " << getId() << endl;
	active = false;
	for (TClusters::iterator iter = beginEdges(); iter!= endEdges(); ++iter) {
	  (*iter)->deactivate();  	
	}
  }
}



Cost Cluster::eval( TAssign* a  ) {
	Cost ubold = wcsp->getUb();
	wcsp->getStore()->store();
	wcsp->setLb(MIN_COST);
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
		  //assert(wcsp->getVar(*itp)->getCluster() == getId());
		}
		++itp;
	} 
	cout << "\b}";

	if(sep) {
		cout << " U sep {";
		TVars::iterator its = beginSep();
		while(its != endSep()) {
			cout << *its;
			++its;
			if(its != endSep()) cout << ",";
		}
		cout << "}";
	}

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

bool TreeDecomposition::isInCurrentClusterSubTree(int idc) 
{
	Cluster* ci = getCurrentCluster();
	Cluster* cj = getCluster(idc);
	assert(ci->isActive());
	return ci->isDescendant(cj);
}

bool TreeDecomposition::isActiveAndInCurrentClusterSubTree(int idc) 
{
	Cluster* ci = getCurrentCluster();
	Cluster* cj = getCluster(idc);
	assert(ci->isActive());
	if(!cj->isActive()) return false;
	else return ci->isDescendant(cj);
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
		if (ToulBar2::verbose >= 3) { cout << "fusion testing "; c->print(); }
		
		TClusters::iterator it =  c->beginEdges();
		while(it != c->endEdges()) {
			Cluster* cj = *it;
			assert(cj == clusters[cj->getId()]);			
			if((c->getId() < cj->getId()) && included(c->getVars(), cj->getVars()) ||
			   (c->getId() < cj->getId()) && (cj->getCtrs().size() == 0) ) {
			   	
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
				if (ToulBar2::verbose >= 1) { cout << "fusion ci " <<  c->getId() << ",  cj " <<  cj->getId() << endl; c->print(); }
				delete cj;
				return true;
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

	if(order.size() != wcsp->numberOfVariables()) {
		cout << "Order file " << ToulBar2::varOrder << " has incorrect number of variables." << endl;
		exit(1);
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
	if (ToulBar2::verbose >= 1) { 
		cout << "----- Before fusions process: " << endl;
		for(unsigned int i=0; i < clusters.size(); i++) {
			if(!clusters[i]) continue;
			Cluster* c = clusters[i];
			c->print(); 
		}
		cout << "----- fusions process starting... " << endl;
	}
	fusions();

	for(unsigned int i = 0; i < clusters.size(); i++) {
		Cluster* c = clusters[i];
		c->getDescendants().clear();
	}

	int h = makeRooted();
	cout << "tree height: " << h << endl;
	print();


    for (unsigned int i=0; i<wcsp->numberOfConstraints(); i++) {
    	Constraint* ctr = wcsp->getCtr(i);
    	if (ctr->connected() && !ctr->isSep()) {
    		if(ctr->arity() == 3) {
    			TernaryConstraint* tctr = (TernaryConstraint*) ctr;
				tctr->xy->setCluster( tctr->getCluster() );
				tctr->xz->setCluster( tctr->getCluster() );
				tctr->yz->setCluster( tctr->getCluster() );
    		} 
    	}
    }
	
    for (unsigned int i=0; i<wcsp->numberOfConstraints(); i++) {
    	Constraint* ctr = wcsp->getCtr(i);
    	if (ctr->connected() && !ctr->isSep()) {
    		if(ctr->arity() == 3) {
    			TernaryConstraint* tctr = (TernaryConstraint*) ctr;
    			tctr->setDuplicates();
    			
    			
				assert(tctr->xy->getCluster() == tctr->getCluster() &&
					   tctr->xz->getCluster() == tctr->getCluster() &&
					   tctr->yz->getCluster() == tctr->getCluster() );
    		} 
    	}
    }
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

		TVars cjsep;
		intersection(c->getVars(), cj->getVars(), cjsep);

		//------- Add the constraint separator
		int i = 0;	
		int arity = cjsep.size();
	    EnumeratedVariable** scopeVars = new EnumeratedVariable* [arity];
		TVars::iterator it = cjsep.begin();
		while(it != cjsep.end()) {
			scopeVars[i] = (EnumeratedVariable *) wcsp->getVar(*it);
			++it;
			i++;
		}
    	cj->setSep( new Separator(wcsp,scopeVars,arity) );  
		delete [] scopeVars;
		//------- 
	
		makeRootedRec( cj, visited );
		clusterSum(c->getDescendants(), cj->getDescendants(), c->getDescendants());
		++itj;	
	}	
}

int TreeDecomposition::makeRooted()
{
	TClusters visited;
	roots.clear();
	Cluster* root;
 
	while(visited.size() < clusters.size()) {
		root = getBiggerCluster(visited);
		roots.push_back(root);
		visited.insert(root);
		makeRootedRec(root, visited);
	}

	if(roots.size() > 1) {
		root = new Cluster( this );
		root->id = clusters.size();
		clusters.push_back( root );

	    for (list<Cluster*>::iterator iter = roots.begin(); iter!= roots.end(); ++iter) {
			Cluster* oneroot = *iter;
			root->addEdge(oneroot);
			oneroot->setParent(root);
			root->getDescendants().insert( root );
			clusterSum(root->getDescendants(), oneroot->getDescendants(), root->getDescendants());
	    }
	    roots.clear();
	    roots.push_back( root );
	}
	
	
	for(unsigned int i = 0; i < clusters.size(); i++) {
		Cluster* c = clusters[i];
	
		c->quickdescendants.clear();
		for(unsigned int j = 0; j < clusters.size(); j++) {
			c->quickdescendants.push_back( c->descendants.find( clusters[j] ) != c->descendants.end() );
		}
		
		TCtrs::iterator itctr = c->beginCtrs();
		while(itctr != c->endCtrs()) {
			Constraint* ctr = *itctr;
			ctr->setCluster(c->getId());
			++itctr;
		}
		if(c->sep) c->sep->setSep();
		int posx = 0;
		TVars::iterator itv = c->beginVars();
		while(itv != c->endVars()) {
			Variable* var = wcsp->getVar(*itv);
			if(!c->isSepVar( var->wcspIndex )) var->setCluster(c->getId());
			else {
				var->setSep();
				var->addCluster(c->getId(), posx++);  // we add the cluster and also the position of the variable for the delta stucture
			}
			++itv;
		}
		c->setup();
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
		if(x->assigned()) continue;
		
		Cluster* ci  = clusters[x->getCluster()];
    	if(!ci->isVar(x->wcspIndex) || ci->isSepVar(x->wcspIndex)) {
			cout << "cluster: " << ci->getId() << " , var " << x->wcspIndex << endl;
    		return false;
    	}
    
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
		if (ToulBar2::verbose >= 1) cout << "add delta " << cost << " to var " << x->wcspIndex << " (cluster " << cx->getId() << ") value " << value << " from subtree " << ck->getId() << " (cluster " << cyid << ")" << endl;
		ck->addDelta(posx, value, cost);
	  }
	}
  }
}



Cluster* TreeDecomposition::getBiggerCluster( TClusters& visited ) {
	Cluster* cmax = NULL;
	int maxsize = 0;
	
	for(unsigned int i = 0; i < clusters.size(); i++) {
		Cluster* c = clusters[i];
		if(visited.find(c) == visited.end()) {
			if(c->getNbVars() > maxsize) {
				maxsize = c->getNbVars();
				cmax = c;
			}
		}
	}
	return cmax;
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
