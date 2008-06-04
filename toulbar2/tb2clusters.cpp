
 
#include "tb2clusters.hpp"
#include "tb2naryconstr.hpp"

#include <list>
#include <algorithm>


/************************************************************************************************/
// NaryNogood

Separator::Separator(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in)
  : AbstractNaryConstraint(wcsp, scope_in, arity_in), 
    nonassigned(arity_in, &wcsp->getStore()->storeValue), 
    isUsed(false, &wcsp->getStore()->storeValue), 
    lbPrevious(MIN_COST, &wcsp->getStore()->storeCost), 
    optPrevious(false, &wcsp->getStore()->storeValue)
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
	if (arity_ == 0) {
	  queueSep();
	} else {
	  for(int i=0;i<arity_;i++) {         
		if (getVar(i)->assigned()) assign(i);
	  }           
	}
}


Separator::Separator(WCSP *wcsp)
			: AbstractNaryConstraint(wcsp), 
			nonassigned(0, &wcsp->getStore()->storeValue), 
			isUsed(false, &wcsp->getStore()->storeValue), 
			lbPrevious(MIN_COST, &wcsp->getStore()->storeCost), 
			optPrevious(false, &wcsp->getStore()->storeValue)
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
		assert(cluster->getParent()->getId() == Constraint::cluster);
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

	int nvars = cluster->getVars().size();
	if(!nvars) return;
	
	char* sbuf = new char [cluster->getVars().size()+1];
	int i = 0;
	int nproper = 0;
	it = cluster->beginVars();
	while(it != cluster->endVars()) {
		if (!cluster->isSepVar(*it)) nproper++;
    	sbuf[i] = CHAR_FIRST;
		++it;	
		i++;
	}
	sbuf[nproper] = '\0';
	s = string(sbuf);
	delete [] sbuf; 
}

void Separator::set( Cost c, bool opt ) { 
	int i = 0;
	WCSP* wcsp = cluster->getWCSP();
	Cost deltares = MIN_COST;
	if (ToulBar2::verbose >= 1) cout << "( ";
	TVars::iterator it = vars.begin();
	while(it != vars.end()) {
	    assert(wcsp->assigned(*it));
	    Value val = wcsp->getValue(*it);
		if (ToulBar2::verbose >= 1) cout << "(" << *it << "," << val << ") ";
		t[i] = val + CHAR_FIRST;
	    deltares += delta[i][val];   // delta structrue
		++it;
		i++;
	}
	assert(!opt || c + deltares >= MIN_COST);
	if (ToulBar2::verbose >= 1) cout << ") Learn nogood " << c << " + delta=" << deltares << "(opt=" << opt << ")" << " on cluster " << cluster->getId() << endl;
	//assert(nogoods.find(string(t)) == nogoods.end() || nogoods[string(t)].second <= max(MIN_COST, c + deltares));
	nogoods[t] = TPairNG(max(MIN_COST, c + deltares), opt); 
}    


Cost Separator::getRDS() {
	int i = 0;
	WCSP* wcsp = cluster->getWCSP();
	Cost sumdelta = MIN_COST;	
	TVars::iterator it = vars.begin();
	while(it != vars.end()) {
		  EnumeratedVariable* x = (EnumeratedVariable*) wcsp->getVar(*it);
		  if (wcsp->td->deltaModified[x->wcspIndex]) {
			Cost del = MIN_COST; 
			Value val;
			if(x->assigned()) { 
			  val = x->getValue();
			  del = delta[i][val];      
			} else {
			  for (EnumeratedVariable::iterator itx = x->begin(); itx != x->end(); ++itx) {
				val = *itx;
				if(del < delta[i][val]) del = delta[i][val];
			  }
			}      
			sumdelta += del;
		  }
		  ++it;
		  i++;
	}
	Cost res = cluster->getOpt() - sumdelta;
  	if(res < MIN_COST) res = MIN_COST;
	return res;
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
		TPairNG p = itng->second;
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

bool Separator::solGet(TAssign& a, string& sol) 
{
	int i = 0;
	TVars::iterator it = vars.begin();
	while(it != vars.end()) {
	  Value val = a[*it];
	  t[i] = val + CHAR_FIRST;	 // build the tuple
	  ++it;
	  i++;
	}
	TPairSol p;
	TSols::iterator itsol = solutions.find(t);
	if(itsol != solutions.end()) {
		p = itsol->second;
		sol = p.second;

		if (ToulBar2::verbose >= 1) cout << "asking  solution  sep:" << t << "  cost: " << p.first << endl; 

		return true;
	}
	return false;
}


void Separator::solRec(Cost ub) 
{
	WCSP* wcsp = cluster->getWCSP();


	int i = 0;
	TVars::iterator it = vars.begin();
	while(it != vars.end()) {
	  assert(wcsp->assigned(*it));
	  Value val = wcsp->getValue(*it);
	  t[i] = val + CHAR_FIRST;	 // build the tuple
	  ++it;
	  i++;
	}
	TPairSol p;
	TSols::iterator itsol = solutions.find(t);
	if(itsol != solutions.end()) {
		p = itsol->second;
	    //assert(p.first < ub);
	}

	wcsp->restoreSolution(cluster);
	
	i = 0;
    it = cluster->beginVars();
	while(it != cluster->endVars()) {
	    assert(wcsp->assigned(*it));
		if (!cluster->isSepVar(*it)) {
			Value val = wcsp->getValue(*it);
			s[i] = val + CHAR_FIRST;
			i++;
		}	
		++it;	
	}

	solutions[t] = TPairSol(ub,s);
	
	if (ToulBar2::verbose >= 1) cout << "recording solution  " << " cost: " << ub << " sol: " << s <<  " sep: " << t << endl; 
}


void Separator::resetOpt() 
{
	TNoGoods::iterator it = nogoods.begin();
	while(it != nogoods.end()) {
		(it->second).second = false;
		++it;
	}
}



void Separator::setCluster(int id)
{
	TVars::iterator it = vars.begin();
	while(it != vars.end()) {
	  Variable* v = wcsp->getVar(*it);
	  v->setCluster(id);
	  ++it;
	}
}

void Separator::print(ostream& os) {
{
	os << this << " nogoods(";
	long totaltuples = 1;
	for(int i = 0; i < arity_;i++) {
		os << scope[i]->wcspIndex;
		if(i < arity_-1) os << ",";
		totaltuples = totaltuples * scope[i]->getDomainInitSize();
	}
	os << ")    ";
	os << " |nogoods| = " << nogoods.size() << " / " << totaltuples;
	if (ToulBar2::verbose >= 4) {
		os << "nogoods: {";
		TNoGoods::iterator  it = nogoods.begin();
		while(it != nogoods.end()) {
			TPairNG p = it->second;
			os << "<" << it->first << "," << p.first << ">";
			if(it != nogoods.end()) os << " "; 
			it++;
		} 
		os << "} " << endl;
	}
	os << endl;
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

TClusters::iterator Cluster::removeEdge( TClusters::iterator it ) { 
	TClusters::iterator itaux = it;
	++it;
	edges.erase(itaux);
	return it;
}

void Cluster::removeEdge( Cluster* c )  {
	TClusters::iterator it = edges.find(c); 
	if(it != edges.end()) edges.erase(it);
}	


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

void Cluster::deconnectSep(bool bassign, int dir) {
	if(!sep) return;
	TVars::iterator its = beginSep();
	while(its != endSep()) {
		Variable *x = wcsp->getVar(*its);
	    ConstraintList* xctrs = x->getConstrs();		
	    for (ConstraintList::iterator it=xctrs->begin(); it != xctrs->end(); ++it) {
            Constraint* ctr = (*it).constr;
            Cluster* ctrc = td->getCluster(ctr->getCluster());
			if (ctr->isSep() && isDescendant(ctrc) && ctrc != this) { continue; }
           	if(bassign || ctr->isSep()) { ctr->deconnect(); continue; }
            if(!dir) { if(isDescendant(ctrc)) ctr->deconnect();  }          
            else     { if(!isDescendant(ctrc)) ctr->deconnect(); }  
	    }
	    if(bassign) x->assign( x->getSupport() );
		++its;
	}
}

void Cluster::resetNGSep() {
	if(!sep) return;
	TVars::iterator its = beginSep();
	while(its != endSep()) {
		Variable *x = wcsp->getVar(*its);
	    ConstraintList* xctrs = x->getConstrs();		
	    for (ConstraintList::iterator it=xctrs->begin(); it != xctrs->end(); ++it) {
            Constraint* ctr = (*it).constr;
           	if(ctr->isSep()) ((Separator*) ctr)->resetOpt();
	    }
		++its;
	}
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

Cost Cluster::getLbRecRDS() {
  assert(isActive());
  Cost res = lb; 
  for (TClusters::iterator iter = beginEdges(); iter!= endEdges(); ++iter) {
	if((*iter)->isActive()) {
		Cost propa = (*iter)->getLbRecRDS();
		Cost rds = (*iter)->sep->getRDS();
		res += max(propa,rds);
	}
  } 
  return res; 
}

Cost Cluster::getLbRecNoGoodsRDS() {
  assert(isActive());
  Cost res = lb;
  for (TClusters::iterator iter = beginEdges(); iter!= endEdges(); ++iter) {
	if(!(*iter)->isActive()) continue;
	Cluster* c = *iter;
	Cost propalb = c->getLbRecRDS();
	Cost rds = c->sep->getRDS();
	res += max(propalb, rds);
  } 
  return max(res,((sep)?sep->getRDS():MIN_COST));	
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


void Cluster::forwardNoGood() 
{
	for (TClusters::iterator itc = beginEdges(); itc!= endEdges(); ++itc) {
	  	Cluster* c = *itc;
		string photo("");
		for (TVars::iterator itv = c->beginSep(); itv!= c->endSep(); ++itv) {
			EnumeratedVariable* x = (EnumeratedVariable*) wcsp->getVar(*itv);
			for(unsigned int v = 0; v < x->getDomainInitSize(); v++) {        
				char buff[10];
				if(x->unassigned() && x->canbe( x->toValue(v) )) {
					//sprintf(buff,"%10d",x->getCost( x->toValue(v) ));
				} else {
					//sprintf(buff,"%10d",-1);
				}
				photo = photo + string(buff);
			}
		}        	
		Separator* s = c->sep;
		TNoGoods::iterator itp = s->forwardNG.find(photo);
		if(itp != s->forwardNG.end()) {
			s->forwardNG[photo] = TPairNG(itp->second.first+1, false);
			cout << "Hit!" << endl;
		}
		else s->forwardNG[photo] = TPairNG(0, false); 		  		
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


Cost Cluster::sampling()
{
	wcsp->getStore()->store();

	TVars::iterator it;
	bool valid = true;
    try { 
		it = beginVars();
		while(it != endVars()) {
			Variable* x = wcsp->getVar(*it);
			if(!isSepVar(*it)) {
  			    wcsp->enforceUb();
				x->assign( x->getSupport() );
		   	    wcsp->propagate();
			}
			++it;
		}
	   for (TClusters::iterator iter = beginEdges(); iter!= endEdges(); ++iter) {
	  	  Cluster* c = *iter;
	  	  if(c->isActive()) {
	  	  	TNoGoods::iterator itng = c->beginNG();
	  	  	while(itng != c->endNG()) {
	  	 
				it = beginSep();
				while(it != endSep()) {
					/*Variable* x = wcsp->getVar(*it);
	  			    wcsp->enforceUb();
					x->assign( it->first[i] - CHAR_FIRST );
			   	    wcsp->propagate();*/
					++it;
				}
	  	 			
	  	  		
	  	  		++itng;
	  	  	}
	  	  }	
	   }
	} 
	catch(Contradiction) { valid = false; }  	 

	Cost res = wcsp->getUb();

	wcsp->getStore()->restore();

	return res;	
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



void Cluster::getSolution( TAssign& sol )
{
	TVars::iterator it;
	if(parent == NULL) {
		if(vars.size() == 0) {
		} else {
		    it = beginVars();
			while(it != endVars()) {
				sol[*it] = wcsp->getValue(*it);
				++it;	
			}			
		}
	}
	string s;
	if(sep) {
		sep->solGet(sol, s);
		int i = 0;
	    it = beginVars();
		while(it != endVars()) {
			if (!isSepVar(*it)) {
				sol[*it] = ((EnumeratedVariable*) wcsp->getVar(*it))->toValue(s[i] - CHAR_FIRST);
				i++;
			}	
			++it;	
		}
	} 
    for (TClusters::iterator iter = beginEdges(); iter!= endEdges(); ++iter) {
		Cluster* cluster = *iter;
		cluster->getSolution(sol);
    } 
}

bool Cluster::isEdge( Cluster* c ) {
    for (TClusters::iterator iter = beginEdges(); iter!= endEdges(); ++iter) {
		Cluster* cluster = *iter;
    	if(c == cluster) return true;
    }
    return false;
}

void Cluster::setup()
{ 
  	if(sep) sep->setup(this); 
  	
    TVars::iterator it = beginVars();
	while(it != endVars()) {
		bool notsep = true;
	    for (TClusters::iterator iter = beginEdges(); iter!= endEdges(); ++iter) {
			Cluster* c = *iter;
			if(c->isSepVar(*it)) notsep = false;
	    }
		if(notsep) varsNotSep.insert(*it);
		++it;	
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
		  //cout << *itp << "C" << wcsp->getVar(*itp)->getCluster() << ",";
		  assert(wcsp->getVar(*itp)->getCluster()==-1 || wcsp->getVar(*itp)->getCluster() == getId());
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

  	/*cout << " ctrs {";
  	TCtrs::iterator itctr = beginCtrs();
  	while(itctr != endCtrs()) {
  	  Constraint* ctr = *itctr;
  	  cout << "( "; 
  	  for(int i=0;i<ctr->arity();i++) cout << ctr->getVar(i)->wcspIndex << " ";
  	  cout << ">C" << ctr->getCluster() << ")"; 
  	  ++itctr;
  	}
  	cout << "}";

	cout << " descendants {";
	TClusters::iterator itd = beginDescendants();
	while(itd != endDescendants()) {
      cout << (*itd)->getId();
      ++itd;
      if(itd != endDescendants()) cout << ",";
	}
	cout << "}";*/
	
	cout << endl;
}

/*****************************************************************************************/
/* ClusteredWCSP																		 */
/*****************************************************************************************/

TreeDecomposition::TreeDecomposition(WCSP* wcsp_in) : wcsp(wcsp_in), 
  currentCluster(-1, &wcsp_in->getStore()->storeValue) {
  deltaModified = vector<StoreInt>(wcsp_in->numberOfVariables(), StoreInt(false, &wcsp_in->getStore()->storeValue));
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

bool TreeDecomposition::isDescendant( Variable* x, Variable* y ) {
	Cluster* cx = getCluster( x->getCluster() );
	Cluster* cy = getCluster( y->getCluster() );
	return cx->isDescendant(cy);
}


void TreeDecomposition::fusions()
{
	while(fusion());

	TClusters visited;
	Cluster* croot = getBiggerCluster(visited);
	fusionRec(croot, croot);


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
	if (ToulBar2::verbose >= 1) cout << "Tree decomposition width: " << treewidth - 1 << endl;
}

// fusion process to minimize tree height when the seprators are included
void TreeDecomposition::fusionRec( Cluster* c, Cluster* noc ) {
	    TClusters::iterator it =  c->beginEdges();
		while(it != c->endEdges()) {
			Cluster* cj = *it;
			++it;
			if(cj == c) continue;
			if(cj == noc) continue;
			fusionRec( cj, c );
		}

		it =  c->beginEdges();
		while(it != c->endEdges()) {
			Cluster* cj = *it;
			++it;
			if(cj == c) continue;
			if(cj == noc) continue;
			TVars varsum;
			TVars varinter;
			sum(c->getVars(), cj->getVars(), varsum);
			intersection(c->getVars(), cj->getVars(), varinter);

			int dif = 2;
			bool bf1 = (varinter.size() > 2) && (varsum.size() <= c->getVars().size() + dif);
			bool bf2 = (varinter.size() > 2) && (varsum.size() <= cj->getVars().size() + dif);
			bool bf3 = (varinter.size() > 100);
			if(bf1 || bf2 || bf3) {
				//fusion(c,cj);
			}
		}
}


void TreeDecomposition::fusion( Cluster* ci, Cluster* cj )
{
	if(!ci) return;
	if(!cj) return;

    if (ToulBar2::verbose >= 1) cout << "fusion: " << ci->getId() << " " << cj->getId() << endl;
			
	ci->addVars(cj->getVars());
	ci->addCtrs(cj->getCtrs());
	ci->addEdges(cj->getEdges());
	TClusters::iterator itk =  cj->beginEdges();
	while(itk != cj->endEdges()) {
		Cluster* ck = *itk;
		ck->removeEdge(cj);
		ck->addEdge(ci);
		++itk;
	}
	ci->removeEdge(ci);
	clusters[ cj->getId() ] = NULL;
	if (ToulBar2::verbose >= 1) { cout << "fusion ci " <<  ci->getId() << ",  cj " <<  cj->getId() << endl; ci->print(); }
	delete cj;
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
			   (c->getId() < cj->getId()) && (cj->getCtrs().size() == 0)) { 
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
        if (ToulBar2::verbose >= 1) cout << "No order file specified... taking index order." << endl;
		for(unsigned int i=0;i<wcsp->numberOfVariables();i++) order.push_back(i);
    } else {    	
	    while(file) {
	    	int ix;
	    	file >> ix;
	    	if(file) order.push_back(ix);
	    } 
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
            	for(int k=0; k < ctr->arity(); k++) if (ctr->getVar(k)->unassigned()) c->addVar( ctr->getVar(k) );
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
	if (ToulBar2::verbose >= 1) cout << "tree height: " << h << endl;

    for (unsigned int i=0; i<wcsp->numberOfConstraints(); i++) {
    	Constraint* ctr = wcsp->getCtr(i);
   		ctr->assignCluster();
    	if (ctr->connected() && !ctr->isSep()) {
    		if(ctr->arity() == 3) {
    			TernaryConstraint* tctr = (TernaryConstraint*) ctr;
				tctr->xy->setCluster( tctr->getCluster() );
				tctr->xz->setCluster( tctr->getCluster() );
				tctr->yz->setCluster( tctr->getCluster() );
    		} 
    	}
    }
    for (int i=0; i<wcsp->elimBinOrder; i++) if (wcsp->elimBinConstrs[i]->connected()) {
    	Constraint* ctr = wcsp->elimBinConstrs[i];
   		ctr->assignCluster();
    	if (ctr->connected() && !ctr->isSep()) {
    		if(ctr->arity() == 3) {
    			TernaryConstraint* tctr = (TernaryConstraint*) ctr;
				tctr->xy->setCluster( tctr->getCluster() );
				tctr->xz->setCluster( tctr->getCluster() );
				tctr->yz->setCluster( tctr->getCluster() );
    		} 
    	}
    }
    for (int i=0; i<wcsp->elimTernOrder; i++) if (wcsp->elimTernConstrs[i]->connected()) {
    	Constraint* ctr = wcsp->elimTernConstrs[i];
   		ctr->assignCluster();
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
    for (int i=0; i<wcsp->elimTernOrder; i++) if (wcsp->elimTernConstrs[i]->connected()) {
    	Constraint* ctr = wcsp->elimTernConstrs[i];
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
    
	if (ToulBar2::verbose >= 1) {
	  print();
	  cout << endl;
	}
	assert(verify());
}

void TreeDecomposition::makeDescendants( Cluster* c )
{
	c->getDescendants().insert(c);
	sum(c->getVarsTree(), c->getVars(), c->getVarsTree());
	TClusters::iterator itj = c->beginEdges();
	while(itj != c->endEdges()) {
		Cluster* cj = *itj;
		++itj;	
		makeDescendants(cj);		
		clusterSum(c->getDescendants(), cj->getDescendants(), c->getDescendants());
		sum(c->getVarsTree(), cj->getVarsTree(), c->getVarsTree());
	}
}



bool TreeDecomposition::reduceHeight( Cluster* c )
{
	bool changed = false;
	Cluster* cparent = c->getParent();
	TClusters::iterator itj;
	if(cparent) {
		itj =  c->beginEdges();
		while(itj != c->endEdges()) {
			Cluster* cj = *itj;
			++itj;	
			TVars cjsep;
   	    	intersection(c->getVars(), cj->getVars(), cjsep);
			if(included(cjsep, cparent->getVars())) {
				// replacing the cluster higher in the tree
				c->removeEdge(cj);
				cparent->addEdge(cj);
				cj->setParent(cparent);
				changed = true;
				return true;
			}
		}
	}
	itj =  c->beginEdges();
	while(itj != c->endEdges()) {
		Cluster* cj = *itj;
		cj->removeEdge(c);
		++itj;	
		changed = changed || reduceHeight(cj); 		
		if(changed) return true;
	}
	return false;
}


void TreeDecomposition::makeRootedRec( Cluster* c,  TClusters& visited )
{
	TClusters::iterator itj =  c->beginEdges();
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
		++itj;	
	}	
}

int TreeDecomposition::makeRooted()
{
	TClusters visited;
	roots.clear();
	Cluster* root;
 
    bool selected = false;
	while(visited.size() < clusters.size()) {
		if(!selected && ToulBar2::btdRootCluster >= 0 && ToulBar2::btdRootCluster < (int)clusters.size()) {
			root = getCluster(ToulBar2::btdRootCluster);
			selected = true;
		} else root = getBiggerCluster(visited);
		
		roots.push_back(root);
		visited.insert(root);
		makeRootedRec(root, visited);
		while(reduceHeight(root));
		makeDescendants(root);
	}

	// if it is a forest
	if(roots.size() > 1) {
		root = new Cluster( this );
		root->id = clusters.size();
		clusters.push_back( root );

	    for (list<Cluster*>::iterator iter = roots.begin(); iter!= roots.end(); ++iter) {
			Cluster* oneroot = *iter;
			
	 	    EnumeratedVariable** scopeVars = new EnumeratedVariable* [1];
	    	oneroot->setSep( new Separator(wcsp,scopeVars,0) );  
			if (oneroot->getVars().size() <= 1 && oneroot->getDescendants().size() == 1) oneroot->sep->unqueueSep();
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
	rdsroot = NULL;
	
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

void TreeDecomposition::newSolution( Cost lb ) 
{
	TAssign a;
		
	Cluster* root = getRoot();
    wcsp->restoreSolution(root);
	root->getSolution( a );	

	if (ToulBar2::elimDegree>0 && root->getVars().size() == 0) {
	  // recorded solutions in clusters containing a single variable eliminated in preprocessing may be wrong due to variable elimination in preprocessing; must be recovered after complete assignment and restoreSolution
	  for (unsigned int i=0; i< wcsp->numberOfVariables(); i++) {
		if (wcsp->enumerated(i)) {
		  EnumeratedVariable *x = (EnumeratedVariable *) wcsp->getVar(i);
		  x->assignWhenEliminated( a[i] );
		}
	  }
	  wcsp->restoreSolution();
	  for (unsigned int i=0; i< wcsp->numberOfVariables(); i++) {
		if (wcsp->enumerated(i)) {
		  a[i] = wcsp->getValue(i);
		}
	  }
	}

    if (ToulBar2::showSolutions) {
		TAssign::iterator it = a.begin();
		while(it != a.end()) {
			Value v = it->second;
			cout << v << " ";
			++it;
		}
		cout << endl;
    }
	    
    if (ToulBar2::writeSolution) {
        ofstream file("sol");
        if (!file) {
          cerr << "Could not write file " << "solution" << endl;
          exit(EXIT_FAILURE);
        }
		TAssign::iterator it = a.begin();
		while(it != a.end()) {
			Value v = it->second;
			file << v << " ";
			++it;
		}
		file << endl;
    }
	
	if(ToulBar2::xmlflag) {
		cout << "o " << lb << endl;
		wcsp->solution_XML();
	} 
	else if(ToulBar2::uai) {
		wcsp->solution_UAI(lb);		
	}
	
    
	
}


void TreeDecomposition::printStats( Cluster* c )
{
	if(!c) return;
	c->printStats();
	TClusters::iterator itj =  c->beginEdges();
	while(itj != c->endEdges()) {
		Cluster* cj = *itj;
		++itj;
		printStats(cj);
	}
}


void TreeDecomposition::resetOptRec( Cluster* c )
{
	TClusters::iterator itj =  c->beginEdges();
	while(itj != c->endEdges()) {
		Cluster* cj = *itj;
		c->resetOpt();
		++itj;
		resetOptRec(cj);
	}
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
	assert(x->clusters.size() > 0);
	if (cost>MIN_COST && !deltaModified[x->wcspIndex]) deltaModified[x->wcspIndex] = true;	  
	x->beginCluster();        	
	while( x->nextCluster(ckid,posx) ) {
	  Cluster* ck = getCluster( ckid );
	  if(ck->isDescendant(cy)) {
		if (ToulBar2::verbose >= 2) cout << "add delta " << cost << " to var " << x->wcspIndex << " (cluster " << cx->getId() << ") value " << value << " from subtree " << ck->getId() << " (cluster " << cyid << ")" << endl;
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
		if(!c) continue;
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



