/*
 * **************** Backtrack with Tree Decomposition *******************
 * 
 */

#include "tb2solver.hpp"
#include "tb2domain.hpp"
#include "tb2pedigree.hpp"
#include "tb2bep.hpp"
#include "tb2clusters.hpp"

/*
 * Variable ordering heuristics
 * 
 */

int Solver::getNextUnassignedVar(Cluster *cluster)
{
  if (unassignedVars->empty()) return -1;
  for (TVars::iterator iter = cluster->beginVars(); iter!= cluster->endVars(); ++iter) {
	if (wcsp->unassigned(*iter)) return *iter;
  }
  return -1;
}

int Solver::getVarMinDomainDivMaxDegree(Cluster *cluster)
{
  if (unassignedVars->empty()) return -1;
  int varIndex = -1;
  Cost worstUnaryCost = MIN_COST;
  double best = MAX_VAL - MIN_VAL;
  
  for (TVars::iterator iter = cluster->beginVars(); iter!= cluster->endVars(); ++iter) {
	if (wcsp->unassigned(*iter)) {
        // remove following "+1" when isolated variables are automatically assigned
        double heuristic = (double) wcsp->getDomainSize(*iter) / (wcsp->getDegree(*iter) + 1);
        if (varIndex < 0 || heuristic < best - 1./100001.
            || (heuristic < best + 1./100001. && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
        }
	}
  }
  return varIndex;   
}

int Solver::getVarSup(Cluster *cluster)
{
  if (unassignedVars->empty()) return -1;
  
  int varIndex = -1;
  double best = wcsp->getUb();    
  for (TVars::iterator iter = cluster->beginVars(); iter!= cluster->endVars(); ++iter) {
	EnumeratedVariable* x = (EnumeratedVariable*) ((WCSP*)wcsp)->getVar(*iter);
	if (wcsp->unassigned(x->wcspIndex)) {
		int dom = wcsp->getDomainSize(x->wcspIndex);
		Cost sumc = MIN_COST; 
		Cost minc = wcsp->getUb();
		int nsup = 0;
		for (EnumeratedVariable::iterator itx = x->begin(); itx != x->end(); ++itx) {        
			Cost cv = x->getCost(*itx);
			sumc += cv;
			if(cv > MIN_COST) { if(minc > cv) minc = cv; }
			else nsup++;
		}        
		double h = (double) dom / (double)(dom - nsup + 1); 
		if(h < best) { best = h; varIndex = x->wcspIndex; }		
    }
  }
  return varIndex;
}


int Solver::getVarMinDomainDivMaxDegreeLastConflict(Cluster *cluster)
{
  if (unassignedVars->empty()) return -1;
	
  if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar) && cluster->isVar(lastConflictVar)) return lastConflictVar;
  
  int varIndex = -1;
  Cost worstUnaryCost = MIN_COST;
  double best = MAX_VAL - MIN_VAL;
    
  for (TVars::iterator iter = cluster->beginVars(); iter!= cluster->endVars(); ++iter) {
	if (wcsp->unassigned(*iter)) {
        // remove following "+1" when isolated variables are automatically assigned
        double heuristic = (double) wcsp->getDomainSize(*iter) / (wcsp->getDegree(*iter) + 1);
        if (varIndex < 0 || heuristic < best - 1./100001.
            || (heuristic < best + 1./100001. && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
        }
    }
  }
  return varIndex;
}



int Solver::getVarFreedom(Cluster *cluster)
{
  if (unassignedVars->empty()) return -1;
	
  if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar) && cluster->isVar(lastConflictVar)) return lastConflictVar;
  
  int varIndex = -1;
  Cost worstUnaryCost = MIN_COST;
  double best = MAX_VAL - MIN_VAL;
    
  for (TVars::iterator iter = cluster->beginVarsTree(); iter!= cluster->endVarsTree(); ++iter) {
	if (wcsp->unassigned(*iter)) {
        // remove following "+1" when isolated variables are automatically assigned
        double heuristic = (double) wcsp->getDomainSize(*iter) / (wcsp->getDegree(*iter) + 1);
        if (varIndex < 0 || heuristic < best - 1./100001.
            || (heuristic < best + 1./100001. && wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
            best = heuristic;
            varIndex = *iter;
            worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
        }
    }
  }
  return varIndex;
}



/*
 * Choice points
 * 
 */

Cost Solver::binaryChoicePoint(Cluster *cluster, Cost lbgood, Cost cub, int varIndex, Value value, Cluster* onlyson, bool freedom)
{
	TreeDecomposition* td = wcsp->getTreeDec();
    assert(wcsp->unassigned(varIndex));
    assert(wcsp->canbe(varIndex,value));
    bool dichotomic = (ToulBar2::dichotomicBranching && ToulBar2::dichotomicBranchingSize < wcsp->getDomainSize(varIndex));
    Value middle = value;
    bool increasing = true;
    if (dichotomic) {
      middle = (wcsp->getInf(varIndex) + wcsp->getSup(varIndex)) / 2;
      if (value <= middle) increasing = true;
      else increasing = false;
    }
    try {
        store->store();
		assert(freedom || td->getCurrentCluster() == cluster);
		assert(freedom || wcsp->getLb() == cluster->getLbRec());
		wcsp->setUb(cub);
		if (CUT(lbgood, cub)) THROWCONTRADICTION;
		if(!freedom && ToulBar2::btdMode == 2) {
			Cost rds = td->getLbRecNoGoodsRDS();
			lbgood = max(rds,lbgood);
			if(CUT(lbgood, cub)) THROWCONTRADICTION;
		}
        lastConflictVar = varIndex;
        if (dichotomic) {
    	  if (increasing) decrease(varIndex, middle); else increase(varIndex, middle+1);
    	} else assign(varIndex, value);
        lastConflictVar = -1;
        Cost res;
        if(!freedom) res = recursiveSolve(cluster, lbgood, cub, onlyson);
        else	     res = recursiveSolveFreedom(cluster, lbgood, cub);
		cub = min(res,cub);
    } catch (Contradiction) {
        wcsp->whenContradiction();
    }
    store->restore();
    nbBacktracks++;
    try {
        store->store();
		assert(freedom || wcsp->getTreeDec()->getCurrentCluster() == cluster);
		assert(freedom || wcsp->getLb() == cluster->getLbRec());
		wcsp->setUb(cub);
		if (CUT(lbgood, cub)) THROWCONTRADICTION;
		if(!freedom && ToulBar2::btdMode == 2) {
			Cost rds = td->getLbRecNoGoodsRDS();
			lbgood = max(rds,lbgood);
			if(CUT(lbgood, cub)) THROWCONTRADICTION;
		}
		if (dichotomic) {
		  if (increasing) increase(varIndex, middle+1); else decrease(varIndex, middle);
		} else remove(varIndex, value);
		Cost res;
		if(!freedom)  res = recursiveSolve(cluster, lbgood, cub, onlyson);
		else		  res = recursiveSolveFreedom(cluster, lbgood, cub);
		cub = min(res,cub);
    } catch (Contradiction) {
	  wcsp->whenContradiction();
    }
    store->restore();
	return cub;
}

/*
 * Backtrack with Tree Decomposition
 * 
 */






Cost Solver::recursiveSolveFreedom(Cluster *cluster, Cost lbgood, Cost cub)
{		
  TreeDecomposition* td = wcsp->getTreeDec();
  if (ToulBar2::verbose >= 1) cout << "[" << store->getDepth() << "] freedom recursive solve     cluster: " << cluster->getId() << "     cub: " << cub << "     clb: " << cluster->getLb() << "     lbgood: " << lbgood << "     wcsp->lb: " << wcsp->getLb() << "     wcsp->ub: " << wcsp->getUb() << endl;
  int varIndex = -1;
  varIndex = getVarFreedom(cluster);
  int nassign = wcsp->numberOfVariables() - wcsp->numberOfUnassignedVariables() - cluster->getNbSepVars();
  if (varIndex < 0 || nassign > 7) {
	Cost lb = wcsp->getLb();
	if (lb < cub) {
		if(cluster == td->getRoot()) {
			//cout << "New solution: " <<  lb << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << store->getDepth() << ")" << endl;
			//td->newSolutionTree();
		}
	}
	return lb;
  } else {
	if (wcsp->enumerated(varIndex)) {
	  assert(wcsp->canbe(varIndex, wcsp->getSupport(varIndex)));
	  cub = binaryChoicePoint(cluster, lbgood, cub, varIndex, wcsp->getSupport(varIndex), NULL, true);
	} else {
	  cub = binaryChoicePoint(cluster, lbgood, cub, varIndex, wcsp->getInf(varIndex), NULL, true);
	}
	if (ToulBar2::verbose >= 1) cout << "[" << store->getDepth() << "] C" << cluster->getId() << " return " << cub << endl;
	return cub;
  }
}









Cost Solver::recursiveSolve(Cluster *cluster, Cost lbgood, Cost cub, Cluster *onlyson)
{		
  TreeDecomposition* td = wcsp->getTreeDec();

  if (ToulBar2::verbose >= 1) cout << "[" << store->getDepth() << "] recursive solve     cluster: " << cluster->getId() << "     cub: " << cub << "     clb: " << cluster->getLb() << "     lbgood: " << lbgood << "     wcsp->lb: " << wcsp->getLb() << "     wcsp->ub: " << wcsp->getUb() << endl;
  int varIndex = -1;
  if (ToulBar2::lastConflict) varIndex = getVarMinDomainDivMaxDegreeLastConflict(cluster);
  else varIndex = getVarMinDomainDivMaxDegree(cluster);
  //else varIndex = getVarSup(cluster);  // experimental
  
  if (varIndex < 0) {
	Cost lb = wcsp->getLb();
	if (ToulBar2::verbose >= 1) cout << "[" << store->getDepth() << "] C" << cluster->getId() << " lb= " << lb << endl;
	for (TClusters::iterator iter = cluster->beginEdges(); lb < cub && iter!= cluster->endEdges(); ++iter) {
	  Cluster* c = *iter;
	  if(onlyson) { if(c != onlyson) continue; }	

	  bool opt = false;
	  Cost lbSon = MIN_COST;
	  bool good = false;
	  if (!c->isActive()) {
		c->reactivate();
		lbSon = c->nogoodGet(opt);
		good = true;
	  } else {
		lbSon = c->getLbRec();
	  }
	  if (!opt) {
		Cost ubSon = cub - lb + lbSon;
		td->setCurrentCluster(c);
		wcsp->setUb(ubSon);
		wcsp->setLb((good)?c->getLbRec():lbSon);
		try {
		  store->store();
		  assert(ubSon > lbSon);
		  wcsp->enforceUb();
		  wcsp->propagate();
		  Cost newlbSon = recursiveSolve(c, lbSon, ubSon, onlyson);
		  c->nogoodRec(newlbSon, (newlbSon < ubSon));
		  assert(newlbSon > lbSon || (newlbSon == lbSon && newlbSon < ubSon));
		  lb += newlbSon - lbSon;
		} catch (Contradiction) {
		  wcsp->whenContradiction();
		  c->nogoodRec(ubSon, false);
		  lb = cub;
		}
		store->restore();
	  }
	}
	if (lb < cub) {
		cluster->solutionRec(lb);
		if(cluster == td->getRoot() || cluster == td->rdsroot) {
			cout << "New solution: " <<  lb << " (" << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << store->getDepth() << ")" << endl;
			if(cluster == td->getRoot())  td->newSolution();
		}
	}
	if (ToulBar2::verbose >= 1) cout << "[" << store->getDepth() << "] C" << cluster->getId() << " return " << lb << endl;
	return lb;
  } 
  else {
	if (wcsp->enumerated(varIndex)) {
	  assert(wcsp->canbe(varIndex, wcsp->getSupport(varIndex)));
	  cub = binaryChoicePoint(cluster, lbgood, cub, varIndex, wcsp->getSupport(varIndex), onlyson);
	} else {
	  cub = binaryChoicePoint(cluster, lbgood, cub, varIndex, wcsp->getInf(varIndex), onlyson);
	}
	if (ToulBar2::verbose >= 1) cout << "[" << store->getDepth() << "] C" << cluster->getId() << " return " << cub << endl;
	return cub;
  }
}




Cost Solver::recursiveSolveRDS(Cluster *cluster)
{	
	TreeDecomposition* td = wcsp->getTreeDec();
	TClusters::iterator it = cluster->beginEdges();
	while(it != cluster->endEdges()) {
		recursiveSolveRDS(*it);
 		++it;
	}
	Cost res; 
    Cost cub = wcsp->getUb();
	try {
	  store->store();
	  cluster->deconnectSep();	
	  td->setCurrentCluster(cluster);
	  cluster->setLb(wcsp->getLb());
	  res = recursiveSolve(cluster, wcsp->getLb(), cub);
	} catch (Contradiction) {
	  wcsp->whenContradiction();
	}
	cub = min(res,cub);
	cluster->setLb_opt(cub);	
	store->restore();	
	
	return cub;
}





void Solver::solveClusters2by2(Cluster *c, Cost cub)
{	
	TreeDecomposition* td = wcsp->getTreeDec();
	TClusters::iterator it = c->beginEdges();
	while(it != c->endEdges()) {
		Cluster* ci = *it;
 		++it;
		if(ci->sepSize() == 0) continue;
	    
		try {
		  store->store();
  	      c->deconnectSep(true);
		 
		  // deconnect brothers
		  TClusters::iterator itj = c->beginEdges();
		  while(itj != c->endEdges()) {
			 Cluster* cj = *itj;
			 if(ci != cj) cj->deconnectSep(false);
		 	 ++itj;
		  }

		  // deconnect sons
		  itj = ci->beginEdges();
		  while(itj != ci->endEdges()) {
			 Cluster* cj = *itj;
			 cj->deconnectSep(false);
		 	 ++itj;
		  }
		  c->setLb(MIN_COST);
		  wcsp->setLb(MIN_COST);
		  wcsp->setUb(cub);
		  td->setCurrentCluster(c);
		  cout << "--- Solving clusters " << c->id << " " << ci->id << " ..." << endl;
		  Cost res = recursiveSolve(c, MIN_COST, cub, ci);
		  ci->printStats();
		  cout << "---  done  cost = " << res << " ("    << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << store->getDepth() << ")" << endl << endl;
		} catch (Contradiction) {
		  wcsp->whenContradiction();
		}
		store->restore();	
		ci->resetOpt();
	}
	
	it = c->beginEdges();
	while(it != c->endEdges()) {
		solveClusters2by2(*it, cub);
		++it;
	}
}



void Solver::solveClusters()
{	
	TreeDecomposition* td = wcsp->getTreeDec();
	int nclusters = td->getNbOfClusters();
	Cost cub = wcsp->getUb();

	for(int i=0; i<nclusters; i++) {
		Cluster* c = td->getCluster(i);
		try {
		  store->store();
  	      c->deconnectSep(true);
	      c->setClusterSep( c->getId() );
		  TClusters::iterator itj = c->beginEdges();
		  while(itj != c->endEdges()) {
			 Cluster* cj = *itj;
			 cj->deconnectSep(false,0);
	         cj->setClusterSep( c->getId() );
		 	 ++itj;
		  }
		  c->setLb(MIN_COST);
		  wcsp->setLb(MIN_COST);
		  wcsp->setUb(cub);
		  td->setCurrentCluster(c);
		  cout << "--- Solve cluster " << c->id << " ..." << endl;
		  Cost res = recursiveSolve(c, MIN_COST, cub, c);
		  cout << "---  done  cost = " << res << " ("    << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << store->getDepth() << ")" << endl << endl;
		  c->setLb_opt(res);
		} catch (Contradiction) {
		  wcsp->whenContradiction();
		}
		store->restore();	
	}
}



void Solver::solveClustersUb()
{	
	TreeDecomposition* td = wcsp->getTreeDec();
	Cost cub = wcsp->getUb();
	Cost cubpart = cub / 100; 
	Cost cubi = cubpart;
	bool nosolution = true;

	while(nosolution) {
		try {
		  store->store();
		  wcsp->setUb(cubi);
		  cout << "--- Solve ub = " << cubi << " ..." << endl;
		  Cost res = recursiveSolve(td->getRoot(), wcsp->getLb(), cubi);
		  td->printStats(td->getRoot());	
		  cout << "---  done  cost = " << res << " ("    << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << store->getDepth() << ")" << endl << endl;
		  nosolution = res >= cubi;
		} catch (Contradiction) {
		  wcsp->whenContradiction();
		}
		store->restore();
		//td->resetOptRec(td->getRoot());	
	    cubi += cubpart;
	}
}


void Solver::solveClustersSubTree(Cluster *c, Cost cub)
{	
	TreeDecomposition* td = wcsp->getTreeDec();

	TClusters::iterator it = c->beginEdges();
	while(it != c->endEdges()) {
		solveClustersSubTree(*it, cub);
		++it;
	}

	try {
	  store->store();
	  
	  if(c != td->getRoot()) {
	      c->deconnectSep(true);
		  c->setLb(MIN_COST);
		  wcsp->setUb(cub - wcsp->getLb());
		  wcsp->setLb(MIN_COST);
	  } else {
		wcsp->setUb(cub);
	  }
	  td->setCurrentCluster(c);
	  td->rdsroot = c;
	  
	  //cub = c->sampling();

	  Cost lbfreedom = 0;
	  
	  /*store->store();
	  cout << "--- Solving freedom cluster subtree " << c->id << " ... "; flush(cout);
	  lbfreedom = recursiveSolveFreedom(c, MIN_COST, cub);	  
  	  store->restore();	
	  cout << "lbfreedom = " << lbfreedom << endl; */

	  cout << "--- Solving cluster subtree " << c->id << " ..." << endl;
	  Cost res = recursiveSolve(c, lbfreedom, cub);
	  c->setLb_opt(res);
	  c->printStatsRec();
	  cout << "---  done  cost = " << res << " ("    << nbBacktracks << " backtracks, " << nbNodes << " nodes, depth " << store->getDepth() << ")" << endl << endl;
	} catch (Contradiction) {
	  wcsp->whenContradiction();
	}
	store->restore();
    c->resetNGSep();
}
