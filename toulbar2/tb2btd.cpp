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

/*
 * Choice points
 * 
 */

Cost Solver::binaryChoicePoint(Cluster *cluster, Cost cub, int varIndex, Value value)
{
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
		assert(wcsp->getLb() == cluster->getLbRec());
		wcsp->setUb(cub);
        lastConflictVar = varIndex;
        if (dichotomic) {
    	  if (increasing) decrease(varIndex, middle); else increase(varIndex, middle+1);
    	} else assign(varIndex, value);
        lastConflictVar = -1;
        Cost res = recursiveSolve(cluster, cub);
		cub = min(res,cub);
    } catch (Contradiction) {
        wcsp->whenContradiction();
    }

    store->restore();
    nbBacktracks++;
    try {
        store->store();
		assert(wcsp->getLb() == cluster->getLbRec());
		wcsp->setUb(cub);
		if (dichotomic) {
		  if (increasing) increase(varIndex, middle+1); else decrease(varIndex, middle);
		} else remove(varIndex, value);
		Cost res = recursiveSolve(cluster, cub);
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

Cost Solver::recursiveSolve(Cluster *cluster, Cost cub)
{		
  if (ToulBar2::verbose >= 1) cout << "[" << store->getDepth() << "] recursive solve     cluster: " << cluster->getId() << "     ub: " << cub << "     clb: " << cluster->getLb() << "     wcsp->lb: " << wcsp->getLb() << "     wcsp->ub: " << wcsp->getUb() << endl;
  int varIndex = -1;
  if (ToulBar2::lastConflict) varIndex = getVarMinDomainDivMaxDegreeLastConflict(cluster);
  else varIndex = getVarMinDomainDivMaxDegree(cluster);
  if (varIndex < 0) {
	Cost lb = cluster->getLbRec();
	if (ToulBar2::verbose >= 1) cout << "[" << store->getDepth() << "] C" << cluster->getId() << " lb= " << lb << endl;
	for (TClusters::iterator iter = cluster->beginEdges(); lb < cub && iter!= cluster->endEdges(); ++iter) {
	  Cluster* c = *iter;
	  Cost lbSon = c->getLbRec();
	  Cost ubSon = cub - lb + lbSon;
	  wcsp->getTreeDec()->setCurrentCluster(c);
	  wcsp->setLb(lbSon);
	  wcsp->setUb(ubSon);
	  assert(ubSon > lbSon);
	  wcsp->enforceUb();
	  wcsp->propagate();
	  lb += recursiveSolve(c, ubSon) - lbSon;
	}
	if (ToulBar2::verbose >= 1) cout << "[" << store->getDepth() << "] C" << cluster->getId() << " return " << lb << endl;
	return lb;
  } 
  else {
	if (wcsp->enumerated(varIndex)) {
	  assert(wcsp->canbe(varIndex, wcsp->getSupport(varIndex)));
	  cub = binaryChoicePoint(cluster, cub, varIndex, wcsp->getSupport(varIndex));
	} else {
	  cub = binaryChoicePoint(cluster, cub, varIndex, wcsp->getInf(varIndex));
	}
	if (ToulBar2::verbose >= 1) cout << "[" << store->getDepth() << "] C" << cluster->getId() << " return " << cub << endl;
	return cub;
  }
}
