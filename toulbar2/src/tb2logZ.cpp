#include "tb2solver.hpp"
#include "tb2trienum.hpp"
#include "tb2domain.hpp"
#include "tb2binconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2clusters.hpp"

TLogProb Solver::Zub()  // Calculate an uper-bound on Z before exploration (step 0)
{

  TLogProb newlogU=-numeric_limits<TLogProb>::infinity();
  vector<Cost> vbinmin(unassignedVars->getSize(), 0);
  Cost newCost = wcsp->getLb() + wcsp->getNegativeLb();
  switch (ToulBar2::isZUB) {
    //Upper bound on Z edition 0
  case 0 :
    newCost += wcsp->LogProb2Cost(unassignedVars->getSize() * Log(wcsp->getMaxDomainSize()));
    newlogU = wcsp->LogSumExp(ToulBar2::logU, newCost);
    break;
    //Upper bound on Z edition 1
  case 1 :
    //cout<<"p0 : "<<wcsp->Cost2Prob(newCost)<<endl;
    for (BTList<Value>::iterator iter_variable = unassignedVars->begin(); iter_variable != unassignedVars->end(); ++iter_variable) { // Loop on the unassigned variables
      EnumeratedVariable *var = (EnumeratedVariable *)((WCSP *) wcsp)->getVar(*iter_variable);
      //cout<<"Variable "<<var->getName()<<endl;
      Cost SumUnaryCost = MAX_COST;
      if (wcsp->enumerated(*iter_variable)) {
	for (EnumeratedVariable::iterator iter_value = var->begin(); iter_value != var->end(); ++iter_value) { // loop over the domain of the variable
	  SumUnaryCost = wcsp->LogSumExp(SumUnaryCost, var->getCost(*iter_value)); // Sum of the exponential of Unary cost over the domain.
	  //cout<<"UNARYZ : "<<wcsp->Cost2Prob(var->getCost(*iter_value))<<endl;
	}
	newCost += SumUnaryCost; //Sum the older cost with the new log(sum(unarycost))
	//cout<<"NewCost : "<<wcsp->Cost2Prob(newCost)<<endl;
      } else {
	newCost += wcsp->LogProb2Cost(Log(wcsp->getDomainSize(*iter_variable)));
      }
    }
    newlogU = wcsp->LogSumExp(ToulBar2::logU, newCost);
    break;
    
    //Upper bound on Z edition 2
  case 2 :
    newlogU = wcsp->LogSumExp(ToulBar2::logU, wcsp->spanningTreeZ(newCost));
    break;
  }
  
  return newlogU;
}


TLogProb Solver::MeanFieldZ() // Compute a lower bound on logZ using MF approximation
{
	//~ cout<<"Cout "<<wcsp->getLb()<< " "<<wcsp->getNegativeLb()<<" "<<ToulBar2::markov_log<<endl;
	//~ cout<<"Log "<<wcsp->Cost2LogProb(wcsp->getLb())<< " "<<wcsp->Cost2LogProb(wcsp->getNegativeLb())<<" "<<ToulBar2::markov_log<<endl;
	//~ cout<<"Prob "<<wcsp->Cost2Prob(wcsp->getLb())<< " "<<wcsp->Cost2Prob(wcsp->getNegativeLb())<<" "<<exp(ToulBar2::markov_log)<<endl;

  TLogProb newLogL=wcsp->Cost2LogProb(wcsp->getLb() + wcsp->getNegativeLb()) ; // New Lower Bound on Z init to c_0
  TLogProb q_i;
  TLogProb q_j;
  unsigned int n_conv=1;
  
  for (unsigned int i = 0; i < wcsp->numberOfVariables() ; i++) { // Loop on the unassigned variables
	if (wcsp->unassigned(i)) {
        EnumeratedVariable *y = (EnumeratedVariable *)((WCSP *) wcsp)->getVar(i);  // get the variable i
        y->initMFdistrib();
      }
  }
   
  for (unsigned int k=0; k<n_conv;k++){
    for (unsigned int i = 0; i < wcsp->numberOfVariables() ; i++) { // Loop on the unassigned variables
      if (wcsp->unassigned(i)) {
        EnumeratedVariable *y = (EnumeratedVariable *)((WCSP *) wcsp)->getVar(i);  // get the variable i
        y->UpdateMFdistrib();
      }
    }
  }

  //~ for (unsigned int i = 0; i < wcsp->numberOfVariables() ; i++) { // Loop on the unassigned variables
	//~ if (wcsp->unassigned(i)) {
        //~ EnumeratedVariable *y = (EnumeratedVariable *)((WCSP *) wcsp)->getVar(i);  // get the variable i
        //~ y->showMFdistrib();
      //~ }
  //~ }

   for (unsigned int i = 0; i < wcsp->numberOfVariables() ; i++) { // Loop on variable i
	 if (wcsp->unassigned(i)) {
		EnumeratedVariable *x = (EnumeratedVariable *)((WCSP *) wcsp)->getVar(i); // variable i = x
		for (EnumeratedVariable::iterator itx = x->begin(); itx != x->end(); ++itx) { // Loop on x values
			q_i=x->MFdistrib[x->toCurrentIndex(*itx)] ; // get q_i(x_i)
			//cout << x->getName() << " " <<x->toCurrentIndex(*itx) << " "<<q_i<<endl;
			if (q_i == 0 ) continue; // if q_i(x_i) = 0 then the whole contribution is equal to 0
			else newLogL += - q_i*Log(q_i); // Subtract the entropy q_i(x_i) * log(q_i(x_i))
			TLogProb newContribution = wcsp->Cost2LogProb(x->getCost(*itx)) ;
			for (unsigned int j = i + 1; j < wcsp->numberOfVariables(); j++) { // Loop on variable i
				if (wcsp->unassigned(j)) {
					EnumeratedVariable *y = (EnumeratedVariable *)((WCSP *) wcsp)->getVar(j); // variable j = y;
					BinaryConstraint *bctr = y->getConstr(x);// get constr that link i to j
					if (bctr != NULL) {
						for (EnumeratedVariable::iterator ity = y->begin(); ity != y->end(); ++ity) { // Loop on y values
							if(wcsp->Cost2Prob(bctr->getCost(x, y, *itx, *ity))==0) continue; // if p=0 then we consider log(p)=0 in order to not take in account impossible proba
							q_j=y->MFdistrib[y->toCurrentIndex(*ity)] ; //get q_j(x_j)
							newContribution += q_j * wcsp->Cost2LogProb(bctr->getCost(x, y, *itx, *ity)); // Contribution q_j(x_j) * c_ij(x_i,x_j)
							//~ cout << x->getName() << " "<< y->getName()<< " ";
							//~ cout << x->toCurrentIndex(*itx) << " "<< y->toCurrentIndex(*ity)<<" ";
							//~ cout << q_i << " "<<q_j<<" ";
						}
					}
				}
			}
			newLogL += q_i * newContribution;
		}
	 } 
   }
   
   return newLogL;
}

///////////////: Gumbel Perturbation (Non concluant) ////////////
TLogProb Solver::GumofThrone()
{
  TLogProb LogZhat = 0;
  for (auto iter : ToulBar2::trieZ->get_sols()) {
    LogZhat += iter;
  }
  LogZhat = LogZhat / ToulBar2::trieZ->get_sols().size();
  return LogZhat;
}


//////////////////// PROD SUM Diffusion ////////////
void Solver::ProdSumDiffusion()
{
	int times = 0;
	bool change = true;
	//cout << "ProdSumDiffusion"<<endl;

	while (change && (times < ToulBar2::prodsumDiffusion)) {
		change = false;
		int nchanged = 0;
		for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) //loop on the variable
			if (wcsp->unassigned(i)) {
				//~ cout<< "Variable "<<i<<endl;
				EnumeratedVariable *var = (EnumeratedVariable *)((WCSP *) wcsp)->getVar(i);  // get the variable i
				if (var->Marginalisation()) {
					change = true;
					nchanged++;
				}
			}
		times++;
	}
}

bool EnumeratedVariable::Marginalisation()
{
	bool change = false;
	EnumeratedVariable *x;
	Constraint *ctr = NULL;

	for (ConstraintList::iterator itc = getConstrs()->begin(); itc != getConstrs()->end(); ++itc) { //constraints linking the "this" variable
		ctr = (*itc).constr;
		if (ctr->arity() == 2 && !ctr->isSep() && ctr != NULL) {
			BinaryConstraint *bctr = (BinaryConstraint *) ctr;
			x = (EnumeratedVariable *) bctr->getVarDiffFrom((Variable *) this); // get variable different from "this"
			//~ cout<<"link to variable "<< x->getName()<<endl;
			for (iterator it = this->begin(); it != this->end(); ++it) {// Loop on "this" values
				Cost csum = MAX_COST;
				//~ cout<<"Binary Cost before : ";
				for (iterator itx = x->begin(); itx != x->end(); ++itx) { // Loop on x values
					Cost cbin = bctr->getCost(this, x, *it, *itx); // get the binary cost linking "this,it" and "x,itx"
					//~ cout<<wcsp->Cost2Prob(cbin)<<' ';
					csum = wcsp->LogSumExp(csum, cbin);
				}
				//~ cout<<endl;
				Double extc = to_double(csum - getCost(*it)) / 2. ;
				Cost cost_extc = (Long) extc;
				//~ cout<<"Sum : "<<wcsp->Cost2Prob(csum - getCost(*it))<<" Delta : "<<wcsp->Cost2Prob(extc)<<endl;
				//~ cout<<"After : ";
				for (iterator itx = x->begin(); itx != x->end(); ++itx) {
					//~ if(LUBTEST(bctr->getCost(this, x, *it, *itx), cost_extc)){
					//~ cout<<"Negative binary cost appear !!!!!!"<<endl;
					//~ }
					//~ cout<<wcsp->Cost2Prob(bctr->getCost(this, x, *it, *itx))<<' ';
					bctr->addcost(this, x, *it, *itx, -cost_extc);
					//~ cout<<wcsp->Cost2Prob(bctr->getCost(this, x, *it, *itx))<<endl;
					//~ }else{
					//~ cout<<"Negative Cost :"<<bctr->getCost(this, x, *it, *itx) + cost_extc <<endl;
					//~ tobechanged=false;
					//~ break;
					//~ }
					//~ cout<<wcsp->Cost2Prob(bctr->getCost(this,x,*it,*itx))<<' ';
				}
				//~ cout<<endl;
				//~ if((cost_extc + getCost(*it)) <= MIN_COST){
				//~ cout<<"Negative unary cost appear !!!!!!"<<endl;
				//~ cout<<"Sum : "<<wcsp->Cost2Prob(cost_extc + getCost(*it))<<' ';
				//~ }
				project(*it, cost_extc);
				//~ cout<<" After :"<<wcsp->Cost2Prob(getCost(*it))<<endl;
				change = true;
			}
		}
	}
	//UnaryNormalization();
	return change;
}

//////NORMALISATION ON BINARY AND UNARY COST EPTs (TRY Non concluant)
void Solver::PropagateNoc()
{
	//cout<<"Propagate NOC "<<endl;
	bool change = true;
	while (change) {
		change = false;
		for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) //loop on the variable
			if (wcsp->unassigned(i)) {
				//~ cout<< "Variable "<<i<<endl;
				EnumeratedVariable *var = (EnumeratedVariable *)((WCSP *) wcsp)->getVar(i);  // get the variable i
				if (var->Normalization()) {
					change = true;
				}
			}
	}
}

void EnumeratedVariable::UnaryNormalization()
{
	Cost csum = MAX_COST;
	for (iterator it = begin(); it != end(); ++it) {
		//cout<<"Before Value "<<toIndex(*it)<<" : "<<wcsp->Cost2Prob(getCost(*it))<<endl;
		csum = wcsp->LogSumExp(csum, getCost(*it));
	}
	if (wcsp->Cost2Prob(csum) != 1.) {
		for (iterator it = begin(); it != end(); ++it) {
			project(*it, -csum);
			//cout<<"After Value "<<toIndex(*it)<<" : "<<wcsp->Cost2Prob(getCost(*it))<<endl;
		}
		if (csum < MIN_COST) {
			//cout<<"OldNegLB : "<<wcsp->Cost2Prob(wcsp->getNegativeLb())<<endl;
			//cout<<"DecreaseLB : "<<wcsp->Cost2Prob(csum)<<endl;
			wcsp->decreaseLb(csum);
			//cout<<"NewNegLB : "<<wcsp->Cost2Prob(wcsp->getNegativeLb())<<endl;
		} else wcsp->increaseLb(csum);
	}
}

bool EnumeratedVariable::Normalization()
{
	//~ vector<Cost> vsum;
	bool change = false;
	Constraint *ctr = NULL;
	EnumeratedVariable *x;

	for (ConstraintList::iterator itc = getConstrs()->begin(); itc != getConstrs()->end(); ++itc) { //constraints linking the "this" variable
		ctr = (*itc).constr;
		if (ctr->arity() == 2 && !ctr->isSep() && ctr != NULL) {
			BinaryConstraint *bctr = (BinaryConstraint *) ctr;
			x = (EnumeratedVariable *) bctr->getVarDiffFrom((Variable *) this); // get variable different from "this"
			//~ cout<<"link to variable "<< x->getName()<<endl;
			for (iterator it = this->begin(); it != this->end(); ++it) { // Loop on "this" values
				Cost csum = MAX_COST;
				for (iterator itx = x->begin(); itx != x->end(); ++itx) { // Loop on x values
					csum = wcsp->LogSumExp(csum, bctr->getCost(this, x, *it, *itx));
				}
				if (csum < 0) {
					change = true;
					for (iterator itx = x->begin(); itx != x->end(); ++itx) {
						bctr->addcost(this, x, *it, *itx, -csum);
					}
					cout << endl;
					project(*it, csum);
				} else {change = false;}
			}
		}
	}
	UnaryNormalization();
	return change;
}

pair<TLogProb,TLogProb> Solver::GetOpen_LB_UB (OpenList &open){

  OpenList copy_open=open;
  TLogProb OpenLB = -numeric_limits<TLogProb>::infinity(); //init
  TLogProb OpenUB = -numeric_limits<TLogProb>::infinity(); //init
  
  while(!copy_open.empty()){
    OpenNode nd = copy_open.top();
    copy_open.pop();
    OpenLB = wcsp->LogSumExp(OpenLB,nd.getZlb());
    OpenUB = wcsp->LogSumExp(OpenUB,nd.getZub());
  }
  pair<TLogProb,TLogProb> LB_UB(OpenLB,OpenUB);
  
  return LB_UB;
}

void Solver::showZGap()
{
	if (Store::getDepth() == initialDepth) {
    //cout << ToulBar2::GlobalLogLbZ +ToulBar2::markov_log << " <= Log(Z) <= ";
    //cout << ToulBar2::GlobalLogUbZ +ToulBar2::markov_log << endl;
    cout << wcsp->LogSumExp(ToulBar2::logZ,ToulBar2::GlobalLogLbZ) + ToulBar2::markov_log<< " <= Log(Z) <= ";
    cout << wcsp->LogSumExp(wcsp->LogSumExp(ToulBar2::logZ ,ToulBar2::GlobalLogUbZ),ToulBar2::logU) + ToulBar2::markov_log << endl;
    //cout<<"Z chapo : " <<ToulBar2::logZ +ToulBar2::markov_log<<endl;
	}
}

void TimeZOut(){	
WeightedCSP *wcsp=new WCSP(0,NULL);
cout << wcsp->LogSumExp(ToulBar2::logZ,ToulBar2::GlobalLogLbZ) + ToulBar2::markov_log<< " <= Log(Z) <= ";
cout << wcsp->LogSumExp(wcsp->LogSumExp(ToulBar2::logZ ,ToulBar2::GlobalLogUbZ),ToulBar2::logU) + ToulBar2::markov_log << endl;
cout <<"Epsilon : "<<Exp(wcsp->LogSumExp(wcsp->LogSumExp(ToulBar2::logZ ,ToulBar2::GlobalLogUbZ),ToulBar2::logU) - wcsp->LogSumExp(ToulBar2::logZ,ToulBar2::GlobalLogLbZ))-1<<endl;
cout <<"Time   : " << cpuTime() - ToulBar2::startCpuTime << " seconds" << endl;

exit(0);
}

void Solver::hybridCounting(TLogProb Zlb, TLogProb Zub)
{
  // Need to adapt this solver to Z mode
  assert(Zlb < Zub);
  if (ToulBar2::hbfs) {
	ToulBar2::timeOut=TimeZOut;
    CPStore *cp_ = NULL;
    OpenList *open_ = NULL;
	
    // normal BFS without BTD, i.e., hybridSolve is not reentrant
    if (cp != NULL) delete cp;
    cp = new CPStore();
    cp_ = cp;
    if (open != NULL) delete open;
    open = new OpenList();
    open_ = open;
		
    cp_->store();
    if (open_->size() == 0 ) { // start a new list of open nodes if needed
      nbHybridNew++;
      // reinitialize current open list and insert empty node
      *open_ = OpenList();
      addOpenNode(*cp_, *open_, wcsp->getLb(),Zlb,Zub);
    } else nbHybridContinue++;
    
    nbHybrid++; // do not count empty root cluster
        
    while (!open_->empty() ){
      hbfsLimit = ((ToulBar2::hbfs > 0) ? (nbBacktracks + ToulBar2::hbfs) : LONGLONG_MAX);
      int storedepthBFS = Store::getDepth();
      try {
	Store::store();
        //cout<<"open size before pop : "<<open_->size()<<endl;
	OpenNode nd = open_->top();
	open_->pop();
        restore(*cp_, nd);
	Cost bestlb = MAX(nd.getCost(), wcsp->getLb());
	//cout <<"Sub LB : "<< nd.getZlb()<< " to : "<< ToulBar2::GlobalLogLbZ<<" "<<(ToulBar2::GlobalLogLbZ>nd.getZlb())<<endl;
	//cout <<"Sub UB : "<< nd.getZub()<< " to : "<< ToulBar2::GlobalLogUbZ<<" "<<(ToulBar2::GlobalLogUbZ>nd.getZub())<<endl;
	//cout <<endl; 
	//ToulBar2::GlobalLogLbZ = wcsp->LogSubExp(ToulBar2::GlobalLogLbZ, (nd.getZlb()) ); // remove pop open node from global Z LB
	//ToulBar2::GlobalLogUbZ = wcsp->LogSubExp(ToulBar2::GlobalLogUbZ, (nd.getZub()) ); // remove pop open node from global Z UB
	//cout<<endl;
	recursiveSolve(bestlb);
      } catch (Contradiction) {
	wcsp->whenContradiction();
      }
        
      Store::restore(storedepthBFS);
      cp_->store();

      if (cp_->size() >= static_cast<std::size_t>(ToulBar2::hbfsCPLimit) || open_->size() >= static_cast<std::size_t>(ToulBar2::hbfsOpenNodeLimit)) {
	ToulBar2::hbfs = 0;
	ToulBar2::hbfsGlobalLimit = 0;
	hbfsLimit = LONGLONG_MAX;
      }
      pair<TLogProb,TLogProb> LB_UB = GetOpen_LB_UB(*open);
      ToulBar2::GlobalLogLbZ = get<0>(LB_UB);
      ToulBar2::GlobalLogUbZ = get<1>(LB_UB);
      if (ToulBar2::verbose >=1) showZGap();
          
      // if LogUb - LogLb < log(1+epsilon) then finish
      //cout <<"Epsilon : "<<exp(wcsp->LogSumExp(wcsp->LogSumExp(ToulBar2::logZ ,ToulBar2::GlobalLogUbZ),ToulBar2::logU) - wcsp->LogSumExp(ToulBar2::logZ,ToulBar2::GlobalLogLbZ))-1<<endl;
      if (wcsp->LogSumExp(ToulBar2::logZ ,ToulBar2::GlobalLogUbZ) - wcsp->LogSumExp(ToulBar2::logZ,ToulBar2::GlobalLogLbZ)<= Log1p(ToulBar2::sigma)) break;
      if (ToulBar2::hbfs && nbRecomputationNodes > 0) { // wait until a nonempty open node is restored (at least after first global solution is found)
	assert(nbNodes > 0);
	if (nbRecomputationNodes > nbNodes / ToulBar2::hbfsBeta && ToulBar2::hbfs <= ToulBar2::hbfsGlobalLimit) ToulBar2::hbfs *= 2;
	else if (nbRecomputationNodes < nbNodes / ToulBar2::hbfsAlpha && ToulBar2::hbfs >= 2) ToulBar2::hbfs /= 2;
	if (ToulBar2::debug >= 2) cout << "HBFS backtrack limit: " << ToulBar2::hbfs << endl;
      }
    }
    //cout<<"open size after finished : "<<open_->size()<<" "<<open_->empty()<<endl;
  } else {
    hbfsLimit = LONGLONG_MAX;
    recursiveSolve();
  }
}


TLogProb Solver::binaryChoicePointBTDZ(Cluster *cluster, int varIndex, Value value)
{
  TLogProb logZcluster = -numeric_limits<TLogProb>::infinity();
	if (ToulBar2::interrupted) throw TimeOut();
  TreeDecomposition *td = wcsp->getTreeDec();
	assert(wcsp->unassigned(varIndex));
	assert(wcsp->canbe(varIndex, value));
	bool dichotomic = (ToulBar2::dichotomicBranching && ToulBar2::dichotomicBranchingSize < wcsp->getDomainSize(varIndex));
	Value middle = value;
	bool increasing = true;
	if (dichotomic) {
		middle = (wcsp->getInf(varIndex) + wcsp->getSup(varIndex)) / 2;
		if (value <= middle) increasing = true;
		else increasing = false;
	}
	try {
    
		Store::store();
		assert(wcsp->getTreeDec()->getCurrentCluster() == cluster);
		//assert(wcsp->getLb() == cluster->getLbRec());
		//wcsp->setUb(cub);
		lastConflictVar = varIndex;
		if (dichotomic) {
			if (increasing) decrease(varIndex, middle);
			else increase(varIndex, middle + 1);
		} 
    else assign(varIndex, value);
		lastConflictVar = -1;
		TLogProb logres = BTD_sharpZ(cluster);
		logZcluster = wcsp->LogSumExp(logZcluster,logres);
	} catch (Contradiction) {
		wcsp->whenContradiction();
	}
	Store::restore();
	nbBacktracks++;
	if (ToulBar2::restart > 0 && nbBacktracks > nbBacktracksLimit) throw NbBacktracksOut();
  cluster->nbBacktracks++;
	try {
		Store::store();
		assert(wcsp->getTreeDec()->getCurrentCluster() == cluster);
		if (dichotomic) {
			if (increasing) increase(varIndex, middle + 1, cluster->nbBacktracks >= cluster->hbfsLimit || nbBacktracks >= cluster->hbfsGlobalLimit);
			else decrease(varIndex, middle, cluster->nbBacktracks >= cluster->hbfsLimit || nbBacktracks >= cluster->hbfsGlobalLimit);
		} else remove(varIndex, value, cluster->nbBacktracks >= cluster->hbfsLimit || nbBacktracks >= cluster->hbfsGlobalLimit);
    
    if (!ToulBar2::hbfs && cluster == td->getRoot() && initialDepth + 1 == Store::getDepth()) {initialDepth++;};

		TLogProb logres = BTD_sharpZ(cluster);
		logZcluster = wcsp->LogSumExp(logZcluster,logres); // Normaly sum but here we are in the logdomain so logsomexp
	} catch (Contradiction) {
		wcsp->whenContradiction();
	}
	Store::restore();
	return logZcluster;
}

TLogProb Solver::BTD_sharpZ(Cluster *cluster)
{

	TreeDecomposition *td = wcsp->getTreeDec();
	TLogProb logZcluster = -numeric_limits<TLogProb>::infinity();
  TLogProb logres=-numeric_limits<TLogProb>::infinity();
	TCtrs totalList;
	//if (ToulBar2::verbose >= 1) 
  cout << "[" << Store::getDepth() << "] recursive solve     cluster: " << cluster->getId() << " **************************************************************" << endl;

	int varIndex = -1;
	if (ToulBar2::Static_variable_ordering) varIndex = getNextUnassignedVar(cluster);
	else if (ToulBar2::weightedDegree && ToulBar2::lastConflict) varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxWeightedDegreeLastConflictRandomized(cluster) : getVarMinDomainDivMaxWeightedDegreeLastConflict(cluster));
	else if (ToulBar2::lastConflict) varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxDegreeLastConflictRandomized(cluster) : getVarMinDomainDivMaxDegreeLastConflict(cluster));
	else if (ToulBar2::weightedDegree) varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxWeightedDegreeRandomized(cluster) : getVarMinDomainDivMaxWeightedDegree(cluster));
	else varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxDegreeRandomized(cluster) : getVarMinDomainDivMaxDegree(cluster));

	if (varIndex < 0) {
		// Current cluster is completely assigned
    logZcluster =  wcsp->Cost2LogProb(cluster->getLb());
    //~ logZcluster =  wcsp->Cost2LogProb(cluster->getLb()+ wcsp->getNegativeLb());
    //cout<<"Add "<<logZcluster<<" on cluster "<<cluster->getId()<<endl;
    //cluster->add2logZ(logZcluster); // Add solution to the logZcluster
		//if (ToulBar2::verbose >= 1) 
    cout << "[" << Store::getDepth() << "] C" << cluster->getId() << " logZ= " << logZcluster << endl;
    
		for (TClusters::iterator iter = cluster->beginSortedEdges(); iter != cluster->endSortedEdges(); ++iter) {
			// Solves each cluster son
			Cluster *c = *iter;
      td->setCurrentCluster(c);
      //~ Cost lbSon = MIN_COST;
      //~ bool good =false; 
      //~ if (!c->isActive()) {
				//~ c->reactivate();
				//~ c->nogoodGet(lbSon, ubSon, &c->open);
				//~ good = true;
			//~ } 
      //~ else {
				//~ lbSon = c->getLbRec();
			//~ }
      if (!c->islogZset()){ // if already computed don't go there BOOL REMINDER
        try {
            Store::store();
            wcsp->propagate();
            logres = BTD_sharpZ(c);
            cout<<"setlogZ : "<<logres<<" ending of  cluster "<<c->getId()<<endl;
            c->setlogZ(logres);
        } catch (Contradiction) {
            cout<<"Catch Contradiction"<<endl;
            wcsp->whenContradiction();
            logres = -numeric_limits<TLogProb>::infinity();
            c->setlogZ(logres);
        }
        Store::restore();
      }
      else{
        logres=c->getlogZ();
        cout<<"getlogZ : "<<logres<<" of cluster  "<<c->getId()<<endl;
      }
      logZcluster += logres; // normally we need to multiply but here we are in the log domain e.g Z = Z*res <=> log(Z) = log(Z) + log(res)
    }
		return logZcluster;
	} else {
		// Enumerates cluster proper variables
		if (wcsp->enumerated(varIndex)) {
			assert(wcsp->canbe(varIndex, wcsp->getSupport(varIndex)));
			// Reuse last solution found if available
			Value bestval = wcsp->getBestValue(varIndex);
			logZcluster = binaryChoicePointBTDZ(cluster, varIndex, (wcsp->canbe(varIndex, bestval)) ? bestval : wcsp->getSupport(varIndex));
		} else {
			logZcluster = binaryChoicePointBTDZ(cluster, varIndex, wcsp->getInf(varIndex));
		}
		//if (ToulBar2::verbose >= 1) 
    cout << "[" << Store::getDepth() << "] C" << cluster->getId() << " return " << logZcluster << endl;
		return logZcluster;
	}
}
