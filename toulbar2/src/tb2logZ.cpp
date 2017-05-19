#include "tb2solver.hpp"
#include "tb2trienum.hpp"
#include "tb2domain.hpp"
#include "tb2binconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2clusters.hpp"

void Solver::enforceZUb(Cluster *cluster)
{
    TLogProb newlogU;
    if (cluster){
      //~ cout<< ToulBar2::logU<<endl;
      //~ newlogU = wcsp->LogSumExp(ToulBar2::logU,Zub(cluster)); // DOES NOT GIVE EPSILON APPROXIMATION
      //~ newlogU = wcsp->LogSumExp(cluster->getlogU(),Zub(cluster)); //Compute Cluster upper bound
      //~ newlogU = wcsp->LogSumExp(ToulBar2::logU,Zub());
      //if (ToulBar2::verbose >= 0)  cout << "ZCUT Using bound " << ToulBar2::isZUB << " U : " << newlogU << " log Z "<< logZcluster.first << " Log(eps x Z) : " << logZcluster.first + ToulBar2::logepsilon << " " << Store::getDepth() << endl;
      //~ if (newlogU < ToulBar2::logepsilon + logZcluster.first) {
          //~ if (ToulBar2::verbose >= 1)  cout << "ZCUT Using bound " << ToulBar2::isZUB << " U : " << newlogU << " log Z "<< logZcluster.first << " Log(eps x Z) : " << logZcluster.first + ToulBar2::logepsilon << " " << Store::getDepth() << endl;
          //~ cluster->setlogU(newlogU);
          //~ ToulBar2::logU = wcsp->LogSumExp(ToulBar2::logU,newlogU); //FAIL
          //~ ToulBar2::logU = newlogU; //FAIL
          //~ THROWCONTRADICTION;
      //~ }
    }
    else{
    newlogU = wcsp->LogSumExp(ToulBar2::logU,Zub());
      if (newlogU < ToulBar2::logepsilon + ToulBar2::logZ) {
          if (ToulBar2::verbose >= 1)  cout << "ZCUT Using bound " << ToulBar2::isZUB << " U : " << newlogU << " log Z "<< ToulBar2::logZ << " Log(eps x Z) : " << ToulBar2::logZ + ToulBar2::logepsilon << " " << Store::getDepth() << endl;
          ToulBar2::logU = newlogU;
          ToulBar2::GlobalLogUbZ=wcsp->LogSumExp(ToulBar2::GlobalLogUbZ,ToulBar2::logU);
          THROWCONTRADICTION;
      }
    }
}

TLogProb Solver::Zub(Cluster *cluster)  // Calculate an upper bound on Z (on a cluster if specified)
{
  Cost newCost = wcsp->getLb() + wcsp->getNegativeLb();
  switch (ToulBar2::isZUB) {
    //Upper bound on Z edition 0
  case 0 :
    if (cluster){
      newCost = cluster->getLb() + wcsp->LogProb2Cost(cluster->getNbVars() * Log(cluster->MaxDomainSize()));
      return wcsp->Cost2LogProb(newCost);
    }
    else{
      newCost += wcsp->LogProb2Cost(unassignedVars->getSize() * Log(wcsp->getMaxDomainSize()));
      return wcsp->Cost2LogProb(newCost);
    }
    break;
    //Upper bound on Z edition 1
  case 1 :
  if (cluster){
    newCost=cluster->getLb();
    for (TVars::iterator iter_variable = cluster->beginVars(); iter_variable != cluster->endVars(); ++iter_variable) { // Loop on the unassigned variables
      EnumeratedVariable *var = (EnumeratedVariable *)((WCSP *) wcsp)->getVar(*iter_variable);
      Cost SumUnaryCost = MAX_COST;
      if (var->unassigned()){
        if (wcsp->enumerated(*iter_variable)) {
          for (EnumeratedVariable::iterator iter_value = var->begin(); iter_value != var->end(); ++iter_value) { // loop over the domain of the variable
            //cout<<"Sum " SumUnaryCost<<
            SumUnaryCost = wcsp->LogSumExp(SumUnaryCost, var->getCost(*iter_value)); // Sum of the exponential of Unary cost over the domain.
          }
          newCost += SumUnaryCost; //Sum the older cost with the new log(sum(unarycost))
        } else {
          newCost += wcsp->LogProb2Cost(Log(wcsp->getDomainSize(*iter_variable)));
        }
      }
    }
    return wcsp->Cost2LogProb(newCost);
  }
  else{
    for (BTList<Value>::iterator iter_variable = unassignedVars->begin(); iter_variable != unassignedVars->end(); ++iter_variable) { // Loop on the unassigned variables
      EnumeratedVariable *var = (EnumeratedVariable *)((WCSP *) wcsp)->getVar(*iter_variable);
      Cost SumUnaryCost = MAX_COST;
      if (wcsp->enumerated(*iter_variable)) {
        for (EnumeratedVariable::iterator iter_value = var->begin(); iter_value != var->end(); ++iter_value) { // loop over the domain of the variable
          SumUnaryCost = wcsp->LogSumExp(SumUnaryCost, var->getCost(*iter_value)); // Sum of the exponential of Unary cost over the domain.
        }
        newCost += SumUnaryCost; //Sum the older cost with the new log(sum(unarycost))
      } else {
        newCost += wcsp->LogProb2Cost(Log(wcsp->getDomainSize(*iter_variable)));
      }
    }
    return wcsp->Cost2LogProb(newCost);
  }
  break;

  case 2 : //Upper bound on Z edition 2
  if (cluster){
    cerr<<"We do not implemented UbT with BTD, yet"<<endl;
    exit(0);
  }
  else{
    return wcsp->Cost2LogProb(wcsp->spanningTreeZ(newCost));
  }
  break;
  }
  return -numeric_limits<TLogProb>::infinity();
}


TLogProb Solver::MeanFieldZ(Cluster *cluster) // Compute a lower bound on logZ using MF approximation
{
	//~ cout<<"Cout "<<wcsp->getLb()<< " "<<wcsp->getNegativeLb()<<" "<<ToulBar2::markov_log<<endl;
	//~ cout<<"Log "<<wcsp->Cost2LogProb(wcsp->getLb())<< " "<<wcsp->Cost2LogProb(wcsp->getNegativeLb())<<" "<<ToulBar2::markov_log<<endl;
	//~ cout<<"Prob "<<wcsp->Cost2Prob(wcsp->getLb())<< " "<<wcsp->Cost2Prob(wcsp->getNegativeLb())<<" "<<exp(ToulBar2::markov_log)<<endl;

  TLogProb newLogL=wcsp->Cost2LogProb(wcsp->getLb() + wcsp->getNegativeLb()) ; // New Lower Bound on Z init to c_0
  TLogProb q_i;
  TLogProb q_j;
  unsigned int n_conv=1;
 
  if (cluster){
    newLogL = wcsp->Cost2LogProb(cluster->getLb());
    for (TVars::iterator i = cluster->beginVars(); i != cluster->endVars() ; i++) { // Loop on the unassigned variables
      EnumeratedVariable *y = (EnumeratedVariable *)((WCSP *) wcsp)->getVar(*i);  // get the variable i
      if(y->unassigned()){
          y->initMFdistrib();
      }
    }
   
    for (unsigned int k=0; k<n_conv;k++){
      for (TVars::iterator i = cluster->beginVars(); i != cluster->endVars() ; i++) { // Loop on the unassigned variables
        EnumeratedVariable *y = (EnumeratedVariable *)((WCSP *) wcsp)->getVar(*i);  // get the variable i
        if(y->unassigned()){
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

    for (TVars::iterator i = cluster->beginVars(); i != cluster->endVars() ; i++) { // Loop on variable i
      EnumeratedVariable *x = (EnumeratedVariable *)((WCSP *) wcsp)->getVar(*i);  // get the variable i
      if (x->unassigned()) {
        for (EnumeratedVariable::iterator itx = x->begin(); itx != x->end(); ++itx) { // Loop on x values
          q_i=x->MFdistrib[x->toCurrentIndex(*itx)] ; // get q_i(x_i)
          //~ cout << x->getName() << " " <<x->toCurrentIndex(*itx) << " "<<q_i<<endl;
          if (q_i == 0 ) continue; // if q_i(x_i) = 0 then the whole contribution is equal to 0
          else newLogL += - q_i*Log(q_i); // Subtract the entropy q_i(x_i) * log(q_i(x_i))
          TLogProb newContribution = wcsp->Cost2LogProb(x->getCost(*itx)) ;
          for (TVars::iterator j = cluster->beginVars(); j != cluster->endVars() ; j++) { // Loop on variable i
            EnumeratedVariable *y = (EnumeratedVariable *)((WCSP *) wcsp)->getVar(*j); // variable j = y;
              if (y->unassigned() && stoi(y->getName()) > stoi(x->getName())) {
                BinaryConstraint *bctr = y->getConstr(x);// get constr that link i to j
                //cout << x->getName() << " "<< y->getName()<< " "<<y->getConstr(x)<<" "<<(bctr != NULL)<<endl;
                  if (bctr != NULL) {
                    for (EnumeratedVariable::iterator ity = y->begin(); ity != y->end(); ++ity) { // Loop on y values
                      if(wcsp->Cost2Prob(bctr->getCost(x, y, *itx, *ity))==0) continue; // if p=0 then we consider log(p)=0 in order to not take in account impossible proba
                        q_j=y->MFdistrib[y->toCurrentIndex(*ity)] ; //get q_j(x_j)
                        newContribution += q_j * wcsp->Cost2LogProb(bctr->getCost(x, y, *itx, *ity)); // Contribution q_j(x_j) * c_ij(x_i,x_j)
                        //~ cout << x->getName() << " "<< y->getName()<< " ";
                        //~ cout << x->toCurrentIndex(*itx) << " "<< y->toCurrentIndex(*ity)<<endl;
                        //~ cout << q_i << " "<<q_j<<" ";
                    }
                  }
              }
          }
          newLogL += q_i * newContribution;
        }
      } 
    }
  }
  else{ 
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


TPairLz Separator::getLz()
{
	int i = 0;
	Cost deltares = MIN_COST; // cost that run away from the cluster
	if (ToulBar2::verbose >= 1) cout << "( ";
	TVars::iterator it = vars.begin();
	while (it != vars.end()) {
		assert(cluster->getWCSP()->assigned(*it));
		Value val = cluster->getWCSP()->getValue(*it);
		if (ToulBar2::verbose >= 1) cout << "(" << *it << "," << val << ") ";
		t[i] = val + CHAR_FIRST;	 // build the tuple
		deltares -= delta[i][val];      // delta structure
		++it;
		i++;
	}
	TLZGoods::iterator itlz = lzgoods.find(t);
  if (ToulBar2::verbose >=1){
  cout << " < C" << cluster->getId() << " , ";
  Cout << t;
  cout << ",delta: "<<deltares <<">";
  }
	if (itlz != lzgoods.end()) {
		TLogProb p = (itlz->second).first;
		if (ToulBar2::verbose >= 1)	cout << ") Use #good with LogZ = " << p - wcsp->Cost2LogProb(deltares)<< " LogZ on cluster " << cluster->getId() << endl;
		return make_pair((p - wcsp->Cost2LogProb(deltares)),true); 
	} else {
		if (ToulBar2::verbose >= 1)	cout << ") NOT FOUND for cluster " <<  cluster->getId() << endl;
    return make_pair(-numeric_limits<TLogProb>::infinity(),false);
	}
}


void Separator::setLz(TLogProb logz)
{
	int i = 0;
	WCSP *wcsp = cluster->getWCSP();
	Cost deltares = MIN_COST;
	if (ToulBar2::verbose >= 1)	cout << "( ";
	TVars::iterator it = vars.begin();
	while (it != vars.end()) {
		assert(wcsp->assigned(*it));
		Value val = wcsp->getValue(*it);
		if (ToulBar2::verbose >= 1)	cout << "(" << *it << "," << val << ") ";
		t[i] = val + CHAR_FIRST;
		deltares += delta[i][val];
		++it;
		i++;
	}
  assert(deltares >= MIN_COST);
  if (ToulBar2::verbose >= 1){
  cout << " < C" << cluster->getId() << " , ";
  Cout << t;
  cout << ",delta: "<<deltares;
  }
  TLZGoods::iterator itlz = lzgoods.find(t);
  if (itlz == lzgoods.end()) {
    if (ToulBar2::verbose >= 1) cout << ", set lzgoods from " << lzgoods[t].first<<" to "<<logz- wcsp->Cost2LogProb(deltares)<<" >"<<endl;
    lzgoods[t] = TPairLz(logz - wcsp->Cost2LogProb(deltares),true);
  }
  else{
    cout<<", already existing !"<<endl;
    exit(0); 
  }
}

void Separator::add2Lz(TLogProb logz) // Add solution to logz of the separator.
{
	int i = 0;
	WCSP *wcsp = cluster->getWCSP();
	Cost deltares = MIN_COST;
	if (ToulBar2::verbose >= 1)	cout << "( ";
	TVars::iterator it = vars.begin();
	while (it != vars.end()) {
		assert(wcsp->assigned(*it));
		Value val = wcsp->getValue(*it);
		if (ToulBar2::verbose >= 1)	cout << "(" << *it << "," << val << ") ";
		t[i] = val + CHAR_FIRST;
		deltares += delta[i][val];
		++it;
		i++;
	}
  TLZGoods::iterator itng = lzgoods.find(t);
  if (ToulBar2::verbose >=1){
  cout << " < C" << cluster->getId() << " , ";
  Cout << t;
  cout << ",delta: "<<deltares << ",";
  }
  if (itng == lzgoods.end()) { // There isn't any logZ count
    if (ToulBar2::verbose >=1) cout <<" No exist: ";
    lzgoods[t].first = logz - wcsp->Cost2LogProb(deltares); // init logZ count
    if (ToulBar2::verbose >=1) cout<<"lzgood = "<<lzgoods[t].first- wcsp->Cost2LogProb(deltares)<<endl;
  }
  else{ // add to current count
    if (ToulBar2::verbose >=1) cout <<" Already exist: ";
    itng->second.first= wcsp->LogSumExp(itng->second.first,logz - wcsp->Cost2LogProb(deltares));
    if (ToulBar2::verbose >=1)cout<<"add logz " << logz - wcsp->Cost2LogProb(deltares)<< " to " << lzgoods[t].first<<" >"<<endl;
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
  //~ if (ToulBar2::logepsilon>-numeric_limits<TLogProb>::infinity()) {
    //~ enforceZUb(cluster);
    //~ enforceZUb();
	//~ }
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
    //~ if (cluster->nbBacktracks >= cluster->hbfsLimit || nbBacktracks >= cluster->hbfsGlobalLimit) { // if backtrack fuel is empty then add to open list
        //~ addOpenNode(*(cluster->cp), *(cluster->open), bestlb, cluster->getCurrentDelta()); 
		//~ } else {
      TLogProb logres = BTD_sharpZ(cluster);
		//~ }
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
	if (ToulBar2::verbose >= 1) cout << "[" << Store::getDepth() << "] recursive solve     cluster: " << cluster->getId() << " **************************************************************" << endl;

	int varIndex = -1;
	if (ToulBar2::Static_variable_ordering) varIndex = getNextUnassignedVar(cluster);
	else if (ToulBar2::weightedDegree && ToulBar2::lastConflict) varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxWeightedDegreeLastConflictRandomized(cluster) : getVarMinDomainDivMaxWeightedDegreeLastConflict(cluster));
	else if (ToulBar2::lastConflict) varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxDegreeLastConflictRandomized(cluster) : getVarMinDomainDivMaxDegreeLastConflict(cluster));
	else if (ToulBar2::weightedDegree) varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxWeightedDegreeRandomized(cluster) : getVarMinDomainDivMaxWeightedDegree(cluster));
	else varIndex = ((ToulBar2::restart > 0) ? getVarMinDomainDivMaxDegreeRandomized(cluster) : getVarMinDomainDivMaxDegree(cluster));

	if (varIndex < 0) {
		// Current cluster is completely assigned
    cout<<"C"<<cluster->getId()<<" lb "<<cluster->getLb()<<" neglb "<<cluster->getNegativeLb()<<endl;
    logZcluster =  wcsp->Cost2LogProb(cluster->getLb()+cluster->getNegativeLb());
    
    if(ToulBar2::hbfs){
      if (ToulBar2::verbose >= 1) cout<<"Add "<<logZcluster<<" on cluster "<<cluster->getId()<<endl;
      cluster->add2logZ(logZcluster); // Add solution to the logZcluster
    }
    
		for (TClusters::iterator iter = cluster->beginSortedEdges(); iter != cluster->endSortedEdges(); ++iter) {
			// Solves each cluster son
			Cluster *c = *iter;
      td->setCurrentCluster(c);
      if (ToulBar2::verbose >= 1) cout << "[" << Store::getDepth() << "] C" << c->getId() << endl;
      //~ TLogProb ZubCluster = Zub(c);
      //~ TLogProb ZlbCluster = MeanFieldZ(c);
      //~ if (ToulBar2::verbose >= 1) cout<<ZlbCluster<<" <= Z(C"<<c->getId()<<") <= "<< ZubCluster <<endl;
      TPairLz ClusterCount = c->getlogZ();
      
      bool isCount = ClusterCount.second;
      if (isCount) { // if already counted 
        logres=ClusterCount.first;
        if (ToulBar2::verbose >= 1) cout<<"getlogZ : "<<logres<<" of cluster  "<<c->getId()<<endl;
      }
      else{ 
        try {
            Store::store();
            wcsp->propagate();
            logres = BTD_sharpZ(c);
            if(!ToulBar2::hbfs){
              if (ToulBar2::verbose >= 1) cout<<"setlogZ : "<<logres<<" ending of  cluster "<<c->getId()<<endl;
              c->setlogZ(logres);
            }
        } catch (Contradiction) {
            if (ToulBar2::verbose >= 1) cout<<"Catch Contradiction"<<endl;
            wcsp->whenContradiction();
            logres = -numeric_limits<TLogProb>::infinity();
            c->setlogZ(logres);
        }
        Store::restore();
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
		if (ToulBar2::verbose >= 1) cout << "[" << Store::getDepth() << "] C" << cluster->getId() << " return " << logZcluster << endl;
		return logZcluster;
	}
}

//~ void Solver::hybridCounting(Cluster *cluster, TLogProb Zlb, TLogProb Zub)
//~ {
  //~ assert(Zlb < Zub);
  //~ if (ToulBar2::verbose >= 1 && cluster) cout << "hybridSolve C" << cluster->getId() << " " << Zlb << " " << Zub << endl;
  //~ if (ToulBar2::hbfs) {
	//~ ToulBar2::timeOut=TimeZOut;
    //~ CPStore *cp_ = NULL;
    //~ OpenList *open_ = NULL;
    //~ Cost delta = MIN_COST;
    //~ // BFS with BTD on current cluster (can be root or not)
    //~ assert(cluster->cp);
    //~ cp_ = cluster->cp;
    //~ if (cluster == wcsp->getTreeDec()->getRoot()) {
      //~ if (!cluster->open) cluster->open = new OpenList();
      //~ cluster->setUb(cub); // global problem upper bound
    //~ } else {
      //~ delta = cluster->getCurrentDelta();
      //~ if (!cluster->open) {
        //~ //cluster->nogoodRec(clb, MAX_COST, &cluster->open); // create an initial empty open list
        //~ //cluster->setUb(MAX_COST); // no initial solution found for this cluster
      //~ }
    //~ }
    //~ assert(cluster->open);
    //~ open_ = cluster->open;

    //~ cp_->store();
    //~ if (open_->size() == 0) { // start a new list of open nodes if needed
      //~ if (cluster->getNbVars() > 0) nbHybridNew++;
      //~ // reinitialize current open list and insert empty node
      //~ *open_ = OpenList();
      //~ addOpenNode(*cp_, *open_, wcsp->getLb(),Zlb,Zub);
    //~ } else if (cluster->getNbVars() > 0) nbHybridContinue++;
    
    //~ if (cluster->getNbVars() > 0) nbHybrid++; // do not count empty root cluster
    //~ cluster->hbfsGlobalLimit = ((ToulBar2::hbfsGlobalLimit > 0) ? (nbBacktracks + ToulBar2::hbfsGlobalLimit) : LONGLONG_MAX);

    //~ while (!open_->empty() && nbBacktracks <= cluster->hbfsGlobalLimit){
      //~ cluster->hbfsLimit = ((ToulBar2::hbfs > 0) ? (cluster->nbBacktracks + ToulBar2::hbfs) : LONGLONG_MAX);
      //~ assert(wcsp->getTreeDec()->getCurrentCluster() == cluster);
      //~ assert(cluster->isActive());
      //~ int storedepthBFS = Store::getDepth();
      //~ try {
        //~ Store::store();
        //~ //cout<<"open size before pop : "<<open_->size()<<endl;
        //~ OpenNode nd = open_->top();
        //~ open_->pop();
        //~ if (wcsp->getTreeDec() && ToulBar2::verbose >= 1 ){
          //~ cout << "[C" << wcsp->getTreeDec()->getCurrentCluster()->getId() << "] ";
          //~ cout << "[ " << nd.getCost(delta) << ", " << cub <<  "] ( " << open_->size() << "+1 still open)" << endl;
        //~ }
        //~ restore(*cp_, nd);
        //~ Cost bestlb = MAX(nd.getCost(), wcsp->getLb());
        //~ //cout <<"Sub LB : "<< nd.getZlb()<< " to : "<< ToulBar2::GlobalLogLbZ<<" "<<(ToulBar2::GlobalLogLbZ>nd.getZlb())<<endl;
        //~ //cout <<"Sub UB : "<< nd.getZub()<< " to : "<< ToulBar2::GlobalLogUbZ<<" "<<(ToulBar2::GlobalLogUbZ>nd.getZub())<<endl;
        //~ //cout <<endl; 
        //~ //ToulBar2::GlobalLogLbZ = wcsp->LogSubExp(ToulBar2::GlobalLogLbZ, (nd.getZlb()) ); // remove pop open node from global Z LB
        //~ //ToulBar2::GlobalLogUbZ = wcsp->LogSubExp(ToulBar2::GlobalLogUbZ, (nd.getZub()) ); // remove pop open node from global Z UB
        //~ //cout<<endl;
        //~ pair<Cost, Cost> res = recursiveSolve(cluster, bestlb, cub);
        //~ open_->updateClosedNodesLb(res.first, delta);
        //~ open_->updateUb(res.second, delta);
				
      //~ } catch (Contradiction) {
        //~ wcsp->whenContradiction();
      //~ }
        
      //~ Store::restore(storedepthBFS);
      //~ cp_->store();

      //~ if (cp_->size() >= static_cast<std::size_t>(ToulBar2::hbfsCPLimit) || open_->size() >= static_cast<std::size_t>(ToulBar2::hbfsOpenNodeLimit)) {
        //~ ToulBar2::hbfs = 0;
        //~ ToulBar2::hbfsGlobalLimit = 0;
        //~ cluster->hbfsGlobalLimit = LONGLONG_MAX;
        //~ cluster->hbfsLimit = LONGLONG_MAX;
      //~ }
      //~ pair<TLogProb,TLogProb> LB_UB = GetOpen_LB_UB(*open);
      //~ ToulBar2::GlobalLogLbZ = get<0>(LB_UB);
      //~ ToulBar2::GlobalLogUbZ = get<1>(LB_UB);
      //~ if (ToulBar2::verbose >=1) showZGap();
          
      //~ // if LogUb - LogLb < log(1+epsilon) then finish
      //~ //cout <<"Epsilon : "<<exp(wcsp->LogSumExp(wcsp->LogSumExp(ToulBar2::logZ ,ToulBar2::GlobalLogUbZ),ToulBar2::logU) - wcsp->LogSumExp(ToulBar2::logZ,ToulBar2::GlobalLogLbZ))-1<<endl;
      //~ if (wcsp->LogSumExp(ToulBar2::logZ ,ToulBar2::GlobalLogUbZ) - wcsp->LogSumExp(ToulBar2::logZ,ToulBar2::GlobalLogLbZ)<= Log1p(ToulBar2::sigma)) break;
      //~ if (ToulBar2::hbfs && nbRecomputationNodes > 0) { // wait until a nonempty open node is restored (at least after first global solution is found)
        //~ assert(nbNodes > 0);
        //~ if (nbRecomputationNodes > nbNodes / ToulBar2::hbfsBeta && ToulBar2::hbfs <= ToulBar2::hbfsGlobalLimit) ToulBar2::hbfs *= 2;
        //~ else if (nbRecomputationNodes < nbNodes / ToulBar2::hbfsAlpha && ToulBar2::hbfs >= 2) ToulBar2::hbfs /= 2;
      //~ if (ToulBar2::debug >= 2) cout << "HBFS backtrack limit: " << ToulBar2::hbfs << endl;
      //~ }
    //~ }
    //~ //cout<<"open size after finished : "<<open_->size()<<" "<<open_->empty()<<endl;
  //~ } else {
			//~ cluster->hbfsGlobalLimit = LONGLONG_MAX;
			//~ cluster->hbfsLimit = LONGLONG_MAX;
			//~ pair<Cost, Cost> res = recursiveSolve(cluster, clb, cub);
  //~ }
//~ }

