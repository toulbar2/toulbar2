#include "tb2solver.hpp"
#include "tb2trienum.hpp"
#include "tb2binconstr.hpp"
#include "tb2enumvar.hpp"


TLogProb Solver::Zub()  // Calculate an uper-bound on Z before exploration (step 0)
{
  //~ if(ToulBar2::prodsumDiffusion>0)
  //~ {
  //~ //PropagateNoc();
  //~ ProdSumDiffusion();
  //~ }
  TLogProb newlogU;
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
  case 3: //test
    
    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {// Loop on the unassigned variables
      if (wcsp->unassigned(i)) {
	EnumeratedVariable *y = (EnumeratedVariable *)((WCSP *) wcsp)->getVar(i);  // get the variable i
	//cout<<"Variable "<<y->getName()<<endl;
	for (unsigned int j = i + 1; j < wcsp->numberOfVariables(); j++) {
	  if (wcsp->unassigned(j)) {
	    EnumeratedVariable *x = (EnumeratedVariable *)((WCSP *) wcsp)->getVar(j); // get the variable j;
	    //cout<<"link to Variable "<<x->getName()<<endl;
	    BinaryConstraint *bctr = x->getConstr(y);// get constr that link i to j
	    if (bctr != NULL) {
	      Cost cbinmin = MAX_COST;
	      for (EnumeratedVariable::iterator ity = y->begin(); ity != y->end(); ++ity) { // Loop on "this" values
		for (EnumeratedVariable::iterator itx = x->begin(); itx != x->end(); ++itx) { // Loop on x values
		  Cost cbin = bctr->getCost(y, x, *ity, *itx); // get the binary cost linking "this,it" and "x,itx"
		  //cout<<cbin<< ' ';
		  if (cbin < cbinmin) cbinmin = cbin;
		}
	      }
	      //cout<<endl;
	      vbinmin[i] = cbinmin;
	    } else {vbinmin[i] = 0;}
	  }
	}
	Cost SumUnaryCost = MAX_COST;
	for (EnumeratedVariable::iterator ity = y->begin(); ity != y->end(); ++ity) { // Loop on "this" values
	  SumUnaryCost = wcsp->LogSumExp(SumUnaryCost, y->getCost(*ity)); // Sum of the exponential of Unary cost over the domain.
	}
	
	//cout << wcsp->Cost2Prob(SumUnaryCost) <<endl;
	//cout<<"minbin : "<< wcsp->Cost2Prob(vbinmin[i])<<endl;
	newCost += (SumUnaryCost + vbinmin[i]); //Sum the older cost with the new log(sum(unarycost))
	//newCost += vbinmin[i];
      }
    }
    newlogU = wcsp->LogSumExp(ToulBar2::logU, newCost);
    break;
  case 4:
    for (BTList<Value>::iterator iter_variable = unassignedVars->begin(); iter_variable != unassignedVars->end(); ++iter_variable) { // Loop on the unassigned variables
      EnumeratedVariable *var = (EnumeratedVariable *)((WCSP *) wcsp)->getVar(*iter_variable);
      //cout<<"Variable "<<var->getName()<<endl;
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
    
    for (unsigned int i = 0; i < wcsp->numberOfVariables() - 1; i++) { // Loop on the unassigned variables
      if (wcsp->unassigned(i)) {
	EnumeratedVariable *y = (EnumeratedVariable *)((WCSP *) wcsp)->getVar(i);  // get the variable i
	for (unsigned int j = i + 1; j < wcsp->numberOfVariables(); j++) {
	  //~ cout<<"Variable "<<i<<" "<<j<<endl;
	  if (wcsp->unassigned(j)) {
	    EnumeratedVariable *x = (EnumeratedVariable *)((WCSP *) wcsp)->getVar(j); // get the variable j;
	    BinaryConstraint *bctr = x->getConstr(y);; // get constr that link i to j
	    if (bctr != NULL) {
	      Cost SumBinaryCost = MAX_COST;
	      //~ cout<<"BinaryZ : "<<endl;
	      for (EnumeratedVariable::iterator ity = y->begin(); ity != y->end(); ++ity) { // Loop on "this" values
		for (EnumeratedVariable::iterator itx = x->begin(); itx != x->end(); ++itx) { // Loop on x values
		  Cost cbin = bctr->getCost(y, x, *ity, *itx); // get the binary cost linking "this,it" and "x,itx"
		  //~ cout<<wcsp->Cost2Prob(cbin)<<' ';
		  SumBinaryCost = wcsp->LogSumExp(SumBinaryCost, cbin);
		}
		//~ cout<<endl;
	      }
	      //~ cout<< "SumBinaryCost : " <<wcsp->Cost2Prob(SumBinaryCost)<<endl;
	      newCost += SumBinaryCost;
	    }
	    //cout<<"Newcost : " <<wcsp->Cost2LogProb(newCost)<<endl;
	  }
	}
      }
    }
    newlogU = wcsp->LogSumExp(ToulBar2::logU, newCost);
    break;
    // Full Partition function
  default :
    newlogU = numeric_limits<TLogProb>::infinity();
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
    Cost delta = MIN_COST;
	
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
	Cost bestlb = MAX(nd.getCost(delta), wcsp->getLb());
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
      if (wcsp->LogSumExp(wcsp->LogSumExp(ToulBar2::logZ ,ToulBar2::GlobalLogUbZ),ToulBar2::logU) - wcsp->LogSumExp(ToulBar2::logZ,ToulBar2::GlobalLogLbZ)<= Log1p(ToulBar2::sigma)) break;
      if (wcsp->LogSumExp(ToulBar2::logZ ,ToulBar2::GlobalLogUbZ) - wcsp->LogSumExp(ToulBar2::logZ,ToulBar2::GlobalLogLbZ) < Log(ToulBar2::sigma)) break;
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

