#include "tb2solver.hpp"
#include "tb2trienum.hpp"
#include "tb2binconstr.hpp"
#include "tb2enumvar.hpp"

TLogProb Solver::Zub(){ // Calculate an uper-bound on Z before exploration (step 0)

    TLogProb newlogU;
    //cout<<wcsp->getNegativeLb()<<endl;
    Cost newCost = wcsp->getLb() + wcsp->getNegativeLb();
        
    switch(ToulBar2::isZUB){
      //Upper bound on Z edition 0
    case 0 :
      newCost += wcsp->LogProb2Cost(unassignedVars->getSize() * Log10(wcsp->getMaxDomainSize()));
      newlogU = wcsp->LogSumExp(ToulBar2::logU, newCost);
      break;
      //Upper bound on Z edition 1
    case 1 :
      for (BTList<Value>::iterator iter_variable = unassignedVars->begin(); iter_variable != unassignedVars->end(); ++iter_variable) { // Loop on the unassigned variables
        EnumeratedVariable *var = (EnumeratedVariable *) ((WCSP *) wcsp)->getVar(*iter_variable);
        Cost SumUnaryCost = MAX_COST;
        if (wcsp->enumerated(*iter_variable)) {
          for (EnumeratedVariable::iterator iter_value = var->begin(); iter_value != var->end(); ++iter_value) { // loop over the domain of the variable
            SumUnaryCost = wcsp->LogSumExp(SumUnaryCost,var->getCost(*iter_value)); // Sum of the exponential of Unary cost over the domain.
          }
        }
        else {
          newCost += wcsp->LogProb2Cost(Log10(wcsp->getDomainSize(*iter_variable)));
        }
        newCost += SumUnaryCost; //Sum the older cost with the new log(sum(unarycost))
      }
      newlogU = wcsp->LogSumExp(ToulBar2::logU, newCost);
      break;

      //Upper bound on Z edition 2
    case 2 :
      newlogU = wcsp->LogSumExp(ToulBar2::logU, wcsp->spanningTreeZ(newCost));
      break;
    // Full Partition function
    default :
      newlogU=numeric_limits<TLogProb>::infinity();
    }
     return newlogU;
}

TLogProb Solver::UnderTheZ(){

  TLogProb LogZhat = -numeric_limits<TProb>::infinity();
  for(auto iter : ToulBar2::trieZ->get_sols() ){
    TLogProb iter_logprob = wcsp->Cost2LogProb((iter + wcsp->getNegativeLb()));
    LogZhat = wcsp->LogSumExp(LogZhat,iter_logprob);
  }
  return LogZhat;
}


TLogProb Solver::GumofThrone(){

  TLogProb LogZhat = 0;
  for(auto iter : ToulBar2::trieZ->get_sols() ){
    LogZhat += iter;
  }
  LogZhat = LogZhat / ToulBar2::trieZ->get_sols().size();
  return LogZhat;
}

void Solver::ProdSumDiffusion()
{
	for (int times = 0; times < 2; times++) {
		bool change = true;
		int maxit = ToulBar2::prodsumDiffusion;
		cout << "ProdSumDiffusion: " << endl;
		//~ cout << "   max iterations " << maxit << endl;
		//cout << "   C0 = " << wcsp->getLb() << endl;
		int ntimes = 0;
		while (change && (ntimes < maxit)) { 
			change = false;
			int nchanged = 0;
			for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) //loop on the variable
				if (wcsp->unassigned(i)) {
					EnumeratedVariable *var = (EnumeratedVariable *) ((WCSP *) wcsp)->getVar(i); // get the variable i
					if (var->Marginalisation()) { // if marginalisation is over 
						change = true;
						nchanged++;
						//var->findSupport();
					}
				}
			ntimes++;
			//~ cout << "iteration " << ntimes << "   changed: " << nchanged << endl;
		}
		//~ cout << "   done iterations: " << ntimes << endl;
	}
}

bool EnumeratedVariable::Marginalisation()
{
	//Cost Top = wcsp->getUb();
	bool change = false;
	EnumeratedVariable* x;
	Constraint* ctr = NULL;
  Cost csum = MAX_COST;
	ConstraintList::iterator itc = getConstrs()->begin();
	if(itc != getConstrs()->end())	ctr = (*itc).constr;
	while(ctr) { // While there is constraint linking the "this" variable
		if(ctr->arity() == 2 && !ctr->isSep()) {
			BinaryConstraint* bctr = (BinaryConstraint*) ctr;
			x = (EnumeratedVariable*) bctr->getVarDiffFrom( (Variable*) this ); // get variable different from "this"
			for (EnumeratedVariable::iterator it = this->begin(); it != this->end(); ++it) { // Loop on "this" values
				Cost cu = getCost(*it); // get the unary cost of the ith value of "this"
				//Cost cmin = Top;
				for (EnumeratedVariable::iterator itx = x->begin(); itx != x->end(); ++itx) { // Loop on x values
					Cost cbin = bctr->getCost(this,x,*it,*itx); // get the binary cost linking "this,it" and "x,itx"
          csum = wcsp->LogSumExp(csum,cbin);
        }
				Double mean = to_double(csum) / 2.;
        //cout<<"BinSum : "<<mean<<' ';
				Double extc = to_double(cu) + mean;
        //cout<<"Delta : "<< extc<<' '<<endl;				 
				//~ if(abs(extc) >= 1) {
				  //~ Cost costi = (Long) extc;
				  //~ for (iterator itx = x->begin(); itx != x->end(); ++itx) {
            //~ bctr->addcost(this,x,*it,*itx,-costi);				
				  //~ }
          //~ project(*it, costi); 
				  //~ change = true;
				//~ }
			}
      cout<<endl;
		}
		++itc;
		if(itc != getConstrs()->end()) ctr = (*itc).constr;
		else ctr = NULL;
	}
	return change;
}
