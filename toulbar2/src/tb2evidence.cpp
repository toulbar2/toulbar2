#include "tb2solver.hpp"
#include "tb2trienum.hpp"
#include "tb2binconstr.hpp"
#include "tb2enumvar.hpp"

TLogProb Solver::Zub(){ // Calculate an uper-bound on Z before exploration (step 0)
    
    //~ if(ToulBar2::prodsumDiffusion>0)
    //~ {
      //~ //PropagateNoc();
      //~ ProdSumDiffusion();
    //~ }
    TLogProb newlogU;
     vector<Cost> vbinmin(unassignedVars->getSize(),0);
    Cost newCost = wcsp->getLb() + wcsp->getNegativeLb();
    switch(ToulBar2::isZUB){
      //Upper bound on Z edition 0
    case 0 :
      newCost += wcsp->LogProb2Cost(unassignedVars->getSize() * Log(wcsp->getMaxDomainSize()));
      newlogU = wcsp->LogSumExp(ToulBar2::logU, newCost);
      break;
      //Upper bound on Z edition 1
    case 1 :
      //cout<<"p0 : "<<wcsp->Cost2Prob(newCost)<<endl;
      for (BTList<Value>::iterator iter_variable = unassignedVars->begin(); iter_variable != unassignedVars->end(); ++iter_variable) { // Loop on the unassigned variables
        EnumeratedVariable *var = (EnumeratedVariable *) ((WCSP *) wcsp)->getVar(*iter_variable);
        //cout<<"Variable "<<var->getName()<<endl;
        Cost SumUnaryCost = MAX_COST;
        if (wcsp->enumerated(*iter_variable)) {
          for (EnumeratedVariable::iterator iter_value = var->begin(); iter_value != var->end(); ++iter_value) { // loop over the domain of the variable
            SumUnaryCost = wcsp->LogSumExp(SumUnaryCost,var->getCost(*iter_value)); // Sum of the exponential of Unary cost over the domain.
            //cout<<"UNARYZ : "<<wcsp->Cost2Prob(var->getCost(*iter_value))<<endl;
          }
          newCost += SumUnaryCost; //Sum the older cost with the new log(sum(unarycost))
          //cout<<"NewCost : "<<wcsp->Cost2Prob(newCost)<<endl;
        }
        else {
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
          EnumeratedVariable *y = (EnumeratedVariable *) ((WCSP *) wcsp)->getVar(i); // get the variable i
          //cout<<"Variable "<<y->getName()<<endl;
          for(unsigned int j = i+1; j < wcsp->numberOfVariables(); j++ ) {
            if (wcsp->unassigned(j)) {
              EnumeratedVariable *x=(EnumeratedVariable *) ((WCSP *) wcsp)->getVar(j); // get the variable j;
              //cout<<"link to Variable "<<x->getName()<<endl;
              BinaryConstraint *bctr = x->getConstr(y);// get constr that link i to j
              
              if(bctr != NULL){
                Cost cbinmin = MAX_COST;
                for (EnumeratedVariable::iterator ity = y->begin(); ity != y->end(); ++ity) { // Loop on "this" values
                  for (EnumeratedVariable::iterator itx = x->begin(); itx != x->end(); ++itx) { // Loop on x values
                    Cost cbin = bctr->getCost(y,x,*ity,*itx); // get the binary cost linking "this,it" and "x,itx"
                    //cout<<cbin<< ' ';
                    if (cbin<cbinmin) cbinmin = cbin;
                  }
                }
                //cout<<endl;
                vbinmin[i]=cbinmin;
              }
              else{vbinmin[i]=0;}
            }
          }
          Cost SumUnaryCost = MAX_COST;
          for (EnumeratedVariable::iterator ity = y->begin(); ity != y->end(); ++ity) { // Loop on "this" values
            SumUnaryCost = wcsp->LogSumExp(SumUnaryCost,y->getCost(*ity)); // Sum of the exponential of Unary cost over the domain.
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
        EnumeratedVariable *var = (EnumeratedVariable *) ((WCSP *) wcsp)->getVar(*iter_variable);
        //cout<<"Variable "<<var->getName()<<endl;
        Cost SumUnaryCost = MAX_COST;
        if (wcsp->enumerated(*iter_variable)) {
          for (EnumeratedVariable::iterator iter_value = var->begin(); iter_value != var->end(); ++iter_value) { // loop over the domain of the variable
            SumUnaryCost = wcsp->LogSumExp(SumUnaryCost,var->getCost(*iter_value)); // Sum of the exponential of Unary cost over the domain.
          }
          newCost += SumUnaryCost; //Sum the older cost with the new log(sum(unarycost))
        }
        else {
          newCost += wcsp->LogProb2Cost(Log(wcsp->getDomainSize(*iter_variable)));
        }
      }
      
      for (unsigned int i = 0; i < wcsp->numberOfVariables()-1; i++) { // Loop on the unassigned variables
        if (wcsp->unassigned(i)) {
        EnumeratedVariable *y = (EnumeratedVariable *) ((WCSP *) wcsp)->getVar(i); // get the variable i
        for(unsigned int j = i+1; j < wcsp->numberOfVariables(); j++ ) {
           //~ cout<<"Variable "<<i<<" "<<j<<endl;
           if (wcsp->unassigned(j)) {
            EnumeratedVariable *x=(EnumeratedVariable *) ((WCSP *) wcsp)->getVar(j); // get the variable j;
            BinaryConstraint *bctr = x->getConstr(y);; // get constr that link i to j
            if(bctr != NULL){
              Cost SumBinaryCost = MAX_COST;
              //~ cout<<"BinaryZ : "<<endl;
              for (EnumeratedVariable::iterator ity = y->begin(); ity != y->end(); ++ity) { // Loop on "this" values
                for (EnumeratedVariable::iterator itx = x->begin(); itx != x->end(); ++itx) { // Loop on x values
                  Cost cbin = bctr->getCost(y,x,*ity,*itx); // get the binary cost linking "this,it" and "x,itx"
                  //~ cout<<wcsp->Cost2Prob(cbin)<<' ';
                  SumBinaryCost = wcsp->LogSumExp(SumBinaryCost,cbin);
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
      newlogU=numeric_limits<TLogProb>::infinity();
    }
     return newlogU;
}
//~ 
//~ // Function that took a node root and a list of sons of this root
//~ // and return a list of (value,cost) 
//~ WCSP::ResultVisitZ WCSP::visitZ(int root, WCSP::Sons &sons)
//~ {
    //~ assert(unassigned(root));
    //~ EnumeratedVariable *x = (EnumeratedVariable *) getVar(root);
    //~ ResultVisitZ res;
    //~ 
    //~ for (EnumeratedVariable::iterator iter = x->begin(); iter != x->end(); ++iter) { // Loop over the domain of the variable root
        //~ res.push_back( pair<Value, Cost>(*iter, x->getCost(*iter)) ); // Fill the res list with the unarycost of root's sons
    //~ }
    //~ 
    //~ for (unsigned int i=0; i<sons[root].size(); i++) { // Loop until there is no root's sons.
		//~ //Recurcive call of visitZ over the root's son i 
        //~ ResultVisitZ resi = visitZ(sons[root][i].first, sons);
//~ 
        //~ // Reuse code developed for EnumeratedVariable::elimVar( BinaryConstraint* ctr )
        //~ EnumeratedVariable *y = (EnumeratedVariable *) getVar(sons[root][i].first); // Take the i-th sons of the root
        //~ BinaryConstraint *ctr = sons[root][i].second; // Take the binary constraint that link the i-th sons and the root
        //~ int pos = 0;
        //~ for (EnumeratedVariable::iterator iter1 = x->begin(); iter1 != x->end(); ++iter1, pos++) { //Loop over the value in root
            //~ Cost mincost = MAX_COST;
            //~ for (unsigned int iter = 0; iter < resi.size(); ++iter) {
				//~ //SUM UNARY of the i-th sons AND BINARY that link the i-th sons to the root
                //~ Cost curcost = resi[iter].second + ctr->getCost(y, x, resi[iter].first, *iter1); 
                //~ mincost = LogSumExp(mincost,curcost);
            //~ }
            //~ res[pos].second += mincost;
        //~ }
    //~ }
    //~ 
    //~ return res;
//~ }
//~ 
//~ TLogProb WCSP::spanningTreeZ(Cost c0)
//~ {
  //~ pair<vector<int>,Sons> RootSon = spanningTree();
  //~ vector<int> roots = RootSon.first;
  //~ Sons sons = RootSon.second;
  //~ TLogProb res = 0 ;
//~ 
  //~ for (int i = roots.size()-1; i >= 0; i--) { // Loop over all the roots nodes
      //~ ResultVisitZ resi = visitZ(roots[i], sons); // Construct a list of (value,cost) for each roots. The value term contain the bynary cost and the unary cost.
      //~ Cost mincost = MAX_COST;
      //~ for (unsigned int iter = 0; iter < resi.size(); ++iter) { // Loop over all the nodes in the root's list
		  //~ //Bring back together all the costs that go from a node to a leaf in the spanning tree : Dynamic programming
          //~ mincost = LogSumExp(mincost,resi[iter].second); 
      //~ }
		//~ // Total Sum on the spanning tree
      //~ res += Cost2LogProb(mincost) ; 
  //~ }
  //~ return res + Cost2LogProb(c0);
//~ }


///////////////: Gumbel Perturbation (Non concluant) ////////////
TLogProb Solver::GumofThrone(){

  TLogProb LogZhat = 0;
  for(auto iter : ToulBar2::trieZ->get_sols() ){
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
					EnumeratedVariable *var = (EnumeratedVariable *) ((WCSP *) wcsp)->getVar(i); // get the variable i
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
	EnumeratedVariable* x;
	Constraint* ctr = NULL;
  
	for(ConstraintList::iterator itc = getConstrs()->begin();itc != getConstrs()->end(); ++itc ) { //constraints linking the "this" variable
		ctr = (*itc).constr;
    if(ctr->arity() == 2 && !ctr->isSep() && ctr != NULL) {
			BinaryConstraint* bctr = (BinaryConstraint*) ctr;
			x = (EnumeratedVariable*) bctr->getVarDiffFrom( (Variable*) this ); // get variable different from "this"
			//~ cout<<"link to variable "<< x->getName()<<endl;
      for (iterator it = this->begin(); it != this->end(); ++it) {// Loop on "this" values
        Cost csum = MAX_COST;
        //~ cout<<"Binary Cost before : ";
				for (iterator itx = x->begin(); itx != x->end(); ++itx) { // Loop on x values
          Cost cbin = bctr->getCost(this,x,*it,*itx); // get the binary cost linking "this,it" and "x,itx"
          //~ cout<<wcsp->Cost2Prob(cbin)<<' ';
          csum = wcsp->LogSumExp(csum,cbin);
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
              bctr->addcost(this,x,*it,*itx,-cost_extc);
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
          project(*it,cost_extc);
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
					EnumeratedVariable *var = (EnumeratedVariable *) ((WCSP *) wcsp)->getVar(i); // get the variable i
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
  if(wcsp->Cost2Prob(csum)!=1.){
    for (iterator it = begin(); it != end(); ++it) {
      project(*it,-csum);
      //cout<<"After Value "<<toIndex(*it)<<" : "<<wcsp->Cost2Prob(getCost(*it))<<endl;
    }
    if (csum < MIN_COST){
       //cout<<"OldNegLB : "<<wcsp->Cost2Prob(wcsp->getNegativeLb())<<endl;
      //cout<<"DecreaseLB : "<<wcsp->Cost2Prob(csum)<<endl;
      wcsp->decreaseLb(csum);
      //cout<<"NewNegLB : "<<wcsp->Cost2Prob(wcsp->getNegativeLb())<<endl;
    }else wcsp->increaseLb(csum);
  }
}

bool EnumeratedVariable::Normalization()
{
  //~ vector<Cost> vsum;
  bool change = false;
  Constraint* ctr = NULL;
  EnumeratedVariable* x;
  
  for(ConstraintList::iterator itc = getConstrs()->begin();itc != getConstrs()->end(); ++itc ) { //constraints linking the "this" variable
		ctr = (*itc).constr;
    if(ctr->arity() == 2 && !ctr->isSep() && ctr != NULL) {
			BinaryConstraint* bctr = (BinaryConstraint*) ctr;
			x = (EnumeratedVariable*) bctr->getVarDiffFrom( (Variable*) this ); // get variable different from "this"
			//~ cout<<"link to variable "<< x->getName()<<endl;
      for (iterator it = this->begin(); it != this->end(); ++it) { // Loop on "this" values
        Cost csum = MAX_COST;
				for (iterator itx = x->begin(); itx != x->end(); ++itx) { // Loop on x values
          csum = wcsp->LogSumExp(csum,bctr->getCost(this,x,*it,*itx));
        }
        if(csum<0){
          change=true;
          for (iterator itx = x->begin(); itx != x->end(); ++itx) {
            bctr->addcost(this,x,*it,*itx,-csum);
          }
          cout<<endl;
          project(*it,csum);
        }
        else{change=false;}
      }
    }
	}
  UnaryNormalization();
  return change;
}

