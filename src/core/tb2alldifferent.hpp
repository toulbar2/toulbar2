/** \file tb2alldifferent.hpp
 *  \brief Propagates an All-Different or a Permutation Constraint.
 *
 *  It maintains different levels of consistency:
 *  - tests if the constraint is always satisfied or violated when enough variables are assigned (NC level)
 *  TODO - enforces domain/bound arc consistency considering the hard constraint only (AC level)
 *  TODO - approximates full zero-inverse soft consistency (FDAC or EDAC or VAC/AC)
 *  TODO - enforces virtual arc consistency (VAC-lin)
 *
 *  TOBEVERIFIED full zero-inverse soft consistency may be not enforced because it requires propagating at every unary cost move and not only when it increases/projects from zero???
 *
 */

#ifndef TB2ALLDIFFERENT_HPP_
#define TB2ALLDIFFERENT_HPP_

#include <utility>
#include <variant>
#include "tb2abstractconstr.hpp"
#include "tb2ternaryconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2wcsp.hpp"
#include "search/tb2clusters.hpp"
#include "utils/tb2lapjv.hpp"

using namespace std;

class AllDifferentConstraint : public AbstractNaryConstraint {
    Cost Original_ub; // initial upper bound when creating the constraint
    StoreCost lb; // projected cost to problem lower bound (if it is zero then all deltaCosts must be zero)
    StoreCost assigneddeltas;
    vector<Long> conflictWeights; // used by weighted degree heuristics
    vector<vector<StoreCost>> deltaCosts; // extended costs from unary costs to values in the cost function
    vector<StoreInt> storeResults; // store the last optimal assignment
    StoreInt isResults;
    int *rowsol = nullptr;
    Cost* ReduceCostRow = nullptr;
    Cost* ReduceCostCol= nullptr;
    Cost * cost_matrix = nullptr;
    int NbValues;
    vector<int> AssignedVar;
    vector<int> AssignedVal;
    vector<int> NoAssignedVar;
    int NbAssigned;
    int NbNoAssigned;
    vector<bool> isAssignedValue;
    vector<bool> varAlreadyProcessed;
    vector<Value> exceptedValues;

    void projectLB(Cost c)
    {
        if (c > MIN_COST) {
            lb += c;
            Constraint::projectLB(c);
        }
    }
	
 
    // returns true if the constraint can be projected to small local cost function in extension
    //bool canbeProjectedInExtension();

    // Depending of the value and the cost, extend or project the cost on the index value of the variable var
    void ExtOrProJ(int var, Value value, Cost C)
    {
        int value_idx = scope[var]->toIndex(value);
        TreeDecomposition* td = wcsp->getTreeDec();
        if (C > MIN_COST) {
                if (td && scope[var]->canbe(value))
                    td->addDelta(cluster, scope[var], value, -C);
                deltaCosts[var][value_idx] += C;
                assert(scope[var]->getCost(value) >= C);
                scope[var]->extend(value, C);
        } else if (C < MIN_COST) {
                if (td && scope[var]->canbe(value))
                    td->addDelta(cluster, scope[var], value, -C);
                deltaCosts[var][value_idx] += C;
                scope[var]->project(value, -C, true);
        }
    }
	

public:
    AllDifferentConstraint(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in)
        : AbstractNaryConstraint(wcsp, scope_in, arity_in)
        , Original_ub(wcsp->getUb())
        , lb(0)
        , assigneddeltas(0)
		, isResults(false)
    {
    	if (arity_in > 0) {
    		conflictWeights.push_back(0);
    		NbValues = (int)scope[0]->getDomainInitSize() ;

    		assert( NbValues >= arity_in); //TODO: do not assume identical initial domain size equal to the arity
    		deltaCosts.emplace_back(NbValues, MIN_COST);
    		int DomainSize;
    		for (int i = 1; i < arity_in; i++) {
    			conflictWeights.push_back(0);
    			DomainSize = (int)scope[i]->getDomainInitSize();
    			assert(DomainSize >= arity_in); //TODO: do not assume identical initial domain size equal to the arity
    			deltaCosts.emplace_back(DomainSize, MIN_COST);
    			if (DomainSize> NbValues) NbValues = DomainSize;
    		}
    		if(NbValues < arity_in) THROWCONTRADICTION;
    		storeResults = vector<StoreInt>(arity_in, StoreInt(-1));
    		NoAssignedVar = vector<int>(arity_in, -1);
    		AssignedVar= vector<int>(arity_in, -1);
    		AssignedVal= vector<int>(arity_in, -1);
    		propagate();
    	} else {
            deconnect();
        }
    }

    virtual ~AllDifferentConstraint() {}

    void read(istream& file) // TODO: add a parameter for controlling the level of propagation if necessary
    {
        int nbExcepted = 0;
        file >> nbExcepted;
        for (int v=0; v < nbExcepted; v++) {
            Value except;
            file >> except;
            exceptedValues.push_back(except);
        }
    }

    bool extension() const FINAL { return false; } // TODO: allows functional variable elimination but no other preprocessing

    Cost getDefCost() FINAL { return MAX_COST; }

    Long getConflictWeight() const override { return Constraint::getConflictWeight(); }
    Long getConflictWeight(int varIndex) const override
    {
        assert(varIndex >= 0);
        assert(varIndex < arity_);
        return conflictWeights[varIndex] + Constraint::getConflictWeight();
    }
    void incConflictWeight(Constraint* from) override
    {
        assert(from != NULL);
        if (from == this) {
            if (getNonAssigned() == arity_ || deconnected()) {
                Constraint::incConflictWeight(1);
            } else {
                for (int i = 0; i < arity_; i++) {
                    if (connected(i)) {
                        conflictWeights[i]++;
                    }
                }
            }
        } else if (deconnected()) {
            for (int i = 0; i < from->arity(); i++) {
                int index = getIndex(from->getVar(i));
                if (index >= 0) { // the last conflict constraint may be derived from two binary constraints (boosting search), each one derived from an n-ary constraint with a scope which does not include parameter constraint from
                    assert(index < arity_);
                    conflictWeights[index]++;
                }
            }
        }
    }
    void resetConflictWeight() override
    {
        conflictWeights.assign(conflictWeights.size(), 0);
        Constraint::resetConflictWeight();
    }

    /// \brief returns true if constraint always satisfied and has (less than) zero cost only
    //bool universal(Cost zero = MIN_COST) override;

    Cost eval(const Tuple& s) override
    {
        // returns the cost of the corresponding assignment s
        Cost res = -lb + assigneddeltas;
        Cost nbsame = 0;
        vector<bool> alreadyUsed(NbValues, false);
        for (int i = 0; i < arity_; i++) {
            assert(s[i] < NbValues);
            res += deltaCosts[i][s[i]];
            if (alreadyUsed[s[i]]) {
                nbsame++;
            } else {
                alreadyUsed[s[i]] = true;
            }
        }
        if (nbsame > 0 || res > wcsp->getUb()) {
            if (nbsame > 0 && Original_ub < wcsp->getUb() && 1.0L * Original_ub * nbsame < wcsp->getUb()) {
                res = Original_ub * nbsame; // VNS-like methods may exploit a relaxation of the constraint
            } else {
                res = wcsp->getUb();
            }
        }
        assert(res <= wcsp->getUb());
        assert(res >= MIN_COST);
        return res;
    }

    Cost getCost(int index, Value val)
    {
        assert(index >= 0 && index < arity_);
        return deltaCosts[index][scope[index]->toIndex(val)];
    }

    double computeTightness() override { return 0; } // TODO: compute factorial(n)/(n**n)?

    // TODO: needed for dominance test by DEE
    // pair<pair<Cost, Cost>, pair<Cost, Cost>> getMaxCost(int index, Value a, Value b)

    Cost getMaxFiniteCost() override ///< \brief returns the maximum finite cost for any valid tuple less than wcsp->getUb()
    {
        Cost sumdelta = - lb + assigneddeltas;
        for (int i = 0; i < arity_; i++) {
            Cost m = *max_element(deltaCosts[i].begin(), deltaCosts[i].end());
            if (m > MIN_COST)
                sumdelta += m;
        }
        if (CUT(sumdelta, wcsp->getUb()))
            return MAX_COST;
        else
            return sumdelta;
    }

    // void setInfiniteCost(Cost ub)
    void setInfiniteCost(Cost ub) override
    {
        Original_ub = min(ub, Original_ub);
    }

    void assign(int varIndex) override
    {
        if (ToulBar2::verbose >= 7) {
            cout << "assign " << scope[varIndex]->getName() << " in " << *this << endl;
        }
        if (connected(varIndex)) {
            deconnect(varIndex);
            assert(getNonAssigned() >= 0);

            if (getNonAssigned() <= NARYPROJECTIONSIZE && (getNonAssigned() <= 1 || prodInitDomSize <= NARYPROJECTIONPRODDOMSIZE || maxInitDomSize <= NARYPROJECTION3MAXDOMSIZE || (getNonAssigned() == 2 && maxInitDomSize <= NARYPROJECTION2MAXDOMSIZE))) {
                deconnect();
                projectNary();
            } else {
                // TODO: incremental bound propagation
                propagate();
                if (ToulBar2::FullEAC)
                    reviseEACGreedySolution();
            }
        }
    }
  bool filtreAndPropagate()
  {  
        isAssignedValue = vector<bool>(NbValues, false);
        varAlreadyProcessed = vector<bool>(arity_, false);
        AssignedVar.clear();
        AssignedVar.reserve(arity_);       
        AssignedVal.clear();
        AssignedVal.reserve(arity_);
        NoAssignedVar.clear();
        NoAssignedVar.reserve(arity_);

        for (int var = 0; var < arity_; ++var) {
            auto* variable = scope[var];
            if (!variable->assigned()) {
                NoAssignedVar.push_back(var);
            } else {
                int index_val = variable->toIndex(variable->getValue());
                if(!isAssignedValue[index_val]){
	                isAssignedValue[index_val] = true;
	                AssignedVar.push_back(var);
	                AssignedVal.push_back(index_val); 
                }
                else {
                    THROWCONTRADICTION;
                }
            }
        }
    
        NbAssigned = AssignedVar.size();
        bool VarAssigned;
        bool NaryPro = false;
    
        do {
            VarAssigned = false;
            vector<int> newlyAssignedVars;
    
            for (int val = 0; val < NbValues; ++val) {
                if (!isAssignedValue[val]) continue;
    
                for (int varIndex : NoAssignedVar) {
                    if (varAlreadyProcessed[varIndex]) continue;
    
                    auto* variable = scope[varIndex];
                    Value value = variable->toValue(val);
    
                    if (variable->canbe(value)) {
                        variable->remove(value);
                        if (variable->assigned()) {
                            newlyAssignedVars.push_back(varIndex);
                            varAlreadyProcessed[varIndex] = true;

                        }
                    }
                }
            }
    
            for (int varIndex : newlyAssignedVars) {
                if (connected(varIndex)) {
                    deconnect(varIndex);
                    NbNoAssigned--;
                    NbAssigned++;   
                    if (getNonAssigned() <= NARYPROJECTIONSIZE && (getNonAssigned() <= 1 || prodInitDomSize <= NARYPROJECTIONPRODDOMSIZE || maxInitDomSize <= NARYPROJECTION3MAXDOMSIZE || (getNonAssigned() == 2 && maxInitDomSize <= NARYPROJECTION2MAXDOMSIZE))) {
                        deconnect();
                        projectNary();
                        NaryPro = true;
                        break;
                    } else {
                        VarAssigned = true;
                        auto* variable = scope[varIndex];
                        int valIndex = variable->toIndex(variable->getValue());
                        isAssignedValue[valIndex] = true;
                        AssignedVar.push_back(varIndex);
                        AssignedVal.push_back(valIndex);
                    }
                }
                else return (false);
            }
    
            if (VarAssigned && !NaryPro) {
                vector<int> updatedNoAssigned;
                for (int varIndex : NoAssignedVar) {
                    if (!scope[varIndex]->assigned()) {
                        updatedNoAssigned.push_back(varIndex);
                    }
                }
                NoAssignedVar.swap(updatedNoAssigned);
            }
    
        } while (VarAssigned);
        
        return (!NaryPro);
   }

    void propagate() override
    {
        if (ToulBar2::dumpWCSP % 2) // do not propagate if problem is dumped before preprocessing
            return;
        if (ToulBar2::interrupted && !ToulBar2::isZ) {
            throw TimeOut();
        }
        // propagates from scratch the constraint
        if (connected()) {
            if (ToulBar2::verbose >= 3) {
                cout << "propagate " << *this << endl;
            }
            bool b = false;
            for (int i = 0; !b && connected() && i < arity_; i++) { // SdG: do not continue to assign if one assign already done
                if (connected(i) && scope[i]->assigned()) {
                    assign(i);
                    b = true;
                }
            }
            if (!b && connected()) {
                // TODO: compute a lower bound and prune forbidden domain values
            
                // Determine whether propagation should be skipped
		bool skipPropagation = false;
		if (isResults) {
			skipPropagation = true;
			for (int var = 0; var < arity_; var++) {
				if (scope[var]->cannotbe(scope[var]->toValue(storeResults[var])) || scope[var]->getCost(scope[var]->toValue(storeResults[var])) > MIN_COST) {
					skipPropagation = false;
					break;
				}
			}
		}
  		if (!skipPropagation) {
                    // Initialize total cost and determine assigned/unassigned variables
                    NbNoAssigned = getNonAssigned();
                    
                    if (NbNoAssigned < arity_) {
               
                        if (filtreAndPropagate()) {
                        
                            vector<int> NoAssignedVal;
			    for (int val = 0; val < NbValues; ++val) {
				if (!isAssignedValue[val]) {
				    NoAssignedVal.push_back(val);
				}
			    }
                            
                            int NbNoAssignedVal =  NoAssignedVal.size();
                           // Initialize cost matrix for the Jonker algorithm
                            Cost curent_ub = wcsp->getUb() - wcsp->getLb() ;
			    delete[] cost_matrix;
		            delete[] ReduceCostRow;
		            delete[] ReduceCostCol;
		            delete[] rowsol;
			    
			    cost_matrix = new Cost[NbNoAssigned*NbNoAssignedVal];			    
			    rowsol = new int[NbNoAssigned];
			    ReduceCostRow = new Cost[NbNoAssigned];
			    ReduceCostCol = new Cost[NbNoAssignedVal];
			  
                            for (int i = 0; i < NbNoAssigned; ++i) {
                                int varIndex = NoAssignedVar[i];
                                auto* variable = scope[varIndex];
                                for (int j = 0; j < NbNoAssignedVal; ++j) {
                                    int valIndex = NoAssignedVal[j];
                                    Value val = variable->toValue(valIndex);
                                    cost_matrix[i * NbNoAssignedVal + j] = variable->canbe(val) ? variable->getCost(val) : curent_ub;    
                                }
                            }
			    Cost TotalCost = lapjv(NbNoAssigned, NbNoAssignedVal, cost_matrix ,  rowsol,  ReduceCostRow, ReduceCostCol, curent_ub );
                    
		            if (TotalCost >= curent_ub ) {
					THROWCONTRADICTION;
                            }else if (TotalCost >=0 ){
								
                                Cost jonker =TotalCost;
                                for (int i = 0; i < NbNoAssigned; ++i) {

                                    storeResults[NoAssignedVar[i]] = NoAssignedVal[rowsol[i]];
                                }

                                for (int i = 0; i < NbAssigned; ++i) {
                                    storeResults[AssignedVar[i]] = AssignedVal[i];
                                }
		                isResults = true;
		                // Update the lower bound with the Jonker algorithm's total cost
		                projectLB(jonker);
		            
		                // Modify unary costs using ExtOrProJ, for unassigned variables only

		                for (int i = 0; i < NbNoAssigned; ++i) {
		                    int var = NoAssignedVar[i];
		                    auto* variable = scope[var];
		            
		                    for (int j = 0; j < NbNoAssignedVal; ++j) {
		                        int index_val = NoAssignedVal[j];
		                        Value val = variable->toValue(index_val);
		            
		                        if (variable->canbe(val)) {
		                            //Adjust the cost using reduced row and column costs

		                            ExtOrProJ(var, val, (ReduceCostRow[i] + ReduceCostCol[j]));
		                        }

		                    }

		                }

	  		        // Update support if needed              
				for (int i = 0; i < NbNoAssigned; ++i) {
				    int var = NoAssignedVar[i];
	                    	    auto* variable = scope[var];		            
				    Value opt_val = variable->toValue(storeResults[var]);
				    if (variable->getSupport() != opt_val)  {
				  	  if (ToulBar2::verbose > 0)
						cout << "CHANGE ALLDIFF SUPPORT " << variable->getName() << " from " << variable->getSupport() << " to " << opt_val << endl;
#ifndef NDEBUG
					variable->queueEAC1(); // EAC support may have been lost
#endif
					variable->setSupport(opt_val);

					assert(variable->getCost(variable->getSupport()) == MIN_COST);
				    }
			    	}
		           }
			}
		    }
 

		    else{

                        // Initialize cost matrix for the Jonker algorithm
				   
			//Cost *cost_matrix;
    		        delete[] cost_matrix;
		        delete[] ReduceCostRow;
		        delete[] ReduceCostCol;
		        delete[] rowsol;

			cost_matrix = new Cost[arity_ * NbValues];
		        
		        Cost curent_ub = wcsp->getUb() - wcsp->getLb() +1;	
		        for (int varIndex = 0; varIndex < arity_; ++varIndex) {
		                auto* variable = scope[varIndex];

		            
		                for (int valIndex = 0; valIndex < NbValues; ++valIndex) {
		                    Value val = variable->toValue(valIndex);
		            
		                    cost_matrix[varIndex * NbValues + valIndex] = variable->canbe(val) ? variable->getCost(val) : curent_ub;
		                }
		        }

		            
	
		            // Solve the Linear Assignment Problem (LAP) using the Jonker algorithm
			    rowsol = new int[arity_];
			    ReduceCostRow = new Cost[arity_];
			    ReduceCostCol = new Cost[NbValues];
			    Cost TotalCost = lapjv(arity_, NbValues , cost_matrix ,  rowsol,  ReduceCostRow, ReduceCostCol, curent_ub );

	
		            if (TotalCost >= curent_ub  ) {
		   
		            
		                THROWCONTRADICTION;
		            } else if (TotalCost >=0) {
	       
				    for (int var = 0; var < arity_; var++){
					storeResults[var] = rowsol[var];

				    }
				    isResults = true;

				    Cost jonker = TotalCost;
					      
				   // Update the lower bound with the Jonker algorithm's total cost
				   projectLB(jonker);
			    
				  // Modify unary costs using ExtOrProJ
				  for (int var = 0; var < arity_; ++var) {
				      auto* variable = scope[var];
			    
				      for (int index_val = 0; index_val < NbValues; ++index_val) {
				          Value val = variable->toValue(index_val);
			    
				          if (variable->canbe(val)) {
				            	ExtOrProJ(var, val, (ReduceCostRow[var] + ReduceCostCol[index_val]));
				          }
				     }
				  }

		 
				// Update support if needed
				for (int var = 0; var < arity_; ++var) {
				    auto* variable = scope[var];		            
				    
				    Value opt_val = variable->toValue(storeResults[var]);
				    if (variable->getSupport() != opt_val)  {
				  	  if (ToulBar2::verbose > 0)
						cout << "CHANGE ALLDIFF SUPPORT " << variable->getName() << " from " << variable->getSupport() << " to " << opt_val << endl;
		#ifndef NDEBUG
					variable->queueEAC1(); // EAC support may have been lost
		#endif

					variable->setSupport(opt_val);
					if (variable->cannotbe(variable->getSupport()))
						cout<< "var "<<endl;
					assert(variable->getCost(variable->getSupport()) == MIN_COST);
				    }
			   	 }		 
		             }		
		        }
		  }	
            }           
        }
    }

  //TODO: checks that the constraint is still satisfiable (called by WCSP::verify in Debug mode at each search node)
 //bool verify() override;
 bool verify() override
    {
        if (!isResults) {
            return false;
        }
        vector<bool> alreadyUsed(NbValues, false);
        for (int i = 0; i < arity_; i++) {
            if (alreadyUsed[storeResults[i]] || scope[i]->cannotbe(scope[i]->toValue(storeResults[i])) || scope[i]->getCost(scope[i]->toValue(storeResults[i])) > MIN_COST) {
                if (alreadyUsed[storeResults[i]]) {
                    cout << "variable " << scope[i]->getName() << " value " << scope[i]->toValue(storeResults[i]) << " used twice!" << endl;
                } else if (scope[i]->cannotbe(scope[i]->toValue(storeResults[i]))) {
                    cout << "variable " << scope[i]->getName() << " value " << scope[i]->toValue(storeResults[i]) << " has been removed!" << endl;
                } else if (scope[i]->getCost(scope[i]->toValue(storeResults[i])) > MIN_COST) {
                    cout << "variable " << scope[i]->getName() << " value " << scope[i]->toValue(storeResults[i]) << " has nonzero cost!" << endl;
                }
                return false;
            } else {
                alreadyUsed[storeResults[i]] = true;
            }
        }
        return true;
    }

    void increase(int index) override
    {
        remove(index);
    }
    void decrease(int index) override
    {
        remove(index);
    }
    void remove(int index) override
    {
        if (scope[index]->unassigned()) {
            bool revise = ToulBar2::FullEAC && getVar(index)->cannotbe(getVar(index)->getSupport());
            if (!isResults || scope[index]->cannotbe(scope[index]->toValue(storeResults[index]))) {
	        propagate();
	    }
            if (revise)
                reviseEACGreedySolution();
        } else {
            assign(index);
        }
    }
    void projectFromZero(int index) override
    {
        // TODO: incremental cost propagation
        bool revise = ToulBar2::FullEAC && (getVar(index)->cannotbe(getVar(index)->getSupport()) || getVar(index)->getCost(getVar(index)->getSupport()) > MIN_COST);
        if (!isResults || scope[index]->cannotbe(scope[index]->toValue(storeResults[index])) || scope[index]->getCost(scope[index]->toValue(storeResults[index])) > MIN_COST) {
	    propagate();
	}
        if (revise)
            reviseEACGreedySolution();
    }

    //bool checkEACGreedySolution(int index = -1, Value supportValue = 0) FINAL; // TODO

    bool reviseEACGreedySolution(int index = -1, Value supportValue = 0) FINAL
    {
        bool result = checkEACGreedySolution(index, supportValue);
        if (!result) {
            if (index >= 0) {
                getVar(index)->unsetFullEAC();
            } else {
                int a = arity();
                for (int i = 0; i < a; i++) {
                    getVar(i)->unsetFullEAC();
                }
            }
        }
        return result;
    }
	
	
    void print(ostream& os) override
    {
        os << this << " alldifferent(";

        int unassigned_ = 0;
        for (int i = 0; i < arity_; i++) {
            if (scope[i]->unassigned())
                unassigned_++;
            os << wcsp->getName(scope[i]->wcspIndex);
            if (i < arity_ - 1)
                os << ",";
        }
        os << ") "
           << " cost: " << -lb << " + " << assigneddeltas << " + (";
        for (int i = 0; i < arity_; i++) {
            for (unsigned int j = 0; j < deltaCosts[i].size(); j++) {
                os << deltaCosts[i][j];
                if (j < deltaCosts[i].size() - 1)
                    os << "|";
            }
            if (i < arity_ - 1)
                os << ",";
        }
        os << ") ";
        os << "/" << getTightness();
        if (ToulBar2::weightedDegree) {
            os << "/" << getConflictWeight();
            for (int i = 0; i < arity_; i++) {
                os << "," << conflictWeights[i];
            }
        }
        os << " exceptedvalues:";
        for (Value v : exceptedValues) {
            os << " " << v;
        }
        os << " arity: " << arity_;
        os << " unassigned: " << getNonAssigned() << "/" << unassigned_ << endl;
     }

    void dump(ostream& os, bool original = true) override
    {
        assert(lb == MIN_COST); // TODO: how to dump with deltaCosts?
        if (original) {
            os << arity_;
            for (int i = 0; i < arity_; i++)
                os << " " << scope[i]->wcspIndex;
            os << " -1 alldiff ";
        } else {
            os << getNonAssigned();
            for (int i = 0; i < arity_; i++)
                if (scope[i]->unassigned())
                    os << " " << scope[i]->getCurrentVarId();
            os << " -1 alldiff ";
        }
        os << exceptedValues.size();
        for (Value v : exceptedValues) {
            os << " " << v;
        }
        os << endl;
    }

    void dump_CFN(ostream& os, bool original = true) override
    {
        assert(lb == MIN_COST); // TODO: how to dump with deltaCosts?
        bool printed = false;
        os << "\"F_";

        if (original) {
            printed = false;
            for (int i = 0; i < arity_; i++) {
                if (printed)
                    os << "_";
                os << scope[i]->wcspIndex;
                printed = true;
            }

            os << "\":{\"scope\":[";
            printed = false;
            for (int i = 0; i < arity_; i++) {
                if (printed)
                    os << ",";
                os << "\"" << name2cfn(scope[i]->getName()) << "\"";
                printed = true;
            }
            os << "],\"type\":\"alldiff\",\"params\":{\"exceptedvalues\":[";
        } else {
            for (int i = 0; i < arity_; i++)
                if (scope[i]->unassigned()) {
                    if (printed)
                        os << "_";
                    os << scope[i]->getCurrentVarId();
                    printed = true;
                }
            os << "\":{\"scope\":[";
            printed = false;
            for (int i = 0; i < arity_; i++)
                if (scope[i]->unassigned()) {
                    if (printed)
                        os << ",";
                    os << "\"" << name2cfn(scope[i]->getName()) << "\"";
                    printed = true;
                }
            os << "],\"type\":\"alldiff\",\"params\":{\"exceptedvalues\":[";
        }
        printed = false;
        for (Value v : exceptedValues) {
            if (printed)
                os << ",";
            os << v;
            printed = true;
        }
        os << "]}},\n";
    }
};
#endif /*TB2ALLDIFFERENT_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
