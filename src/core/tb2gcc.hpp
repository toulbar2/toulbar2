/** \file tb2gcc.hpp
 *  \brief Propagates a Global Cardinality Constraint.
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

#ifndef TB2GCC_HPP_
#define TB2GCC_HPP_

#include <utility>
#include <variant>
#include "tb2abstractconstr.hpp"
#include "tb2ternaryconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2wcsp.hpp"
#include "search/tb2clusters.hpp"
#include "utils/tb2lapjv.hpp"
#include <unordered_set>
#include <unordered_map>

using namespace std;

class GlobalCardinalityConstraint : public AbstractNaryConstraint {
    Cost Original_ub; // initial upper bound when creating the constraint
    StoreCost lb; // Projected cost to problem lower bound (if zero, all deltaCosts must be zero)
    StoreCost assigneddeltas; // Accumulated deltas from assigned values (used in cost propagation)
    vector<Long> conflictWeights; // Used by weighted degree heuristics to prioritize variables
    vector<vector<StoreCost>> deltaCosts; // Extended unary costs to all values (2D cost matrix)
    vector<StoreValue> storeLastAssignment; // Stores the last optimal assignment of variables
    StoreInt storeAssignment; // True if store optimal assignment
    int *rowSol = nullptr; // Row solution pointer for assignment
    Cost* ReduceCostRow = nullptr; // Row reduction costs 
    Cost* ReduceCostCol = nullptr; // Column reduction costs 
    vector<Cost> costMatrix; // Flattened cost matrix 
    int NbValues; // Total number of values in the domain
    vector<int> AssignedVar; // List of currently assigned variables
    vector<int> AssignedVal; // List of values assigned to variables
    vector<int> NoAssignedVar; // List of unassigned variables
    int NbAssigned; // Number of assigned variables
    int NbNoAssigned; // Number of unassigned variables
    vector<bool> isAssignedValue; // Flags for whether each value is assigned
    vector<bool> varAlreadyProcessed; // Flags for whether a variable has already been processed
    vector<Value> exceptedValues; // Values to be excluded 
    vector<int> exceptedValIndex; // Indices of excluded values
    bool excepted; // Flag indicating if no excluded values
    bool isSquare; // Indicates if the cost matrix is square 
    vector<string> UnionVarDomain; // Combined domain of all variables 
    vector<int> VarDomainSize; // Domain size per variable
    unordered_map<string, int> mapDomainValToIndex; // Maps domain values to indices in UnionVarDomain
    bool SameDomain; // True if all variables have the same domain
    map<Value, pair<int, int>> bounds; // lower and upper bound capacities for every value

    void projectLB(Cost c)
    {
        if (c > MIN_COST) {
            if(!isSquare) lb += c;
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
        	if(!isSquare)
        	{
		        if (td && scope[var]->canbe(value))
		            td->addDelta(cluster, scope[var], value, -C);
		        deltaCosts[var][value_idx] += C;
                }
                assert(scope[var]->getCost(value) >= C);
                scope[var]->extend(value, C);
        } else if (C < MIN_COST) {
        	if(!isSquare)
        	{
		        if (td && scope[var]->canbe(value))
		            td->addDelta(cluster, scope[var], value, -C);
		        deltaCosts[var][value_idx] += C;
                }
                scope[var]->project(value, -C, true);
        }
    }

public:
    GlobalCardinalityConstraint(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in)
        : AbstractNaryConstraint(wcsp, scope_in, arity_in)
        , Original_ub(wcsp->getUb())
        , lb(0) 
        , assigneddeltas(0) 
        , storeAssignment(false) 
        , NbValues(0)
        , NbAssigned(0)
        , NbNoAssigned(arity_in)
        , excepted(false)
        , isSquare(false)
        , SameDomain(true) 
        {
        if (arity_in > 0) { 
            unordered_set<string> seen; 
            isSquare = true; 
            NbValues = 0; 
            int DomainSize;
            int varIndex = 0;
            // add dummy value names needed by AllDifferent
            for (int var_ind = 0; var_ind < arity_in; var_ind++) {
                auto* variable = scope[var_ind];
                if(variable->isValueNames()) continue;
                for (int val_ind = 0; val_ind < (int)variable->getDomainInitSize() ; val_ind++) {
                    variable->addValueName(to_string("v") + to_string(variable->toValue(val_ind)));
                }
            }
            conflictWeights.push_back(0); 

            // Get the domain size of the first variable
            auto* variable = scope[varIndex];
            DomainSize = (int)variable->getDomainInitSize();
            VarDomainSize.push_back(DomainSize);
           // cout<<"variable : x"<<varIndex<<" domain size : "<<DomainSize<<endl;

            // Collect values from the first variable's domain
            for (int valIndex= 0; valIndex< DomainSize; valIndex++) {
               // cout<<variable->getValueName(valIndex)<<" ";
                if (seen.insert(variable->getValueName(valIndex)).second) {
                    UnionVarDomain.push_back(variable->getValueName(valIndex));
                    NbValues++;
                }
                
            }
           // cout<<endl;

            // Repeat for the rest of the variables
            for (int varIndex = 1; varIndex< arity_in; varIndex++) {
                conflictWeights.push_back(0); 
                auto* variable = scope[varIndex];
                DomainSize = (int)variable->getDomainInitSize();
               // cout<<"variable : x"<<varIndex<<" domain size : "<<DomainSize<<endl;
                VarDomainSize.push_back(DomainSize);
                if(VarDomainSize[varIndex -1] != VarDomainSize[varIndex]) SameDomain = false;

                // Add only new unique values to the union
                for (int valIndex = 0; valIndex < DomainSize; valIndex++) {
                    //cout<<variable->getValueName(valIndex)<<" ";
                    if (seen.insert(variable->getValueName(valIndex)).second) {
                        UnionVarDomain.push_back(variable->getValueName(valIndex));
                        NbValues++;
                        SameDomain = false; // Different domains detected
                    }
                }

            }

            // Map each unique value to its index
            for (int valIndex = 0; valIndex < NbValues; valIndex++) {
                mapDomainValToIndex[UnionVarDomain[valIndex]] = valIndex;
            }

            // If more values than variables, cost matrix is not square
            if (NbValues > arity_in) {
                isSquare = false;
                for (int varIndex = 0; varIndex < arity_in; varIndex++) {
                    // Initialize extended cost vector for each variable
                    deltaCosts.emplace_back(VarDomainSize[varIndex], MIN_COST);
                }
            }
            // Initialize 
            storeLastAssignment = vector<StoreValue>(arity_in, StoreValue(WRONG_VAL));
            NoAssignedVar = vector<int>(arity_in, -1); 
            AssignedVar = vector<int>(arity_in, -1);
            AssignedVal = vector<int>(arity_in, -1); 
            costMatrix = vector<Cost>(arity_ * NbValues, MAX_COST);

            // Allocate memory 
            rowSol = new int[arity_]; 
            ReduceCostRow = new Cost[arity_]; 
            ReduceCostCol = new Cost[NbValues]; 
        } else {
            deconnect();
        }
    }

    virtual ~GlobalCardinalityConstraint() {}

    void read(istream& file) // TODO: add a parameter for controlling the level of propagation if necessary
    {
        int nbValues = 0;
        int sumlb = 0;
        int sumub = 0;
        int lower = 0;
        int upper = 0;
        file >> nbValues;
        for (int v = 0; v < nbValues; v++) {
            Value value;
            file >> value;
            file >> lower;
            sumlb += lower;
            file >> upper;
            sumub += upper;
            bounds[value] = {lower, upper};
        }
        if (sumlb > arity_) THROWCONTRADICTION;
        if (sumub < arity_) THROWCONTRADICTION;
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
        if (isSquare)
        {
            Cost res = MIN_COST;
            Cost nbsame = 0;
            vector<bool> alreadyUsed(arity_, false);
            for (int varIndex = 0; varIndex < arity_; varIndex++) {
                if (alreadyUsed[s[varIndex]]) {
                    auto it = find(exceptedValues.begin(), exceptedValues.end(), scope[varIndex]->toValue(s[varIndex]));
                    nbsame += (it == exceptedValues.end());
                } else {
                    alreadyUsed[s[varIndex]] = true;
                }
            }
            if (nbsame > 0) {
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

        Cost res = -lb + assigneddeltas;
        Cost nbsame = 0;
        vector<bool> alreadyUsed(NbValues, false);
        for (int varIndex = 0; varIndex < arity_; varIndex++) {
            res += deltaCosts[varIndex][s[varIndex]];
            int valIndex = mapDomainValToIndex[scope[varIndex]->getValueName(s[varIndex])];
            if (alreadyUsed[valIndex]) {
                auto it = find(exceptedValIndex.begin(), exceptedValIndex.end(), valIndex);
                nbsame += (it == exceptedValIndex.end());
            } else {
                alreadyUsed[valIndex] = true;
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

    Cost getMaxFiniteCost() override ///< \brief returns the maximum finite cost for any valid tuple less than wcsp->getUb()
    {
    	if(isSquare && !excepted) return MIN_COST;
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

    Cost getCost(int index, Value val)
    {
        assert(index >= 0 && index < arity_);
        return deltaCosts[index][scope[index]->toIndex(val)];
    }

    double computeTightness() override { return 0; } // TODO: compute factorial(n)/(n**n)?

    // TODO: needed for dominance test by DEE
    // pair<pair<Cost, Cost>, pair<Cost, Cost>> getMaxCost(int index, Value a, Value b)

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


    bool RemoveAssignVar()
    {  

        /**
        * @brief Filters and removes assigned values from unassigned variables in the AllDifferent constraint.
        * 
        * This function enforces the AllDifferent constraint by:
        * - Verifying that assigned variables do not share the same value.
        * - Removing values from unassigned variable domains if already assigned elsewhere.
        * - Propagating newly assigned variables and updating internal state.
        * - Triggering full n-ary projection if thresholds on domain sizes or variable count are met.
        * 
        * @return true if constraint propagation was successful without needing full projection.
        * @return false if the constraint switched to full n-ary projection.
        * 
        * @throws Throws a contradiction exception (THROWCONTRADICTION) if two variables are assigned the same value.
        */

        // Initialize tracking vectors for assigned values and processed variables
        isAssignedValue = vector<bool>(NbValues, false);
        varAlreadyProcessed = vector<bool>(arity_, false);

        // Clear and reserve space for assigned and unassigned variables containers
        AssignedVar.clear();
        AssignedVar.reserve(arity_);       
        AssignedVal.clear();
        AssignedVal.reserve(arity_);
        NoAssignedVar.clear();
        NoAssignedVar.reserve(arity_);

        // Identify assigned and unassigned variables
        for (int varIndex = 0; varIndex < arity_; ++varIndex) {
            auto* variable = scope[varIndex];

            if (!variable->assigned()) {
                // Variable not assigned yet: add to unassigned list
                NoAssignedVar.push_back(varIndex);
            } else {
                // Variable assigned: find the index of its assigned value
                int valIndex = mapDomainValToIndex[variable->getValueName(variable->toIndex(variable->getValue()))];
                auto it = find(exceptedValIndex.begin(), exceptedValIndex.end(), valIndex);
                if (!isAssignedValue[valIndex]) {
                    AssignedVar.push_back(varIndex);
                    if(it == exceptedValIndex.end()){
                        // Mark this value as assigned and track the variable
                        isAssignedValue[valIndex] = true;
                        AssignedVal.push_back(valIndex); 
                    } 
                       
                } else {
                    // Contradiction: same value assigned to more than one variable
                    if(it == exceptedValIndex.end()) THROWCONTRADICTION;
                }
            }
        }

        NbAssigned = AssignedVar.size();
        bool VarAssigned;
        bool NaryPro = false;

        // Propagation loop: remove assigned values from unassigned variables' domains
        do {
            VarAssigned = false;
            vector<int> newlyAssignedVars;

            // For each assigned value, check all unassigned variables
            for (int valIndex = 0; valIndex < NbValues; ++valIndex) {
                if (!isAssignedValue[valIndex]) continue; // Skip values not assigned

                for (int varIndex : NoAssignedVar) {
                    if (varAlreadyProcessed[varIndex]) continue; // Skip already processed vars

                    auto* variable = scope[varIndex];
                    Value value = variable->toValue(variable->toIndex(UnionVarDomain[valIndex]));

                    if (variable->canbe(value) && (find(exceptedValIndex.begin(), exceptedValIndex.end(), valIndex) == exceptedValIndex.end())) {
                        // Remove the assigned value from the domain of unassigned variable
                        variable->remove(value);

                        if (variable->assigned()) {
                            // If variable becomes assigned after removal, mark for propagation
                            newlyAssignedVars.push_back(varIndex);
                            varAlreadyProcessed[varIndex] = true;
                        }
                    }
                }
            }

            // Process newly assigned variables after domain filtering
            for (int varIndex : newlyAssignedVars) {
                if (connected(varIndex)) {
                    // Disconnect variable from constraint and update counts
                    deconnect(varIndex);
                    NbNoAssigned--;
                    NbAssigned++;

                    // Check whether to switch to full n-ary projection based on thresholds
                    if (getNonAssigned() <= NARYPROJECTIONSIZE &&
                        (getNonAssigned() <= 1 ||
                         prodInitDomSize <= NARYPROJECTIONPRODDOMSIZE ||
                         maxInitDomSize <= NARYPROJECTION3MAXDOMSIZE ||
                         (getNonAssigned() == 2 && maxInitDomSize <= NARYPROJECTION2MAXDOMSIZE))) {
                        deconnect();
                        projectNary();
                        NaryPro = true;
                        break;
                    } else {
                        // Update assigned values and variables tracking
                        VarAssigned = true;
                        auto* variable = scope[varIndex];
                        int valIndex = mapDomainValToIndex[variable->getValueName(variable->toIndex(variable->getValue()))];
                        if(find(exceptedValIndex.begin(), exceptedValIndex.end(), valIndex) == exceptedValIndex.end()){
                            isAssignedValue[valIndex] = true;
                            AssignedVal.push_back(valIndex);
                        }
                        AssignedVar.push_back(varIndex);
                        
                    }
                } else {
                    // Variable not connected to constraint: stop propagation
                    return false;
                }
            }

            // Update list of unassigned variables after propagation
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
    /**
     * @brief Propagate the AllDifferent constraint to enforce consistency and improve problem lower bound.
     * 
     * This method performs constraint propagation by checking assigned variables,
     * removing inconsistent values, computing lower bounds using the Jonker algorithm,
     * and updating unary costs and supports accordingly.
     * 
     * It handles different cases depending on whether variables are assigned or not,
     * if the domains are the same, and whether "excepted" values are considered.
     * 
     * @throws TimeOut if the propagation is interrupted.
     * @throws ContradictionException if the constraint becomes unsatisfiable.
     */

        if (ToulBar2::dumpWCSP % 2) // skip propagation if the problem is dumped before preprocessing
            return;

        if (ToulBar2::interrupted && !ToulBar2::isZ) {
            throw TimeOut();
        }

        // Propagate only if the constraint is connected
        if (connected()) {
            if (ToulBar2::verbose >= 3) {
                cout << "propagate " << *this << endl;
            }

            bool b = false;

            // Attempt to assign variables already connected and assigned
            for (int varIndex = 0; !b && connected() && varIndex < arity_; varIndex++) {
                if (connected(varIndex) && scope[varIndex]->assigned()) {
                    assign(varIndex);
                    b = true;
                }
            }

            if (!b && connected()) {
                // Check whether propagation should be skipped based on stored results
                bool skipPropagation = false;
                if (storeAssignment) {
                    skipPropagation = true;
                    for (int varIndex = 0; varIndex < arity_; varIndex++) {
                        if (scope[varIndex]->cannotbe(storeLastAssignment[varIndex]) || scope[varIndex]->getCost(storeLastAssignment[varIndex]) > MIN_COST) {
                            skipPropagation = false;
                            break;
                        }
                    }
                }

                if (!skipPropagation) {
                    // Initialize number of unassigned variables
                    NbNoAssigned = getNonAssigned();
                    bool filtreExcepted = false;

                    // If there are assigned variables, handle accordingly
                    if (NbNoAssigned < arity_) {
                        if (excepted) {
                            filtreExcepted = RemoveAssignVar();
                        } else if (SameDomain && RemoveAssignVar()) {
                            // Collect unassigned values
                            vector<int> NoAssignedVal;
                            for (int valIndex = 0; valIndex < NbValues; ++valIndex) {
                                if (!isAssignedValue[valIndex]) {
                                    NoAssignedVal.push_back(valIndex);
                                }
                            }

                            int NbNoAssignedVal = NoAssignedVal.size();

                            // Initialize cost matrix for the Jonker algorithm
                            Cost current_ub = wcsp->getUb() - wcsp->getLb();

                            for (int varInd = 0; varInd < NbNoAssigned; ++varInd) {
                                int varIndex = NoAssignedVar[varInd];
                                auto* variable = scope[varIndex];

                                for (int valInd = 0; valInd < NbNoAssignedVal; ++valInd) {
                                    int valIndex = NoAssignedVal[valInd];
                                    Value val = variable->toValue(valIndex);
                                    costMatrix[varInd * NbNoAssignedVal + valInd] = variable->canbe(val) ? variable->getCost(val) : current_ub;
                                }
                            }

                            // Solve assignment problem using Jonker algorithm
                            Cost TotalCost = lapjv(NbNoAssigned, NbNoAssignedVal, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub);

                            if (TotalCost >= current_ub) {
                                THROWCONTRADICTION;
                            } else if (TotalCost >= 0) {
                                Cost jonker = TotalCost;

                                // Store results from Jonker algorithm
                                for (int varIndex = 0; varIndex < NbNoAssigned; ++varIndex) {
                                    auto* variable = scope[NoAssignedVar[varIndex]];
                                    storeLastAssignment[NoAssignedVar[varIndex]] = variable->toValue(variable->toIndex(UnionVarDomain[NoAssignedVal[rowSol[varIndex]]]));
                                }

                                for (int varIndex  = 0; varIndex  < NbAssigned; ++varIndex) {
                                    auto* variable = scope[AssignedVar[varIndex]];
                                    storeLastAssignment[AssignedVar[varIndex]] = variable->toValue(variable->toIndex(UnionVarDomain[AssignedVal[varIndex]]));
                                }

                                storeAssignment = true;

                                // Update the lower bound with the Jonker algorithm's total cost
                                projectLB(jonker);

                                // Adjust unary costs for unassigned variables based on reduced costs
                                for (int varInd = 0; varInd < NbNoAssigned; ++varInd) {
                                    int varIndex = NoAssignedVar[varInd];
                                    auto* variable = scope[varIndex];

                                    for (int valInd = 0; valInd < NbNoAssignedVal; ++valInd) {
                                        int valIndex = NoAssignedVal[valInd];
                                        Value value = variable->toValue(valIndex);

                                        if (variable->canbe(value)) {
                                            ExtOrProJ(varIndex, value, (ReduceCostRow[varInd] + ReduceCostCol[valInd]));
                                        }
                                    }
                                }

                                // Update support values if needed for unassigned variables
                                for (int varInd  = 0; varInd  < NbNoAssigned; ++varInd ) {
                                    int varIndex = NoAssignedVar[varInd];
                                    auto* variable = scope[varIndex];
                                    Value optimalValue = storeLastAssignment[varIndex];
                                    if (variable->getSupport() != optimalValue) {
                                        if (ToulBar2::verbose > 0)
                                            cout << "CHANGE GCC SUPPORT " << variable->getName() << " from " << variable->getSupport() << " to " << optimalValue << endl;
#ifndef NDEBUG
                                        variable->queueEAC1(); // EAC support may have been lost
#endif
                                        variable->setSupport(optimalValue);
                                        assert(variable->getCost(variable->getSupport()) == MIN_COST);
                                    }
                                }
                            }
                        }
                    }

                    // Case when all variables are unassigned or filtering exception or domains differ
                    if (NbNoAssigned == arity_ || filtreExcepted || !SameDomain) {
                        Cost current_ub = wcsp->getUb() - wcsp->getLb();

                        // Initialize the cost matrix for all variables and their domain values
                        for (int varIndex = 0; varIndex < arity_; ++varIndex) {
                            auto* variable = scope[varIndex];

                            for (int valIndex = 0; valIndex < VarDomainSize[varIndex]; ++valIndex) {
                                Value value = variable->toValue(valIndex);
                                string valName = variable->getValueName(valIndex);
                                costMatrix[varIndex * NbValues + mapDomainValToIndex[valName]] = variable->canbe(value) ? variable->getCost(value) : current_ub;
                            }
                        }

                        // Solve the Linear Assignment Problem (LAP) using the Jonker algorithm
                        Cost TotalCost;
                        if (excepted)
                            TotalCost = lapjv(arity_, NbValues, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub, exceptedValIndex);
                        else
                            TotalCost = lapjv(arity_, NbValues, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub);

                        if (TotalCost >= current_ub) {
                            THROWCONTRADICTION;
                        } else if (TotalCost >= 0) {
                            // Store results from Jonker algorithm
                            for (int varIndex = 0; varIndex < arity_; varIndex++) {
                                auto* variable = scope[varIndex];
                                storeLastAssignment[varIndex] = variable->toValue(variable->toIndex(UnionVarDomain[rowSol[varIndex]]));
                            }

                            storeAssignment = true;

                            // Update lower bound with the total cost found
                            projectLB(TotalCost);

                            // Adjust unary costs using reduced row and column costs
                            for (int varIndex = 0; varIndex < arity_; ++varIndex) {
                                auto* variable = scope[varIndex];

                                for (int valIndex = 0; valIndex < VarDomainSize[varIndex]; ++valIndex) {
                                    Value value = variable->toValue(valIndex);
                                    string valName = variable->getValueName(valIndex);

                                    if (variable->canbe(value)) {
                                        ExtOrProJ(varIndex, value, (ReduceCostRow[varIndex] + ReduceCostCol[mapDomainValToIndex[valName]]));
                                    }
                                }
                            }

                            // Update support if needed for all variables
                            for (int varIndex = 0; varIndex < arity_; ++varIndex) {
                                auto* variable = scope[varIndex];
                                Value optimalValue = storeLastAssignment[varIndex];
                                if (variable->getSupport() != optimalValue) {
                                    if (ToulBar2::verbose > 0)
                                        cout << "CHANGE GCC SUPPORT " << variable->getName() << " from " << variable->getSupport() << " to " << optimalValue << endl;
#ifndef NDEBUG
                                    variable->queueEAC1(); // EAC support may have been lost
#endif
                                    variable->setSupport(optimalValue);
                                    assert(variable->getCost(variable->getSupport()) == MIN_COST);
                                }
                            }
                        }
                    }
                }
            }
        }
    }



    //checks that the constraint is still satisfiable (called by WCSP::verify in Debug mode at each search node)
    bool verify() override
    {
        if (!storeAssignment) {
            return false;
        }
        vector<bool> alreadyUsed(NbValues, false);
        for (int i = 0; i < arity_; i++) {
            int valIndex = mapDomainValToIndex[scope[i]->getValueName(scope[i]->toIndex(storeLastAssignment[i]))];
            if ((alreadyUsed[valIndex] && (!excepted || find(exceptedValIndex.begin(), exceptedValIndex.end(), valIndex) == exceptedValIndex.end())) || scope[i]->cannotbe(storeLastAssignment[i]) || scope[i]->getCost(storeLastAssignment[i]) > MIN_COST) {
                if (alreadyUsed[valIndex]) {
                    cout << "variable " << scope[i]->getName() << " value " << storeLastAssignment[i] << " used twice!" << endl;
                } else if (scope[i]->cannotbe(storeLastAssignment[i])) {
                    cout << "variable " << scope[i]->getName() << " value " << storeLastAssignment[i] << " has been removed!" << endl;
                } else if (scope[i]->getCost(storeLastAssignment[i]) > MIN_COST) {
                    cout << "variable " << scope[i]->getName() << " value " << storeLastAssignment[i] << " has nonzero cost!" << endl;
                }
                return false;
            } else {
                alreadyUsed[valIndex] = true;
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
            if (!storeAssignment || scope[index]->cannotbe(storeLastAssignment[index])) {
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
        if (!storeAssignment || scope[index]->cannotbe(storeLastAssignment[index]) || scope[index]->getCost(storeLastAssignment[index]) > MIN_COST) {
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
        os << this << " gcc(";

        int unassigned_ = 0;
        for (int i = 0; i < arity_; i++) {
            if (scope[i]->unassigned())
                unassigned_++;
            os << wcsp->getName(scope[i]->wcspIndex);
            if (i < arity_ - 1)
                os << ",";
        }
		if(!isSquare)
		{
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
		}
        os << "/" << getTightness();
        if (ToulBar2::weightedDegree) {
            os << "/" << getConflictWeight();
            for (int i = 0; i < arity_; i++) {
                os << "," << conflictWeights[i];
            }
        }
        os << " " << bounds.size() << ":[";
        for (const auto& [key, bound] : bounds) {
            os << " (" << key;
            os << "," << bound.first;
            os << "," << bound.second;
            os << " )";
        }
        os << "] arity: " << arity_;
        os << " unassigned: " << getNonAssigned() << "/" << unassigned_ << endl;
    }

    void dump(ostream& os, bool original = true) override
    {
        assert(lb == MIN_COST); // TODO: how to dump with deltaCosts?
        if (original) {
            os << arity_;
            for (int i = 0; i < arity_; i++)
                os << " " << scope[i]->wcspIndex;
            os << " -1 gcc ";
        } else {
            os << getNonAssigned();
            for (int i = 0; i < arity_; i++)
                if (scope[i]->unassigned())
                    os << " " << scope[i]->getCurrentVarId();
            os << " -1 gcc ";
        }
        os << bounds.size();
        for (const auto& [key, bound] : bounds) {
            os << " " << key;
            os << " " << bound.first;
            os << " " << bound.second;
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
            os << "],\"type\":\"gcc\",\"params\":{\"bounds\":[";
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
            os << "],\"type\":\"gcc\",\"params\":{\"bounds\":[";
        }
        printed = false;
        for (const auto& [key, bound] : bounds) {
            if (printed)
                os << ",";
            os << "[" << key;
            os << "," << bound.first;
            os << "," << bound.second;
            os << "]";
            printed = true;
        }
        os << "]}},\n";
    }

};
#endif /*TB2GCC_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

