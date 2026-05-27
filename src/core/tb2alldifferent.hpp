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

#pragma once
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
#include <boost/heap/pairing_heap.hpp>
#include <random>
#include <algorithm>

using namespace std;
using PQ = boost::heap::pairing_heap<pair<Cost, int>, boost::heap::compare<greater<>>>;

class AllDifferentConstraint : public AbstractNaryConstraint {
    Cost Original_ub;           // Initial upper bound when creating the constraint
    StoreCost lb;               // Projected cost to problem lower bound (if zero, all deltaCosts must be zero)
    StoreCost assigneddeltas;   // Accumulated deltas from assigned values (used in cost propagation)
    vector<Long> conflictWeights;           // Used by weighted degree heuristics to prioritize variables
    vector<StoreCost> deltaCosts;   // Extended unary costs to all values 
    vector<int> ValuesCapacity; 
    vector<StoreValue> storeLastAssignment; // Stores the last optimal assignment of variables
    StoreInt storeAssignment;   // True if an optimal assignment is currently stored
    int* rowSol = nullptr;      // Row solution array: rowSol[var] = value index assigned to variable var by the LAP solver
    Cost* ReduceCostRow = nullptr; // Reduced costs per row (per variable) after solving the LAP
    Cost* ReduceCostCol = nullptr; // Reduced costs per column (per value) after solving the LAP
    vector<Cost> costMatrix;    // Flattened row-major cost matrix passed to the Jonker LAP solver
    int NbValues;               // Total number of distinct values across all variable domains (size of UnionVarDomain)
    vector<int> AssignedVar;    // Indices (in scope) of currently assigned variables
    vector<int> AssignedVal;    // Value indices (in UnionVarDomain) assigned to each variable in AssignedVar
    vector<int> NoAssignedVar;  // Indices (in scope) of currently unassigned variables
    int NbAssigned;             // Number of currently assigned variables
    int NbNoAssigned;           // Number of currently unassigned variables
    int NbNoAssignedVal;        // Number of values not yet taken by any assigned variable (columns in the reduced cost matrix)
    vector<bool> isAssignedValue;       // isAssignedValue[v] = true if value v (index in UnionVarDomain) is already taken by an assigned variable
    vector<bool> varAlreadyProcessed;   // varAlreadyProcessed[i] = true if variable i has already been processed during domain filtering
    vector<Value> exceptedValues;       // Values exempt from the AllDifferent constraint (may be shared by multiple variables)
    vector<int> exceptedValIndex;       // Indices in UnionVarDomain of each excepted value
    vector<uint8_t> isExceptedVal;      // isExceptedVal[v] = 1 if value v (index in UnionVarDomain) is an excepted value
    bool excepted;              // True if at least one excepted value exists, false otherwise
    bool isSquare;              // True if the cost matrix is square (number of variables == number of values)
    vector<string> UnionVarDomain;      // Sorted union of all variable domains (unique value names)
    vector<int> VarDomainSize;          // Initial domain size for each variable (indexed by variable position in scope)
    unordered_map<string, int> mapDomainValToIndex; // Maps each value name to its index in UnionVarDomain
    bool SameDomain;            // True if all variables share the same domain
    vector<Cost> ReduceCostMatrix;      // Reduced cost matrix after LAP solving: ReduceCostMatrix[var * NbNoAssignedVal + val] = reduced cost of assigning value val to variable var
    vector<uint8_t> visited;    // visited[v] = 1 if variable/value node v has been finalized by Dijkstra
    vector<uint8_t> inHeap;     // inHeap[v] = 1 if variable/value node v is currently in the priority queue
    vector<Cost> distanceToVar; // Shortest distances from the Dijkstra source variable to each other variable node in the residual graph
    vector<Cost> distanceToVal; // Shortest distances from the Dijkstra source variable to each value node in the residual graph (reverse graph)
    int Q;                      // Number of Dijkstra runs to perform during reduced-cost filtering (derived from FiltLevel and NbNoAssigned)
    double FiltLevel;           // Fraction of variables used as Dijkstra sources during filtering (0 = disabled, 1 = all variables)
    unordered_set<int> trackingList;    // Set of variable/value nodes not yet finalized by Dijkstra, used by the bimodal traversal strategy
    Cost MaxReducedCost;        // Maximum reduced cost observed over all arcs in the current residual graph (used as filtering threshold)
    Cost ReducedCost;           // Temporary reduced cost of a single arc computed during cost projection
    Cost VarMaxReducedCost;     // Maximum reduced cost among all arcs incident to a given variable (used to guide source selection)
    int findConflict;           // Flag returned by the LAP solver: non-zero if a conflict (infeasible sub-problem) was detected, value encodes the number of conflicting variables
    int lower, upper;           // Lower and upper thresholds on NbNoAssigned controlling when reduced-cost filtering is activated (lower) and when Dijkstra-based filtering is applied (upper)


    void projectLB(Cost c)
    {
        if (c > MIN_COST) {
            //if (!isSquare)
               // lb += c;
            Constraint::projectLB(c);
        }
    }

    // returns true if the constraint can be projected to small local cost function in extension
    // bool canbeProjectedInExtension();

    // Depending of the value and the cost, extend or project the cost on the index value of the variable var
    void ExtOrProJ(int var, Value value, Cost C)
    {
        //int value_idx = scope[var]->toIndex(value);
        TreeDecomposition* td = wcsp->getTreeDec();
        if (C > MIN_COST) {
            if (!isSquare) {
                if (td && scope[var]->canbe(value))
                    td->addDelta(cluster, scope[var], value, -C);

            }
            assert(scope[var]->getCost(value) >= C);
            scope[var]->extend(value, C);
        } else if (C < MIN_COST) {
            if (!isSquare) {
                if (td && scope[var]->canbe(value))
                    td->addDelta(cluster, scope[var], value, -C);
            }
            scope[var]->project(value, -C, true);
        }
    }

public:
    // Constructor for the AllDifferentConstraint class
    AllDifferentConstraint(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in)
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
        , FiltLevel(0.0) 

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
                if (variable->isValueNames())
                    continue;
                for (int val_ind = 0; val_ind < (int)variable->getDomainInitSize(); val_ind++) {
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
            for (int valIndex = 0; valIndex < DomainSize; valIndex++) {
                // cout<<variable->getValueName(valIndex)<<" ";
                if (seen.insert(variable->getValueName(valIndex)).second) {
                    UnionVarDomain.push_back(variable->getValueName(valIndex));
                    NbValues++;
                }
            }
            // cout<<endl;

            // Repeat for the rest of the variables
            for (int varIndex = 1; varIndex < arity_in; varIndex++) {
                conflictWeights.push_back(0);
                auto* variable = scope[varIndex];
                DomainSize = (int)variable->getDomainInitSize();
                // cout<<"variable : x"<<varIndex<<" domain size : "<<DomainSize<<endl;
                VarDomainSize.push_back(DomainSize);
                if (VarDomainSize[varIndex - 1] != VarDomainSize[varIndex])
                    SameDomain = false;

                // Add only new unique values to the union
                for (int valIndex = 0; valIndex < DomainSize; valIndex++) {
                    // cout<<variable->getValueName(valIndex)<<" ";
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
                deltaCosts.assign(NbValues, MIN_COST);
            }

            int rate = ToulBar2::ReducedCostsFiltering;
        
            if (rate > 0){
                FiltLevel = rate / 100.0;
            }
           

            // Test value symmetries
            //            for (unsigned int a = 0; a < NbValues; ++a) {
            //                for (unsigned int b = a+1; b < NbValues; ++b) {
            //                    if (valueSymmetry(UnionVarDomain[a], UnionVarDomain[b])) {
            //                        if (ToulBar2::verbose >= 1) {
            //                            cout << "detect value symmetry between " <<  UnionVarDomain[a] << " and " << UnionVarDomain[b] << endl;
            //                        }
            //                    }
            //                }
            //            }

            // Initialize

           // lower = static_cast<int>(arity_in * 0.5);
            upper = static_cast<int>(arity_in * 0.8);
            storeLastAssignment = vector<StoreValue>(arity_in, StoreValue(WRONG_VAL));
            NoAssignedVar = vector<int>(arity_in, -1);
            isExceptedVal = vector<uint8_t>(NbValues, 0);
            ValuesCapacity =  vector<int>(NbValues, 1);
            AssignedVar = vector<int>(arity_in, -1);
            AssignedVal = vector<int>(arity_in, -1);
            costMatrix = vector<Cost>((arity_+1) * NbValues, MAX_COST);
            ReduceCostMatrix = vector<Cost>(arity_ * NbValues, MAX_COST);

            // Allocate memory
            rowSol = new int[arity_+1];
            ReduceCostRow = new Cost[arity_+1];
            ReduceCostCol = new Cost[NbValues];
        } else {
            deconnect();
        }
    }

    virtual ~AllDifferentConstraint()
    {
        delete[] rowSol;
        delete[] ReduceCostRow;
        delete[] ReduceCostCol;
    }
    void read(istream& files) // TODO: add a parameter for controlling the level of propagation if necessary
    {
        string line;
        getline(files, line);
        stringstream file(line);
        int nbExcepted = 0;
        file >> nbExcepted;
        excepted = nbExcepted > 0;
        for (int v = 0; v < nbExcepted; v++) {
            Value except;
            file >> except;
            exceptedValues.push_back(except);
            for (int varIndex = 0; varIndex < arity_; ++varIndex) {
                auto* variable = scope[varIndex];
                if (variable->canbe(except)) {
                    exceptedValIndex.push_back(mapDomainValToIndex[variable->getValueName(variable->toIndex(except))]);
                    isExceptedVal[mapDomainValToIndex[variable->getValueName(variable->toIndex(except))]] = 1;
                    ValuesCapacity[mapDomainValToIndex[variable->getValueName(variable->toIndex(except))]] = arity_;       
                    break;
                }
            }
        }
        if (excepted) {
            isSquare = false;
            deltaCosts.assign(NbValues, MIN_COST);
            /*for (int varIndex = 0; varIndex < arity_; varIndex++) {
                deltaCosts.emplace_back(NbValues, MIN_COST);
            }*/
        } else {
            // contradiction
            if (NbValues < arity_)
                THROWCONTRADICTION;
        }
        
        if(!isSquare){
            int nbDelta;
            file >> nbDelta;
            if(nbDelta){
                for (int v = 0; v < nbDelta; v++) {
                    Value value;
                    Cost delta;
                    file >> value;
                    file >> delta;
                    for (int varIndex = 0; varIndex < arity_; ++varIndex) {
                        auto* variable = scope[varIndex];
                        if (variable->canbe(value)) {
                            deltaCosts[mapDomainValToIndex[variable->getValueName(variable->toIndex(value))]] = delta;
                            break;
                        }
                    }
                }
            }
        }
    }

    bool extension() const FINAL { return false; } // TODO: allows functional variable elimination but no other preprocessing
    bool isAllDiff() const FINAL { return true; }
    bool isAllDiffSquare() const FINAL { return isSquare; }

    vector<Value> getExceptedValues() const FINAL { return exceptedValues; }

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

                if (arity_ < (int)wcsp->numberOfVariables()) { // SdG: no need to increase weightedDegree of all problem variables
                    Constraint::incConflictWeight(1);
                }
            } else {

                for (int i = 0; i < arity_; i++) {
                //for (int i = 0; i < numConflictVars; i++) {
                    auto* variable = scope[i];

                    if (!variable->unassigned() ){
                    //if (connected(i)) {
                        //conflictWeights[lastConflictVars[i]]++;
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
    // bool universal(Cost zero = MIN_COST) override;

    bool implies(Constraint* ctr) FINAL
    {
        if (scopeIncluded(ctr) && ctr->isBinary()) {
            BinaryConstraint* bctr = (BinaryConstraint*)ctr;
            EnumeratedVariable* x = (EnumeratedVariable*)(bctr->getVar(0));
            EnumeratedVariable* y = (EnumeratedVariable*)(bctr->getVar(1));
            if (x->isValueNames() && y->isValueNames()) {
                for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
                    for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
                        if ((x->getValueName(x->toIndex(*iterx)) != y->getValueName(y->toIndex(*itery)) || (excepted && (find(exceptedValues.begin(), exceptedValues.end(), *iterx) != exceptedValues.end() || find(exceptedValues.begin(), exceptedValues.end(), *itery) != exceptedValues.end()))) && bctr->getCost(*iterx, *itery) > MIN_COST) {
                            return false;
                        }
                    }
                }
            } else {
                for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
                    for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
                        if ((*iterx != *itery || (excepted && find(exceptedValues.begin(), exceptedValues.end(), *iterx) != exceptedValues.end())) && bctr->getCost(*iterx, *itery) > MIN_COST) {
                            return false;
                        }
                    }
                }
            }
            return true;
        }
        return false;
    }

    void projects(Constraint* ctr) FINAL
    {
        if (scopeIncluded(ctr) && ctr->isBinary()) {
            BinaryConstraint* bctr = (BinaryConstraint*)ctr;
            EnumeratedVariable* x = (EnumeratedVariable*)(bctr->getVar(0));
            EnumeratedVariable* y = (EnumeratedVariable*)(bctr->getVar(1));
            Cost mult_ub = (wcsp->getUb() < (MAX_COST / MEDIUM_COST)) ? (max(LARGE_COST, wcsp->getUb() * MEDIUM_COST)) : wcsp->getUb();
            if (x->isValueNames() && y->isValueNames()) {
                for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
                    string s = x->getValueName(x->toIndex(*iterx));
                    unsigned int yindex = y->toIndex(s);
                    Value yval = y->toValue(yindex);
                    if ((!excepted || (find(exceptedValues.begin(), exceptedValues.end(), *iterx) == exceptedValues.end() && find(exceptedValues.begin(), exceptedValues.end(), yval) == exceptedValues.end())) && y->canbe(yval) && !CUT(bctr->getCost(*iterx, yval), wcsp->getUb())) {
                        bctr->addcost(*iterx, yval, mult_ub - bctr->getCost(*iterx, yval));
                    }
                }
            } else {
                for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
                    if ((!excepted || find(exceptedValues.begin(), exceptedValues.end(), *iterx) == exceptedValues.end()) && y->canbe(*iterx) && !CUT(bctr->getCost(*iterx, *iterx), wcsp->getUb())) {
                        bctr->addcost(*iterx, *iterx, mult_ub - bctr->getCost(*iterx, *iterx));
                    }
                }
            }
        }
    }

    bool valueSymmetry(string stra, string strb)
    {
        if (!isSquare) {
            return false;
        }
        for (int j = 0; j < arity_; j++) {
            EnumeratedVariable* xj = (EnumeratedVariable*)getVar(j);
            unsigned int ida = xj->toIndex(stra);
            Value xja = xj->toValue(ida);
            unsigned int idb = xj->toIndex(strb);
            Value xjb = xj->toValue(idb);
            if (xj->cannotbe(xja) || xj->cannotbe(xjb) || xj->getCost(xja) != xj->getCost(xjb)) {
                return false;
            }
            if (xj->unassigned() && xj->getDegree() > 1) {
                ConstraintList* constrsj = xj->getConstrs();
                for (ConstraintList::iterator it = constrsj->begin(); it != constrsj->end(); ++it) {
                    Constraint* ctr = (*it).constr;
                    if (ctr->isBinary() && !ctr->isSep() && scopeIncluded(ctr)) {
                        BinaryConstraint* cjk = (BinaryConstraint*)ctr;
                        EnumeratedVariable* xk = (EnumeratedVariable*)((cjk->getVar(0) == xj) ? cjk->getVar(1) : cjk->getVar(0));
                        unsigned int ida = xk->toIndex(stra);
                        Value xka = xk->toValue(ida);
                        unsigned int idb = xk->toIndex(strb);
                        Value xkb = xk->toValue(idb);
                        for (EnumeratedVariable::iterator iterk = xk->begin(); iterk != xk->end(); ++iterk) {
                            if (*iterk != xka && *iterk != xkb && cjk->getCost(xj, xk, xja, *iterk) != cjk->getCost(xj, xk, xjb, *iterk)) {
                                return false;
                            }
                        }
                    } else if (!ctr->isAllDiff()) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    Cost eval(const Tuple& s) override
    {
        // returns the cost of the corresponding assignment s

        if (isSquare) {
            Cost res = MIN_COST;
            Cost nbsame = 0;
            
            vector<uint8_t> alreadyUsed(arity_, 0);
            for (int varIndex = 0; varIndex < arity_; varIndex++) {
                auto val = s[varIndex];

                if (alreadyUsed[val]) {
                    nbsame += 1;

                } else {
                    alreadyUsed[val] = 1;
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

        //Cost res = -lb + assigneddeltas;
        Cost res = 0;
        Cost nbsame = 0;
        vector<int> alreadyUsed(NbValues, 0);
        for (int varIndex = 0; varIndex < arity_; varIndex++) {
            int valIndex = mapDomainValToIndex[scope[varIndex]->getValueName(s[varIndex])];

            if (alreadyUsed[valIndex]) {
                int it = excepted == 0 ? 0 : isExceptedVal[valIndex];
                if(it == 0){ 
                    nbsame += 1;
                }
                
            } 
            alreadyUsed[valIndex] += 1 ;
        
        }
        if(nbsame == 0){
            if(excepted){
                for (int valIndex = 0; valIndex < NbValues; valIndex++) {
                    res += (ValuesCapacity[valIndex] - alreadyUsed[valIndex]) * deltaCosts[valIndex];
                }
            }
            else{
                for (int valIndex = 0; valIndex < NbValues; valIndex++) {
                    if(alreadyUsed[valIndex] == 0)
                         res += deltaCosts[valIndex];
                 }
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
        if (isSquare && !excepted)
            return MIN_COST;
        Cost sumdelta = -lb + assigneddeltas;

        Cost m = *max_element(deltaCosts.begin(), deltaCosts.end());
        if (m > MIN_COST)
                sumdelta += m;
  
        if (CUT(sumdelta, wcsp->getUb()))
            return MAX_COST;
        else
            return sumdelta;
    }

    Cost getCost(int index, Value val)
    {
        assert(index >= 0 && index < arity_);
        return deltaCosts[scope[index]->toIndex(val)];
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

    /**
     * @brief Remove assigned values from unassigned
     *        variable domains, detecting contradictions, and triggering n-ary
     *        projection when thresholds are met.
     * @return true if propagation succeeded without n-ary projection, false otherwise.
     * @throws THROWCONTRADICTION if two variables share the same assigned value.
     */
    bool RemoveAssignVar()
    {
        // Reset tracking vectors for this propagation pass
        isAssignedValue = vector<bool>(NbValues, false);
        varAlreadyProcessed = vector<bool>(arity_, false);

        AssignedVar.clear();
        AssignedVar.reserve(arity_);
        AssignedVal.clear();
        AssignedVal.reserve(arity_);
        NoAssignedVar.clear();
        NoAssignedVar.reserve(arity_);
        unordered_map<int, int> conflictvar; // Maps value index -> variable index for contradiction detection

        // Partition variables into assigned / unassigned
        for (int varIndex = 0; varIndex < arity_; ++varIndex) {
            auto* variable = scope[varIndex];

            if (variable->unassigned()) {
                NoAssignedVar.push_back(varIndex);
            } else {
                int valIndex = mapDomainValToIndex[variable->getValueName(variable->toIndex(variable->getValue()))];
                int it = excepted == 0 ? 0 : isExceptedVal[valIndex];
                if (!isAssignedValue[valIndex]) {
                    AssignedVar.push_back(varIndex);
                    if (it == 0) {
                        isAssignedValue[valIndex] = true;   // Mark value as taken
                        AssignedVal.push_back(valIndex);
                        conflictvar[valIndex] = varIndex;   // Record owner for conflict detection
                    }
                } else {
                    // Two variables share the same non-excepted value: contradiction
                    if (it == 0){
                        int varIndex2 = conflictvar[valIndex];
                        conflictWeights[varIndex2]++;
                        conflictWeights[varIndex]++;
                       
                        THROWCONTRADICTION;
                    }
                }
            }
        }

        NbAssigned = AssignedVar.size();
        bool VarAssigned;
        bool NaryPro = false;

        // Iteratively remove assigned values from unassigned variable domains
        do {
            VarAssigned = false;
            vector<int> newlyAssignedVars; // Variables that become assigned during this iteration

            for (int valIndex = 0; valIndex < NbValues; ++valIndex) {
                if (!isAssignedValue[valIndex])
                    continue; // Only process taken values

                for (int varIndex : NoAssignedVar) {
                    if (varAlreadyProcessed[varIndex])
                        continue;

                    auto* variable = scope[varIndex];
                    Value value = variable->toValue(variable->toIndex(UnionVarDomain[valIndex]));

                    // Remove value if present in domain and not excepted
                    if (variable->canbe(value) && (find(exceptedValIndex.begin(), exceptedValIndex.end(), valIndex) == exceptedValIndex.end())) {
                        variable->remove(value);

                        if (variable->assigned()) {
                            // Variable became singleton: schedule for propagation
                            newlyAssignedVars.push_back(varIndex);
                            varAlreadyProcessed[varIndex] = true;
                        }
                    }
                }
            }

            // Integrate newly assigned variables into the constraint state
            for (int varIndex : newlyAssignedVars) {
                if (connected(varIndex)) {
                    deconnect(varIndex);
                    NbNoAssigned--;
                    NbAssigned++;

                    // Switch to full n-ary projection if size thresholds are reached
                    if (getNonAssigned() <= NARYPROJECTIONSIZE && (getNonAssigned() <= 1 || prodInitDomSize <= NARYPROJECTIONPRODDOMSIZE || maxInitDomSize <= NARYPROJECTION3MAXDOMSIZE || (getNonAssigned() == 2 && maxInitDomSize <= NARYPROJECTION2MAXDOMSIZE))) {
                        deconnect();
                        projectNary();
                        NaryPro = true;
                        break;
                    } else {
                        VarAssigned = true; // Signal another propagation iteration
                        auto* variable = scope[varIndex];
                        int valIndex = mapDomainValToIndex[variable->getValueName(variable->toIndex(variable->getValue()))];
                        int it = excepted == 0 ? 0 : isExceptedVal[valIndex];
                        if (it == 0) {
                            isAssignedValue[valIndex] = true;
                            AssignedVal.push_back(valIndex);
                        }
                        AssignedVar.push_back(varIndex);
                    }
                } else {
                    return false; // Variable disconnected externally: abort
                }
            }

            // Rebuild unassigned variable list after each iteration
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

 

        
    /**
     * @brief Bimodal Dijkstra on the forward residual graph.
     *        Computes shortest distances from a source variable to all other
     *        variable nodes. At each settled node, chooses between iterating
     *        over the (sparse) adjacency list or the (small) unsettled set,
     *        whichever is cheaper.
     * @param source    Index of the source variable node.
     * @param VarList   VarList[val] = list of variable indices that have val as a non-matching arc.
     * @param ValList   ValList[var][val] = 1 if arc (var,val) exists in the residual graph.
     * 
     *      Source : Bimodal Depth-First Search for Scalable GAC for AllDifferent.
     *      Sulian Le Bozec Chiffoleau; Nicolas Beldiceanu; Charles Prud'homme; 
     *      Gilles Simonin; and Xavier Lorca, roceedings of the Thirty-Fourth 
     *      International Joint Conference on Artificial Intelligence, IJCAI 2025, 
     *      Montreal, Canada, August 16-22, 2025. 
     */
    void BimodalDijkstra(int source, vector<vector<int>>& VarList, vector<vector<uint8_t>>& ValList)
    {
        vector<PQ::handle_type> handles(NbNoAssigned);
        PQ pq;

        distanceToVar.assign(NbNoAssigned, MAX_COST); // Initialize all distances to infinity
        visited.assign(NbNoAssigned, 0);

        // trackingList holds unsettled variable nodes for the dense-side traversal
        trackingList.clear();
        trackingList.reserve(NbNoAssigned);
        for (int var = 0; var < NbNoAssigned; ++var)
            trackingList.insert(var);

        inHeap.assign(NbNoAssigned, 0);
        distanceToVar[source] = 0;
        handles[source] = pq.push({0, source});
        inHeap[source] = 1;

        while (!pq.empty()) {
            auto [d, var] = pq.top();
            pq.pop();
            visited[var] = 1;
            trackingList.erase(var); // Remove from unsettled set

            int val = rowSol[var];           // Matching arc: var -> its assigned value
            auto& varlist = VarList[val];    // Variables reachable via value val (non-matching arcs)

            // Bimodal choice: sparse side (adjacency list) vs dense side (unsettled set)
            if ((int)varlist.size() < (int)trackingList.size()) {
                // Sparse: iterate over adjacency list
                for (int nextVar : varlist) {
                    if (visited[nextVar]) continue;
                    Cost alt = d + ReduceCostMatrix[nextVar * NbNoAssignedVal + val];
                    if (alt < distanceToVar[nextVar]) {
                        distanceToVar[nextVar] = alt;
                        if (!inHeap[nextVar]) {
                            handles[nextVar] = pq.push({alt, nextVar});
                            inHeap[nextVar] = 1;
                        } else {
                            pq.decrease(handles[nextVar], {alt, nextVar});
                        }
                        if (alt == d) { // Zero-cost arc: settle immediately
                            trackingList.erase(nextVar);
                            visited[nextVar] = 1;
                        }
                    }
                }
            } else {
                // Dense: iterate over unsettled set and probe ValList
                auto it = trackingList.begin();
                while (it != trackingList.end()) {
                    auto nextVar = *it;
                    if (ValList[nextVar][val]) { // Arc exists in residual graph
                        Cost alt = d + ReduceCostMatrix[nextVar * NbNoAssignedVal + val];
                        if (alt < distanceToVar[nextVar]) {
                            distanceToVar[nextVar] = alt;
                            if (!inHeap[nextVar]) {
                                handles[nextVar] = pq.push({alt, nextVar});
                                inHeap[nextVar] = 1;
                            } else {
                                pq.decrease(handles[nextVar], {alt, nextVar});
                            }
                            if (alt == d) { // Zero-cost arc: settle immediately
                                it = trackingList.erase(it);
                                visited[nextVar] = 1;
                                continue;
                            }
                        }
                    }
                    ++it;
                }
            }
        }
    }


    /**
     * @brief Bimodal Dijkstra on the reverse residual graph.
     *        Computes shortest distances from the value matched to the source
     *        variable back to all other value nodes, used by the Régin landmark
     *        upper-bound test.
     * @param source  Index of the source variable (its matched value is the start node).
     * @param ValList ValList[var][val] = 1 if arc (var,val) exists in the residual graph.
     */
    void BimodalDijkstraReverseGraph(int source, vector<vector<uint8_t>>& ValList)
    {
        vector<PQ::handle_type> handles(NbNoAssignedVal);
        PQ pq;

        distanceToVal.assign(NbNoAssignedVal, 0); // Initialize all value distances to infinity
        visited.assign(NbNoAssignedVal, 0);
        vector<int> colSol(NbNoAssignedVal); // Inverse of rowSol: colSol[val] = var matched to val

        // Build inverse matching and initialize unsettled value set
        trackingList.clear();
        trackingList.reserve(NbNoAssigned);
        for (int var = 0; var < NbNoAssigned; ++var) {
            int val = rowSol[var];
            colSol[val] = var;
            trackingList.insert(val); // Only matched values are nodes in the reverse graph
            distanceToVal[val] = MAX_COST;
        }

        int val = rowSol[source]; // Start from the value matched to source
        inHeap.assign(NbNoAssignedVal, 0);
        distanceToVal[val] = 0;
        handles[val] = pq.push({0, val});
        inHeap[val] = 1;

        while (!pq.empty()) {
            auto [d, val] = pq.top();
            pq.pop();
            visited[val] = 1;
            trackingList.erase(val);

            int var = colSol[val]; // Variable matched to current value node

            // Dense traversal over unsettled value nodes
            auto it = trackingList.begin();
            while (it != trackingList.end()) {
                auto nextVal = *it;
                if (ValList[var][nextVal]) { // Reverse arc exists
                    Cost alt = d + ReduceCostMatrix[var * NbNoAssignedVal + nextVal];
                    if (alt < distanceToVal[nextVal]) {
                        distanceToVal[nextVal] = alt;
                        if (!inHeap[nextVal]) {
                            handles[nextVal] = pq.push({alt, nextVal});
                            inHeap[nextVal] = 1;
                        } else {
                            pq.decrease(handles[nextVal], {alt, nextVal});
                        }
                        if (alt == d) { // Zero-cost arc: settle immediately
                            it = trackingList.erase(it);
                            visited[nextVal] = 1;
                            continue;
                        }
                    }
                }
                ++it;
            }
        }
    }

    /**
     * @brief Standard Dijkstra on the forward residual graph (sparse variant).
     *        Used when the bimodal strategy is not needed (e.g. dense graphs
     *        where VarList is always smaller than trackingList).
     * @param source  Index of the source variable node.
     * @param VarList VarList[val] = list of variable indices reachable via value val.
     */
    void Dijkstra(int source, vector<vector<int>>& VarList)
    {
        vector<PQ::handle_type> handles(NbNoAssigned);
        PQ pq;

        distanceToVar.assign(NbNoAssigned, MAX_COST);
        visited.assign(NbNoAssigned, 0);
        inHeap.assign(NbNoAssigned, 0);
        distanceToVar[source] = 0;
        handles[source] = pq.push({0, source});
        inHeap[source] = 1;

        while (!pq.empty()) {
            auto [d, var] = pq.top();
            pq.pop();
            visited[var] = 1;

            int val = rowSol[var];        // Matching arc: traverse to matched value
            auto& varlist = VarList[val]; // Variables reachable from that value (non-matching arcs)

            for (int nextVar : varlist) {
                if (visited[nextVar]) continue;
                Cost alt = d + ReduceCostMatrix[nextVar * NbNoAssignedVal + val];
                if (alt < distanceToVar[nextVar]) {
                    distanceToVar[nextVar] = alt;
                    if (!inHeap[nextVar]) {
                        handles[nextVar] = pq.push({alt, nextVar});
                        inHeap[nextVar] = 1;
                    } else {
                        pq.decrease(handles[nextVar], {alt, nextVar});
                    }
                    if (alt == d) // Zero-cost arc: settle immediately without re-queuing
                        visited[nextVar] = 1;
                }
            }
        }
    }

    /**
     * @brief Régin landmark upper-bound test.
     *        Runs Bimodal Dijkstra in both directions from the source variable,
     *        then checks whether the sum of the two maximum distances is within
     *        the remaining cost budget H. If so, no arc in the residual graph
     *        can be pruned and filtering can stop early.
     * @param source  Source variable index used as the landmark.
     * @param VarList VarList[val] = variables reachable via value val (forward graph).
     * @param ValList ValList[var][val] = 1 if arc (var,val) exists in residual graph.
     * @param H       Remaining cost budget (wcsp upper bound minus current lower bound).
     * @return true if all values are consistent (no pruning possible), false otherwise.
     */
    bool ReginLandmarkUpperBound(int source, vector<vector<int>>& VarList, vector<vector<uint8_t>>& ValList, Cost H)
    {
        // Forward pass: distances from source to all variable nodes
        BimodalDijkstra(source, VarList, ValList);
        // Backward pass: distances from source's matched value to all value nodes
        BimodalDijkstraReverseGraph(source, ValList);

        Cost maxDistanceFromVar = *max_element(distanceToVar.begin(), distanceToVar.end());
        Cost maxDistanceToVar   = *max_element(distanceToVal.begin(), distanceToVal.end());

        // If max_forward + max_backward <= H - MaxReducedCost, no arc exceeds the budget
        bool AllValConsistent = maxDistanceFromVar + maxDistanceToVar <= H - MaxReducedCost;
        return AllValConsistent;
    }

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
    void propagate() override
    {

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
                    b=true;
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
                        //    cout<<"done\n";
                            vector<int> NoAssignedVal;
                            for (int valIndex = 0; valIndex < NbValues; ++valIndex) {
                                if (!isAssignedValue[valIndex]) {
                                    NoAssignedVal.push_back(valIndex);
                                }
                            }

                            NbNoAssignedVal = NoAssignedVal.size();

                            // Initialize cost matrix for the Jonker algorithm
                            Cost current_ub = wcsp->getUb() - wcsp->getLb();

                            for (int varInd = 0; varInd < NbNoAssigned; ++varInd) {
                                int varIndex = NoAssignedVar[varInd];
                                auto* variable = scope[varIndex];

                                for (int valInd = 0; valInd < NbNoAssignedVal; ++valInd) {
                                    int valIndex = NoAssignedVal[valInd];
                                    Value val = variable->toValue(valIndex);
                                    if(variable->canbe(val)){
                                        auto valCost = variable->getCost(val);
                                        costMatrix[varInd * NbNoAssignedVal + valInd] =  valCost < current_ub ? valCost : current_ub;
                                    }
                                    else{
                                        costMatrix[varInd * NbNoAssignedVal + valInd] = current_ub;
                                    }
                                }
                            }

                            // Solve assignment problem using Jonker algorithm
                            Cost TotalCost;
                            if(!isSquare){
                                for (int valInd = 0; valInd < NbNoAssignedVal; ++valInd) {
                                    costMatrix[NbNoAssigned * NbNoAssignedVal + valInd] = deltaCosts[NoAssignedVal[valInd]];
                                }
                                TotalCost = lapjv(NbNoAssigned+1, NbNoAssignedVal, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub, findConflict);

                            }
                            else{

                            
                              TotalCost = lapjv(NbNoAssigned, NbNoAssignedVal, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub, findConflict);
                            }
                            if (TotalCost >= current_ub) {
                                
                                if(findConflict){
                                    int varIndex;
                                    for(int var=0; var< findConflict; var++){
                                        if(rowSol[var] >= NbNoAssigned) continue;
                                        varIndex = NoAssignedVar[rowSol[var]];
                                        conflictWeights[varIndex]++;
                                    } 
                                } 
                                wcsp->revise(this);
                                THROWCONTRADICTION;

                                
                            } else if (TotalCost >= 0) {
                                Cost jonker = TotalCost;

                                // Store results from Jonker algorithm
                                for (int varIndex = 0; varIndex < NbNoAssigned; ++varIndex) {
                                    auto* variable = scope[NoAssignedVar[varIndex]];
                                    storeLastAssignment[NoAssignedVar[varIndex]] = variable->toValue(variable->toIndex(UnionVarDomain[NoAssignedVal[rowSol[varIndex]]]));
                                }

                                for (int varIndex = 0; varIndex < NbAssigned; ++varIndex) {
                                    auto* variable = scope[AssignedVar[varIndex]];
                                    storeLastAssignment[AssignedVar[varIndex]] = variable->toValue(variable->toIndex(UnionVarDomain[AssignedVal[varIndex]]));
                                }

                                storeAssignment = true;
                                                              
                               if(!isSquare){
                                 
                                    Cost mindelta  = ReduceCostRow[NbNoAssigned];
                                    int valIndex;
                                    for (int valInd = 0; valInd < NbNoAssignedVal; ++valInd) {
                                        valIndex =  NoAssignedVal[valInd]; 
                                        deltaCosts[valIndex] -= (ReduceCostCol[valInd] + mindelta);
                                    }
                                    if (mindelta > 0){

                                        jonker += (mindelta* (NbValues - arity_)) - mindelta;
                                        
                                        if(jonker >= current_ub){
                                            wcsp->revise(this);
                                            THROWCONTRADICTION;
                                        }
                                    }

                                }

                                // Update the lower bound with the Jonker algorithm's total cost
                                projectLB(jonker);
                                Cost newCurrentUb = current_ub - jonker;

                                // Adjust unary costs for unassigned variables based on reduced costs

                                vector<vector<int>> VarList(NbNoAssignedVal);
                                vector<vector<uint8_t>> ValList(NbNoAssigned);
                                MaxReducedCost = 0;

                                for (int varInd = 0; varInd < NbNoAssigned; ++varInd) {
                                    int varIndex = NoAssignedVar[varInd];
                                    auto* variable = scope[varIndex];
                                    ValList[varInd] = vector<uint8_t>(NbNoAssignedVal, 0);
                                    for (int valInd = 0; valInd < NbNoAssignedVal; ++valInd) {
                                        if(costMatrix[varInd * NbNoAssignedVal + valInd] < current_ub ){
                                            int valIndex = NoAssignedVal[valInd];
                                            Value value = variable->toValue(valIndex);  
                                            ExtOrProJ(varIndex, value, (ReduceCostRow[varInd] + ReduceCostCol[valInd]));

                                            if (FiltLevel > 0 && (NbNoAssigned <= upper)) {
                                                ReducedCost = costMatrix[varInd * NbNoAssignedVal + valInd] - (ReduceCostRow[varInd] + ReduceCostCol[valInd]);
                                                if(ReducedCost < newCurrentUb){
                                                    ReduceCostMatrix[varInd * NbNoAssignedVal + valInd] = ReducedCost;
                                                    if(ReducedCost > MaxReducedCost) MaxReducedCost = ReducedCost;
                                                    if (rowSol[varInd] != valInd){
                                                        VarList[valInd].push_back(varInd);
                                                        ValList[varInd][valInd] = 1;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }



                                // Update support values if needed for unassigned variables
                                for (int varInd = 0; varInd < NbNoAssigned; ++varInd) {
                                    int varIndex = NoAssignedVar[varInd];
                                    auto* variable = scope[varIndex];
                                    Value optimalValue = storeLastAssignment[varIndex];
                                    if (variable->getSupport() != optimalValue) {
                                        if (ToulBar2::verbose > 0)
                                            cout << "CHANGE ALLDIFF SUPPORT " << variable->getName() << " from " << variable->getSupport() << " to " << optimalValue << endl;
#ifndef NDEBUG
                                        variable->queueEAC1(); // EAC support may have been lost
#endif
                                        variable->setSupport(optimalValue);
                                        assert(variable->getCost(variable->getSupport()) == MIN_COST);
                                    }
                                }


                                   /* (BEGIN) : Filtering of variables domains with Sellmann or Cambazard method 

                                   Source : Claus, G.; Cambazard, H.; and Jost, V. 2020. Analysis of Re-
                                            duced Costs Filtering for Alldifferent and Minimum Weight
                                            Alldifferent Global Constraints. In Giacomo, G. D.; Catal´a,
                                            A.; Dilkina, B.; Milano, M.; Barro, S.; Bugar´ın, A.; and
                                            Lang, J., eds., ECAI 2020 - 24th European Conference on
                                            Artificial Intelligence, volume 325, 323–330. Santiago de
                                            Compostela, Spain: IOS Press.  */

                                 if (FiltLevel > 0 && (NbNoAssigned <= upper) && (MaxReducedCost >= newCurrentUb/NbNoAssigned) ) {
                                    
                                    int position; 
                                    int varInde; 
                                    vector<int> VariableList(NbNoAssigned);     
                                    iota(VariableList.begin(), VariableList.end(), 0);                           
                                    
                                    if (FiltLevel >= 1. - (double)ToulBar2::epsilon) {
                                        Q = NbNoAssigned;    
                                    } else {
                                        Q = min(NbNoAssigned, 1 + static_cast<int>(NbNoAssigned * FiltLevel));  
                                        myrearrange(VariableList); 
                                    }
                                    
                                    int comp = 0;
                                    int numzeroremoved = 0;
                                    uint8_t valremoved;

                                    int stop = Q/2;
                                    bool AllValConsistent = false;

                                    for (int ind = 0; ind < Q; ++ind) {
                                        if(numzeroremoved >= stop){
                                            break;
                                        }
                                        varInde = VariableList[ind];
                                        numzeroremoved++;
                                        
                                        if (comp == 0){
                                        // Apply Régin Upper Bounds of Shortest Paths algorithm with landmark 
                                        

                                            AllValConsistent = ReginLandmarkUpperBound(varInde , VarList, ValList, newCurrentUb);
                                        }
                                        else{

                                            //Bimodal Dijkstra’s shortest path algorithm 
                                            BimodalDijkstra(varInde, VarList, ValList ); 
                                            
                                        }
                                        if(AllValConsistent){
                                            
                                            break;
                                        }
                                        comp++;
                                        int valInd;
                                        for (int row = 0; row < NbNoAssigned; ++row) {
                                            if (distanceToVar[row] >= MAX_COST)
                                                continue;
                                            valremoved = 0;
                                            valInd = rowSol[row];
                                            position = -1;
                                            auto& varlist = VarList[valInd];
                                            for (int varInd : varlist) {
                                                position++;
                                                if (distanceToVar[varInd] >= MAX_COST)
                                                    continue;
                                                Cost reducedCost = ReduceCostMatrix[varInd * NbNoAssignedVal + valInd] - distanceToVar[varInd] + distanceToVar[row];

                                                if (reducedCost >= newCurrentUb) {
                                                    int varIndex = NoAssignedVar[varInd];
                                                    int valIndex = NoAssignedVal[valInd];
                                                    auto* variable = scope[varIndex];
                                                    varlist[position] = -1;
                                                    ValList[varInd][valInd] = 0;
                                                    valremoved = 1;

                                                    Value value = variable->toValue(variable->toIndex(UnionVarDomain[valIndex]));

                                                    if (variable->canbe(value)) {
                                                        if (ToulBar2::verbose > 0)
                                                            cout << "REMOVE VALUE " << value << " from " << variable->getName() << endl;
                                                        ExtOrProJ(varIndex, value, -newCurrentUb); //SdG: project infinite cost on this value and avoid to skip and reenter the AllDiff constraint without finishing the current filtering
                                                    }
                                                }
                                            }
                                            if(valremoved){
                                                varlist.erase(std::remove(varlist.begin(), varlist.end(), -1), varlist.end());
                                                numzeroremoved = 0;
                                            }

                                        }

                                    }
                                }
                                /* (END) : Filtering of variables domains with Sellmann or Cambazard method */
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
                                
                                if(variable->canbe(value)){
                                    auto valCost = variable->getCost(value);
                                    costMatrix[varIndex * NbValues + mapDomainValToIndex[valName]] =  valCost < current_ub ? valCost : current_ub;
                                }
                                else{
                                    costMatrix[varIndex * NbValues + mapDomainValToIndex[valName]] = current_ub;
                                }
                            }
                        }

                        // Solve the Linear Assignment Problem (LAP) using the Jonker algorithm
                        Cost TotalCost;

                        if(!isSquare){
                            for (int valInd = 0; valInd < NbValues; ++valInd) {
                                costMatrix[arity_ * NbValues + valInd] = deltaCosts[valInd];
                            }
                            if (excepted) {
                                TotalCost = lapjv(arity_+1, NbValues, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub, exceptedValIndex, findConflict);
                            }
                            else{
                                TotalCost = lapjv(arity_ +1, NbValues, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub, findConflict);
                            }
                        
                        }
                        else {
                            TotalCost = lapjv(arity_, NbValues, costMatrix, rowSol, ReduceCostRow, ReduceCostCol, current_ub, findConflict);
                        
                        }
                        
                        if (TotalCost >= current_ub) {
                           
                               if(findConflict){
                                    int varIndex; 
                                    for(int var=0; var < findConflict; var++){
                                        varIndex = rowSol[var];
                                        if(varIndex >= arity_) continue;
                                        conflictWeights[varIndex]++;
                                    } 
                                } 
                                wcsp->revise(this);
                                THROWCONTRADICTION;
                                
                            }else if (TotalCost >= 0) {
                            // Store results from Jonker algorithm
                            for (int varIndex = 0; varIndex < arity_; varIndex++) {
                                auto* variable = scope[varIndex];
                                storeLastAssignment[varIndex] = variable->toValue(variable->toIndex(UnionVarDomain[rowSol[varIndex]]));
                            }

                            storeAssignment = true;

                            if(!isSquare){
                                Cost mindelta = ReduceCostRow[arity_];
                                for (int valInd = 0; valInd < NbValues; ++valInd) {
                                    deltaCosts[valInd] -= (ReduceCostCol[valInd] + mindelta);
                                }
                                
                                if (mindelta > 0){
                                    int nbexcep = exceptedValIndex.size();
                                    if(excepted){
                                        TotalCost += (mindelta *(NbValues - arity_ - nbexcep + (nbexcep*arity_)) - mindelta);

                                    }
                                    else{
                                        TotalCost += (mindelta *(NbValues - arity_) - mindelta);
                                    }
                                    if(TotalCost >= current_ub){
                                         wcsp->revise(this);
                                         THROWCONTRADICTION;
                                    }
                                }
                            }

                            // Update lower bound with the total cost found
                            projectLB(TotalCost);

                            // Adjust unary costs using reduced row and column costs
                            for (int varIndex = 0; varIndex < arity_; ++varIndex) {
                                auto* variable = scope[varIndex];

                                for (int valIndex = 0; valIndex < VarDomainSize[varIndex]; ++valIndex) {
                                    string valName = variable->getValueName(valIndex);
                                    if(costMatrix[varIndex * NbValues + mapDomainValToIndex[valName]] < current_ub){ 
                                        Value value = variable->toValue(valIndex);
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
                                        cout << "CHANGE ALLDIFF SUPPORT " << variable->getName() << " from " << variable->getSupport() << " to " << optimalValue << endl;
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

    // checks that the constraint is still satisfiable (called by WCSP::verify in Debug mode at each search node)
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

    // bool checkEACGreedySolution(int index = -1, Value supportValue = 0) FINAL; // TODO

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
        if (!isSquare) {
            os << ") "
               << " cost: " << -lb << " + " << assigneddeltas << " + (";
            //for (int i = 0; i < arity_; i++) {
                for (unsigned int j = 0; j < deltaCosts.size(); j++) {
                    os << deltaCosts[j];
                    if (j < deltaCosts.size() - 1)
                        os << "|";
                }
                //if (i < arity_ - 1)
                    //os << ",";
            //}
            os << ") ";
        }
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
        
        if(!isSquare){
            Cost maxdelta = *max_element(deltaCosts.begin(), deltaCosts.end());
            if (maxdelta > 0){
                int current_val = 0;
                unordered_map<Value, Cost> mapValuesDeltaCosts;
                for(int valInd = 0; valInd < NbValues ; valInd++){
                    int valIndex = NbValues - valInd;
                    if(deltaCosts[valIndex] == 0) continue;
                    Value value;
                    for(int varIndex =  0; varIndex < arity_ ; varIndex++){
                        auto* variable = scope[varIndex];
                        value = variable->toValue(variable->toIndex(UnionVarDomain[valIndex]));
                        if(variable->canbe(value)){
                            mapValuesDeltaCosts[value] = deltaCosts[valIndex];
                            current_val++;
                            break;
                        }     
                    }
                }
                os <<" " <<current_val;
                for (auto& [key, delta] : mapValuesDeltaCosts) {               
                    os << " " << key;
                    os << " " << delta;
                }
            }
        }
        os << endl;
    }

    void dump_CFN(ostream& os, bool original = true) override
    {
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
        os << "],\"deltacosts\":[";
        if(isSquare){
            os << "]}},\n";
        }
        else{
            Cost maxdelta = *max_element(deltaCosts.begin(), deltaCosts.end());
            if(maxdelta > 0){
                int current_val = 0;
                unordered_map<Value, Cost> mapValuesDeltaCosts;
                for(int valInd = 0; valInd < NbValues ; valInd++){
                    int valIndex = NbValues - valInd;
                    if(deltaCosts[valIndex] == 0) continue;
                    Value value;
                    for(int varIndex =  0; varIndex < arity_ ; varIndex++){
                        auto* variable = scope[varIndex];
                        value = variable->toValue(variable->toIndex(UnionVarDomain[valIndex]));
                        if(variable->canbe(value)){
                            mapValuesDeltaCosts[value] = deltaCosts[valIndex];
                            current_val++;
                            break;
                        }     
                    }
                }
                printed = false;
                for (auto& [key, delta] : mapValuesDeltaCosts) {
                    if (printed)
                        os << ",";
                    os << "[" << key;
                    os << "," << delta;
                    os << "]";
                    printed = true;
                }
            }
            os << "]}},\n";
        }
    }
};

#endif /*TB2ALLDIFFERENT_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
