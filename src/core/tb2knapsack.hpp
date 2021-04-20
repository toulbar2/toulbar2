#ifndef TB2KNAPSACK_HPP_
#define TB2KNAPSACK_HPP_
#include <utility>
#include <variant>
#include "tb2abstractconstr.hpp"
#include "tb2ternaryconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2wcsp.hpp"
#include "../utils/tb2store.hpp"
#include <numeric>
#include <chrono>
#include <thread>
class KnapsackConstraint : public AbstractNaryConstraint {
    StoreLong capacity; // knapsack capacity
    vector<vector<Long>> weights; // knapsack linear positive integer coefficients
    vector<StoreInt> GreatestWeightIdx;
    vector<StoreInt> LowestWeightIdx;
    StoreCost lb; // projected cost to problem lower bound (if it is zero then all deltaCosts must be zero)
    vector<vector<StoreCost>> deltaCosts; // extended costs from unary costs for value 1 to the cost function
    StoreCost assigneddeltas;
    StoreInt nonassigned; // number of non-assigned variables during search, must be backtrackable!
    vector<Long> conflictWeights; // used by weighted degree heuristics
    StoreLong  MaxWeight;
    StoreLong MinWeight;
    Long Original_capacity;
    vector<StoreCost> LastCostforKnapsack; // Store the cost used in the last propagation
    vector<Cost> CostforKnapsack; // temporary data structure for propagate
    vector<vector<Cost>> Profit; // temporary data structure for propagate
    vector<Cost> UnaryCost1;
    vector<Cost> UnaryCost0;
    vector<vector<Double>> OptSol;
    vector<int> arrvar; // temporary data structure for propagate
    vector<vector<Value>> VarVal;
    vector<vector<Value>> NotVarVal;
    vector<Value> lastval1;
    vector<Value> lastval0;
    vector<StoreInt> lastval0ok;
    vector<StoreInt> assigned;
    vector<int> current_scope_idx;
    vector<vector<int>> current_val_idx;
    vector<vector<Double>> Slopes;
    StoreInt carity;
    vector<int> nbValue;
    vector<int> tempveStoreint;
    int LastDomSize;

    void projectLB(Cost c)
    {
        if(c>0){
            lb += c;
            Constraint::projectLB(c);
        }

    }

    static Double Ceil(Double v){
        const Double epsilon = 1e-7;

        if(floorl(v)+epsilon>v)
            return floorl(v);
        else
            return ceill(v);
    }
    void Updatelastval0(int idx){
        if (!lastval0ok[idx] && !scope[idx]->canbe(lastval0[idx]) ) {
            int last = lastval0[idx];
            unsigned int j = 0;
            while (j < NotVarVal[idx].size() && last == lastval0[idx]) {
                if (scope[idx]->canbe(NotVarVal[idx][j])) {
                    lastval0[idx] = NotVarVal[idx][j];
                    VarVal[idx].back() = lastval0[idx];
                } else
                    j++;
            }
            if(last==lastval0[idx])
                lastval0ok[idx]=true;
        }
    }
    //Return if the variable scope[idx] is unassigned
    bool isunassigned(int idx) {
        if (assigned[idx] == 0 && scope[idx]->getDomainSize() > 1) {
            Updatelastval0(idx);
            if (!scope[idx]->canbe(lastval1[idx])) {
                int last = lastval1[idx];
                unsigned int j = 0;
                while (j < VarVal[idx].size() - 1 && last == lastval1[idx]) {
                    if (scope[idx]->canbe(VarVal[idx][j]))
                        lastval1[idx] = VarVal[idx][j];
                    else
                        j++;
                }
                if (last == lastval1[idx]) {
                    return false;
                }
            }
            assert(scope[idx]->canbe(VarVal[idx].back()) || lastval0ok[idx]);
            assert(lastval1[idx] != lastval0[idx]);
            return true;
        } else
            return false;
    }

    void Group_extendNVV(int var, Cost C){
        for(unsigned int i=0;i<NotVarVal[var].size();i++){
            if(scope[var]->canbe(NotVarVal[var][i])){
                assert(scope[var]->getCost(NotVarVal[var][i])>=C);
                scope[var]->extend(NotVarVal[var][i],C);}
        }
    }
    void Group_ProjectNVV(int var, Cost C){
        for(unsigned int i=0;i<NotVarVal[var].size();i++){
            if(scope[var]->canbe(NotVarVal[var][i]))
                scope[var]->project(NotVarVal[var][i],C,true);
        }
    }

    //Depending of the value and the cost, extend or project the cost on the value 1 or 0 of the variable var
    void ExtOrProJ(int var,int value,Cost C )
    {
        if(C>0){
            if(value<(int)VarVal[var].size()-1){
                assert(scope[var]->getCost(VarVal[var][value])>=C);
                scope[var]->extend(VarVal[var][value],C);
                deltaCosts[var][value]+=C;
            }
            else{
                Group_extendNVV(var, C);
                deltaCosts[var].back()+=C;
            }
        } else {
            if(value<(int)VarVal[var].size()-1)
            {
                scope[var]->project(VarVal[var][value],-C,true);
                deltaCosts[var][value]+=C;
            }
            else{
                Group_ProjectNVV(var, -C);
                deltaCosts[var].back()+=C;
            }
        }
    }
    void get_current_scope()
    {
        // recover current scope
        bool greatok;
        bool lowok;
        Long w1;
        unsigned int k=0;
        for (int i = 0; i < arity_; i++) {
            greatok= false;
            lowok=false;
            if (assigned[i]==0 ) {
                Updatelastval0(i);
                tempveStoreint.clear();
                    for (int j = 0; j < (int)VarVal[i].size(); ++j) {
                        if(scope[i]->canbe(VarVal[i][j])) {
                            if(GreatestWeightIdx[i]==j)
                                greatok=true;
                            if(LowestWeightIdx[i]==j)
                                lowok=true;
                            tempveStoreint.push_back(j);
                        }
                    }
                    if(!greatok){
                        MaxWeight-=weights[i][GreatestWeightIdx[i]];
                        GreatestWeightIdx[i]=LowestWeightIdx[i];
                        w1=weights[i][LowestWeightIdx[i]];
                        for (int j = 0; j < (int)VarVal[i].size(); ++j) {
                            if(scope[i]->canbe(VarVal[i][j]) && weights[i][j]>=w1) {
                                w1=weights[i][j];
                                GreatestWeightIdx[i]=j;
                            }
                        }
                        MaxWeight+=weights[i][GreatestWeightIdx[i]];
                    }
                    if(!lowok){
                        MinWeight-=weights[i][LowestWeightIdx[i]];
                        LowestWeightIdx[i]=GreatestWeightIdx[i];
                        w1=weights[i][GreatestWeightIdx[i]];
                        for (int j = 0; j < (int)VarVal[i].size(); ++j) {
                            if(scope[i]->canbe(VarVal[i][j]) && weights[i][j]<=w1) {
                                w1=weights[i][j];
                                LowestWeightIdx[i]=j;
                            }
                        }
                        MinWeight+=weights[i][LowestWeightIdx[i]];
                    }
                    if(k>=current_scope_idx.size()) {
                        current_scope_idx.push_back(i);
                        current_val_idx.push_back(tempveStoreint);
                        nbValue.push_back(current_val_idx.back().size());
                    }else {
                        current_scope_idx[k] = i;
                        current_val_idx[k] = tempveStoreint;
                        nbValue[k] = current_val_idx[k].size();
                    }
                k++;
                assert(scope[i]->canbe(VarVal[i][LowestWeightIdx[i]]));
                assert(scope[i]->canbe(VarVal[i][GreatestWeightIdx[i]]));
                LastDomSize+=scope[i]->getDomainSize();
            }
        }
        if(k<current_val_idx.size()){
            current_scope_idx.resize(k,0);
            current_val_idx.resize(k,current_scope_idx);
            nbValue.resize(k,0);
        }
        assert(current_scope_idx.size()==current_val_idx.size());
        assert(current_scope_idx.size()==k);
        assert(current_val_idx.size()==k);
        assert(nbValue.size()==k);
        carity = current_scope_idx.size();

        assert(carity==current_scope_idx.size());
    }

public:
    KnapsackConstraint(WCSP *wcsp, EnumeratedVariable **scope_in, int arity_in, Long capacity_in,
                       vector<vector<Long>> weights_in,Long MaxWeight_in, vector<vector<int>> VarVal_in,vector<vector<int>> NotVarVal_in)
            : AbstractNaryConstraint(wcsp, scope_in, arity_in)
            , capacity(capacity_in)
            , weights(std::move(weights_in))
            , lb(MIN_COST)
            , assigneddeltas(MIN_COST)
            , nonassigned(arity_in)
            , MaxWeight(MaxWeight_in)
            , MinWeight(0)
            , Original_capacity(capacity_in)
            , VarVal(std::move(VarVal_in))
            , NotVarVal(std::move(NotVarVal_in))
            , carity(arity_in)
    {
        for (int i = 0; i < arity_in; i++) {
            assert(VarVal[i].size()>1);
            lastval0ok.emplace_back(false);
            OptSol.emplace_back(weights[i].size(),MIN_COST);
            Profit.emplace_back(weights[i].size(),MIN_COST);
            current_val_idx.emplace_back(VarVal[i].size(),MIN_COST);
            deltaCosts.emplace_back(weights[i].size(),MIN_COST);
            conflictWeights.push_back(0);
            CostforKnapsack.push_back(MIN_COST);
            LastCostforKnapsack.push_back(MIN_COST);
            assigned.emplace_back(0);
            UnaryCost0.push_back(MIN_COST);
            current_scope_idx.emplace_back(0);
            tempveStoreint.emplace_back(0);
            nbValue.emplace_back(0);
            GreatestWeightIdx.emplace_back(max_element(weights[i].begin(), weights[i].end())-weights[i].begin());
            LowestWeightIdx.emplace_back(min_element(weights[i].begin(), weights[i].end())-weights[i].begin());
            if (NotVarVal[i].empty()) {
                lastval0.push_back(scope[i]->getInf()-1);
            } else {
                lastval0.push_back(NotVarVal[i][0]);
            }
            lastval1.push_back(VarVal[i][0]);
            assert(VarVal[i].size()==weights[i].size());
        }
#ifndef NDEBUG
        Long Sumw=0;
        for (int i = 0; i < arity_; ++i) {
            Sumw+=weights[i][GreatestWeightIdx[i]];
        }
        assert(MaxWeight==Sumw);
#endif
        if(universal()) {
            deconnected();
        }
        get_current_scope();
    }

    virtual ~KnapsackConstraint() {}

    bool extension() const FINAL { return false; } // TODO: allows functional variable elimination but not other preprocessing

    void reconnect()
    {
        if (deconnected()) {
            nonassigned = arity_;
            AbstractNaryConstraint::reconnect();
        }
    }
    int getNonAssigned() const { return nonassigned; }

    Long getConflictWeight() const { return Constraint::getConflictWeight(); }
    Long getConflictWeight(int varIndex) const
    {
        assert(varIndex >= 0);
        assert(varIndex < arity_);
        return conflictWeights[varIndex] + Constraint::getConflictWeight();
    }
    void incConflictWeight(Constraint* from)
    {
        assert(from!=NULL);
        if (from == this) {
            if (deconnected() || nonassigned==arity_) {
                Constraint::incConflictWeight(1);
            } else {
                for (int i = 0; i < arity_; i++) {
                    if (connected(i)) {
                        conflictWeights[i]++; //TODO: increase only variables participating to the relaxed optimal solution cost
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
    void resetConflictWeight()
    {
        conflictWeights.assign(conflictWeights.size(), 0);
        Constraint::resetConflictWeight();
    }

    bool universal()
    {
        // returns true if constraint always satisfied
        if(capacity<=0)
            return true;
        Long minweight=0;
        for (int i = 0; i < carity; ++i) {
            minweight+=weights[current_scope_idx[i]][LowestWeightIdx[current_scope_idx[i]]];
        }
        if (minweight >= capacity)
            return true;
        else
            return false;
    }

    Cost eval(const Tuple& s)
    {
        // returns the cost of the corresponding assignment s
        Long W=0;
        Cost res = -lb+assigneddeltas;
        for (int i = 0; i < arity_; i++) {
            if (ToulBar2::verbose >= 1)
                cout<<s[i]<<" ";
            auto* var = (EnumeratedVariable*)getVar(i);
            auto it=find(VarVal[i].begin(), VarVal[i].end(), var->toValue(s[i]));
            if( it == VarVal[i].end()) {
                res += deltaCosts[i].back();
                W+=weights[i].back();
            }
            else
            {
                W+=weights[i][distance(VarVal[i].begin(),it)];
                res+=deltaCosts[i][distance(VarVal[i].begin(),it)];
            }
        }
        if(W < Original_capacity || res > wcsp->getUb())
            res=wcsp->getUb();
        assert(res <= wcsp->getUb());
        assert(res>=MIN_COST);
        return res;
    }
    Cost evalsubstr(const Tuple& s, Constraint* ctr) FINAL { return evalsubstrAny(s, ctr); }
    Cost evalsubstr(const Tuple& s, NaryConstraint* ctr) FINAL { return evalsubstrAny(s, ctr); }
    template <class T>
    Cost evalsubstrAny(const Tuple& s, T* ctr)
    {
        int count = 0;

        for (int i = 0; i < arity_; i++) {
            int ind = ctr->getIndex(getVar(i));
            if (ind >= 0) {
                evalTuple[i] = s[ind];
                count++;
            }
        }
        assert(count <= arity_);

        Cost cost;
        if (count == arity_)
            cost = eval(evalTuple);
        else
            cost = MIN_COST;

        return cost;
    }
    Cost getCost() FINAL
    {
        for (int i = 0; i < arity_; i++) {
            EnumeratedVariable* var = (EnumeratedVariable*)getVar(i);
            evalTuple[i] = var->toIndex(var->getValue());
        }
        return eval(evalTuple);
    }

    double computeTightness() { return MIN_COST; } //TODO: compute a ratio of feasible tuples divided by getDomainSizeProduct()

    //TODO: needed for dominance test by DEE
    //pair<pair<Cost, Cost>, pair<Cost, Cost>> getMaxCost(int index, Value a, Value b)

    //Cost getMaxFiniteCost() //TODO: return the maximum finite cost for any valid tuple less than wcsp->getUb()
    Cost getMaxFiniteCost()
    {
        Cost delta=0;
        for(int i=0; i<arity_;i++)
        {
            Cost m = *max_element(deltaCosts[i].begin(),deltaCosts[i].end());
            if(m>MIN_COST)
                delta+=m;
        }
        Cost sumdelta = ((lb-assigneddeltas > MIN_COST) ?  delta-lb+assigneddeltas : MIN_COST);
        if (CUT(sumdelta, wcsp->getUb()))
            return MAX_COST;
        else
            return sumdelta;
    }

    //void setInfiniteCost(Cost ub)
    void assign(int varIndex)
    {
        if (assigned[varIndex]==0) {
            assigned[varIndex]=1;
            if(scope[varIndex]->getDomainSize()==1){
                nonassigned=nonassigned-1;
                assigned[varIndex]=2;
                deconnect(varIndex);
            }
            assert(nonassigned >= 0);
            //Update the problem
            auto it=find(VarVal[varIndex].begin(), VarVal[varIndex].end(), scope[varIndex]->getInf());
            if( it == VarVal[varIndex].end()) {
                capacity -= weights[varIndex].back();
                assigneddeltas+=deltaCosts[varIndex].back();
            }else
            {
                capacity -= weights[varIndex][distance(VarVal[varIndex].begin(),it)];
                assigneddeltas+=deltaCosts[varIndex][distance(VarVal[varIndex].begin(),it)];
            }
            fill(deltaCosts[varIndex].begin(),deltaCosts[varIndex].end(),0);
            MaxWeight -= weights[varIndex][GreatestWeightIdx[varIndex]];
            MinWeight-=weights[varIndex][LowestWeightIdx[varIndex]];
            get_current_scope();
            if (universal()) {
                deconnect();
                Cost TobeProjected = -lb + assigneddeltas;
                lb = MIN_COST;
                Cost mindelta;
                for (int i = 0; i < carity; i++) {
                    mindelta=MAX_COST;
                    for (int j = 0; j < nbValue[i]; ++j) {
                        if (mindelta > deltaCosts[current_scope_idx[i]][current_val_idx[i][j]])
                            mindelta = deltaCosts[current_scope_idx[i]][current_val_idx[i][j]];
                    }
                    TobeProjected += mindelta;
                    for (int j = 0; j < nbValue[i] ; ++j) {
                        assert(mindelta <= deltaCosts[current_scope_idx[i]][current_val_idx[i][j]]);
                        if(current_val_idx[i][j]==(int)VarVal[current_scope_idx[i]].size()-1){
                            Group_ProjectNVV(current_scope_idx[i],deltaCosts[current_scope_idx[i]][current_val_idx[i][j]] - mindelta);
                        }
                        else
                            scope[current_scope_idx[i]]->project(VarVal[current_scope_idx[i]][current_val_idx[i][j]],deltaCosts[current_scope_idx[i]][current_val_idx[i][j]] - mindelta);
                    }
                    scope[current_scope_idx[i]]->findSupport();
                }
                assert(TobeProjected >= MIN_COST);
                Constraint::projectLB(TobeProjected);
            } else if (nonassigned <= 3 ) {
            deconnect(); // this constraint is removed from the current WCSP problem
            if (ToulBar2::verbose >= 2)
                cout << "EVALUATION " << endl;
            projectNary();
            } else {
                //TODO: incremental bound propagation

                propagate();
                if (ToulBar2::FullEAC)
                    reviseEACGreedySolution();
            }
        }else{
            if(assigned[varIndex]==1 && scope[varIndex]->getDomainSize()==1){
                nonassigned=nonassigned-1;
                assigned[varIndex]=2;
                if(nonassigned <= 3){
                    deconnect(); // this constraint is removed from the current WCSP problem
                    if (ToulBar2::verbose >= 2)
                        cout<<"EVALUATION "<<endl;
                    projectNary();
                }
            }
        }
    }

    //Return True if a value has been deleted else return False
    bool BoundConsistency()
    {
        int k = 0, k2;
        bool b = false;
        int currentvar;
        while (k < carity && !b) {
            currentvar=current_scope_idx[k];
            //Determine if at least one value of the current variable can be erased.
            if (MaxWeight - weights[currentvar][GreatestWeightIdx[currentvar]] + weights[currentvar][LowestWeightIdx[currentvar]]<capacity) {
                k2 = 0;
                if(assigned[currentvar]==0 && isunassigned(currentvar)) {
                    while (k2 < nbValue[k] && !b) {
                        if (MaxWeight - weights[currentvar][GreatestWeightIdx[currentvar]] + weights[currentvar][current_val_idx[k][k2]] < capacity) {
                            if (current_val_idx[k][k2] == (int)VarVal[currentvar].size() - 1) {
                                for (unsigned int i = 0; i < NotVarVal[currentvar].size(); ++i) {
                                    if (scope[currentvar]->canbe(NotVarVal[currentvar][i])) {
                                        scope[currentvar]->remove(NotVarVal[currentvar][i]);
                                    }
                                }
                            } else {
                                scope[currentvar]->remove(VarVal[currentvar][current_val_idx[k][k2]]);
                            }
                            if(!connected() || assigned[currentvar]>0)
                                b= true;
                        }
                        k2++;
                    }
                    b=true;
                }
            }
            k++;
        }
        return b;
    }

    bool ComputeProfit()
    {
        Cost verifopt=-lb+assigneddeltas;   //Used to check if the last optimal solution has still a cost of 0
        Long verifweight=0;        //Used to check if the last optimal solution has still a cost of 0
        Double storec,storew;   //Used to check if the last optimal solution has still a cost of 0
        int sumOPT=0;
        Cost UN;
        bool b1=true;
        for (int i = 0; i < carity; i++) {
            assert(assigned[current_scope_idx[i]]==0 && scope[current_scope_idx[i]]->getDomainSize()>1);
            storec = 0;
            storew = 0;
            for (int j = 0; j < nbValue[i] ; j++) {
                assert(scope[current_scope_idx[i]]->getCost(VarVal[current_scope_idx[i]][current_val_idx[i][j]])>= MIN_COST);
                if(current_val_idx[i][j]==(int)VarVal[current_scope_idx[i]].size()-1){
                    UnaryCost0[current_scope_idx[i]]=MAX_COST;
                    for (unsigned int l = 0; l < NotVarVal[current_scope_idx[i]].size(); l++) {
                        if (scope[current_scope_idx[i]]->canbe(NotVarVal[current_scope_idx[i]][l])) {
                            assert(scope[current_scope_idx[i]]->getCost(NotVarVal[current_scope_idx[i]][l])>=MIN_COST);
                            UN = scope[current_scope_idx[i]]->getCost(NotVarVal[current_scope_idx[i]][l]);
                            if (UN < UnaryCost0[current_scope_idx[i]]) {
                                UnaryCost0[current_scope_idx[i]] = UN;

                            }
                        }
                    }
                    Profit[current_scope_idx[i]][current_val_idx[i][j]] =  UnaryCost0[current_scope_idx[i]] + deltaCosts[current_scope_idx[i]][current_val_idx[i][j]];
                }
                else{
                    Profit[current_scope_idx[i]][current_val_idx[i][j]]=scope[current_scope_idx[i]]->getCost(VarVal[current_scope_idx[i]][current_val_idx[i][j]]) + deltaCosts[current_scope_idx[i]][current_val_idx[i][j]];
                }
                storew += weights[current_scope_idx[i]][current_val_idx[i][j]] * OptSol[current_scope_idx[i]][current_val_idx[i][j]];
                storec += Profit[current_scope_idx[i]][current_val_idx[i][j]] * OptSol[current_scope_idx[i]][current_val_idx[i][j]];
                sumOPT+=OptSol[current_scope_idx[i]][current_val_idx[i][j]];
            }
            //Compute the cost of the last optimal solution
            verifopt += Ceil(storec);
            verifweight += Ceil(storew);
            if(b1 && (*max_element(Profit[i].begin(),Profit[i].end())!=MIN_COST || *min_element(Profit[i].begin(),Profit[i].end())!=MIN_COST))
                b1= false;
            assert(*max_element(Profit[i].begin(),Profit[i].end())<MAX_COST);
        }
        assert(verifopt>= 0 || verifweight<capacity || sumOPT<carity);
        if(b1)
            return false;
        if(verifopt>0 || verifweight<capacity || sumOPT<carity)
            return true;
        else
            return false;
    }

    //Return a vector containing the slopes, a slope is defined between 2 values of a variable : <Var, Val1, Val2, slopes>
    void ComputeSlopes(Long *W,Cost *c)
    {
        int item1 = 0;
        int k = 0;
        Slopes.clear();
        for (int i = 0; i < carity; i++) {
            fill(OptSol[current_scope_idx[i]].begin(),OptSol[current_scope_idx[i]].end(),0);
            //Sort the value in ascending weight
            arrvar.clear();
            arrvar=current_val_idx[i];
            sort(arrvar.begin(), arrvar.end(),
                    [&](int x, int y) {if(weights[current_scope_idx[i]][x] == weights[current_scope_idx[i]][y]){
                        if(Profit[current_scope_idx[i]][x] == Profit[current_scope_idx[i]][y]){
                            if(scope[current_scope_idx[i]]->getSupport()==VarVal[current_scope_idx[i]][x]){
                                return true;
                            }
                            else
                                return false;
                        }else
                            return Profit[current_scope_idx[i]][x] > Profit[current_scope_idx[i]][y];
                    }else
                        return weights[current_scope_idx[i]][x] < weights[current_scope_idx[i]][y];});
            //Find the value with the heaviest weight
            k = arrvar.size() - 1;
            item1 = arrvar[k];
            while (k > 0) {
                k--;
                //We don't consider dominated items : p_i < p_j && w_i > w_j   (i dominate j)
                if (Profit[current_scope_idx[i]][arrvar[k]] < Profit[current_scope_idx[i]][item1] &&
                    weights[current_scope_idx[i]][arrvar[k]] < weights[current_scope_idx[i]][item1]) {
                    //If it is the first slope for the current variable, we directy add it : <Var, val1, val2, slope>
                    //Else we compare the new slope with the precedent one to verify if the precedent value isn't dominated.
                    //If it is dominated, we replace the last slope else we add a new one
                    if (Slopes.empty() || Slopes.back()[0] != current_scope_idx[i] || Slopes.back()[3] >=Double((Profit[current_scope_idx[i]][item1] - Profit[current_scope_idx[i]][arrvar[k]])) /(weights[current_scope_idx[i]][item1] -weights[current_scope_idx[i]][arrvar[k]]))
                        Slopes.push_back(
                                      {Double(current_scope_idx[i]), Double(arrvar[k]), Double(item1),
                                       Double((Profit[current_scope_idx[i]][item1] - Profit[current_scope_idx[i]][arrvar[k]])) /
                                       (weights[current_scope_idx[i]][item1] - weights[current_scope_idx[i]][arrvar[k]])});
                    else{
                        Slopes.back() = {Double(current_scope_idx[i]), Double(arrvar[k]), Slopes.back()[2],
                                     Double((Profit[current_scope_idx[i]][Slopes.back()[2]] - Profit[current_scope_idx[i]][arrvar[k]])) /
                                     (weights[current_scope_idx[i]][Slopes.back()[2]] - weights[current_scope_idx[i]][arrvar[k]])};
                        while (Slopes.size()>1 && Slopes.end()[-2][0]==Slopes.back()[0] && Slopes.back()[3]>=Slopes.end()[-2][3]){
                            Slopes.end()[-2]= {Double(current_scope_idx[i]), Slopes.back()[1], Slopes.end()[-2][2],
                                             Double((Profit[current_scope_idx[i]][Slopes.end()[-2][2]] - Profit[current_scope_idx[i]][Slopes.back()[1]])) /
                                             (weights[current_scope_idx[i]][Slopes.end()[-2][2]] - weights[current_scope_idx[i]][Slopes.back()[1]])};
                            Slopes.pop_back();
                        }

                    }
                    item1 = arrvar[k];
                }
            }
            //Compute a first Solution
            if (Slopes.size() == 0 || Slopes.back()[0] != current_scope_idx[i]) {
                OptSol[current_scope_idx[i]][item1] = 1;
                *W += weights[current_scope_idx[i]][item1];
                *c += Profit[current_scope_idx[i]][item1];
            } else {
                *W += weights[current_scope_idx[i]][Slopes.back()[1]];
                OptSol[current_scope_idx[i]][Slopes.back()[1]] = 1;
                *c += Profit[current_scope_idx[i]][Slopes.back()[1]];
            }
        }
    }

    //Find the optimal solution. Follow the order of the slopes and modify OptSol until we fill the constraint
    void FindOpt(vector<vector<Double>> Slopes, Long *W,Cost *c,Double *xk, int *iter)
    {
        int currentVar;
        Long capacityLeft;
        while (*W < capacity) {
            currentVar = int(Slopes[*iter][0]);
            if (*W + weights[currentVar][Slopes[*iter][2]] - weights[currentVar][Slopes[*iter][1]] >=capacity) {
                capacityLeft = capacity - *W;
                *xk = Double(capacityLeft) / (weights[currentVar][Slopes[*iter][2]] - weights[currentVar][Slopes[*iter][1]]);
                OptSol[currentVar][Slopes[*iter][2]] = *xk;
                OptSol[currentVar][Slopes[*iter][1]] = 1 - *xk;
                *W += weights[currentVar][Slopes[*iter][2]] - weights[currentVar][Slopes[*iter][1]];
                *c += Ceil(*xk *(Profit[currentVar][Slopes[*iter][2]] - Profit[currentVar][Slopes[*iter][1]]));
                assert(capacityLeft > 0);
            } else {
                assert(OptSol[currentVar][Slopes[*iter][1]]==1);
                OptSol[currentVar][Slopes[*iter][1]] = 0;
                OptSol[currentVar][Slopes[*iter][2]] = 1;
                *W += weights[currentVar][Slopes[*iter][2]] - weights[currentVar][Slopes[*iter][1]];
                *c += Profit[currentVar][Slopes[*iter][2]] - Profit[currentVar][Slopes[*iter][1]];
                *iter=*iter+1;
            }
        }
    }

    //Do the Extension/Projection
    void ExtensionProjection(vector<Double> y_i, Double y_cc)
    {
        for (int i = 0; i < carity; i++) {
            for (int j = 0; j < nbValue[i]; ++j) {
                    if (OptSol[current_scope_idx[i]][current_val_idx[i][j]] == 1) {
                        if (current_val_idx[i][j] == (int)VarVal[current_scope_idx[i]].size() - 1) {
                            Group_extendNVV(current_scope_idx[i], UnaryCost0[current_scope_idx[i]]);
                            deltaCosts[current_scope_idx[i]][current_val_idx[i][j]] += UnaryCost0[current_scope_idx[i]];
                        } else {
                            deltaCosts[current_scope_idx[i]][current_val_idx[i][j]] += scope[current_scope_idx[i]]->getCost(VarVal[current_scope_idx[i]][current_val_idx[i][j]]);
                            scope[current_scope_idx[i]]->extend(VarVal[current_scope_idx[i]][current_val_idx[i][j]], scope[current_scope_idx[i]]->getCost(VarVal[current_scope_idx[i]][current_val_idx[i][j]]));
                        }
                    } else if (OptSol[current_scope_idx[i]][current_val_idx[i][j]] == 0) {
                        if (Ceil(-deltaCosts[current_scope_idx[i]][current_val_idx[i][j]] + y_i[i] + y_cc * weights[current_scope_idx[i]][current_val_idx[i][j]])!=0)
                            ExtOrProJ(current_scope_idx[i], current_val_idx[i][j],Ceil(-deltaCosts[current_scope_idx[i]][current_val_idx[i][j]] + y_i[i] + y_cc * weights[current_scope_idx[i]][current_val_idx[i][j]]));
                    } else {
                        if (current_val_idx[i][j] == (int)VarVal[current_scope_idx[i]].size() - 1) {
                            Group_extendNVV(current_scope_idx[i], UnaryCost0[current_scope_idx[i]]);
                            deltaCosts[current_scope_idx[i]][current_val_idx[i][j]] += UnaryCost0[current_scope_idx[i]];
                        } else {
                            deltaCosts[current_scope_idx[i]][current_val_idx[i][j]] += scope[current_scope_idx[i]]->getCost(VarVal[current_scope_idx[i]][current_val_idx[i][j]]);
                            scope[current_scope_idx[i]]->extend(VarVal[current_scope_idx[i]][current_val_idx[i][j]], scope[current_scope_idx[i]]->getCost(VarVal[current_scope_idx[i]][current_val_idx[i][j]]));
                        }
                    }
            }
            scope[current_scope_idx[i]]->findSupport();
        }
    }

    void propagate()
    {
        // propagates from scratch the constraint
        //auto start0 = std::chrono::system_clock::now();
        if (connected()) {
            bool b=false;
            for (int i = 0; connected() && i < arity_; i++) {
                if (assigned[i]==0 && !isunassigned(i)){
                    assign(i);}  /* indique à la contrainte que cette  variable est affectée (donc met à jour nonassigned..) */
                else
                    assert(assigned[i]>0 || scope[i]->getDomainSize()>1);
            }
            if (connected() ) {
                if (!verify()) {
                    THROWCONTRADICTION;
                } else if (nonassigned > 3 && ToulBar2::LcLevel >= LC_AC) {
                    get_current_scope();
#ifndef NDEBUG
                    for (int i = 0; i < carity; ++i) {
                        assert(scope[current_scope_idx[i]]->canbe(VarVal[current_scope_idx[i]].back()) || lastval0ok[current_scope_idx[i]]);
                        assert(scope[current_scope_idx[i]]->getDomainSize()>1);
                        assert(assigned[current_scope_idx[i]]==0);
                        for (int j = 0; j < nbValue[i]; ++j) {
                            assert(scope[current_scope_idx[i]]->canbe(VarVal[current_scope_idx[i]][current_val_idx[i][j]]));
                            assert(scope[current_scope_idx[i]]->canbe(VarVal[current_scope_idx[i]][GreatestWeightIdx[current_scope_idx[i]]]));
                            assert(scope[current_scope_idx[i]]->canbe(VarVal[current_scope_idx[i]][LowestWeightIdx[current_scope_idx[i]]]));
                            assert(weights[current_scope_idx[i]][GreatestWeightIdx[current_scope_idx[i]]]>=weights[current_scope_idx[i]][current_val_idx[i][j]]);
                            assert(weights[current_scope_idx[i]][LowestWeightIdx[current_scope_idx[i]]]<=weights[current_scope_idx[i]][current_val_idx[i][j]]);
                        }

                        assert(nbValue[i]>1);
                    }
#endif

                    //Bound propagation, return true if a variable has been assigned
                    b=BoundConsistency();
                    if(!b){
#ifndef NDEBUG
                        for (int i = 0; i < carity; ++i) {
                            assert(scope[current_scope_idx[i]]->canbe(VarVal[current_scope_idx[i]].back()) || lastval0ok[current_scope_idx[i]]);
                            assert(scope[current_scope_idx[i]]->getDomainSize()>1);
                            assert(assigned[current_scope_idx[i]]==0);
                            for (int j = 0; j < nbValue[i]; ++j) {
                                assert(MaxWeight-weights[current_scope_idx[i]][GreatestWeightIdx[current_scope_idx[i]]]+weights[current_scope_idx[i]][current_val_idx[i][j]]>=capacity);
                                assert(scope[current_scope_idx[i]]->canbe(VarVal[current_scope_idx[i]][current_val_idx[i][j]]));
                                assert(scope[current_scope_idx[i]]->canbe(VarVal[current_scope_idx[i]][GreatestWeightIdx[current_scope_idx[i]]]));
                                assert(scope[current_scope_idx[i]]->canbe(VarVal[current_scope_idx[i]][LowestWeightIdx[current_scope_idx[i]]]));
                                assert(weights[current_scope_idx[i]][GreatestWeightIdx[current_scope_idx[i]]]>=weights[current_scope_idx[i]][current_val_idx[i][j]]);
                                assert(weights[current_scope_idx[i]][LowestWeightIdx[current_scope_idx[i]]]<=weights[current_scope_idx[i]][current_val_idx[i][j]]);
                            }
                            assert(nbValue[i]>1);
                        }
#endif
                        // Return True if the last optimal solution has a cost of 0
                        b=ComputeProfit();

                        if (connected() && b) {
                            Long W = 0;
                            Cost c = -lb + assigneddeltas;
                            ComputeSlopes(&W,&c); // temporary data structure for propagate

                            if(ToulBar2::verbose >= 4) {
                                cout << "cap is " << capacity << endl;
                                for (int i = 0; i < carity; ++i) {
                                        cout << scope[current_scope_idx[i]]->getName();
                                        for (unsigned int j = 0; j < VarVal[current_scope_idx[i]].size(); ++j) {
                                            cout << " " << j << " : " << Profit[current_scope_idx[i]][j] << "/" << weights[current_scope_idx[i]][j];
                                        }
                                        cout << endl;
                                }
                            }
#ifndef NDEBUG
                            for (unsigned int i = 0; i < Slopes.size(); ++i) {
                                assert(Slopes[i].size() == 4);
                                assert(Slopes[i][3] >= MIN_COST);
                                assert(Slopes[i][3] < MAX_COST);
                            }
                            Long Sumw=0;
                            for (int i = 0; i < carity; ++i) {
                                Sumw+=weights[current_scope_idx[i]][GreatestWeightIdx[current_scope_idx[i]]];
                            }
                            assert(Sumw==MaxWeight);
                            assert(c < MAX_COST);
                            assert(W >= 0);
#endif
                            //Sort the Slopes
                            sort(Slopes.begin(), Slopes.end(),
                                 [&](vector<Double> x, vector<Double> y) {
                                     if (x[3] == y[3]) {
                                         if (x[0]==y[0])
                                             return weights[int(x[0])][int(x[1])]<=weights[int(y[0])][int(y[1])];
                                         else
                                             return scope[int(x[0])]->getDACOrder() < scope[int(y[0])]->getDACOrder();
                                     }else
                                         return x[3] < y[3];
                                 });
                            int iter = 0;
                            Double xk = 0;

                            //Find the optimal solution
                            FindOpt(Slopes,&W,&c,&xk,&iter);
                            assert(W >= capacity);
                            assert(xk >= 0 && xk <= 1);
                            assert(iter <= Slopes.size());
                            assert(c > -1);
#ifndef NDEBUG
                            Double t = 0;
                            double a = 0;
                            double epsi = 0.00001;
                            for (int i = 0; i < carity; ++i) {
                                a = 0;
                                for (int j = 0; j < OptSol[current_scope_idx[i]].size(); ++j) {
                                    a += OptSol[current_scope_idx[i]][j];
                                }
                                assert(a < 1 + epsi && a > 1 - epsi);
                                for (int j = 0; j < VarVal[current_scope_idx[i]].size(); ++j) {
                                    t += OptSol[current_scope_idx[i]][j] * weights[current_scope_idx[i]][j];
                                }
                             }
                            assert(Ceil(t) > capacity - epsi && Ceil(t) < capacity + epsi || iter == 0);
#endif

                            if (c > 0) {
                                //Compute the dual variable y_cc and y_i using the optimal primal solution
                                Double y_cc = 0;
                                if(!Slopes.empty())
                                     y_cc = Slopes[iter][3];
                                if (xk == 0)
                                    y_cc = 0;
                                vector<Double> y_i;
                                assert(y_cc >= MIN_COST);
                                int k=0;
                                for (int i = 0; i < carity; i++) {
                                    k = 0;
                                    while (OptSol[current_scope_idx[i]][current_val_idx[i][k]] == 0)
                                        k++;
                                    y_i.push_back(Profit[current_scope_idx[i]][current_val_idx[i][k]] - y_cc * weights[current_scope_idx[i]][current_val_idx[i][k]]);
                                    assert(y_i[i] < MAX_COST);
                                }
                                //Use y_cc and y_i to extend/project the right cost
                                ExtensionProjection(y_i,y_cc);

                                if(ToulBar2::verbose >=4) {
                                    cout<<"projected cost "<<c<<" LB : "<<lb-assigneddeltas<<endl;
                                    for (int i = 0; i < carity; ++i) {
                                        cout << " Delta " << scope[current_scope_idx[i]]->getName() << " : ";
                                        for (unsigned int j = 0; j < VarVal[current_scope_idx[i]].size(); ++j) {
                                            cout << " " << deltaCosts[current_scope_idx[i]][j];
                                        }
                                        cout << endl;
                                    }
                                }
                                assert(c>0);
                                projectLB(c);
                                assert(getMaxFiniteCost()>=0);
                            }
                        }
                    }
                } else {
                    get_current_scope();
                    if (universal()) {
                        deconnect();
                        Cost TobeProjected = -lb + assigneddeltas;
                        lb = MIN_COST;
                        Cost mindelta;
                        for (int i = 0; i < carity; i++) {
                            mindelta=MAX_COST;
                            if(scope[current_scope_idx[i]]->getDomainSize()>1) {
                                for (int j = 0; j < nbValue[i]; ++j) {
                                    if (mindelta > deltaCosts[current_scope_idx[i]][current_val_idx[i][j]])
                                        mindelta = deltaCosts[current_scope_idx[i]][current_val_idx[i][j]];
                                }
                                TobeProjected += mindelta;
                                for (int j = 0; j < nbValue[i] ; ++j) {
                                    assert(mindelta <= deltaCosts[current_scope_idx[i]][j]);
                                    if(current_val_idx[i][j]==(int)VarVal[current_scope_idx[i]].size() - 1)
                                        Group_ProjectNVV(current_scope_idx[i],deltaCosts[current_scope_idx[i]][current_val_idx[i][j]] - mindelta);
                                    else
                                        scope[current_scope_idx[i]]->project(VarVal[current_scope_idx[i]][j],deltaCosts[current_scope_idx[i]][current_val_idx[i][j]] - mindelta);
                                }
                                scope[current_scope_idx[i]]->findSupport();
                            }else{
                                auto it1=find(VarVal[current_scope_idx[i]].begin(), VarVal[current_scope_idx[i]].end(), scope[current_scope_idx[i]]->getInf());
                                if( it1 == VarVal[current_scope_idx[i]].end()) {
                                    TobeProjected+=deltaCosts[current_scope_idx[i]].back();
                                }else
                                {
                                    TobeProjected+=deltaCosts[current_scope_idx[i]][distance(VarVal[current_scope_idx[i]].begin(),it1)];
                                }
                            }
                        }
                        assert(TobeProjected >= MIN_COST && TobeProjected < MAX_COST);
                        Constraint::projectLB(TobeProjected);
                    }else if (nonassigned <= 3) {
                        deconnect(); // this constraint is removed from the current WCSP problem
                        if (ToulBar2::verbose >= 2)
                            cout << "EVALUATION " << endl;
                        projectNary(); // and replaced by a ternary constraint in extension
                     }
                }
            }
        }
    }

    bool verify()
    {
        // checks that propagation has been done correctly such that at least there exists one valid tuple with zero cost (called by WCSP::verify in Debug mode at each search node)
        if(capacity<=MaxWeight) {
            return true;
        }
        else
            return false;
    }
    void increase(int index) override {
        remove(index);
    }
    void decrease(int index) override {
        remove(index);
    }
    void remove(int idx) override
    {
        if (isunassigned(idx)) {
            get_current_scope();
            propagate();
        } else if (assigned[idx] < 2) {
            assign(idx);
        }
    }
    void projectFromZero(int index)
    {
        //TODO: incremental cost propagation
        if((int)current_scope_idx.size()!=carity )
            get_current_scope();
        else{
            int sumVa=0;
            for (int i = 0; i < carity; ++i) {
                sumVa+=scope[current_scope_idx[i]]->getDomainSize();
            }
            if(sumVa!=LastDomSize)
                get_current_scope();
        }

        propagate();
    }

    bool checkEACGreedySolution(int index = -1, Value supportValue = 0) FINAL
    {
        Long W=0;
        Cost res = -lb+assigneddeltas;
        for (int i = 0; i < arity_; i++) {
            Value support = ((i == index) ? supportValue : scope[i]->getSupport());
            auto it=find(VarVal[i].begin(), VarVal[i].end(), support);
            if( it == VarVal[i].end()) {
                res += deltaCosts[i].back();
                W+=weights[i].back();
            }
            else
            {
                W+=weights[i][distance(VarVal[i].begin(),it)];
                res+=deltaCosts[i][distance(VarVal[i].begin(),it)];
            }
        }
        if(W < Original_capacity)
            res=wcsp->getUb();
        return (res==MIN_COST);
    }

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

    void dump(ostream& os, bool original = true)
    {
        bool iszerodeltas = (lb == MIN_COST);
        for (int i = 0; i < arity_; ++i) {
            for (auto it = deltaCosts[i].begin(); it != deltaCosts[i].end(); ++it) {
                Cost d = (*it);
                if (d != MIN_COST) {
                    iszerodeltas = false;
                    break;
                }
            }
        }
        if (original) {
            os << arity_;
            for (int i = 0; i < arity_; i++)
                os << " " << scope[i]->wcspIndex;
            if (iszerodeltas) {
                os << " " << -1 << " knapsackp " << Original_capacity;
                for (int i = 0; i < arity_; i++) {
                    os << " " << VarVal[i].size();
                    for (unsigned int j = 0; j < VarVal[i].size(); ++j) {
                        os << " "<< VarVal[i][j];
                        os << " " << weights[i][j];
                    }
                }
                os << endl;
            } else {
                os << " " << 0 << " " << getDomainSizeProduct() << endl;
                Tuple t;
                Cost c;
                firstlex();
                while (nextlex(t, c)) {
                    for (int i = 0; i < arity_; i++) {
                        os << t[i] << " ";
                    }
                    os << c << endl;
                }
            }
        } else {
            os << nonassigned;
            for (int i = 0; i < arity_; i++)
                if (scope[i]->unassigned())
                    os << " " << scope[i]->getCurrentVarId();
            if (iszerodeltas) {
                os << " " << -1 << " knapsack " << capacity;
                for (int i = 0; i < arity_; i++) {
                    if (scope[i]->unassigned()) {
                        os << " " << VarVal[i].size();
                        for (unsigned int j = 0; j < VarVal[i].size(); ++j) {
                            os << " " << scope[i]->toCurrentIndex(VarVal[i][j]);
                            os << " " << weights[i][0];
                        }
                    }
                }
                os << endl;
            } else {
                os << " " << 0 << " " << getDomainSizeProduct() << endl;
                Tuple t;
                Cost c;
                firstlex();
                while (nextlex(t, c)) {
                    for (int i = 0; i < arity_; i++) {
                        if (scope[i]->unassigned())
                            os << scope[i]->toCurrentIndex(scope[i]->toValue(t[i])) << " ";
                    }
                    os << min(wcsp->getUb(), c) << endl;
                }
            }
        }
    }

    void dump_CFN(ostream& os, bool original = true)
    {
        bool printed = false;
        os << "\"F_";

        bool iszerodeltas = (lb == MIN_COST);
        for (int i = 0; i < arity_; ++i) {
            for (auto it = deltaCosts[i].begin(); it != deltaCosts[i].end(); ++it) {
                Cost d = (*it);
                if (d != MIN_COST) {
                    iszerodeltas = false;
                    break;
                }
            }
        }

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
                os << scope[i]->getName();
                printed = true;
            }
            os << "],\n\"type\":\"knapsackp\",\n\"params\":{\"capacity\":" << Original_capacity << ",\n\t\"weightedvalues\":[";

            if (iszerodeltas) {
                printed = false;
                for (int i = 0; i < arity_; i++) {
                    if (printed)
                        os << ",";
                    os << "[";
                    bool printedbis = false;
                    for (unsigned int j = 0; j < VarVal[i].size(); ++j) {
                        if (printedbis) os << ",";
                        os << VarVal[i][j] << "," << weights[i][j];
                        printedbis=true;
                    }
                    printed = true;
                    os << "]";
                }
            } else {
                Tuple t;
                Cost c;
                printed = false;
                firstlex();
                while (nextlex(t, c)) {
                    os << endl;
                    for (int i = 0; i < arity_; i++) {
                        if (printed)
                            os << ",";
                        os << ((scope[i]->isValueNames()) ? scope[i]->getValueName(t[i]) : std::to_string(t[i]));
                        printed = true;
                    }
                    os << "," << wcsp->Cost2RDCost(c);
                }
            }
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
                    os << scope[i]->getName();
                    printed = true;
                }
            os << "],\n\"type\":\"knapsackp\",\n\"params\":{\"capacity\":" << capacity << ",\n\t\"weightedvalues\":[";

            if (iszerodeltas) {
                printed = false;
                for (int i = 0; i < arity_; i++) {
                    if (scope[i]->unassigned()) {
                        if (printed)
                            os << ",";
                        os << "[";
                        bool printedbis = false;
                        for (unsigned int j = 0; j < VarVal[i].size(); ++j) {
                            if (printedbis) os << ",";
                            os << scope[i]->toCurrentIndex(VarVal[i][j]) << "," << weights[i][j];
                            printedbis=true;
                        }
                        printed = true;
                        os << "]";
                    }
                }
            } else {
                Tuple t;
                Cost c;
                printed = false;
                firstlex();
                while (nextlex(t, c)) {
                    os << endl;
                    for (int i = 0; i < arity_; i++) {
                        if (scope[i]->unassigned()) {
                            if (printed)
                                os << ",";
                            os << scope[i]->toCurrentIndex(scope[i]->toValue(t[i]));
                            printed = true;
                        }
                    }
                    os << "," << wcsp->Cost2RDCost(min(wcsp->getUb(), c));
                }
            }
        }
        os << "]}},\n";
    }
};
#endif /*TB2KNAPSACK_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
