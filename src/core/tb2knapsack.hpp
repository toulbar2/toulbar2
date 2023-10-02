#ifndef TB2KNAPSACK_HPP_
#define TB2KNAPSACK_HPP_
#include <utility>
#include <variant>
#include "tb2abstractconstr.hpp"
#include "tb2ternaryconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2wcsp.hpp"
#include "../utils/tb2store.hpp"

class KnapsackConstraint : public AbstractNaryConstraint {
    int carity;
    Long Original_capacity;
    Cost Original_ub;
    StoreInt nonassigned; // number of non-assigned variables during search, must be backtrackable!
    StoreLong capacity; // knapsack capacity
    StoreLong MaxWeight;
    StoreCost lb; // projected cost to problem lower bound (if it is zero then all deltaCosts must be zero)
    StoreCost assigneddeltas;
    vector<Long> conflictWeights; // used by weighted degree heuristics
    vector<StoreInt> LowestWeightIdx;
    vector<Cost> InitLargestWeight;
    vector<StoreInt> GreatestWeightIdx;
    vector<StoreInt> assigned;
    vector<StoreInt> lastval0ok;
    vector<int> nbValue;
    vector<int> current_scope_idx;
    vector<Value> lastval1;
    vector<Value> lastval0;
    vector<Cost> UnaryCost0;
    vector<vector<int>> current_val_idx;
    vector<vector<Long>> weights; // knapsack linear positive integer coefficients
    vector<vector<Value>> VarVal;
    vector<vector<Value>> NotVarVal;
    vector<vector<StoreCost>> deltaCosts; // extended costs from unary costs for value 1 to the cost function
    vector<vector<Double>> OptSol;
    vector<int> arrvar; // temporary data structure for propagate
    vector<vector<Cost>> Profit; // temporary data structure for propagate
    vector<Double> y_i, tempAMOy_i;
    vector<vector<Double>> yAMO_i;
    vector<std::array<Double, 4>> Slopes;
    int nbRealVar;
    vector<vector<pair<int, int>>> AMO; //Represent AMO constraints
    vector<vector<Long>> Original_weigths;
    vector<int> VirtualVar, CorrAMO; //VirtualVar has size Number of AMO constraint + number of variable in no AMO constraint, it returns if a given index represents a virtual var (an AMO constraint) or a real one.
    //CorrAMO has size arity and gives on which AMO constraint the given variable belongs (0 if it is in no AMO constraint).
    vector<StoreInt> nbVirtualVar;
    StoreInt AlwaysSatisfied;
    bool DoDyn = false;
    Long NbNodes = -1;
    bool sameweight;
    vector<vector<Cost>> tempdeltaCosts;
    vector<bool> Booleanvar;

    void projectLB(Cost c)
    {
        if (c > 0) {
            lb += c;
            Constraint::projectLB(c);
        }
    }

    static Double Ceil(Double v)
    {
        const Double epsilon = 1e-7;

        if (floorl(v) + epsilon > v)
            return floorl(v);
        else
            return ceill(v);
    }
    void Updatelastval0(int idx)
    {
        if (!lastval0ok[idx] && !scope[idx]->canbe(lastval0[idx])) {
            int last = lastval0[idx];
            unsigned int j = 0;
            while (j < NotVarVal[idx].size() && last == lastval0[idx]) {
                if (scope[idx]->canbe(NotVarVal[idx][j])) {
                    lastval0[idx] = NotVarVal[idx][j];
                    VarVal[idx].back() = lastval0[idx];
                } else
                    j++;
            }
            if (last == lastval0[idx])
                lastval0ok[idx] = true;
        }
    }
    //Return if the variable scope[idx] is unassigned
    bool isunassigned(int idx)
    {

        if (CorrAMO[idx] == 0) {
            if (assigned[idx] == 0 && scope[idx]->unassigned()) {
                if (Booleanvar[idx])
                    return true;
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
        } else {
            if (scope[idx]->assigned())
                return false;
            else
                return true;
        }
    }

    void Group_extendNVV(int var, Cost C)
    {
        for (unsigned int i = 0; i < NotVarVal[var].size(); i++) {
            if (scope[var]->canbe(NotVarVal[var][i])) {
                assert(scope[var]->getCost(NotVarVal[var][i]) >= C);
                scope[var]->extend(NotVarVal[var][i], C);
            }
        }
    }
    void Group_ProjectNVV(int var, Cost C)
    {
        for (unsigned int i = 0; i < NotVarVal[var].size(); i++) {
            if (scope[var]->canbe(NotVarVal[var][i]))
                scope[var]->project(NotVarVal[var][i], C, true);
        }
    }

    //Depending of the value and the cost, extend or project the cost on the value 1 or 0 of the variable var
    void ExtOrProJ(int var, int value, Cost C)
    {
        if (C > 0) {
            if (value < (int)VarVal[var].size() - 1) {
                assert(scope[var]->getCost(VarVal[var][value]) >= C);
                scope[var]->extend(VarVal[var][value], C);
                deltaCosts[var][value] += C;
            } else {
                Group_extendNVV(var, C);
                deltaCosts[var].back() += C;
            }
        } else {
            if (value < (int)VarVal[var].size() - 1) {
                scope[var]->project(VarVal[var][value], -C, true);
                deltaCosts[var][value] += C;
            } else {
                Group_ProjectNVV(var, -C);
                deltaCosts[var].back() += C;
            }
        }
    }
    void get_current_scope()
    {
        // recover current scope
        bool greatok;
        bool lowok;
        Long w1;
        int k1 = 0;
        carity = 0;
        for (int i = 0; i < nbRealVar; i++) {
            greatok = false;
            lowok = false;
            if (assigned[i] == 0) {
                nbValue[k1] = 0;
                if (Booleanvar[i] && scope[i]->unassigned()) {
                    current_scope_idx[k1] = i;
                    carity++;
                    current_val_idx[k1][0] = 0;
                    current_val_idx[k1][1] = 1;
                    nbValue[k1] = 2;
                    k1++;
                } else {
                    Updatelastval0(i);
                    for (int j = 0; j < (int)VarVal[i].size(); ++j) {
                        if (scope[i]->canbe(VarVal[i][j])) {
                            if (GreatestWeightIdx[i] == j)
                                greatok = true;
                            if (LowestWeightIdx[i] == j)
                                lowok = true;
                            current_val_idx[k1][nbValue[k1]] = j;
                            nbValue[k1]++;
                        }
                    }
                    if (!greatok) {
                        MaxWeight -= weights[i][GreatestWeightIdx[i]];
                        GreatestWeightIdx[i] = LowestWeightIdx[i];
                        w1 = weights[i][LowestWeightIdx[i]];
                        for (int j = 0; j < nbValue[k1]; ++j) {
                            if (weights[i][current_val_idx[k1][j]] >= w1) {
                                w1 = weights[i][current_val_idx[k1][j]];
                                GreatestWeightIdx[i] = current_val_idx[k1][j];
                            }
                        }
                        MaxWeight += weights[i][GreatestWeightIdx[i]];
                    }
                    if (!lowok) {
                        LowestWeightIdx[i] = GreatestWeightIdx[i];
                        w1 = weights[i][GreatestWeightIdx[i]];
                        for (int j = 0; j < nbValue[k1]; ++j) {
                            if (weights[i][current_val_idx[k1][j]] <= w1) {
                                w1 = weights[i][current_val_idx[k1][j]];
                                LowestWeightIdx[i] = current_val_idx[k1][j];
                            }
                        }
                    }
                    current_scope_idx[k1] = i;
                    k1++;
                    carity++;
                    assert(scope[i]->canbe(VarVal[i][LowestWeightIdx[i]]));
                    assert(scope[i]->canbe(VarVal[i][GreatestWeightIdx[i]]));
                }
            }
        }

        for (unsigned int i = 0; i < AMO.size(); ++i) {
            greatok = false;
            lowok = false;
            nbValue[k1] = 0;
            if (nbVirtualVar[i] >= 1) {
                for (int j = 0; j < (int)AMO[i].size(); ++j) {
                    if (assigned[AMO[i][j].first] == 0) {
                        if (GreatestWeightIdx[nbRealVar + i] == j)
                            greatok = true;
                        if (LowestWeightIdx[nbRealVar + i] == j)
                            lowok = true;
                        current_val_idx[k1][nbValue[k1]] = j;
                        nbValue[k1] = nbValue[k1] + 1;
                    } else {
                        assert(scope[AMO[i][j].first]->assigned());
                    }
                }
                if (GreatestWeightIdx[nbRealVar + i] == (int)AMO[i].size() && !lastval0ok[nbRealVar + i]) {
                    greatok = true;
                }
                if (LowestWeightIdx[nbRealVar + i] == (int)AMO[i].size() && !lastval0ok[nbRealVar + i])
                    lowok = true;
                if (nbValue[k1] > 0) {
                    if (lastval0ok[i + nbRealVar] == false) {
                        current_val_idx[k1][nbValue[k1]] = (int)AMO[i].size();
                        nbValue[k1]++;
                    }
                    if (!greatok) {
                        MaxWeight -= weights[nbRealVar + i][GreatestWeightIdx[nbRealVar + i]];
                        GreatestWeightIdx[nbRealVar + i] = LowestWeightIdx[nbRealVar + i];
                        w1 = weights[nbRealVar + i][LowestWeightIdx[nbRealVar + i]];
                        for (int j = 0; j < nbValue[k1]; ++j) {
                            if (weights[nbRealVar + i][current_val_idx[k1][j]] >= w1) {
                                w1 = weights[nbRealVar + i][current_val_idx[k1][j]];
                                GreatestWeightIdx[nbRealVar + i] = current_val_idx[k1][j];
                            }
                        }
                        MaxWeight += weights[nbRealVar + i][GreatestWeightIdx[nbRealVar + i]];
                    }

                    if (!lowok) {
                        LowestWeightIdx[nbRealVar + i] = GreatestWeightIdx[nbRealVar + i];
                        w1 = weights[nbRealVar + i][GreatestWeightIdx[nbRealVar + i]];
                        for (int j = 0; j < nbValue[k1]; ++j) {
                            if (weights[nbRealVar + i][current_val_idx[k1][j]] <= w1) {
                                w1 = weights[nbRealVar + i][current_val_idx[k1][j]];
                                LowestWeightIdx[nbRealVar + i] = current_val_idx[k1][j];
                            }
                        }
                    }
                    if (nbValue[k1] > 1) {
                        carity++;
                        current_scope_idx[k1] = i + nbRealVar;
                        k1++;
                    }
                }
            }
        }
        assert(current_scope_idx.size() == current_val_idx.size());
    }

public:
    KnapsackConstraint(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in, Long capacity_in,
        vector<vector<Long>> weights_in, Long MaxWeight_in, vector<vector<Value>> VarVal_in, vector<vector<Value>> NotVarVal_in, vector<vector<pair<int, int>>> AMO_in, vector<vector<Long>> Original_weigths_in, vector<int> CorrAMO_in, vector<int> VirtualVar_in, int nonassinged_in, const vector<vector<StoreCost>> InitDel = vector<vector<StoreCost>>(), const StoreCost lb_in = MIN_COST, const StoreCost assigneddelta_in = MIN_COST, const Long Original_capacity_in = 0)
        : AbstractNaryConstraint(wcsp, scope_in, arity_in)
        , carity(arity_in)
        , Original_capacity(capacity_in)
        , Original_ub(wcsp->getUb())
        , nonassigned(nonassinged_in)
        , capacity(capacity_in)
        , MaxWeight(MaxWeight_in)
        , lb(MIN_COST)
        , assigneddeltas(MIN_COST)
        , weights(std::move(weights_in))
        , VarVal(std::move(VarVal_in))
        , NotVarVal(std::move(NotVarVal_in))
        , AMO(std::move(AMO_in))
        , Original_weigths(std::move(Original_weigths_in))
        , VirtualVar(std::move(VirtualVar_in))
        , CorrAMO(std::move(CorrAMO_in))
        , AlwaysSatisfied(0)
    {
        if (weights.size() > 0) {
            if (InitDel.size() != 0)
                deltaCosts = InitDel;
            if (Original_capacity_in != 0)
                Original_capacity = Original_capacity_in;
            unsigned int maxdom = VarVal[0].size();
            sameweight = true;
            for (int i = 0; i < (int)weights.size(); i++) {
                assert(VarVal[i].size() > 1);
                OptSol.emplace_back(weights[i].size(), 0.);
                Profit.emplace_back(weights[i].size(), MIN_COST);
                if (InitDel.size() == 0)
                    deltaCosts.emplace_back(weights[i].size(), MIN_COST);
                conflictWeights.push_back(0);
                assigned.emplace_back(0);
                UnaryCost0.push_back(MIN_COST);
                nbValue.emplace_back(0);
                GreatestWeightIdx.emplace_back(max_element(weights[i].begin(), weights[i].end()) - weights[i].begin());
                LowestWeightIdx.emplace_back(min_element(weights[i].begin(), weights[i].end()) - weights[i].begin());
                InitLargestWeight.emplace_back(weights[i][GreatestWeightIdx[i]]);
                if (NotVarVal[i].empty()) {
                    NotVarVal[i].push_back(VarVal[i].back());
                }
                lastval0.push_back(NotVarVal[i][0]);
                lastval0ok.emplace_back(false);
                lastval1.push_back(VarVal[i][0]);
                assert(VarVal[i].size() == weights[i].size());
                if (maxdom < VarVal[i].size())
                    maxdom = VarVal[i].size();
                if (sameweight && weights[0] != weights[i])
                    sameweight = false;
                if (VarVal[i].size() + NotVarVal[i].size() <= 3) {
                    Booleanvar.push_back(true);
                } else
                    Booleanvar.push_back(false);
            }
            for (int j = 0; j < (int)weights.size(); ++j) {
                current_val_idx.emplace_back(maxdom, MIN_COST);
                current_scope_idx.emplace_back(0);
            }
            nbRealVar = (int)weights.size() - (int)AMO.size();
            if (!AMO.empty()) {
                if (InitDel.size() == 0)
                    deltaCosts.clear();
                assigned.clear();
                conflictWeights.clear();
                for (int i = 0; i < arity_; i++) {
                    conflictWeights.push_back(0);
                    if (scope[i]->unassigned())
                        assigned.emplace_back(0);
                    else
                        assigned.emplace_back(2);
                    if (InitDel.size() == 0)
                        deltaCosts.emplace_back(2, MIN_COST);
                }
                for (int i = 0; i < (int)AMO.size(); ++i) {
                    nbVirtualVar.emplace_back((int)AMO[i].size());
                }
            }
            get_current_scope();
#ifndef NDEBUG
            Long Sumw = 0;
            for (int i = 0; i < carity; ++i) {
                Sumw += weights[current_scope_idx[i]][GreatestWeightIdx[current_scope_idx[i]]];
            }
            assert(MaxWeight == Sumw);
#endif
            if (universal()) {
                deconnect();
            }
        } else {
            if (!verify()) {
                THROWCONTRADICTION;
            }
            deconnect();
        }
    }

    virtual ~KnapsackConstraint() {}

    bool extension() const FINAL { return false; } // TODO: allows functional variable elimination but not other preprocessing
    bool isKnapsack() const FINAL { return true; }

    void reconnect() override
    {
        if (deconnected()) {
            nonassigned = arity_;
            AbstractNaryConstraint::reconnect();
        }
    }
    int getNonAssigned() const { return nonassigned; }

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
        if (!sameweight) { //TODO: else Constraint::incConflictWeight(1)
            if (from == this) {
                if (deconnected() || nonassigned == arity_) {
                    Constraint::incConflictWeight(1);
                } else {
                    get_current_scope();
                    if (!verify()) {
                        Long SumW = 0, MaxW;
                        vector<int> EDACORDER;
                        for (int i = 0; i < nbRealVar; ++i) {
                            EDACORDER.push_back(i);
                        }
                        sort(EDACORDER.begin(), EDACORDER.end(),
                            [&](int& x, int& y) {
                                return scope[x]->getDACOrder() < scope[y]->getDACOrder(); //TODO: checks if it favors aborbing more unary costs for the last variables in the DAC order or the opposite?!
                            });
                        for (int i = 0; i < nbRealVar; ++i) {
                            MaxW = -1;
                            if (!scope[i]->canbe(VarVal[i][GreatestWeightIdx[i]])) {
                                for (unsigned int j = 0; j < VarVal[i].size() - 1; ++j) {
                                    if (scope[i]->canbe(VarVal[i][j]) && MaxW < weights[i][j]) {
                                        GreatestWeightIdx[i] = j;
                                        MaxW = weights[i][j];
                                    }
                                }
                                if (!scope[i]->canbe(VarVal[i].back())) {
                                    unsigned int k = 0;
                                    while (k < NotVarVal[i].size() && !scope[i]->canbe(NotVarVal[i][k])) {
                                        k = k + 1;
                                    }
                                    if (k < NotVarVal[i].size() && MaxW < weights[i].back()) {
                                        GreatestWeightIdx[i] = VarVal[i].size() - 1;
                                        MaxW = weights[i].back();
                                    }
                                } else if (MaxW < weights[i].back()) {
                                    GreatestWeightIdx[i] = VarVal[i].size() - 1;
                                    MaxW = weights[i].back();
                                }
                            }
                            SumW += weights[i][GreatestWeightIdx[i]];
                        }
                        for (unsigned int i = 0; i < AMO.size(); ++i) {
                            if (lastval0ok[nbRealVar + i])
                                MaxW = weights[nbRealVar + i].back();
                            else
                                MaxW = -1;
                            for (unsigned int j = 0; j < AMO[i].size(); ++j) {
                                if (scope[AMO[i][j].first]->canbe(AMO[i][j].second) && weights[nbRealVar + i][j] > MaxW) {
                                    GreatestWeightIdx[nbRealVar + i] = j;
                                    MaxW = weights[nbRealVar + i][j];
                                }
                            }
                            SumW += MaxW;
                        }
                        assert(SumW < Original_capacity);
                        for (unsigned int i = 0; i < weights.size(); ++i) {
                            if (VirtualVar[i] == 0) {
                                if (SumW - weights[EDACORDER[i]][GreatestWeightIdx[EDACORDER[i]]] + InitLargestWeight[EDACORDER[i]] >= Original_capacity)
                                    conflictWeights[EDACORDER[i]]++;
                                else
                                    SumW += -weights[EDACORDER[i]][GreatestWeightIdx[EDACORDER[i]]] + InitLargestWeight[EDACORDER[i]];
                            } else {
                                if (SumW - weights[i][GreatestWeightIdx[i]] + InitLargestWeight[i] >= Original_capacity) {
                                    for (unsigned int j = 0; j < AMO[VirtualVar[i] - 1].size(); ++j) {
                                        if (scope[conflictWeights[AMO[VirtualVar[i] - 1][j].first]]->assigned())
                                            conflictWeights[AMO[VirtualVar[i] - 1][j].first]++;
                                    }
                                } else
                                    SumW += -weights[i][GreatestWeightIdx[i]] + InitLargestWeight[i];
                            }
                        }
                    } else {
                        Constraint::incConflictWeight(1); //Increase ConflictWeight of the constraint or the conflict weight of each assigned variables ?
                    }
                }
            } else if (deconnected()) {
                for (int i = 0; i < from->arity(); i++) {
                    if (from->getVar(i)->unassigned()) {
                        int index = getIndex(from->getVar(i));
                        if (index >= 0) { // the last conflict constraint may be derived from two binary constraints (boosting search), each one derived from an n-ary constraint with a scope which does not include parameter constraint from
                            assert(index < arity_);
                            conflictWeights[index]++;
                        }
                    }
                }
            }
        }
    }
    void resetConflictWeight() override
    {
        conflictWeights.assign(conflictWeights.size(), 0);
        Constraint::resetConflictWeight();
    }
    vector<int> GetOrder()
    {
        vector<Double> WeightedProfit;
        vector<int> OutScope;
        vector<int> ROutScope;
        for (int i = 0; i < arity_; ++i) {
            WeightedProfit.push_back(0);
            if (scope[i]->unassigned()) {
                WeightedProfit.back() = Double(scope[i]->getCost(1) - scope[i]->getCost(0) + deltaCosts[i][1] - deltaCosts[i][0]) / (weights[i][1] - weights[i][0]);
                OutScope.push_back(i);
            }
        }
        sort(OutScope.begin(), OutScope.end(),
            [&](int& x, int& y) {
                if (WeightedProfit[x] == WeightedProfit[y]) {
                    return weights[x] < weights[y];
                } else
                    return WeightedProfit[x] < WeightedProfit[y];
            });
        for (int i = 0; i < (int)OutScope.size(); ++i) {
            ROutScope.push_back(scope[OutScope[i]]->wcspIndex);
        }
        return ROutScope;
    }
    Double getLag()
    {
        get_current_scope();
        ComputeProfit();
        Long W = 0;
        Cost c = -lb + assigneddeltas;
        ComputeSlopes(&W, &c); // temporary data structure for propagate
        int iter = 0;
        Double xk = 0;
        if (W < capacity) {
            //Sort the Slopes
            sort(Slopes.begin(), Slopes.end(),
                [&](std::array<Double, 4>& x, std::array<Double, 4>& y) {
                    if (x[3] == y[3]) {
                        if (x[0] == y[0])
                            return weights[int(x[0])][int(x[1])] <= weights[int(y[0])][int(y[1])];
                        else
                            return scope[int(x[0])]->getDACOrder() < scope[int(y[0])]->getDACOrder(); //TODO: checks if it favors aborbing more unary costs for the last variables in the DAC order or the opposite?!
                    } else
                        return x[3] < y[3];
                });
            //Find the optimal solution
            FindOpt(Slopes, &W, &c, &xk, &iter);
        }
        //Compute the dual variable y_cc and y_i using the optimal primal solution
        Double y_cc = 0;
        if (!Slopes.empty())
            y_cc = Slopes[iter][3];
        if (xk == 0)
            y_cc = 0;
        return y_cc;
    }

    bool universal() override
    {
        // returns true if constraint always satisfied
        if (AMO.empty() && capacity <= 0)
            return true;
        Long minweight = 0;
        for (int i = 0; i < carity; ++i) {
            minweight += weights[current_scope_idx[i]][LowestWeightIdx[current_scope_idx[i]]];
        }
        if (!AMO.empty() && minweight >= capacity)
            AlwaysSatisfied = 1;
        else
            AlwaysSatisfied = 0;

        if (!AMO.empty() && *max_element(nbVirtualVar.begin(), nbVirtualVar.end()) > 1)
            return false;

        if (minweight >= capacity)
            return true;
        else
            return false;
    }

    Cost eval(const Tuple& s) override
    {
        // returns the cost of the corresponding assignment s
        Long W = 0;
        int k, k1;
        vector<int> breakAMO;
        bool ok;
        for (unsigned int i = 0; i < AMO.size(); ++i) {
            breakAMO.push_back(0);
        }
        Cost res = -lb + assigneddeltas;
        for (int i = 0; i < arity_; i++) {
            auto* var = (EnumeratedVariable*)getVar(i);
            if (ToulBar2::verbose >= 4)
                cout << var->getName() << " " << s[i] << " ";
            if (CorrAMO[i] == 0) {
                auto it = find(VarVal[i].begin(), VarVal[i].end(), var->toValue(s[i]));
                if (it == VarVal[i].end()) {
                    res += deltaCosts[i].back();
                    W += weights[i].back();
                } else {
                    W += weights[i][distance(VarVal[i].begin(), it)];
                    res += deltaCosts[i][distance(VarVal[i].begin(), it)];
                }
            } else {
                k = nbRealVar;
                ok = false;
                while (!ok) {
                    if (scope[k] == var) {
                        k1 = 0;
                        while (!ok) {
                            if (AMO[CorrAMO[i] - 1][k1].first == k) {
                                ok = true;
                                if (AMO[CorrAMO[i] - 1][k1].second == var->toValue(s[i])) {
                                    W += weights[nbRealVar + CorrAMO[i] - 1][k1];
                                    breakAMO[CorrAMO[i] - 1] = breakAMO[CorrAMO[i] - 1] + 1;
                                }
                            }
                            k1++;
                        }
                    }
                    k++;
                }
                if (breakAMO[CorrAMO[i] - 1] > 1)
                    res = wcsp->getUb();
                res += deltaCosts[i][var->toValue(s[i])];
            }
        }
        for (unsigned int i = 0; i < AMO.size(); ++i) {
            if (breakAMO[i] == 0)
                W += weights[nbRealVar + i].back();
        }
        if (W < Original_capacity || res > wcsp->getUb()) {
            if (W < Original_capacity && Original_ub < wcsp->getUb() && 1.0 * Original_ub * (Original_capacity - W) < wcsp->getUb()) {
                res = Original_ub * (Original_capacity - W); // VNS-like methods may exploit a relaxation of the constraint
            } else {
                res = wcsp->getUb();
            }
        }
        assert(res <= wcsp->getUb());
        assert(res >= MIN_COST);
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
            assert(var->assigned());
            evalTuple[i] = var->toIndex(var->getValue());
        }
        return eval(evalTuple);
    }

    double computeTightness() override
    {
        //TODO: check if multiplying by the sum of median/mean unary costs (only on VarVal?) improve the results
        //TODO: see if arity plays a role (small arity first?)
        //        double ucost = UNIT_COST;
        //        for (int i = 0; i < arity_; i++) {
        //            EnumeratedVariable* var = (EnumeratedVariable*)getVar(i);
        //            int domsize = var->getDomainSize();
        //            ValueCost array[domsize];
        //            wcsp->getEnumDomainAndCost(var->wcspIndex, array);
        //            if (ToulBar2::weightedTightness == 2) {
        //                Cost unarymediancost = stochastic_selection<ValueCost>(array, 0, domsize - 1, domsize / 2).cost;
        //                ucost += (double)unarymediancost;
        //            } else {
        //                Cost unarytotalcost = MIN_COST;
        //                for (auto& elt : array) {
        //                    unarytotalcost += elt.cost;
        //                }
        //                ucost += (double)unarytotalcost / (double)domsize;
        //            }
        //        }
        assert(capacity > 0);
        assert(MaxWeight > 0);
        return (double)((Double)capacity / (Double)MaxWeight); // * ucost ???
    }

    //TODO: needed for dominance test by DEE
    //pair<pair<Cost, Cost>, pair<Cost, Cost>> getMaxCost(int index, Value a, Value b)

    //Cost getMaxFiniteCost() //TODO: return the maximum finite cost for any valid tuple less than wcsp->getUb()
    Cost getMaxFiniteCost() override
    {
        Cost delta = 0;
        for (int i = 0; i < arity_; i++) {
            Cost m = *max_element(deltaCosts[i].begin(), deltaCosts[i].end());
            if (m > MIN_COST)
                delta += m;
        }
        Cost sumdelta = ((lb - assigneddeltas > MIN_COST) ? delta - lb + assigneddeltas : MIN_COST);
        if (CUT(sumdelta, wcsp->getUb()))
            return MAX_COST;
        else
            return sumdelta;
    }

    void addAMOConstraints(vector<vector<pair<int, int>>> clq, vector<Variable*> vars, WCSP* wcsp)
    {
        vector<Long> TempWeight;
        vector<int> SortedVec, Temp;
        vector<Value> TempVal;
        MaxWeight = 0;
        int nbvarinAMO = 0;
        Original_weigths = weights;
        for (unsigned int i = 0; i < VarVal.size(); ++i) {
            if (VarVal[i][0] == 1) {
                reverse(VarVal[i].begin(), VarVal[i].end());
                reverse(Original_weigths[i].begin(), Original_weigths[i].end());
                reverse(deltaCosts[i].begin(), deltaCosts[i].end());
                if ((int)NotVarVal[i].size() > 0)
                    NotVarVal[i][0] = 1;
                else
                    NotVarVal[i].push_back(1);
            }
        }
        vector<int> ScopeWcspIdx;
        for (int i = clq.size() - 1; i > -1; --i) {
            nbvarinAMO += clq[i].size();
            for (int j = clq[i].size() - 1; j > -1; --j) {
                ScopeWcspIdx.push_back(clq[i][j].first);
            }
        }
        for (int i = 0; i < arity_; ++i) {
            SortedVec.push_back(i);
        }

        sort(SortedVec.begin(), SortedVec.end(),
            [&](int x, int y) {auto it = find(ScopeWcspIdx.begin(), ScopeWcspIdx.end(), scope[x]->wcspIndex);
                 auto it2 = find(ScopeWcspIdx.begin(), ScopeWcspIdx.end(), scope[y]->wcspIndex);
                 return(it> it2); });

        vector<EnumeratedVariable*> scopeVars(arity_);
        for (int i = 0; i < arity_; i++) {
            scopeVars[i] = (EnumeratedVariable*)vars[scope[SortedVec[i]]->wcspIndex];
        }

        int m = 0;
        for (unsigned int i = 0; i < clq.size(); ++i) {
            AMO.push_back(clq[i]);
            for (unsigned int j = 0; j < clq[i].size(); ++j) {
                AMO[i][j].first = arity_ - nbvarinAMO + m;
                m++;
            }
        }
        for (int i = 0; i < (int)AMO.size(); ++i) {
            for (int j = 0; j < (int)AMO[i].size(); ++j) {
                assert(scopeVars[AMO[i][j].first]->wcspIndex == clq[i][j].first);
                assert(AMO[i][j].second == clq[i][j].second);
            }
        }

        vector<vector<StoreCost>> sortdel = deltaCosts;
        vector<vector<Long>> sortOW = Original_weigths;
        vector<Long> sortconfW = conflictWeights;
        vector<StoreInt> sortass = assigned;
        assert((int)assigned.size() == arity_);
        assert((int)conflictWeights.size() == arity_);
        assert((int)Original_weigths.size() == arity_);
        assert((int)deltaCosts.size() == arity_);
        for (int i = 0; i < arity_; ++i) {
            deltaCosts[i] = sortdel[SortedVec[i]];
            Original_weigths[i] = sortOW[SortedVec[i]];
            conflictWeights[i] = sortconfW[SortedVec[i]];
            assigned[i] = sortass[SortedVec[i]];
            if (assigned[i] != 0) {
                assert(scopeVars[i]->assigned());
            } else
                assert(scopeVars[i]->unassigned());
        }
        weights.clear();
        VarVal.clear();
        NotVarVal.clear();
        CorrAMO.clear();
        VirtualVar.clear();
        for (int i = 0; i < arity_ - nbvarinAMO; ++i) {
            VirtualVar.push_back(0);
            CorrAMO.push_back(0);
            VarVal.push_back(vector<int>{ 0, 1 });
            NotVarVal.push_back(vector<int>{ 1 });
            weights.push_back(Original_weigths[i]);
            if (scopeVars[i]->unassigned())
                MaxWeight += *max_element(weights[i].begin(), weights[i].end());
        }
        int n = 0;
        for (int i = 0; i < (int)clq.size(); ++i) {
            VirtualVar.push_back(i + 1);
            TempWeight.clear();
            TempVal.clear();
            for (int j = 0; j < (int)clq[i].size(); ++j) {
                CorrAMO.push_back(i + 1);
                TempWeight.push_back(0);
                for (int k = 0; k < (int)clq[i].size(); ++k) {
                    if (k == j) {
                        TempWeight.back() += Original_weigths[arity_ - nbvarinAMO + k + n][clq[i][k].second];
                    } else {
                        TempWeight.back() += Original_weigths[arity_ - nbvarinAMO + k + n][1 - clq[i][k].second];
                    }
                }
            }
            TempWeight.push_back(0);
            for (int j = 0; j < (int)clq[i].size(); ++j) {
                TempVal.push_back(j);
                TempWeight.back() += Original_weigths[arity_ - nbvarinAMO + j + n][1 - clq[i][j].second];
            }
            MaxWeight += *max_element(TempWeight.begin(), TempWeight.end());
            TempVal.push_back(clq[i].size());
            VarVal.push_back(TempVal);
            TempVal.clear();
            TempVal.push_back(clq[i].size());
            NotVarVal.push_back(TempVal);
            weights.push_back(TempWeight);
            n += clq[i].size();
        }
        //TODO: reinitialize current knapsack constraint instead of creating a new constraint
        deconnect();
        new KnapsackConstraint(wcsp, scopeVars.data(), arity_, capacity, weights, MaxWeight, VarVal, NotVarVal, AMO, Original_weigths, CorrAMO, VirtualVar, nonassigned, deltaCosts, lb, assigneddeltas, Original_capacity);
    }

    void setInfiniteCost(Cost ub) override
    {
        Original_ub = min(ub, Original_ub);
    }

    void assign(int varIndex) override
    {
        if (!verify()) {
            THROWCONTRADICTION;
        } else {
            if (assigned[varIndex] == 0) {
                assigned[varIndex] = 1;
                if (scope[varIndex]->assigned()) {
                    nonassigned = nonassigned - 1;
                    assigned[varIndex] = 2;
                    deconnect(varIndex);
                }
                //Update the problem
                if (CorrAMO[varIndex] == 0) {
                    auto it = find(VarVal[varIndex].begin(), VarVal[varIndex].end(), scope[varIndex]->getInf());
                    if (it == VarVal[varIndex].end()) {
                        capacity -= weights[varIndex].back();
                        assigneddeltas += deltaCosts[varIndex].back();
                    } else {
                        capacity -= weights[varIndex][distance(VarVal[varIndex].begin(), it)];
                        assigneddeltas += deltaCosts[varIndex][distance(VarVal[varIndex].begin(), it)];
                    }
                    fill(deltaCosts[varIndex].begin(), deltaCosts[varIndex].end(), MIN_COST);
                    MaxWeight -= weights[varIndex][GreatestWeightIdx[varIndex]];
                } else {
                    bool ok = false;
                    int k = 0;
                    vector<int> toassign;
                    if (nbVirtualVar[CorrAMO[varIndex] - 1] > 0) {
                        toassign.clear();
                        while (!ok) {
                            if (AMO[CorrAMO[varIndex] - 1][k].first == varIndex) {
                                ok = true;
                                if (AMO[CorrAMO[varIndex] - 1][k].second == scope[varIndex]->getInf()) {
                                    capacity -= weights[nbRealVar + CorrAMO[varIndex] - 1][k];
                                    assigneddeltas += deltaCosts[varIndex][scope[varIndex]->toIndex(scope[varIndex]->getInf())];
                                    fill(deltaCosts[varIndex].begin(), deltaCosts[varIndex].end(), MIN_COST);
                                    MaxWeight -= weights[nbRealVar + CorrAMO[varIndex] - 1][GreatestWeightIdx[nbRealVar + CorrAMO[varIndex] - 1]];
                                    nbVirtualVar[CorrAMO[varIndex] - 1] = 0;
                                    for (unsigned int i = 0; i < AMO[CorrAMO[varIndex] - 1].size(); ++i) {
                                        if (assigned[AMO[CorrAMO[varIndex] - 1][i].first] == 0) {
                                            assigneddeltas += deltaCosts[AMO[CorrAMO[varIndex] - 1][i].first][1 - AMO[CorrAMO[varIndex] - 1][i].second];
                                            assigned[AMO[CorrAMO[varIndex] - 1][i].first] = 2;
                                            toassign.push_back(i);
                                            fill(deltaCosts[AMO[CorrAMO[varIndex] - 1][i].first].begin(), deltaCosts[AMO[CorrAMO[varIndex] - 1][i].first].end(), MIN_COST);
                                        }
                                    }
                                    for (unsigned int i = 0; i < toassign.size(); ++i) {
                                        if (scope[AMO[CorrAMO[varIndex] - 1][toassign[i]].first]->unassigned()) {
                                            nonassigned = nonassigned - 1;
                                            scope[AMO[CorrAMO[varIndex] - 1][toassign[i]].first]->remove(AMO[CorrAMO[varIndex] - 1][toassign[i]].second);
                                        }
                                    }
                                } else {
                                    assigneddeltas += deltaCosts[varIndex][scope[varIndex]->toIndex(scope[varIndex]->getInf())];
                                    fill(deltaCosts[varIndex].begin(), deltaCosts[varIndex].end(), MIN_COST);
                                    nbVirtualVar[CorrAMO[varIndex] - 1] = nbVirtualVar[CorrAMO[varIndex] - 1] - 1;
                                    if (nbVirtualVar[CorrAMO[varIndex] - 1] == 0) {
                                        capacity -= weights[nbRealVar + CorrAMO[varIndex] - 1].back();
                                        MaxWeight -= weights[nbRealVar + CorrAMO[varIndex] - 1][GreatestWeightIdx[nbRealVar + CorrAMO[varIndex] - 1]];
                                        nbVirtualVar[CorrAMO[varIndex] - 1] = 0;
                                    } else if (lastval0ok[nbRealVar + CorrAMO[varIndex] - 1] && nbVirtualVar[CorrAMO[varIndex] - 1] == 1) {
                                        ok = false;
                                        int k1 = 0;
                                        while (!ok) {
                                            if (assigned[AMO[CorrAMO[varIndex] - 1][k1].first] == 0) {
                                                ok = true;
                                                capacity -= weights[nbRealVar + CorrAMO[varIndex] - 1][k1];
                                                assigneddeltas += deltaCosts[AMO[CorrAMO[varIndex] - 1][k1].first][AMO[CorrAMO[varIndex] - 1][k1].second];
                                                fill(deltaCosts[AMO[CorrAMO[varIndex] - 1][k1].first].begin(), deltaCosts[AMO[CorrAMO[varIndex] - 1][k1].first].end(), MIN_COST);
                                                MaxWeight -= weights[nbRealVar + CorrAMO[varIndex] - 1][GreatestWeightIdx[nbRealVar + CorrAMO[varIndex] - 1]];
                                                assigned[AMO[CorrAMO[varIndex] - 1][k1].first] = 2;
                                                nonassigned = nonassigned - 1;
                                                nbVirtualVar[CorrAMO[varIndex] - 1] = 0;
                                                if (scope[AMO[CorrAMO[varIndex] - 1][k1].first]->unassigned()) {
                                                    scope[AMO[CorrAMO[varIndex] - 1][k1].first]->remove(1 - AMO[CorrAMO[varIndex] - 1][k1].second);
                                                }
                                            }
                                            k1++;
                                        }
                                    }
                                }
                            }
                            k++;
                        }
                    }
                }
                get_current_scope();
                if (!verify()) {
                    THROWCONTRADICTION;
                }
                if (connected()) {
                    if (universal()) {
                        deconnect();
                        wcsp->revise(this);
                        Cost TobeProjected = -lb + assigneddeltas;
                        lb = MIN_COST;
                        Cost mindelta;
                        for (int i = 0; i < carity; i++) {
                            mindelta = MAX_COST;
                            if (VirtualVar[current_scope_idx[i]] == 0) {
                                for (int j = 0; j < nbValue[i]; ++j) {
                                    if (mindelta > deltaCosts[current_scope_idx[i]][current_val_idx[i][j]])
                                        mindelta = deltaCosts[current_scope_idx[i]][current_val_idx[i][j]];
                                }
                                TobeProjected += mindelta;
                                for (int j = 0; j < nbValue[i]; ++j) {
                                    assert(mindelta <= deltaCosts[current_scope_idx[i]][current_val_idx[i][j]]);
                                    if (current_val_idx[i][j] == (int)VarVal[current_scope_idx[i]].size() - 1) {
                                        Group_ProjectNVV(current_scope_idx[i], deltaCosts[current_scope_idx[i]][current_val_idx[i][j]] - mindelta);
                                    } else
                                        scope[current_scope_idx[i]]->project(VarVal[current_scope_idx[i]][current_val_idx[i][j]], deltaCosts[current_scope_idx[i]][current_val_idx[i][j]] - mindelta);
                                }
                                scope[current_scope_idx[i]]->findSupport();
                            } else {
                                assert(nbValue[i] <= 2);
                                mindelta = min(deltaCosts[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][0]].first][0], deltaCosts[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][0]].first][1]);
                                TobeProjected += mindelta;
                                scope[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][0]].first]->project(0, deltaCosts[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][0]].first][0] - mindelta);
                                scope[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][0]].first]->project(1, deltaCosts[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][0]].first][1] - mindelta);
                                scope[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][0]].first]->findSupport();
                            }
                        }
                        assert(TobeProjected >= MIN_COST);
                        Constraint::projectLB(TobeProjected);

                    } else if (nonassigned <= NARYPROJECTIONSIZE && (nonassigned < 3 || maxInitDomSize <= NARYPROJECTION3MAXDOMSIZE)) {
                        if (connected()) {
                            deconnect(); // this constraint is removed from the current WCSP problem
                            projectNary();
                        }
                    } else if (AlwaysSatisfied == 1) {
                        Cost mindelta, temp0, temp1;
                        for (int i = 0; i < carity; i++) {
                            mindelta = MAX_COST;
                            if (VirtualVar[current_scope_idx[i]] == 0) {
                                if (scope[current_scope_idx[i]]->unassigned()) {
                                    mindelta = min(deltaCosts[current_scope_idx[i]][0], deltaCosts[current_scope_idx[i]][1]);
                                    if (mindelta != 0) {
                                        lb -= mindelta;
                                        temp0 = deltaCosts[current_scope_idx[i]][0];
                                        temp1 = deltaCosts[current_scope_idx[i]][1];
                                        deltaCosts[current_scope_idx[i]][1] = 0;
                                        deltaCosts[current_scope_idx[i]][0] = 0;
                                        scope[current_scope_idx[i]]->project(0, temp0 - mindelta, true);
                                        scope[current_scope_idx[i]]->project(1, temp1 - mindelta, true);
                                        scope[current_scope_idx[i]]->findSupport();
                                    }
                                } else {
                                    mindelta = deltaCosts[current_scope_idx[i]][scope[current_scope_idx[i]]->toIndex(scope[current_scope_idx[i]]->getValue())];
                                    if (mindelta != 0) {
                                        lb -= mindelta;
                                        deltaCosts[current_scope_idx[i]][1] = 0;
                                        deltaCosts[current_scope_idx[i]][0] = 0;
                                    }
                                }
                            } else {
                                for (int j = 0; j < nbValue[i]; ++j) {
                                    if (current_val_idx[i][j] != (int)AMO[VirtualVar[current_scope_idx[i]] - 1].size()) {
                                        assert(current_val_idx[i][j] != (int)AMO[VirtualVar[current_scope_idx[i]] - 1].size());
                                        if (scope[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][j]].first]->unassigned()) {
                                            mindelta = min(deltaCosts[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][j]].first][0], deltaCosts[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][j]].first][1]);
                                            if (mindelta != 0 && mindelta == deltaCosts[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][j]].first][1 - AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][j]].second]) {
                                                lb -= mindelta;
                                                temp0 = deltaCosts[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][j]].first][1 - AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][j]].second];
                                                deltaCosts[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][j]].first][1] = 0;
                                                deltaCosts[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][j]].first][0] = 0;
                                                scope[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][j]].first]->project(AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][j]].second, temp0 - mindelta, true);
                                                scope[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][j]].first]->findSupport();
                                            }
                                        } else {
                                            int varidx = AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][j]].first;
                                            mindelta = deltaCosts[varidx][scope[varidx]->toIndex(scope[varidx]->getValue())];
                                            if (mindelta != 0) {
                                                lb -= mindelta;
                                                deltaCosts[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][j]].first][1] = 0;
                                                deltaCosts[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][j]].first][0] = 0;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        //TODO: incremental bound propagation
                        propagate();
                        if (ToulBar2::FullEAC)
                            reviseEACGreedySolution();
                    }
                }
            } else {
                if (assigned[varIndex] == 1 && scope[varIndex]->assigned()) {
                    nonassigned = nonassigned - 1;
                    assigned[varIndex] = 2;
                    if (nonassigned <= NARYPROJECTIONSIZE && (nonassigned < 3 || maxInitDomSize <= NARYPROJECTION3MAXDOMSIZE)) {
                        if (connected()) {
                            deconnect(); // this constraint is removed from the current WCSP problem
                            projectNary();
                        }
                    }
                }
            }
        }
    }
    //Return True if a value has been deleted else return False
    bool BoundConsistency()
    {
        int k = 0, k2;
        bool b = false;
        int currentvar, currentval;
        while (k < carity && !b) {
            currentvar = current_scope_idx[k];
            //Determine if at least one value of the current variable can be erased.
            if (MaxWeight - weights[currentvar][GreatestWeightIdx[currentvar]] + weights[currentvar][LowestWeightIdx[currentvar]] < capacity) {
                k2 = 0;
                while (k2 < nbValue[k] && !b) {
                    currentval = current_val_idx[k][k2];
                    if (MaxWeight - weights[currentvar][GreatestWeightIdx[currentvar]] + weights[currentvar][currentval] < capacity) {
                        if (currentval == (int)VarVal[currentvar].size() - 1) {
                            if (VirtualVar[currentvar] == 0) {
                                for (unsigned int i = 0; i < NotVarVal[currentvar].size(); ++i) {
                                    if (scope[currentvar]->canbe(NotVarVal[currentvar][i]))
                                        scope[currentvar]->remove(NotVarVal[currentvar][i]);
                                }
                            } else {
                                lastval0[currentvar] = 0;
                                lastval0ok[currentvar] = true;
                                if (nbVirtualVar[VirtualVar[currentvar] - 1] == 1) {
                                    if (scope[AMO[VirtualVar[currentvar] - 1][current_val_idx[k][0]].first]->canbe(1 - AMO[VirtualVar[currentvar] - 1][current_val_idx[k][0]].second))
                                        scope[AMO[VirtualVar[currentvar] - 1][current_val_idx[k][0]].first]->remove(1 - AMO[VirtualVar[currentvar] - 1][current_val_idx[k][0]].second);
                                    b = true;
                                }
                            }
                        } else {
                            if (VirtualVar[currentvar] == 0) {
                                scope[currentvar]->remove(VarVal[currentvar][currentval]);
                            } else {
                                if (scope[AMO[VirtualVar[currentvar] - 1][currentval].first]->canbe(AMO[VirtualVar[currentvar] - 1][currentval].second))
                                    scope[AMO[VirtualVar[currentvar] - 1][currentval].first]->remove(AMO[VirtualVar[currentvar] - 1][currentval].second);
                                b = true;
                            }
                        }
                        if (!b) {
                            if (!connected())
                                b = true;
                            else if (VirtualVar[currentvar] == 0) {
                                if (!isunassigned(currentvar))
                                    b = true;
                            } else if (nbVirtualVar[VirtualVar[currentvar]] == 0)
                                b = true;
                        }
                    }
                    k2++;
                }
                b = true;
            }
            k++;
        }
        return b;
    }

    void ObjConsistency()
    {
        int k = 0, k2;
        int currentvar;
        Cost SumMin = MIN_COST;
        int OP = 0;
        if (AMO.empty()) {
            for (int i = 0; i < carity; ++i) {
                SumMin += deltaCosts[current_scope_idx[i]][LowestWeightIdx[current_scope_idx[i]]];
            }
            bool b = false;
            int size = (int)tempdeltaCosts.size();
            while (k < carity) {
                currentvar = current_scope_idx[k];
                assert(scope[currentvar]->canbe(VarVal[currentvar][LowestWeightIdx[currentvar]]));
                //Determine if at least one value of the current variable can be erased.
                k2 = 0;
                if (SumMin - deltaCosts[current_scope_idx[k]][LowestWeightIdx[current_scope_idx[k]]] + deltaCosts[current_scope_idx[k]][GreatestWeightIdx[current_scope_idx[k]]] > lb - assigneddeltas) {
                    while (k2 < nbValue[k]) {
                        if (OP < size && currentvar == tempdeltaCosts[OP][0] && current_val_idx[k][k2] == tempdeltaCosts[OP][1]) {
                            b = true;
                            OP++;
                        } else {
                            b = false;
                        }
                        if (SumMin - deltaCosts[current_scope_idx[k]][LowestWeightIdx[current_scope_idx[k]]] + deltaCosts[currentvar][current_val_idx[k][k2]] > lb - assigneddeltas) {
                            assert(SumMin - deltaCosts[current_scope_idx[k]][LowestWeightIdx[current_scope_idx[k]]] + deltaCosts[currentvar][current_val_idx[k][k2]] - lb + assigneddeltas > 0);
                            if (b) {
                                tempdeltaCosts[OP - 1][2] -= SumMin - deltaCosts[current_scope_idx[k]][LowestWeightIdx[current_scope_idx[k]]] + deltaCosts[currentvar][current_val_idx[k][k2]] - lb + assigneddeltas;
                            } else {
                                tempdeltaCosts.push_back({ current_scope_idx[k], current_val_idx[k][k2], -(SumMin - deltaCosts[current_scope_idx[k]][LowestWeightIdx[current_scope_idx[k]]] + deltaCosts[currentvar][current_val_idx[k][k2]] - lb + assigneddeltas) });
                            }
                            if (current_val_idx[k][k2] < (int)VarVal[currentvar].size() - 1) {
                                //Group_ProjectNVV(currentvar, SumMin - deltaCosts[current_scope_idx[k]][LowestWeightIdx[current_scope_idx[k]]] + deltaCosts[currentvar][current_val_idx[k][k2]]-lb+assigneddeltas);
                                deltaCosts[currentvar][current_val_idx[k][k2]] -= SumMin - deltaCosts[current_scope_idx[k]][LowestWeightIdx[current_scope_idx[k]]] + deltaCosts[currentvar][current_val_idx[k][k2]] - lb + assigneddeltas;
                            } else {
                                //scope[currentvar]->project(VarVal[currentvar][current_val_idx[k][k2]], SumMin - deltaCosts[current_scope_idx[k]][LowestWeightIdx[current_scope_idx[k]]] + deltaCosts[currentvar][current_val_idx[k][k2]]-lb+assigneddeltas, true);
                                deltaCosts[currentvar][current_val_idx[k][k2]] -= SumMin - deltaCosts[current_scope_idx[k]][LowestWeightIdx[current_scope_idx[k]]] + deltaCosts[currentvar][current_val_idx[k][k2]] - lb + assigneddeltas;
                            }
                            assert(SumMin - deltaCosts[current_scope_idx[k]][LowestWeightIdx[current_scope_idx[k]]] + deltaCosts[currentvar][current_val_idx[k][k2]] <= lb - assigneddeltas);
                        }
                        assert(SumMin - deltaCosts[current_scope_idx[k]][LowestWeightIdx[current_scope_idx[k]]] + deltaCosts[currentvar][current_val_idx[k][k2]] <= lb - assigneddeltas);
                        k2++;
                    }
                } else {
                    while (OP < size && tempdeltaCosts[OP][0] == currentvar)
                        OP++;
                }
                k++;
            }
        }
        for (int i = 0; i < (int)tempdeltaCosts.size(); ++i) {
            if (tempdeltaCosts[i][2] > 0) {
                if (CorrAMO[tempdeltaCosts[i][0]] > 0)
                    scope[tempdeltaCosts[i][0]]->extend(tempdeltaCosts[i][1], tempdeltaCosts[i][2]);
                else {
                    if (tempdeltaCosts[i][1] != (int)VarVal[tempdeltaCosts[i][0]].size() - 1)
                        scope[tempdeltaCosts[i][0]]->extend(VarVal[tempdeltaCosts[i][0]][tempdeltaCosts[i][1]], tempdeltaCosts[i][2]);
                    else
                        Group_extendNVV(tempdeltaCosts[i][0], tempdeltaCosts[i][2]);
                }
            } else {
                if (CorrAMO[tempdeltaCosts[i][0]] > 0)
                    scope[tempdeltaCosts[i][0]]->project(tempdeltaCosts[i][1], -tempdeltaCosts[i][2], true);
                else {
                    if (tempdeltaCosts[i][1] != (int)VarVal[tempdeltaCosts[i][0]].size() - 1)
                        scope[tempdeltaCosts[i][0]]->project(VarVal[tempdeltaCosts[i][0]][tempdeltaCosts[i][1]], -tempdeltaCosts[i][2], true);
                    else
                        Group_ProjectNVV(tempdeltaCosts[i][0], -tempdeltaCosts[i][2]);
                }
            }
        }
        for (int i = 0; i < carity; ++i) {
            if (VirtualVar[current_scope_idx[i]] == 0)
                scope[current_scope_idx[i]]->findSupport();
            else {
                for (int j = 0; j < nbValue[i]; ++j) {
                    if (current_val_idx[i][j] != (int)AMO[VirtualVar[current_scope_idx[i]] - 1].size())
                        scope[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][j]].first]->findSupport();
                }
            }
        }
    }

    bool ComputeProfit()
    {
        Cost verifopt = -lb + assigneddeltas; //Used to check if the last optimal solution has still a cost of 0
        Long verifweight = 0; //Used to check if the last optimal solution has still a cost of 0
        Double storec, storew; //Used to check if the last optimal solution has still a cost of 0
        Cost UN;
        bool diff0;
        int currentvar, currentval;
        for (int i = 0; i < carity; i++) {
            storec = 0;
            storew = 0;
            currentvar = current_scope_idx[i];
            UnaryCost0[currentvar] = MIN_COST;
            diff0 = true;
            if (Booleanvar[currentvar]) {
                Profit[currentvar][0] = scope[currentvar]->getCost(VarVal[currentvar][0]) + deltaCosts[currentvar][0];
                Profit[currentvar][1] = scope[currentvar]->getCost(VarVal[currentvar][1]) + deltaCosts[currentvar][1];
                UnaryCost0[currentvar] = Profit[currentvar][1] - deltaCosts[currentvar][1];
                storew += weights[currentvar][0] * OptSol[currentvar][0] + weights[currentvar][1] * OptSol[currentvar][1];
                storec += Profit[currentvar][1] * OptSol[currentvar][1] + Profit[currentvar][0] * OptSol[currentvar][0];
            } else {
                for (int j = 0; j < nbValue[i]; j++) {
                    currentval = current_val_idx[i][j];
                    if (VirtualVar[currentvar] == 0) {
                        assert(scope[currentvar]->getCost(VarVal[currentvar][currentval]) >= MIN_COST);
                        if (currentval == (int)VarVal[currentvar].size() - 1) {
                            if (!diff0 && !lastval0ok[currentvar] && scope[currentvar]->getCost(VarVal[currentvar].back()) > MIN_COST) {
                                UnaryCost0[currentvar] = MAX_COST;
                                for (unsigned int l = 0; l < NotVarVal[currentvar].size(); l++) {
                                    if (scope[currentvar]->canbe(NotVarVal[currentvar][l])) {
                                        assert(scope[currentvar]->getCost(NotVarVal[currentvar][l]) >= MIN_COST);
                                        UN = scope[currentvar]->getCost(NotVarVal[currentvar][l]);
                                        if (UN < UnaryCost0[currentvar]) {
                                            UnaryCost0[currentvar] = UN;
                                            if (UN == MIN_COST) {
                                                lastval0[currentvar] = NotVarVal[currentvar][l];
                                                VarVal[currentvar].back() = lastval0[currentvar];
                                                break;
                                            }
                                        }
                                    }
                                }
                            } else {
                                UnaryCost0[currentvar] = MIN_COST;
                            }
                            Profit[currentvar][currentval] = UnaryCost0[currentvar] + deltaCosts[currentvar][currentval];
                        } else {
                            Profit[currentvar][currentval] = scope[currentvar]->getCost(VarVal[currentvar][currentval]) + deltaCosts[currentvar][currentval];
                            if (Profit[currentvar][currentval] == deltaCosts[currentvar][currentval])
                                diff0 = false;
                        }
                    } else {
                        Profit[currentvar][currentval] = 0;
                        for (int k = 0; k < (int)AMO[VirtualVar[currentvar] - 1].size(); ++k) {
                            if (k == currentval) {
                                Profit[currentvar][currentval] += scope[AMO[VirtualVar[currentvar] - 1][k].first]->getCost(AMO[VirtualVar[currentvar] - 1][k].second) + deltaCosts[AMO[VirtualVar[currentvar] - 1][k].first][AMO[VirtualVar[currentvar] - 1][k].second];
                            } else {
                                if (scope[AMO[VirtualVar[currentvar] - 1][k].first]->unassigned()) {
                                    Profit[currentvar][currentval] += scope[AMO[VirtualVar[currentvar] - 1][k].first]->getCost(1 - AMO[VirtualVar[currentvar] - 1][k].second) + deltaCosts[AMO[VirtualVar[currentvar] - 1][k].first][1 - AMO[VirtualVar[currentvar] - 1][k].second];
                                }
                            }
                        }
                    }
                    storew += weights[currentvar][currentval] * OptSol[currentvar][currentval];
                    storec += Profit[currentvar][currentval] * OptSol[currentvar][currentval];
                }
            }
            //Compute the cost of the last optimal solution
            verifopt += Ceil(storec);
            verifweight += Ceil(storew);
            assert(*max_element(Profit[i].begin(), Profit[i].end()) < MAX_COST);
        }
        //assert(verifopt >= 0 || verifweight < capacity);
        if (verifopt > 0 || verifweight < capacity)
            return true;
        else
            return false;
    }

    //Return a vector containing the slopes, a slope is defined between 2 values of a variable : <Var, Val1, Val2, slopes>
    void ComputeSlopes(Long* W, Cost* c)
    {
        int item1 = 0;
        int k = 0, currentvar;
        Slopes.clear();
        for (int i = 0; i < carity; i++) {
            currentvar = current_scope_idx[i];
            fill(OptSol[currentvar].begin(), OptSol[currentvar].end(), 0);
            //Sort the value in ascending weight
            arrvar.clear();
            arrvar = current_val_idx[i];
            arrvar.resize(nbValue[i]);
            sort(arrvar.begin(), arrvar.end(),
                [&](int x, int y) {
                    if (weights[currentvar][x] == weights[currentvar][y]) {
                        if (Profit[currentvar][x] == Profit[currentvar][y]) {
                            if (VirtualVar[currentvar] == 0) {
                                if (scope[currentvar]->getSupport() == VarVal[currentvar][x]) {
                                    return true;
                                } else
                                    return false;
                            } else {
                                return false;
                            }
                        } else
                            return Profit[currentvar][x] > Profit[currentvar][y];
                    } else
                        return weights[currentvar][x] < weights[currentvar][y];
                });
            //Find the value with the heaviest weight
            k = (int)arrvar.size() - 1;
            item1 = arrvar[k];
            while (k > 0) {
                k--;
                //We don't consider dominated items : p_i < p_j && w_i > w_j   (i dominate j)
                if (Profit[currentvar][arrvar[k]] < Profit[currentvar][item1] && weights[currentvar][arrvar[k]] < weights[currentvar][item1]) {
                    //If it is the first slope for the current variable, we directy add it : <Var, val1, val2, slope>
                    //Else we compare the new slope with the precedent one to verify if the precedent value isn't dominated.
                    //If it is dominated, we replace the last slope else we add a new one
                    if (Slopes.empty() || Slopes.back()[0] != currentvar || Slopes.back()[3] >= Double((Profit[currentvar][item1] - Profit[currentvar][arrvar[k]])) / (weights[currentvar][item1] - weights[currentvar][arrvar[k]]))
                        Slopes.push_back({ Double(currentvar), Double(arrvar[k]), Double(item1), Double((Profit[currentvar][item1] - Profit[currentvar][arrvar[k]])) / (weights[currentvar][item1] - weights[currentvar][arrvar[k]]) });
                    else {
                        Slopes.back() = { Double(currentvar), Double(arrvar[k]), Slopes.back()[2], Double((Profit[currentvar][Slopes.back()[2]] - Profit[currentvar][arrvar[k]])) / (weights[currentvar][Slopes.back()[2]] - weights[currentvar][arrvar[k]]) };
                        while (Slopes.size() > 1 && Slopes.end()[-2][0] == Slopes.back()[0] && Slopes.back()[3] >= Slopes.end()[-2][3]) {
                            Slopes.end()[-2] = { Double(currentvar), Slopes.back()[1], Slopes.end()[-2][2], Double((Profit[currentvar][Slopes.end()[-2][2]] - Profit[currentvar][Slopes.back()[1]])) / (weights[currentvar][Slopes.end()[-2][2]] - weights[currentvar][Slopes.back()[1]]) };
                            Slopes.pop_back();
                        }
                    }
                    item1 = arrvar[k];
                }
            }
            //Compute a first Solution
            if (Slopes.size() == 0 || Slopes.back()[0] != currentvar) {
                OptSol[currentvar][item1] = 1;
                *W += weights[currentvar][item1];
                *c += Profit[currentvar][item1];
            } else {
                *W += weights[currentvar][Slopes.back()[1]];
                OptSol[currentvar][Slopes.back()[1]] = 1;
                *c += Profit[currentvar][Slopes.back()[1]];
            }
        }
    }

    //Find the optimal solution. Follow the order of the slopes and modify OptSol until we fill the constraint
    void FindOpt(vector<std::array<Double, 4>>& Slopes, Long* W, Cost* c, Double* xk, int* iter)
    {
        int currentVar;
        Long capacityLeft;
        while (*W < capacity) {
            currentVar = int(Slopes[*iter][0]);
            if (*W + weights[currentVar][Slopes[*iter][2]] - weights[currentVar][Slopes[*iter][1]] >= capacity) {
                capacityLeft = capacity - *W;
                *xk = Double(capacityLeft) / (weights[currentVar][Slopes[*iter][2]] - weights[currentVar][Slopes[*iter][1]]);
                OptSol[currentVar][Slopes[*iter][2]] = *xk;
                OptSol[currentVar][Slopes[*iter][1]] = 1 - *xk;
                *W += weights[currentVar][Slopes[*iter][2]] - weights[currentVar][Slopes[*iter][1]];
                *c += Ceil(*xk * (Profit[currentVar][Slopes[*iter][2]] - Profit[currentVar][Slopes[*iter][1]]));
                assert(capacityLeft > 0);
            } else {
                assert(OptSol[currentVar][Slopes[*iter][1]] == 1);
                OptSol[currentVar][Slopes[*iter][1]] = 0;
                OptSol[currentVar][Slopes[*iter][2]] = 1;
                *W += weights[currentVar][Slopes[*iter][2]] - weights[currentVar][Slopes[*iter][1]];
                *c += Profit[currentVar][Slopes[*iter][2]] - Profit[currentVar][Slopes[*iter][1]];
                *iter = *iter + 1;
            }
        }
    }

    //Do the Extension/Projection
    void ExtensionProjection(vector<Double>& y_i, Double y_cc, vector<vector<Double>> yAMO_i, vector<Double> clq)
    {
        int n = 0, currentvar, currentval;
        Cost C;
        tempdeltaCosts.clear();
        for (int i = 0; i < carity; i++) {
            currentvar = current_scope_idx[i];
            if (VirtualVar[currentvar] == 0) {
                for (int j = 0; j < nbValue[i]; ++j) {
                    currentval = current_val_idx[i][j];
                    if (OptSol[currentvar][currentval] > 0) {
                        if (currentval == (int)VarVal[currentvar].size() - 1) {
                            if (UnaryCost0[currentvar] > MIN_COST) {
                                tempdeltaCosts.push_back({ current_scope_idx[i], current_val_idx[i][j], UnaryCost0[currentvar] });
                                //Group_extendNVV(currentvar, UnaryCost0[currentvar]);
                                deltaCosts[currentvar][currentval] += UnaryCost0[currentvar];
                            }
                        } else {
                            C = scope[currentvar]->getCost(VarVal[currentvar][currentval]);
                            if (C > MIN_COST) {
                                tempdeltaCosts.push_back({ current_scope_idx[i], current_val_idx[i][j], C });
                                deltaCosts[currentvar][currentval] += C;
                                //scope[currentvar]->extend(VarVal[currentvar][currentval],C);
                            }
                        }
                    } else {
                        C = Ceil(-deltaCosts[currentvar][currentval] + y_i[i] + y_cc * weights[currentvar][currentval]);
                        if (C != MIN_COST) {
                            tempdeltaCosts.push_back({ current_scope_idx[i], current_val_idx[i][j], C });
                            deltaCosts[currentvar][currentval] += C;
                            //   ExtOrProJ(currentvar, currentval,C);
                        }
                    }
                }
                //scope[currentvar]->findSupport();
            } else {
                for (int j = 0; j < nbValue[i]; ++j) {
                    currentval = current_val_idx[i][j];
                    if (currentval != (int)AMO[VirtualVar[currentvar] - 1].size()) {
                        if (OptSol[currentvar][currentval] == 1) {
                            C = scope[AMO[VirtualVar[currentvar] - 1][currentval].first]->getCost(AMO[VirtualVar[currentvar] - 1][currentval].second);
                            if (C > MIN_COST) {
                                tempdeltaCosts.push_back({ AMO[VirtualVar[currentvar] - 1][currentval].first, AMO[VirtualVar[currentvar] - 1][currentval].second, C });
                                deltaCosts[AMO[VirtualVar[currentvar] - 1][currentval].first][AMO[VirtualVar[currentvar] - 1][currentval].second] += C;
                                // scope[AMO[VirtualVar[currentvar] - 1][currentval].first]->extend(AMO[VirtualVar[currentvar] - 1][currentval].second,C);
                            }
                            C = Ceil(-deltaCosts[AMO[VirtualVar[currentvar] - 1][currentval].first][1 - AMO[VirtualVar[currentvar] - 1][currentval].second] + yAMO_i[n][j] + y_cc * Original_weigths[AMO[VirtualVar[currentvar] - 1][currentval].first][1 - AMO[VirtualVar[currentvar] - 1][currentval].second]);
                            deltaCosts[AMO[VirtualVar[currentvar] - 1][currentval].first][1 - AMO[VirtualVar[currentvar] - 1][currentval].second] += C;
                            if (C > MIN_COST) {
                                assert(scope[AMO[VirtualVar[currentvar] - 1][currentval].first]->getCost(1 - AMO[VirtualVar[currentvar] - 1][currentval].second) >= C);
                                tempdeltaCosts.push_back({ AMO[VirtualVar[currentvar] - 1][currentval].first, 1 - AMO[VirtualVar[currentvar] - 1][currentval].second, C });
                                // scope[AMO[VirtualVar[currentvar] - 1][currentval].first]->extend(1 - AMO[VirtualVar[currentvar] - 1][currentval].second, C);
                            } else if (C < MIN_COST) {
                                tempdeltaCosts.push_back({ AMO[VirtualVar[currentvar] - 1][currentval].first, 1 - AMO[VirtualVar[currentvar] - 1][currentval].second, C });
                                //scope[AMO[VirtualVar[currentvar] - 1][currentval].first]->project(1 - AMO[VirtualVar[currentvar] - 1][currentval].second, -C,true);
                            }
                        } else if (OptSol[currentvar][currentval] == 0) {
                            C = scope[AMO[VirtualVar[currentvar] - 1][currentval].first]->getCost(1 - AMO[VirtualVar[currentvar] - 1][currentval].second);
                            if (C > MIN_COST) {
                                tempdeltaCosts.push_back({ AMO[VirtualVar[currentvar] - 1][currentval].first, 1 - AMO[VirtualVar[currentvar] - 1][currentval].second, C });
                                deltaCosts[AMO[VirtualVar[currentvar] - 1][currentval].first][1 - AMO[VirtualVar[currentvar] - 1][currentval].second] += C;
                                //scope[AMO[VirtualVar[currentvar] - 1][currentval].first]->extend(1 - AMO[VirtualVar[currentvar] - 1][currentval].second,C);
                            }
                            C = Ceil(-deltaCosts[AMO[VirtualVar[currentvar] - 1][currentval].first][AMO[VirtualVar[currentvar] - 1][currentval].second] + yAMO_i[n][j] + y_cc * Original_weigths[AMO[VirtualVar[currentvar] - 1][currentval].first][AMO[VirtualVar[currentvar] - 1][currentval].second] + clq[n]);
                            deltaCosts[AMO[VirtualVar[currentvar] - 1][currentval].first][AMO[VirtualVar[currentvar] - 1][currentval].second] += C;
                            if (C > MIN_COST) {
                                tempdeltaCosts.push_back({ AMO[VirtualVar[currentvar] - 1][currentval].first, AMO[VirtualVar[currentvar] - 1][currentval].second, C });
                                assert(scope[AMO[VirtualVar[currentvar] - 1][currentval].first]->getCost(AMO[VirtualVar[currentvar] - 1][currentval].second) >= C);
                                //scope[AMO[VirtualVar[currentvar] - 1][currentval].first]->extend(AMO[VirtualVar[currentvar] - 1][currentval].second, C);
                            } else if (C < MIN_COST) {
                                tempdeltaCosts.push_back({ AMO[VirtualVar[currentvar] - 1][currentval].first, AMO[VirtualVar[currentvar] - 1][currentval].second, C });
                                // scope[AMO[VirtualVar[currentvar] - 1][currentval].first]->project(AMO[VirtualVar[currentvar] - 1][currentval].second, -C,true);
                            }
                        } else {
                            C = scope[AMO[VirtualVar[currentvar] - 1][currentval].first]->getCost(AMO[VirtualVar[currentvar] - 1][currentval].second);
                            if (C > MIN_COST) {
                                deltaCosts[AMO[VirtualVar[currentvar] - 1][currentval].first][AMO[VirtualVar[currentvar] - 1][currentval].second] += C;
                                tempdeltaCosts.push_back({ AMO[VirtualVar[currentvar] - 1][currentval].first, AMO[VirtualVar[currentvar] - 1][currentval].second, C });
                                // scope[AMO[VirtualVar[currentvar] - 1][currentval].first]->extend(AMO[VirtualVar[currentvar] - 1][currentval].second,C);
                            }
                            C = scope[AMO[VirtualVar[currentvar] - 1][currentval].first]->getCost(1 - AMO[VirtualVar[currentvar] - 1][currentval].second);
                            if (C > MIN_COST) {
                                tempdeltaCosts.push_back({ AMO[VirtualVar[currentvar] - 1][currentval].first, 1 - AMO[VirtualVar[currentvar] - 1][currentval].second, C });
                                deltaCosts[AMO[VirtualVar[currentvar] - 1][currentval].first][1 - AMO[VirtualVar[currentvar] - 1][currentval].second] += C;
                                // scope[AMO[VirtualVar[currentvar] - 1][currentval].first]->extend(1 - AMO[VirtualVar[currentvar] - 1][currentval].second,C);
                            }
                        }
                        //scope[AMO[VirtualVar[currentvar] - 1][currentval].first]->findSupport();
                    }
                }
                n++;
            }
        }
    }

    Cost max_value(const vector<vector<Long>>& W,
        const vector<vector<Cost>>& P,
        Long Cap, Long MW, int removed)
    {
        if (Cap > 0) {
            vector<int> currOPT;
            if (W.empty())
                return 0;
            assert(MW - Cap + 1 > 0);
            vector<pair<Cost, vector<int>>> last(MW - Cap + 1, make_pair(MAX_COST, currOPT));
            last[0].first = 0;
            for (unsigned int i = 0; i < W.size(); ++i) {
                currOPT.push_back(distance(W[i].begin(), max_element(W[i].begin(), W[i].end())));
                last[0].first += P[i][currOPT[i]];
            }
            last[0].second = currOPT;
            vector<pair<Cost, vector<int>>> current = last;
            for (unsigned int i = 0; i < W.size(); ++i) {
                for (unsigned int j = 0; j < W[i].size(); ++j) {
                    for (int k = MW - Cap; k >= 0; --k) {
                        assert(k < (int)last.size());
                        if (last[k].first < MAX_COST && W[i][j] <= W[i][last[k].second[i]] && Cap <= MW - k + W[i][j] - W[i][last[k].second[i]] && P[i][j] < P[i][last[k].second[i]]) {
                            assert(k - W[i][j] + W[i][last[k].second[i]] < (int)current.size());
                            if (current[k - W[i][j] + W[i][last[k].second[i]]].first > last[k].first - P[i][last[k].second[i]] + P[i][j]) {
                                current[k - W[i][j] + W[i][last[k].second[i]]].first = last[k].first - P[i][last[k].second[i]] + P[i][j];
                                current[k - W[i][j] + W[i][last[k].second[i]]].second = last[k].second;
                                current[k - W[i][j] + W[i][last[k].second[i]]].second[i] = j;
                            }
                        }
                    }
                }
                last = current;
            }
            Cost optc = MAX_COST;
            int lp = 0;
            for (unsigned int i = 0; i < last.size(); ++i) {
                if (optc > last[i].first) {
                    optc = last[i].first;
                    lp = i;
                }
            }
            for (int i = 0; i < carity; ++i) {
                if (i != removed) {
                    fill(OptSol[current_scope_idx[i]].begin(), OptSol[current_scope_idx[i]].end(), 0);
                    if (i < removed) {
                        OptSol[current_scope_idx[i]][current_val_idx[i][last[lp].second[i]]] = 1;
                    } else {
                        OptSol[current_scope_idx[i]][current_val_idx[i][last[lp].second[i - 1]]] = 1;
                    }
                }
            }
            return optc;
        } else {
            Cost NegDel = 0;
            for (unsigned int i = 0; i < P.size(); ++i) {
                NegDel += *min_element(P[i].begin(), P[i].end());
            }
            return NegDel;
        }
    }
    void propagate() override
    {
        if (ToulBar2::dumpWCSP % 2) // do not propagate if problem is dumped before preprocessing
            return;
        if (ToulBar2::interrupted) {
            throw TimeOut();
        }
        // propagates from scratch the constraint
        if (connected()) {
            bool b = false;
            for (int i = 0; connected() && i < arity_; i++) {
                if (CorrAMO[i] == 0) {
                    if (assigned[i] == 0 && !isunassigned(i)) {
                        assign(i);
                        b = true;
                    } /* indique à la contrainte que cette  variable est affectée (donc met à jour nonassigned..) */
                    else
                        assert(assigned[i] > 0 || scope[i]->unassigned());
                } else {
                    if (assigned[i] == 0 && scope[i]->assigned()) {
                        assign(i);
                        b = true;
                    }
                }
            }
            if (connected() && !b) {
                if (!verify()) {
                    THROWCONTRADICTION;
                } else if (nonassigned > 3 && ToulBar2::LcLevel >= LC_AC) {
                    get_current_scope();
#ifndef NDEBUG
                    for (int i = 0; i < carity; ++i) {
                        if (VirtualVar[current_scope_idx[i]] == 0) {
                            assert(scope[current_scope_idx[i]]->canbe(VarVal[current_scope_idx[i]].back()) || lastval0ok[current_scope_idx[i]]);
                            assert(scope[current_scope_idx[i]]->unassigned());
                            assert(assigned[current_scope_idx[i]] == 0);
                        }
                        for (int j = 0; j < nbValue[i]; ++j) {
                            if (VirtualVar[current_scope_idx[i]] == 0) {
                                assert(scope[current_scope_idx[i]]->canbe(
                                    VarVal[current_scope_idx[i]][current_val_idx[i][j]]));
                                assert(scope[current_scope_idx[i]]->canbe(
                                    VarVal[current_scope_idx[i]][GreatestWeightIdx[current_scope_idx[i]]]));
                                assert(scope[current_scope_idx[i]]->canbe(
                                    VarVal[current_scope_idx[i]][LowestWeightIdx[current_scope_idx[i]]]));
                            }
                            assert(weights[current_scope_idx[i]][GreatestWeightIdx[current_scope_idx[i]]] >= weights[current_scope_idx[i]][current_val_idx[i][j]]);
                            assert(weights[current_scope_idx[i]][LowestWeightIdx[current_scope_idx[i]]] <= weights[current_scope_idx[i]][current_val_idx[i][j]]);
                        }
                        assert(nbValue[i] > 1);
                    }
#endif
                    //Bound propagation, return true if a variable has been assigned
                    b = BoundConsistency();
                    if (!b && !ToulBar2::addAMOConstraints_ && ToulBar2::LcLevel == LC_EDAC) {
#ifndef NDEBUG
                        for (int i = 0; i < carity; ++i) {
                            if (VirtualVar[current_scope_idx[i]] == 0) {
                                assert(scope[current_scope_idx[i]]->canbe(VarVal[current_scope_idx[i]].back()) || lastval0ok[current_scope_idx[i]]);
                                assert(scope[current_scope_idx[i]]->unassigned());
                                assert(assigned[current_scope_idx[i]] == 0);
                            }
                            for (int j = 0; j < nbValue[i]; ++j) {
                                assert(MaxWeight - weights[current_scope_idx[i]][GreatestWeightIdx[current_scope_idx[i]]] + weights[current_scope_idx[i]][current_val_idx[i][j]] >= capacity);
                                if (VirtualVar[current_scope_idx[i]] == 0) {
                                    assert(scope[current_scope_idx[i]]->canbe(
                                        VarVal[current_scope_idx[i]][current_val_idx[i][j]]));
                                    assert(scope[current_scope_idx[i]]->canbe(
                                        VarVal[current_scope_idx[i]][GreatestWeightIdx[current_scope_idx[i]]]));
                                    assert(scope[current_scope_idx[i]]->canbe(
                                        VarVal[current_scope_idx[i]][LowestWeightIdx[current_scope_idx[i]]]));
                                }
                                assert(weights[current_scope_idx[i]][GreatestWeightIdx[current_scope_idx[i]]] >= weights[current_scope_idx[i]][current_val_idx[i][j]]);
                                assert(weights[current_scope_idx[i]][LowestWeightIdx[current_scope_idx[i]]] <= weights[current_scope_idx[i]][current_val_idx[i][j]]);
                            }
                            assert(nbValue[i] > 1);
                        }
#endif
                        // Return True if the last optimal solution has a cost of 0
                        b = ComputeProfit();

                        if (connected() && b) {
                            Long W = 0;
                            Cost c = -lb + assigneddeltas;
                            ComputeSlopes(&W, &c); // temporary data structure for propagate

                            if (ToulBar2::verbose >= 4) {
                                cout << "cap is " << capacity << endl;
                                for (int i = 0; i < carity; ++i) {
                                    cout << scope[current_scope_idx[i]]->getName();
                                    for (unsigned int j = 0; j < VarVal[current_scope_idx[i]].size(); ++j) {
                                        cout << " " << j << " : " << Profit[current_scope_idx[i]][current_val_idx[i][j]] << "/"
                                             << weights[current_scope_idx[i]][current_val_idx[i][j]];
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
                            Long Sumw = 0;
                            for (int i = 0; i < carity; ++i) {
                                Sumw += weights[current_scope_idx[i]][GreatestWeightIdx[current_scope_idx[i]]];
                            }
                            assert(Sumw == MaxWeight);
                            assert(c < MAX_COST);
                            assert(W >= 0);
#endif
                            int iter = 0;
                            Double xk = 0;
                            if (ToulBar2::knapsackDP > -2 && NbNodes != wcsp->getNbNodes()) {
                                NbNodes = wcsp->getNbNodes();
                                if (NbNodes % ToulBar2::knapsackDP == 0)
                                    DoDyn = true;
                            }
                            if (ToulBar2::knapsackDP > -2 && DoDyn && W < capacity) {
                                Cost c2;
                                vector<vector<Long>> WeightforDyn, NewWeightforDyn;
                                vector<vector<Cost>> ProfforDyn, NewProfforDyn;
                                vector<Long> TmpWeightforDyn;
                                vector<Cost> TmpProfforDyn;
                                for (int i = 0; i < carity; ++i) {
                                    TmpWeightforDyn.clear();
                                    TmpProfforDyn.clear();
                                    for (int j = 0; j < nbValue[i]; ++j) {
                                        TmpWeightforDyn.push_back(weights[current_scope_idx[i]][current_val_idx[i][j]]);
                                        TmpProfforDyn.push_back(Profit[current_scope_idx[i]][current_val_idx[i][j]]);
                                    }
                                    WeightforDyn.push_back(TmpWeightforDyn);
                                    ProfforDyn.push_back(TmpProfforDyn);
                                }
                                c = max_value(WeightforDyn, ProfforDyn, capacity, MaxWeight, carity + 1) - lb + assigneddeltas;
                                assert(c > -1);
                                if (c > 0) {
                                    vector<vector<Double>> RealOptSol = OptSol;
                                    vector<int> EDACORDER;
                                    for (int i = 0; i < carity; ++i) {
                                        EDACORDER.push_back(i);
                                    }
                                    sort(EDACORDER.begin(), EDACORDER.end(),
                                        [&](int& x, int& y) {
                                            return scope[current_scope_idx[x]]->getDACOrder() < scope[current_scope_idx[y]]->getDACOrder(); //TODO: checks if it favors aborbing more unary costs for the last variables in the DAC order or the opposite?!
                                        });

                                    for (int i = 0; i < carity; ++i) {
                                        if (ToulBar2::interrupted) {
                                            throw TimeOut();
                                        }
                                        NewProfforDyn = ProfforDyn;
                                        NewWeightforDyn = WeightforDyn;
                                        NewWeightforDyn.erase(NewWeightforDyn.begin() + EDACORDER[i]);
                                        NewProfforDyn.erase(NewProfforDyn.begin() + EDACORDER[i]);
                                        int j = 0;
                                        while (RealOptSol[current_scope_idx[EDACORDER[i]]][current_val_idx[EDACORDER[i]][j]] == 0)
                                            j++;
                                        OptSol[current_scope_idx[EDACORDER[i]]][current_val_idx[EDACORDER[i]][j]] = 0;
                                        for (int k = 0; k < nbValue[EDACORDER[i]]; ++k) {
                                            if (j != k) {
                                                OptSol[current_scope_idx[EDACORDER[i]]][current_val_idx[EDACORDER[i]][k]] = 1;
                                                c2 = max_value(NewWeightforDyn, NewProfforDyn, capacity - weights[current_scope_idx[EDACORDER[i]]][current_val_idx[EDACORDER[i]][k]], MaxWeight - weights[current_scope_idx[EDACORDER[i]]][GreatestWeightIdx[current_scope_idx[EDACORDER[i]]]], EDACORDER[i]) - lb + assigneddeltas + Profit[current_scope_idx[EDACORDER[i]]][current_val_idx[EDACORDER[i]][k]];
                                                assert(c2 >= c);
                                                Cost GAP = c2 - c;
                                                assert(GAP >= 0);
                                                Cost totrans;
                                                if (current_val_idx[EDACORDER[i]][k] < (int)VarVal[current_scope_idx[EDACORDER[i]]].size() - 1) {
                                                    totrans = scope[current_scope_idx[EDACORDER[i]]]->getCost(VarVal[current_scope_idx[EDACORDER[i]]][current_val_idx[EDACORDER[i]][k]]) - GAP;
                                                    if (totrans < MIN_COST) {
                                                        scope[current_scope_idx[EDACORDER[i]]]->project(VarVal[current_scope_idx[EDACORDER[i]]][current_val_idx[EDACORDER[i]][k]], -totrans, true);
                                                    } else if (totrans > MIN_COST) {
                                                        scope[current_scope_idx[EDACORDER[i]]]->extend(VarVal[current_scope_idx[EDACORDER[i]]][current_val_idx[EDACORDER[i]][k]], totrans);
                                                    }
                                                    deltaCosts[current_scope_idx[EDACORDER[i]]][current_val_idx[EDACORDER[i]][k]] += totrans;
                                                    ProfforDyn[EDACORDER[i]][k] = deltaCosts[current_scope_idx[EDACORDER[i]]][current_val_idx[EDACORDER[i]][k]];
                                                } else {
                                                    totrans = UnaryCost0[current_scope_idx[EDACORDER[i]]] - GAP;
                                                    if (totrans < MIN_COST) {
                                                        Group_ProjectNVV(current_scope_idx[EDACORDER[i]], -totrans);
                                                    } else if (totrans > MIN_COST) {
                                                        Group_extendNVV(current_scope_idx[EDACORDER[i]], totrans);
                                                    }
                                                    deltaCosts[current_scope_idx[EDACORDER[i]]].back() += totrans;
                                                    ProfforDyn[EDACORDER[i]].back() = deltaCosts[current_scope_idx[EDACORDER[i]]].back();
                                                }
                                                OptSol[current_scope_idx[EDACORDER[i]]][current_val_idx[EDACORDER[i]][k]] = 0;
                                            } else {
                                                if (current_val_idx[EDACORDER[i]][k] < (int)VarVal[current_scope_idx[EDACORDER[i]]].size() - 1) {
                                                    Cost C = scope[current_scope_idx[EDACORDER[i]]]->getCost(VarVal[current_scope_idx[EDACORDER[i]]][current_val_idx[EDACORDER[i]][k]]);
                                                    if (C > MIN_COST) {
                                                        deltaCosts[current_scope_idx[EDACORDER[i]]][current_val_idx[EDACORDER[i]][k]] += C;
                                                        scope[current_scope_idx[EDACORDER[i]]]->extend(VarVal[current_scope_idx[EDACORDER[i]]][current_val_idx[EDACORDER[i]][k]], C);
                                                    }
                                                } else {
                                                    if (UnaryCost0[current_scope_idx[EDACORDER[i]]] > 0) {
                                                        Group_extendNVV(current_scope_idx[EDACORDER[i]], UnaryCost0[current_scope_idx[EDACORDER[i]]]);
                                                        deltaCosts[current_scope_idx[EDACORDER[i]]].back() += UnaryCost0[current_scope_idx[EDACORDER[i]]];
                                                    }
                                                }
                                            }
                                        }
                                        OptSol[current_scope_idx[EDACORDER[i]]][current_val_idx[EDACORDER[i]][j]] = 1;
                                        scope[current_scope_idx[EDACORDER[i]]]->findSupport();
                                    }
                                    OptSol = RealOptSol;
                                    projectLB(c);
#ifndef NDEBUG
                                    for (int i = 0; i < carity; ++i) {
                                        for (int j = 0; j < nbValue[i]; ++j) {
                                            ProfforDyn[i][j] = deltaCosts[current_scope_idx[i]][current_val_idx[i][j]];
                                        }
                                    }

                                    c2 = max_value(WeightforDyn, ProfforDyn, capacity, MaxWeight, carity + 1);
                                    assert(c2 == lb - assigneddeltas);
                                    for (int i = 0; i < carity; ++i) {
                                        NewProfforDyn = ProfforDyn;
                                        NewWeightforDyn = WeightforDyn;
                                        NewWeightforDyn.erase(NewWeightforDyn.begin() + i);
                                        NewProfforDyn.erase(NewProfforDyn.begin() + i);
                                        for (int j = 0; j < nbValue[i]; ++j) {
                                            c2 = max_value(NewWeightforDyn, NewProfforDyn, capacity - weights[current_scope_idx[i]][current_val_idx[i][j]], MaxWeight - weights[current_scope_idx[i]][GreatestWeightIdx[current_scope_idx[i]]], i) + deltaCosts[current_scope_idx[i]][current_val_idx[i][j]];
                                            assert(c2 == lb - assigneddeltas);
                                        }
                                    }
#endif
                                }
                                DoDyn = false;
                            } else {
                                if (W < capacity) {
                                    //Sort the Slopes
                                    sort(Slopes.begin(), Slopes.end(),
                                        [&](std::array<Double, 4>& x, std::array<Double, 4>& y) {
                                            if (x[3] == y[3]) {
                                                if (x[0] == y[0])
                                                    return weights[int(x[0])][int(x[1])] <= weights[int(y[0])][int(y[1])];
                                                else
                                                    return scope[int(x[0])]->getDACOrder() < scope[int(
                                                                                                       y[0])]
                                                                                                 ->getDACOrder(); //TODO: checks if it favors aborbing more unary costs for the last variables in the DAC order or the opposite?!
                                            } else
                                                return x[3] < y[3];
                                        });
                                    //Find the optimal solution
                                    FindOpt(Slopes, &W, &c, &xk, &iter);
                                }
                                if (c > 0) {
                                    assert(W >= capacity);
                                    assert(xk >= 0 && xk <= 1);
                                    assert(iter <= (int)Slopes.size());
                                    assert(c > -1);
                                    //Compute the dual variable y_cc and y_i using the optimal primal solution
                                    Double y_cc = 0;
                                    if (!Slopes.empty())
                                        y_cc = Slopes[iter][3];
                                    if (xk == 0)
                                        y_cc = 0;
                                    y_i.clear();
                                    yAMO_i.clear();
                                    assert(y_cc >= MIN_COST);
                                    int k = 0, currentvar;
                                    vector<Double> clq;
                                    int optva = -1;
                                    for (int i = 0; i < carity; i++) {
                                        k = 0;
                                        currentvar = current_scope_idx[i];
                                        if (VirtualVar[currentvar] == 0) {
                                            while (OptSol[currentvar][current_val_idx[i][k]] == 0)
                                                k++;
                                            y_i.push_back(Profit[currentvar][current_val_idx[i][k]] - y_cc * weights[currentvar][current_val_idx[i][k]]);
                                            assert(y_i[i] < MAX_COST);
                                        } else {
                                            tempAMOy_i.clear();
                                            optva = -1;
                                            clq.push_back(0);
                                            for (int j = 0; j < nbValue[i]; ++j) {
                                                if (current_val_idx[i][j] != (int)AMO[VirtualVar[currentvar] - 1].size()) {
                                                    if (OptSol[currentvar][current_val_idx[i][j]] == 1) {
                                                        tempAMOy_i.push_back(0);
                                                        optva = j;
                                                    } else {
                                                        tempAMOy_i.push_back(scope[AMO[VirtualVar[currentvar] - 1][current_val_idx[i][j]].first]->getCost(1 - AMO[VirtualVar[currentvar] - 1][current_val_idx[i][j]].second) + deltaCosts[AMO[VirtualVar[currentvar] - 1][current_val_idx[i][j]].first][1 - AMO[VirtualVar[currentvar] - 1][current_val_idx[i][j]].second] - y_cc * Original_weigths[AMO[VirtualVar[currentvar] - 1][current_val_idx[i][j]].first][1 - AMO[VirtualVar[currentvar] - 1][current_val_idx[i][j]].second]);
                                                        if (clq.back() > scope[AMO[VirtualVar[currentvar] - 1][current_val_idx[i][j]].first]->getCost(AMO[VirtualVar[currentvar] - 1][current_val_idx[i][j]].second) + deltaCosts[AMO[VirtualVar[currentvar] - 1][current_val_idx[i][j]].first][AMO[VirtualVar[currentvar] - 1][current_val_idx[i][j]].second] - y_cc * Original_weigths[AMO[VirtualVar[currentvar] - 1][current_val_idx[i][j]].first][AMO[VirtualVar[currentvar] - 1][current_val_idx[i][j]].second] - tempAMOy_i.back())
                                                            clq.back() = scope[AMO[VirtualVar[currentvar] - 1][current_val_idx[i][j]].first]->getCost(AMO[VirtualVar[currentvar] - 1][current_val_idx[i][j]].second) + deltaCosts[AMO[VirtualVar[currentvar] - 1][current_val_idx[i][j]].first][AMO[VirtualVar[currentvar] - 1][current_val_idx[i][j]].second] - y_cc * Original_weigths[AMO[VirtualVar[currentvar] - 1][current_val_idx[i][j]].first][AMO[VirtualVar[currentvar] - 1][current_val_idx[i][j]].second] - tempAMOy_i.back();
                                                    }
                                                }
                                            }
                                            if (optva > -1) {
                                                tempAMOy_i[optva] = scope[AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].first]->getCost(AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].second) + deltaCosts[AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].first][AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].second] - y_cc * Original_weigths[AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].first][AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].second] - clq.back();
                                                if (tempAMOy_i[optva] + y_cc * Original_weigths[AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].first][1 - AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].second] - scope[AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].first]->getCost(1 - AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].second) - deltaCosts[AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].first][1 - AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].second] > epsilon) {
                                                    tempAMOy_i[optva] = scope[AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].first]->getCost(1 - AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].second) + deltaCosts[AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].first][1 - AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].second] - y_cc * Original_weigths[AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].first][1 - AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].second];
                                                    clq.back() = scope[AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].first]->getCost(AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].second) + deltaCosts[AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].first][AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].second] - y_cc * Original_weigths[AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].first][AMO[VirtualVar[currentvar] - 1][current_val_idx[i][optva]].second] - tempAMOy_i[optva];
                                                }
                                            }
                                            yAMO_i.push_back(tempAMOy_i);
                                        }
                                    }
                                    //Use y_cc and y_i to extend/project the right cost
                                    ExtensionProjection(y_i, y_cc, yAMO_i, clq);
                                    assert(c > 0);
                                    projectLB(c);
                                    ObjConsistency();
                                }
                            }
                            if (ToulBar2::verbose >= 4) {
                                cout << "projected cost " << c << " LB : " << lb - assigneddeltas << endl;
                                for (int i = 0; i < arity_; ++i) {
                                    if (scope[i]->unassigned()) {
                                        cout << " DeltaCosts " << scope[i]->getName() << " : ";
                                        for (unsigned int j = 0; j < deltaCosts[i].size(); ++j) {
                                            cout << " " << deltaCosts[i][j];
                                        }
                                        cout << endl;
                                    }
                                }
                            }
                            assert(getMaxFiniteCost() >= 0);
                        }
                    }
                } else {
                    get_current_scope();
                    if (universal() && connected()) {
                        deconnect();
                        wcsp->revise(this);
                        Cost TobeProjected = -lb + assigneddeltas;
                        lb = MIN_COST;
                        Cost mindelta;
                        for (int i = 0; i < carity; i++) {
                            mindelta = MAX_COST;
                            if (VirtualVar[current_scope_idx[i]] == 0) {
                                for (int j = 0; j < nbValue[i]; ++j) {
                                    if (mindelta > deltaCosts[current_scope_idx[i]][current_val_idx[i][j]])
                                        mindelta = deltaCosts[current_scope_idx[i]][current_val_idx[i][j]];
                                }
                                TobeProjected += mindelta;
                                for (int j = 0; j < nbValue[i]; ++j) {
                                    assert(mindelta <= deltaCosts[current_scope_idx[i]][current_val_idx[i][j]]);
                                    if (current_val_idx[i][j] == (int)VarVal[current_scope_idx[i]].size() - 1) {
                                        Group_ProjectNVV(current_scope_idx[i], deltaCosts[current_scope_idx[i]][current_val_idx[i][j]] - mindelta);
                                    } else
                                        scope[current_scope_idx[i]]->project(VarVal[current_scope_idx[i]][current_val_idx[i][j]], deltaCosts[current_scope_idx[i]][current_val_idx[i][j]] - mindelta);
                                }
                                scope[current_scope_idx[i]]->findSupport();
                            } else {
                                assert(nbValue[i] <= 2);
                                mindelta = min(deltaCosts[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][0]].first][0], deltaCosts[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][0]].first][1]);
                                TobeProjected += mindelta;
                                scope[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][0]].first]->project(0, deltaCosts[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][0]].first][0] - mindelta);
                                scope[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][0]].first]->project(1, deltaCosts[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][0]].first][1] - mindelta);
                                scope[AMO[VirtualVar[current_scope_idx[i]] - 1][current_val_idx[i][0]].first]->findSupport();
                            }
                        }
                        assert(TobeProjected >= MIN_COST);
                        Constraint::projectLB(TobeProjected);
                    } else if (nonassigned <= NARYPROJECTIONSIZE && (nonassigned < 3 || maxInitDomSize <= NARYPROJECTION3MAXDOMSIZE)) {
                        if (connected()) {
                            deconnect(); // this constraint is removed from the current WCSP problem
                            projectNary();
                        }
                    }
                }
            }
        }
    }

    bool verify() override
    {
        // checks that propagation has been done correctly such that  there exists at least one valid tuple with zero cost (called by WCSP::verify in Debug mode at each search node)
        wcsp->revise(this);
        int breakamo;
        for (unsigned int i = 0; i < AMO.size(); ++i) {
            breakamo = 0;
            for (unsigned int j = 0; j < AMO[i].size(); ++j) {
                if (scope[AMO[i][j].first]->assigned() && scope[AMO[i][j].first]->getValue() == AMO[i][j].second)
                    breakamo++;
            }
            if (breakamo > 1) {
                return false;
            }
        }
        if (capacity <= MaxWeight)
            return true;
        else
            return false;
    }
    void increase(int index) override
    {
        remove(index);
    }
    void decrease(int index) override
    {
        remove(index);
    }
    void remove(int idx) override
    {
        if (isunassigned(idx)) {
            propagate();
        } else if (assigned[idx] < 2) {
            assign(idx);
        }
    }
    void projectFromZero(int index) override
    {
        //TODO: incremental cost propagation
        propagate();
    }

    bool checkEACGreedySolution(int index = -1, Value supportValue = 0) FINAL
    {
        Long W = 0;
        Cost res = -lb + assigneddeltas;
        for (int i = 0; i < arity_; i++) {
            if (CorrAMO[i] == 0) {
                Value support = ((i == index) ? supportValue : scope[i]->getSupport());
                auto it = find(VarVal[i].begin(), VarVal[i].end(), support);
                if (it == VarVal[i].end()) {
                    res += deltaCosts[i].back();
                    W += weights[i].back();
                } else {
                    W += weights[i][distance(VarVal[i].begin(), it)];
                    res += deltaCosts[i][distance(VarVal[i].begin(), it)];
                }
            }
        }
        if (W < Original_capacity)
            res = wcsp->getUb();
        return (res == MIN_COST);
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

    void print(ostream& os) override
    {
        os << this << " knapsackp(";

        int unassigned_ = 0;
        int unassignedAMO = 0;
        for (int i = 0; i < arity_; i++) {
            if (scope[i]->unassigned())
                unassigned_++;
            assert(assigned[i] == 0 || !isunassigned(i));
            if (assigned[i] == 0)
                unassignedAMO++;
            os << wcsp->getName(scope[i]->wcspIndex);
            if (i < arity_ - 1)
                os << ",";
        }
        os << ") "
           << " >= " << capacity << " <= " << MaxWeight << " \\ " << lb << " - " << assigneddeltas << " (";
        for (int i = 0; i < arity_; i++) {
            if (AMO.empty()) {
                for (unsigned int j = 0; j < deltaCosts[i].size(); j++) {
                    os << VarVal[i][j];
                    os << ":";
                    os << weights[i][j];
                    os << ":";
                    os << deltaCosts[i][j];
                    if (j < deltaCosts[i].size() - 1)
                        os << "|";
                }
                if (i < arity_ - 1)
                    os << ",";
            } else {
                for (unsigned int j = 0; j < deltaCosts[i].size(); j++) {
                    os << j;
                    os << ":";
                    os << Original_weigths[i][j];
                    os << ":";
                    os << deltaCosts[i][0];
                    if (j < deltaCosts[i].size() - 1)
                        os << "|";
                }
                if (i < arity_ - 1)
                    os << ",";
            }
        }
        os << ") ";
        os << "/" << getTightness();
        if (ToulBar2::weightedDegree) {
            os << "/" << getConflictWeight();
            for (int i = 0; i < arity_; i++) {
                os << "," << conflictWeights[i];
            }
        }
        os << " arity: " << arity_;
        os << " unassigned: " << unassignedAMO << "/" << nonassigned << "/" << unassigned_ << endl;
    }

    void dump(ostream& os, bool original = true) override
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
                if (!AMO.empty()) {
                    os << " " << -1 << " knapsackc " << Original_capacity;
                    for (int i = 0; i < arity_; i++) {
                        os << " 2 "
                           << "0 " << Original_weigths[i][0] << " 1 " << Original_weigths[i][1];
                    }
                    os << " " << AMO.size();
                    for (unsigned int i = 0; i < AMO.size(); ++i) {
                        os << " " << AMO[i].size();
                        for (unsigned int j = 0; j < AMO[i].size(); ++j) {
                            os << " " << scope[AMO[i][j].first]->getCurrentVarId() << " " << AMO[i][j].second;
                        }
                    }
                    os << endl;
                } else {
                    Long wnot = 0;
                    for (int i = 0; i < arity_; i++) {
                        if (!NotVarVal[i].empty())
                            wnot += weights[i].back();
                    }
                    os << " " << -1 << " knapsackp " << Original_capacity - wnot;
                    for (int i = 0; i < arity_; i++) {
                        if (NotVarVal[i].empty()) {
                            os << " " << VarVal[i].size();
                            for (unsigned int j = 0; j < VarVal[i].size(); ++j) {
                                os << " " << VarVal[i][j];
                                os << " " << weights[i][j];
                            }
                        } else {
                            os << " " << VarVal[i].size() - 1;
                            for (unsigned int j = 0; j < VarVal[i].size() - 1; ++j) {
                                os << " " << VarVal[i][j];
                                os << " " << weights[i][j] - weights[i].back();
                            }
                        }
                    }
                    os << endl;
                }
            } else {
                os << " " << wcsp->getUb() << " " << getDomainSizeProduct() << endl;
                Tuple t;
                Cost c;
                firstlex();
                while (nextlex(t, c)) {
                    for (int i = 0; i < arity_; i++) {
                        os << scope[i]->toValue(t[i]) << " ";
                    }
                    os << c << endl;
                }
            }
        } else {
            os << nonassigned;
            for (int i = 0; i < arity_; i++) {
                if (scope[i]->unassigned())
                    os << " " << scope[i]->getCurrentVarId();
            }
            if (iszerodeltas) {
                if (!AMO.empty()) {
                    os << " " << -1 << " knapsackc " << Original_capacity;
                    for (int i = 0; i < arity_; i++) {
                        if (scope[i]->unassigned())
                            os << " 2 "
                               << "0 " << Original_weigths[i][0] << " 1 " << Original_weigths[i][1];
                    }
                    os << " " << AMO.size();
                    for (unsigned int i = 0; i < AMO.size(); ++i) {
                        os << " " << AMO[i].size();
                        for (unsigned int j = 0; j < AMO[i].size(); ++j) {
                            assert(scope[AMO[i][j].first]->unassigned());
                            os << " " << scope[AMO[i][j].first]->getCurrentVarId() << " " << AMO[i][j].second;
                        }
                    }
                    os << endl;
                } else {
                    Long wnot = 0;
                    for (int i = 0; i < arity_; i++) {
                        if (!NotVarVal[i].empty())
                            wnot += weights[i].back();
                    }
                    os << " " << -1 << " knapsackp " << Original_capacity - wnot;
                    for (int i = 0; i < arity_; i++) {
                        if (scope[i]->unassigned()) {
                            if (NotVarVal[i].empty()) {
                                os << " " << VarVal[i].size();
                                for (unsigned int j = 0; j < VarVal[i].size(); ++j) {
                                    assert(scope[i]->canbe(VarVal[i][j]));
                                    os << " " << VarVal[i][j];
                                    os << " " << weights[i][j];
                                }
                            } else {
                                os << " " << VarVal[i].size() - 1;
                                for (unsigned int j = 0; j < VarVal[i].size() - 1; ++j) {
                                    assert(scope[i]->canbe(VarVal[i][j]));
                                    os << " " << VarVal[i][j];
                                    os << " " << weights[i][j] - weights[i].back();
                                }
                            }
                        }
                    }
                    os << endl;
                }
            } else {
                os << " " << wcsp->getUb() << " " << getDomainSizeProduct() << endl;
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

    void dump_CFN(ostream& os, bool original = true) override
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
        if (!AMO.empty()) { //TODO: cfn reader for knapsackc
            cerr << "Error: cannot save file in cfn format with knapsack constraints including at-most-one constraints (knapsackc)!" << endl;
            throw WrongFileFormat();
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
                os << "\"" << scope[i]->getName() << "\"";
                printed = true;
            }

            if (iszerodeltas) {
                Long wnot = 0;
                for (int i = 0; i < arity_; i++) {
                    if (!NotVarVal[i].empty())
                        wnot += weights[i].back();
                }
                os << "],\n\"type\":\"knapsackv\",\n\"params\":{\"capacity\":" << Original_capacity - wnot << ",\n\t\"weightedvalues\":[";
                printed = false;
                for (int i = 0; i < arity_; i++) {
                    if (NotVarVal[i].empty()) {
                        for (unsigned int j = 0; j < VarVal[i].size(); ++j) {
                            if (printed)
                                os << ",";
                            os << "[" << scope[i]->wcspIndex << "," << scope[i]->toIndex(VarVal[i][j]) << "," << weights[i][j] << "]";
                            printed = true;
                        }
                    } else {
                        for (unsigned int j = 0; j < VarVal[i].size() - 1; ++j) {
                            if (printed)
                                os << ",";
                            os << "[" << scope[i]->wcspIndex << "," << scope[i]->toIndex(VarVal[i][j]) << "," << weights[i][j] - weights[i].back() << "]";
                            printed = true;
                        }
                    }
                }
                os << "]}},\n";
            } else {
                os << "],\n\"defaultcost\":" << wcsp->DCost2Decimal(wcsp->Cost2RDCost(wcsp->getUb())) << ",\n\"costs\":[";
                Tuple t;
                Cost c;
                printed = false;
                firstlex();
                while (nextlex(t, c)) {
                    if (c < wcsp->getUb()) {
                        os << endl;
                        for (int i = 0; i < arity_; i++) {
                            if (printed)
                                os << ",";
                            os << t[i];
                            printed = true;
                        }
                        os << "," << wcsp->DCost2Decimal(wcsp->Cost2RDCost(c));
                    }
                }
                os << "]},\n";
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
                    os << "\"" << scope[i]->getName() << "\"";
                    printed = true;
                }

            if (iszerodeltas) {
                Long wnot = 0;
                for (int i = 0; i < arity_; i++) {
                    if (!NotVarVal[i].empty())
                        wnot += weights[i].back();
                }
                os << "],\n\"type\":\"knapsackv\",\n\"params\":{\"capacity\":" << Original_capacity - wnot << ",\n\t\"weightedvalues\":[";
                printed = false;
                for (int i = 0; i < arity_; i++) {
                    if (scope[i]->unassigned()) {
                        if (NotVarVal[i].empty()) {
                            for (unsigned int j = 0; j < VarVal[i].size(); ++j) {
                                if (scope[i]->canbe(VarVal[i][j])) {
                                    if (printed)
                                        os << ",";
                                    os << "[" << scope[i]->getName() << "," << scope[i]->toCurrentIndex(VarVal[i][j]) << "," << weights[i][j] << "]";
                                    printed = true;
                                }
                            }
                        } else {
                            for (unsigned int j = 0; j < VarVal[i].size() - 1; ++j) {
                                if (scope[i]->canbe(VarVal[i][j])) {
                                    if (printed)
                                        os << ",";
                                    os << "[" << scope[i]->getName() << "," << scope[i]->toCurrentIndex(VarVal[i][j]) << "," << weights[i][j] - weights[i].back() << "]";
                                    printed = true;
                                }
                            }
                        }
                    }
                }
                os << "]}},\n";
            } else {
                os << "],\n\"defaultcost\":inf,\n\"costs\":[";
                Tuple t;
                Cost c;
                printed = false;
                firstlex();
                while (nextlex(t, c)) {
                    if (c < wcsp->getUb()) {
                        os << endl;
                        for (int i = 0; i < arity_; i++) {
                            if (scope[i]->unassigned()) {
                                if (printed)
                                    os << ",";
                                os << scope[i]->toCurrentIndex(scope[i]->toValue(t[i]));
                                printed = true;
                            }
                        }
                        os << "," << wcsp->DCost2Decimal(wcsp->Cost2RDCost(c));
                    }
                }
                os << "]},\n";
            }
        }
    }
};
#endif /*TB2KNAPSACK_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
