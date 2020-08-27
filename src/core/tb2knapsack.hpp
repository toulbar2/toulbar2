#ifndef TB2KNAPSACK_HPP_
#define TB2KNAPSACK_HPP_

#include "tb2abstractconstr.hpp"
#include "tb2ternaryconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2wcsp.hpp"
#include "../utils/tb2store.hpp"
#include <numeric>

// warning! we assume binary variables
class KnapsackConstraint : public AbstractNaryConstraint {
    StoreLong capacity; // knapsack capacity
    vector<Long> weights; // knapsack linear positive integer coefficients
    StoreCost lb; // projected cost to problem lower bound (if it is zero then all deltaCosts must be zero)
    vector<StoreCost> deltaCosts0; // extended costs from unary costs for value 0 to the cost function
    vector<StoreCost> deltaCosts1; // extended costs from unary costs for value 1 to the cost function
    StoreInt nonassigned; // number of non-assigned variables during search, must be backtrackable!
    vector<Long> conflictWeights; // used by weighted degree heuristics
    StoreLong  MaxWeight;
    Long Original_capacity;
    vector<Cost> CostforKnapsack; // temporary data structure for propagate
    vector<int> arrvar; // temporary data structure for propagate
    vector<Double> Weightedtprofit; // temporary data structure for propagate
    StoreLong NegCapacity;

    void projectLB(Cost c)
    {
        if(c>0){
            lb += c;
            Constraint::projectLB(c);
        }

    }

    Double Ceil(Double v){
        const Double epsilon = 1e-7;

        if(floorl(v)+epsilon>v)
            return floorl(v);
        else
            return ceill(v);
    }
    Double Trunc(Double v){
        return truncl(v);
    }

public:
    KnapsackConstraint(WCSP *wcsp, EnumeratedVariable **scope_in, int arity_in, Long capacity_in,
                       vector<Long> weights_in, Long MaxWeigth_in, Long NegCapacity_in)
            : AbstractNaryConstraint(wcsp, scope_in, arity_in)
            , capacity(capacity_in)
            , weights(weights_in)
            , lb(MIN_COST)
            , nonassigned(arity_in)
            , MaxWeight(MaxWeigth_in)
            , Original_capacity(capacity_in)
            , NegCapacity(NegCapacity_in)
            {
        deltaCosts0 = vector<StoreCost>(arity_in, StoreCost(MIN_COST));
        deltaCosts1 = vector<StoreCost>(arity_in, StoreCost(MIN_COST));
        for (int i = 0; i < arity_in; i++) {
            assert(scope_in[i]->getDomainInitSize() == 2 && scope_in[i]->toValue(0) == 0 && scope_in[i]->toValue(1) == 1);
            conflictWeights.push_back(0);
            CostforKnapsack.push_back(MIN_COST);
            arrvar.push_back(i);
            Weightedtprofit.push_back(0.);
        }

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
        //assert(fromElim1==NULL);
        //assert(fromElim2==NULL);
        if (from == this) {
            Constraint::incConflictWeight(1);
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

    bool universal() {
        // returns true if constraint always satisfied
        if (capacity + NegCapacity <= 0)
            return true;
        else
            return false;
    }

    Cost eval(const Tuple& s)
    {
        // returns the cost of the corresponding assignment s
        Long W=0;
        Cost res = -lb;
        for (int i = 0; i < arity_; i++) {
            if (ToulBar2::verbose >= 2)
                cout<<s[i]<<" ";
            EnumeratedVariable* var = (EnumeratedVariable*)getVar(i);
            if(var->toValue(s[i])==0)
                res+=deltaCosts0[i];
            else
            {
                W+=weights[i];
                res+=deltaCosts1[i];
            }
        }
        if(W < Original_capacity || res > wcsp->getUb())
            res=wcsp->getUb();
        assert(res <= wcsp->getUb());
        if (ToulBar2::verbose >= 2)
            cout<<"   "<<res<<endl;
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
    //void setInfiniteCost(Cost ub)

    void assign(int varIndex)
    {
        if (connected(varIndex)) {
            if (ToulBar2::verbose >= 2)
                cout<<" var assigned: " <<scope[varIndex]->getName() << " weight: "<<weights[varIndex]<<endl;
            deconnect(varIndex);
            nonassigned = nonassigned - 1;
            assert(nonassigned >= 0);
            if (scope[varIndex]->getValue() == 1) {
                if (weights[varIndex] > 0) {
                    MaxWeight -= weights[varIndex];
                    capacity -= weights[varIndex];
                } else if (weights[varIndex] < 0) {
                    capacity -= weights[varIndex];
                    NegCapacity += weights[varIndex];
                }
            }
            else if(weights[varIndex]>0)
                MaxWeight -= weights[varIndex];
            else
                NegCapacity+=weights[varIndex];
            if (universal()) {
                deconnect();
                if (ToulBar2::verbose >= 2)
                    cout<<"lb : "<<lb<<endl;
                Cost TobeProjected = -lb;
                lb=0;
                for(int i = 0; i < arity_; i++){
                    if (ToulBar2::verbose >= 2)
                        cout<<scope[i]->getName()<<" deltacosts0 : " <<deltaCosts0[i]<<" deltacosts1 :" << deltaCosts1[i]<<endl;
                    if(scope[i]->unassigned()){
                        scope[i]->project(1,deltaCosts1[i],true);
                        deltaCosts1[i]=0;
                        scope[i]->project(0,deltaCosts0[i],true);
                        deltaCosts0[i]=0;
                        scope[i]->findSupport();
                    }
                    else if(scope[i]->canbe(1)) {
                        TobeProjected+=deltaCosts1[i];
                        deltaCosts1[i]=0;
                    }
                    else{
                        TobeProjected+=deltaCosts0[i];
                        deltaCosts0[i]=0;
                    }
                }
                assert(TobeProjected >= MIN_COST);
                Constraint::projectLB(TobeProjected);
            } else if (nonassigned <= 3) {
                deconnect(); // this constraint is removed from the current WCSP problem
                if (ToulBar2::verbose >= 2)
                    cout<<"EVALUATION "<<endl;
                projectNary(); // and replaced by a ternary constraint in extension
            } else {
                //TODO: incremental bound propagation
                propagate();
                if (ToulBar2::FullEAC)
                    reviseEACGreedySolution();
            }
        }
    }

    void propagate()
    {
        // propagates from scratch the constraint
        if (connected()) {
            for(int i = 0; connected() && i < arity_; i++){
                if (connected(i) && scope[i]->assigned()) {
                    assign(i);
                }
            }
            if (connected()) {
                if (!verify()) {
                    THROWCONTRADICTION;
                } else if (nonassigned > 3 && ToulBar2::LcLevel >= LC_AC) {
                    if (ToulBar2::verbose >= 2)
                        cout << " BOUND  PROPAGATION" << endl;
                    //Bound propagation : we verify that each variable can be both assigned to 1 or 0 without breaking the constraint.
                    int k = 0;
                    bool b = false;
                    while (k < arity_ && b == false) {
                        if(weights[k]>0) {
                            if (MaxWeight - weights[k] < capacity && scope[k]->unassigned()) {
                                scope[k]->assign(1);
                                if (ToulBar2::verbose >= 2)
                                    cout << scope[k]->getName() << " has been assigned" << endl;
                                b = true;
                            } else
                                k++;
                        }
                        else {
                                if (MaxWeight + weights[k] < capacity && scope[k]->unassigned()) {
                                    scope[k]->assign(0);
                                    if (ToulBar2::verbose >= 2)
                                        cout << scope[k]->getName() << " has been assigned" << endl;
                                    b = true;
                                } else
                                    k++;
                            }
                    }
                    if (connected() && ToulBar2::LcLevel >= LC_DAC && b == false) {
                        if (ToulBar2::verbose >= 2)
                            cout << "REDUCED COST PROJECTION" << endl;
                        for (int i = 0; i < arity_; i++) {
                            CostforKnapsack[i] = scope[i]->getCost(1) - scope[i]->getCost(0);
                        }
                        if (ToulBar2::verbose >= 2)
                            cout << "capacity is : " << capacity << endl;
                        Cost NegweightNegprofit=0; //Used in the case of cap>0
                        Cost PosweightPosprofit=0; //Used in the case of cap<=0
                        int nbmincost=0;
                        // Compute weighted profit : p_i / w_i. MAX_COST and MIN_COST are used when the weighted profit is no relevant and we need
                        // to impose some variables to be at the end or the beginning of the sorting.
                        if(capacity>=0) {
                            for (int i = 0; i < arity_; i++) {
                                if (scope[i]->unassigned()) {
                                    if (ToulBar2::verbose >= 2)
                                        cout << scope[i]->getName() << " : " << CostforKnapsack[i] << " / "
                                             << weights[i] << endl;
                                    if ((weights[i] > 0 && CostforKnapsack[i] >= MIN_COST) ||
                                        (weights[i] > 0 && CostforKnapsack[i] < MIN_COST)) {
                                        Weightedtprofit[i] = Double(CostforKnapsack[i]) / weights[i];
                                    } else if (weights[i] < 0 && CostforKnapsack[i] >= MIN_COST) {
                                        Weightedtprofit[i] = MAX_COST;
                                    } else if (weights[i] < 0 && CostforKnapsack[i] < MIN_COST) {
                                        Weightedtprofit[i] = Double(CostforKnapsack[i]) / weights[i];
                                        NegweightNegprofit -= weights[i];
                                    }
                                } else {
                                    Weightedtprofit[i] = MAX_COST;
                                }
                            }
                        }
                        else{
                            for (int i = 0; i < arity_; i++) {
                                if (scope[i]->unassigned()) {
                                    if (ToulBar2::verbose >= 2)
                                        cout << scope[i]->getName() << " : " << CostforKnapsack[i] << " / "
                                             << weights[i] << endl;
                                    if (weights[i] < 0 && CostforKnapsack[i] < MIN_COST) {
                                        Weightedtprofit[i] = Double(CostforKnapsack[i]) / weights[i];
                                    } else if (weights[i] < 0 && CostforKnapsack[i] >= MIN_COST) {
                                        Weightedtprofit[i] = MIN_COST;
                                        nbmincost++;
                                    } else if (weights[i] > 0 && CostforKnapsack[i] > MIN_COST) {
                                        Weightedtprofit[i] = Double(CostforKnapsack[i]) / weights[i];
                                        PosweightPosprofit += weights[i];
                                    } else if (weights[i] > 0 && CostforKnapsack[i] <= MIN_COST){
                                        Weightedtprofit[i] = MAX_COST;
                                    }
                                } else {
                                    Weightedtprofit[i] = MIN_COST;
                                    nbmincost++;
                                }
                            }
                        }
                        if (ToulBar2::verbose >= 2){
                            cout << "Unsorted variable by Weighted profit : ";
                            for (int i = 0; i < int(arrvar.size()); i++) {
                                cout << " " << arrvar[i];
                                cout << "-" << Weightedtprofit[arrvar[i]];
                            }
                            cout << endl;
                        }
                        //Sort variables in ascendant or descendant order depending of the sign of the capacity
                        //Use stable sort for reproducibility of the results
                        if(capacity>=0)
                            stable_sort(arrvar.begin(), arrvar.end(),[&](int x, int y) { return Weightedtprofit[x] < Weightedtprofit[y]; });
                        else
                            stable_sort(arrvar.begin(), arrvar.end(),[&](int x, int y) { return Weightedtprofit[x] > Weightedtprofit[y]; });

                        if (ToulBar2::verbose >= 2){
                            cout << "Sorted variable by weighted profit: ";
                            for (int i = 0; i < int(arrvar.size()); i++) {
                                cout << " " << arrvar[i];
                                cout << "-" << Weightedtprofit[arrvar[i]];
                            }
                            cout << endl;}

                        //Find splitting variable x_k 
                        Long W = 0;
                        int splitvar = -1;
                        Cost c = 0; //Cost we will project on c_0, if the capacity is positive it is the profit sum of the variables before x_k
                        // if the capacity is negative it is the profit sum of the variables after x_k 
                        if(capacity >=0) {
                            while (W < capacity + NegweightNegprofit) {
                                splitvar = splitvar + 1;
                                if (weights[arrvar[splitvar]] > 0)
                                    W += weights[arrvar[splitvar]];
                                else
                                    W -= weights[arrvar[splitvar]];
                                if (CostforKnapsack[arrvar[splitvar]] > MIN_COST)
                                    c += CostforKnapsack[arrvar[splitvar]];
                                else if (weights[arrvar[splitvar]] < 0)
                                    c -= CostforKnapsack[arrvar[splitvar]];
                            }
                        }
                        else{
                            while (W >= capacity - PosweightPosprofit && splitvar != arity_-nbmincost-1) {
                                splitvar = splitvar + 1;
                                if (weights[arrvar[splitvar]] < 0)
                                    W += weights[arrvar[splitvar]];
                                else if (CostforKnapsack[arrvar[splitvar]] > MIN_COST)
                                    W -= weights[arrvar[splitvar]];
                                else
                                    W += weights[arrvar[splitvar]];
                            }
                            for (int i=splitvar+1; i < arity_-nbmincost-1; i++) {
                                if(CostforKnapsack[arrvar[i]] > MIN_COST && weights[arrvar[i]]>0)
                                    c+=CostforKnapsack[arrvar[i]];
                                else if(CostforKnapsack[arrvar[i]] < MIN_COST && weights[arrvar[i]]<0)
                                    c-=CostforKnapsack[arrvar[i]];
                            }
                        }
                        if (ToulBar2::verbose >= 2)
                            cout << "splitvar : " << splitvar << endl;
                        Double xk=0;
                        // we add the profit of x_k 
                        if(splitvar>-1) {
                            Long capacityLeft;
                            if(capacity >=0) {
                                if (weights[arrvar[splitvar]] > 0)
                                    capacityLeft = capacity + NegweightNegprofit - W + weights[arrvar[splitvar]];
                                else
                                    capacityLeft = capacity + NegweightNegprofit - W - weights[arrvar[splitvar]];

                                xk = Double(capacityLeft) / weights[arrvar[splitvar]];
                            }
                            else if(W <= capacity - PosweightPosprofit){
                                if (weights[arrvar[splitvar]] > 0)
                                    capacityLeft = capacity - PosweightPosprofit - W - weights[arrvar[splitvar]];
                                else
                                    capacityLeft = capacity - PosweightPosprofit - W + weights[arrvar[splitvar]];

                                xk = Double(capacityLeft) / weights[arrvar[splitvar]];
                            }
                            if (xk < 0)
                                xk = -xk;
                            assert(xk <= 1);
                            assert(xk >= 0);
                            if (capacity>=0) {
                                if (CostforKnapsack[arrvar[splitvar]] > MIN_COST)
                                    c = c - CostforKnapsack[arrvar[splitvar]] +
                                        Ceil(CostforKnapsack[arrvar[splitvar]] * xk);
                                else if (weights[arrvar[splitvar]] < 0)
                                    c = c + CostforKnapsack[arrvar[splitvar]] +
                                        Ceil(-CostforKnapsack[arrvar[splitvar]] * xk);
                            }
                            else if(W <= capacity - PosweightPosprofit){
                                if (CostforKnapsack[arrvar[splitvar]] > MIN_COST)
                                    c = c + Ceil(CostforKnapsack[arrvar[splitvar]] * (1-xk));
                                else if (weights[arrvar[splitvar]] < 0)
                                    c = c + Ceil(-CostforKnapsack[arrvar[splitvar]] * (1-xk));
                            }
                        }
                        //If c<= 0 it means that we use only negative cost, it means we use only variables with cost on value 0
                        if (c > 0) {
                            //New value for p_k to obtain integer cost (p_k' might be decimal)
                            Double Newkcost = 0.;
                            Double epsi = 1e-5;
                            if(splitvar != 0) {
                                if (capacity >= 0)
                                    Newkcost = min(Double(CostforKnapsack[arrvar[splitvar]]),
                                                   max(Double(CostforKnapsack[arrvar[splitvar - 1]]) /
                                                       weights[arrvar[splitvar - 1]] * weights[arrvar[splitvar]],
                                                       Double((Trunc(xk * CostforKnapsack[arrvar[splitvar]]) + epsi) /
                                                              xk)));
                                else
                                    Newkcost = min(Double(CostforKnapsack[arrvar[splitvar]]),
                                                   max(Double(CostforKnapsack[arrvar[splitvar - 1]]) /
                                                       weights[arrvar[splitvar - 1]] * weights[arrvar[splitvar]],
                                                       Double((Trunc((1 - xk) * CostforKnapsack[arrvar[splitvar]]) +
                                                               epsi) / (1 - xk))));

                            }
                            else {
                                    if(capacity>=0)
                                        Newkcost = min(Double(CostforKnapsack[arrvar[splitvar]]), Double((Trunc(xk * CostforKnapsack[arrvar[splitvar]]) + epsi) / xk));
                                    else
                                        Newkcost = min(Double(CostforKnapsack[arrvar[splitvar]]), Double((Trunc((1-xk) * CostforKnapsack[arrvar[splitvar]]) + epsi) / (1-xk)));
                           }
                            //--------------Test-------------
#ifndef NDEBUG
                            Double Testprofit = 0;
                            if(capacity>=0) {
                                for (int i = 0; i < splitvar; i++) {
                                    if (CostforKnapsack[arrvar[i]] > MIN_COST)
                                        Testprofit += CostforKnapsack[arrvar[i]];
                                    else if (weights[arrvar[i]] < 0)
                                        Testprofit -= CostforKnapsack[arrvar[i]];
                                }
                                if (CostforKnapsack[arrvar[splitvar]] > MIN_COST)
                                    Testprofit += xk * Newkcost;
                                else
                                    Testprofit -= xk * Newkcost;
                            }
                            else{
                                for (int i = splitvar+1; i < arity_-nbmincost-1; i++) {
                                    if (CostforKnapsack[arrvar[i]] > MIN_COST && weights[arrvar[i]]>0)
                                        Testprofit += CostforKnapsack[arrvar[i]];
                                    else if (weights[arrvar[i]] < 0 && weights[arrvar[i]]<0)
                                        Testprofit -= CostforKnapsack[arrvar[i]];
                                }
                                if(W <= capacity - PosweightPosprofit){
                                    if (CostforKnapsack[arrvar[splitvar]] > MIN_COST)
                                        Testprofit += (1-xk) * Newkcost;
                                    else
                                        Testprofit -= (1-xk) * Newkcost;
                                }
                            }
                            assert(Ceil(Testprofit) == Ceil(c));
#endif
                            //--------------------------
                            if (ToulBar2::verbose >= 2)
                                cout << "deltaCost : ";
                            // Compute the reduced cost for each variable and add it to deltacost, depending of the sign of the capacity we proceed differently
                            if(capacity>=0) {
                                for (int i = 0; i < arity_; i++) {
                                    if (scope[arrvar[i]]->unassigned()) {
                                        if (i < splitvar) {
                                            if (CostforKnapsack[arrvar[i]] > MIN_COST) {
                                                deltaCosts1[arrvar[i]] += CostforKnapsack[arrvar[i]];
                                                assert(CostforKnapsack[arrvar[i]] <= scope[arrvar[i]]->getCost(1));
                                                scope[arrvar[i]]->extend(1, CostforKnapsack[arrvar[i]]);
                                            } else if (CostforKnapsack[arrvar[i]] < MIN_COST) {
                                                deltaCosts0[arrvar[i]] -= CostforKnapsack[arrvar[i]];
                                                assert(-CostforKnapsack[arrvar[i]] <= scope[arrvar[i]]->getCost(0));
                                                scope[arrvar[i]]->extend(0, -CostforKnapsack[arrvar[i]]);
                                            }
                                        } else if (CostforKnapsack[arrvar[i]] > MIN_COST && weights[arrvar[i]] > 0) {
                                            deltaCosts1[arrvar[i]] += min( Ceil(Newkcost / Double(weights[arrvar[splitvar]]) * weights[arrvar[i]]), Double(c));
                                            assert(min(Ceil(Newkcost / weights[arrvar[splitvar]] * weights[arrvar[i]]), Double(c)) <= scope[arrvar[i]]->getCost(1));
                                            scope[arrvar[i]]->extend(1, min(Ceil(Newkcost / Double(weights[arrvar[splitvar]]) * weights[arrvar[i]]), Double(c)));
                                        } else if (CostforKnapsack[arrvar[i]] < MIN_COST && weights[arrvar[i]] < 0) {
                                            deltaCosts0[arrvar[i]] += min( Ceil(Newkcost / Double(weights[arrvar[splitvar]]) * -weights[arrvar[i]]), Double(c));
                                            assert(min(Ceil(Newkcost / weights[arrvar[splitvar]] * -weights[arrvar[i]]), Double(c)) <= scope[arrvar[i]]->getCost(0));
                                            scope[arrvar[i]]->extend(0, min(Ceil(Newkcost / Double(weights[arrvar[splitvar]]) * -weights[arrvar[i]]), Double(c)));
                                        }
                                    }
                                }
                            }
                            else{
                                for (int i = 0; i < arity_; i++) {
                                    if (scope[arrvar[i]]->unassigned()) {
                                        if (i > splitvar) {
                                            if (CostforKnapsack[arrvar[i]] > MIN_COST) {
                                                deltaCosts1[arrvar[i]] += CostforKnapsack[arrvar[i]];
                                                assert(CostforKnapsack[arrvar[i]] <= scope[arrvar[i]]->getCost(1));
                                                scope[arrvar[i]]->extend(1, CostforKnapsack[arrvar[i]]);
                                            } else if (CostforKnapsack[arrvar[i]] < MIN_COST) {
                                                deltaCosts0[arrvar[i]] -= CostforKnapsack[arrvar[i]];
                                                assert(-CostforKnapsack[arrvar[i]] <= scope[arrvar[i]]->getCost(0));
                                                scope[arrvar[i]]->extend(0, -CostforKnapsack[arrvar[i]]);
                                            }
                                        } else if (CostforKnapsack[arrvar[i]] > MIN_COST && weights[arrvar[i]] > 0) {
                                            deltaCosts1[arrvar[i]] += min( Ceil(Newkcost / Double(weights[arrvar[splitvar]]) * weights[arrvar[i]]), Double(c));
                                            assert(min(Ceil(Newkcost / weights[arrvar[splitvar]] * weights[arrvar[i]]), Double(c)) <=scope[arrvar[i]]->getCost(1));
                                            scope[arrvar[i]]->extend(1, min(Ceil(Newkcost / Double(weights[arrvar[splitvar]]) * weights[arrvar[i]]), Double(c)));
                                        } else if (CostforKnapsack[arrvar[i]] < MIN_COST && weights[arrvar[i]] < 0) {
                                            deltaCosts0[arrvar[i]] += min(Ceil(Newkcost / Double(weights[arrvar[splitvar]]) * -weights[arrvar[i]]), Double(c));
                                            assert(min(Ceil(Newkcost / weights[arrvar[splitvar]] * -weights[arrvar[i]]), Double(c)) <=scope[arrvar[i]]->getCost(0));
                                            scope[arrvar[i]]->extend(0, min(Ceil(Newkcost / Double(weights[arrvar[splitvar]]) * -weights[arrvar[i]]), Double(c)));
                                        }
                                    }
                                }
                            }
                            if (ToulBar2::verbose >= 2) {
                                cout << endl << "c :" << c << " lb : " << lb << endl;
                            }
                            projectLB(c);
                        }
                    }
                } else if (nonassigned <= 3) {
                    assert(connected());
                    deconnect(); // this constraint is removed from the current WCSP problem
                    if (ToulBar2::verbose >= 2)
                        cout<<"EVALUATION 2"<<endl;
                    projectNary();  // and replaced by a ternary constraint in extension
                }
            }
        }
    }


    bool verify()
    {
        // checks that propagation has been done correctly such that at least there exists one valid tuple with zero cost (called by WCSP::verify in Debug mode at each search node)
        if(capacity<=MaxWeight)
            return true;
        else
            return false;
    }
    void increase(int index) {}
    void decrease(int index) {}
    void remove(int index) {}
    void projectFromZero(int index)
    {
        //TODO: incremental cost propagation
        propagate();
    }

    //bool checkEACGreedySolution(int index = -1, Value supportValue = 0) FINAL //TODO: checks if current EAC support has zero cost
    //bool reviseEACGreedySolution(int index = -1, Value supportValue = 0) FINAL

    void print(ostream& os)
    {
        os << endl
           << this << " knapsack(";

        int unassigned_ = 0;
        for (int i = 0; i < arity_; i++) {
            if (scope[i]->unassigned())
                unassigned_++;
            os << weights[i];
            os << "*";
            os << scope[i]->wcspIndex;
            if (i < arity_ - 1)
                os << " + ";
        }
        os << " >= " << capacity << ") / " << lb << " (";
        for (int i = 0; i < arity_; i++) {
            os << deltaCosts0[i];
            if (i < arity_ - 1)
                os << ",";
        }
        os << ") (";
        for (int i = 0; i < arity_; i++) {
            os << deltaCosts1[i];
            if (i < arity_ - 1)
                os << ",";
        }
        os << ") ";
        if (ToulBar2::weightedDegree) {
            os << "/" << getConflictWeight();
            for (int i = 0; i < arity_; i++) {
                os << "," << conflictWeights[i];
            }
        }
        os << " arity: " << arity_;
        os << " unassigned: " << (int)nonassigned << "/" << unassigned_ << endl;
    }

    void dump(ostream& os, bool original = true)
    {
        bool iszerodeltas = (lb == MIN_COST);
        for (vector<StoreCost>::iterator it = deltaCosts0.begin(); it != deltaCosts0.end(); ++it) {
            Cost d = (*it);
            if (d != MIN_COST) {
                iszerodeltas = false;
                break;
            }
        }
        if (iszerodeltas) {
            for (vector<StoreCost>::iterator it = deltaCosts1.begin(); it != deltaCosts1.end(); ++it) {
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
                os << " " << -1 << " knapsack " << capacity;
                for (int i = 0; i < arity_; i++) {
                    os << " " << weights[i];
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
                    if (scope[i]->unassigned())
                        os << " " << weights[i];
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
        for (vector<StoreCost>::iterator it = deltaCosts0.begin(); it != deltaCosts0.end(); ++it) {
            Cost d = (*it);
            if (d != MIN_COST) {
                iszerodeltas = false;
                break;
            }
        }
        if (iszerodeltas) {
            for (vector<StoreCost>::iterator it = deltaCosts1.begin(); it != deltaCosts1.end(); ++it) {
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
            os << "],\n\"type\":\"knapsack\",\n\"params\":{\"capacity\":" << capacity << ",\n\t\"weights\":[";

            if (iszerodeltas) {
                printed = false;
                for (int i = 0; i < arity_; i++) {
                    if (printed)
                        os << ",";
                    os << weights[i];
                    printed = true;
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
            os << "],\n\"type\":\"knapsack\",\n\"params\":{\"capacity\":" << capacity << ",\n\t\"weights\":[";

            if (iszerodeltas) {
                printed = false;
                for (int i = 0; i < arity_; i++) {
                    if (printed)
                        os << ",";
                    os << weights[i];
                    printed = true;
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
