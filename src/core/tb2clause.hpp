#ifndef TB2WCLAUSE_HPP_
#define TB2WCLAUSE_HPP_

#include "tb2abstractconstr.hpp"
#include "tb2ternaryconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2wcsp.hpp"
#include <numeric>

// warning! we assume binary variables
class WeightedClause : public AbstractNaryConstraint {
    Cost cost; // clause weight
    String tuple; // forbidden assignment corresponding to the negation of the clause
    StoreCost lb; // projected cost to problem lower bound
    vector<StoreCost> deltaCosts; // extended costs from unary costs to the cost function
    int support; // index of a variable in the scope with a zero unary cost on its value which satisfies the clause
    StoreInt nonassigned; // number of non-assigned variables during search, must be backtrackable!
    String evalTuple; // temporary data structure
    vector<Long> conflictWeights; // used by weighted degree heuristics

    Value getTuple(int i) { return scope[i]->toValue(tuple[i] - CHAR_FIRST); }
    Value getClause(int i) { return scope[i]->toValue(!(tuple[i] - CHAR_FIRST)); }
    void projectLB(Cost c)
    {
        lb += c;
        assert(lb <= cost);
        Constraint::projectLB(c);
    }
    void extend(Cost c)
    {
        for (int i = 0; i < arity_; i++) {
            EnumeratedVariable* x = scope[i];
            if (x->unassigned()) {
                deltaCosts[i] += c;
                Value v = getClause(i);
                TreeDecomposition* td = wcsp->getTreeDec();
                if (td)
                    td->addDelta(cluster, x, v, -c);
                x->extend(v, c);
            } else
                assert(x->getValue() == getTuple(i));
        }
        projectLB(c);
    }
    void satisfied(int varIndex)
    {
        nonassigned = 0;
        assert(scope[varIndex]->assigned());
        assert(scope[varIndex]->getValue() == getClause(varIndex));
        assert(deltaCosts[varIndex] == lb);
        for (int i = 0; i < arity_; i++) {
            EnumeratedVariable* x = scope[i];
            assert(deconnected(i));
            if (i != varIndex) {
                Cost c = deltaCosts[i];
                Value v = getClause(i);
                if (c > MIN_COST) {
                    deltaCosts[i] = MIN_COST;
                    if (x->unassigned()) {
                        if (!CUT(c + wcsp->getLb(), wcsp->getUb())) {
                            TreeDecomposition* td = wcsp->getTreeDec();
                            if (td)
                                td->addDelta(cluster, x, v, c);
                        }
                        x->project(v, c, true);
                        x->findSupport();
                    } else {
                        if (x->canbe(v)) {
                            Constraint::projectLB(c);
                        }
                    }
                }
            }
        }
    }

public:
    // warning! give the negation of the clause as input
    WeightedClause(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in, Cost cost_in = MIN_COST, String tuple_in = String())
        : AbstractNaryConstraint(wcsp, scope_in, arity_in)
        , cost(cost_in)
        , tuple(tuple_in)
        , lb(MIN_COST)
        , support(0)
        , nonassigned(arity_in)
    {
        if (tuple_in.empty() && arity_in > 0)
            tuple = String(arity_in, CHAR_FIRST);
        deltaCosts = vector<StoreCost>(arity_in, StoreCost(MIN_COST));
        Char* tbuf = new Char[arity_in + 1];
        tbuf[arity_in] = '\0';
        for (int i = 0; i < arity_in; i++) {
            assert(scope_in[i]->getDomainInitSize() == 2);
            tbuf[i] = CHAR_FIRST;
            conflictWeights.push_back(0);
        }
        evalTuple = String(tbuf);
        delete[] tbuf;
    }

    virtual ~WeightedClause() {}

    void setTuple(const String& tin, Cost c) FINAL
    {
        cost = c;
        tuple = tin;
    }

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

    bool universal()
    {
        if (cost != MIN_COST || lb != MIN_COST)
            return false;
        for (int i = 0; i < arity_; i++)
            if (deltaCosts[i] != MIN_COST)
                return false;
        return true;
    }

    Cost eval(const String& s)
    {
        Cost res = -lb;
        bool istuple = true;
        for (int i = 0; i < arity_; i++) {
            if (tuple[i] != s[i]) {
                res += deltaCosts[i];
                istuple = false;
            }
        }
        if (istuple)
            res += cost;
        assert(res >= MIN_COST);
        return res;
    }
    Cost evalsubstr(const String& s, Constraint* ctr) FINAL { return evalsubstrAny(s, ctr); }
    Cost evalsubstr(const String& s, NaryConstraint* ctr) FINAL { return evalsubstrAny(s, ctr); }
    template <class T>
    Cost evalsubstrAny(const String& s, T* ctr)
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
            evalTuple[i] = var->toIndex(var->getValue()) + CHAR_FIRST;
        }
        return eval(evalTuple);
    }

    double computeTightness() { return 1.0 * cost / getDomainSizeProduct(); }

    //    pair< pair<Cost,Cost>, pair<Cost,Cost> > getMaxCost(int index, Value a, Value b) { return make_pair(make_pair(MAX_COST,MAX_COST),make_pair(MAX_COST,MAX_COST)); }

    Cost getMaxFiniteCost()
    {
        Cost sumdelta = accumulate(deltaCosts.begin(), deltaCosts.end(), -lb);
        if (CUT(sumdelta, wcsp->getUb()))
            return MAX_COST;
        if (CUT(cost, wcsp->getUb()))
            return sumdelta;
        else
            return max(sumdelta, cost - lb);
    }
    void setInfiniteCost(Cost ub)
    {
        Cost mult_ub = ((ub < (MAX_COST / MEDIUM_COST)) ? (max(LARGE_COST, ub * MEDIUM_COST)) : ub);
        if (CUT(cost, ub))
            cost = mult_ub;
    }

    void assign(int varIndex)
    {
        if (connected(varIndex)) {
            deconnect(varIndex);
            nonassigned = nonassigned - 1;
            assert(nonassigned >= 0);

            if (scope[varIndex]->getValue() == getClause(varIndex)) {
                deconnect();
                satisfied(varIndex);
                return;
            }

            if (nonassigned <= 3) {
                deconnect();
                projectNary();
            }
        }
    }

    // propagates the minimum between the remaining clause weight and unary costs of all literals to the problem lower bound
    void propagate()
    {
        Cost mincost = (connected() && scope[support]->unassigned()) ? scope[support]->getCost(getClause(support)) : MAX_COST;
        for (int i = 0; connected() && i < arity_; i++) {
            EnumeratedVariable* x = scope[i];
            if (x->assigned()) {
                assign(i);
            } else if (mincost > MIN_COST) {
                Cost ucost = x->getCost(getClause(i));
                if (ucost < mincost) {
                    mincost = ucost;
                    support = i;
                }
            }
        }
        if (connected() && mincost < MAX_COST && mincost > MIN_COST && cost > lb) {
            extend(min(cost - lb, mincost));
        }
    };

    bool verify()
    {
        String t;
        Cost c;
        firstlex();
        while (nextlex(t, c)) {
            if (c == MIN_COST)
                return true;
        }
        return false;
    }
    void increase(int index) {}
    void decrease(int index) {}
    void remove(int index) {}
    void projectFromZero(int index)
    {
        if (index == support && cost > lb)
            propagate();
    }

    void print(ostream& os)
    {
        os << endl
           << this << " clause(";

        int unassigned_ = 0;
        for (int i = 0; i < arity_; i++) {
            if (scope[i]->unassigned())
                unassigned_++;
            if (getClause(i) == 0)
                os << "-";
            os << scope[i]->wcspIndex;
            if (i < arity_ - 1)
                os << ",";
        }
        os << ") s:" << support << " / " << cost << " - " << lb << " (";
        for (int i = 0; i < arity_; i++) {
            os << deltaCosts[i];
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
        Cost maxdelta = MIN_COST;
        for (vector<StoreCost>::iterator it = deltaCosts.begin(); it != deltaCosts.end(); ++it) {
            Cost d = (*it);
            if (d > maxdelta)
                maxdelta = d;
        }
        if (original) {
            os << arity_;
            for (int i = 0; i < arity_; i++)
                os << " " << scope[i]->wcspIndex;
            if (maxdelta == MIN_COST) {
                os << " " << 0 << " " << 1 << endl;
                for (int i = 0; i < arity_; i++) {
                    os << scope[i]->toValue(tuple[i] - CHAR_FIRST) << " ";
                }
                os << cost << endl;
            } else {
                os << " " << 0 << " " << getDomainSizeProduct() << endl;
                String t;
                Cost c;
                firstlex();
                while (nextlex(t, c)) {
                    for (int i = 0; i < arity_; i++) {
                        os << scope[i]->toValue(t[i] - CHAR_FIRST) << " ";
                    }
                    os << c << endl;
                }
            }
        } else {
            os << nonassigned;
            for (int i = 0; i < arity_; i++)
                if (scope[i]->unassigned())
                    os << " " << scope[i]->getCurrentVarId();
            if (maxdelta == MIN_COST) {
                os << " " << 0 << " " << 1 << endl;
                for (int i = 0; i < arity_; i++) {
                    if (scope[i]->unassigned())
                        os << scope[i]->toCurrentIndex(scope[i]->toValue(tuple[i] - CHAR_FIRST)) << " ";
                }
                os << min(wcsp->getUb(), cost) << endl;
            } else {
                os << " " << 0 << " " << getDomainSizeProduct() << endl;
                String t;
                Cost c;
                firstlex();
                while (nextlex(t, c)) {
                    for (int i = 0; i < arity_; i++) {
                        if (scope[i]->unassigned())
                            os << scope[i]->toCurrentIndex(scope[i]->toValue(t[i] - CHAR_FIRST)) << " ";
                    }
                    os << min(wcsp->getUb(), c) << endl;
                }
            }
        }
    }
};
#endif /*TB2WCLAUSE_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
