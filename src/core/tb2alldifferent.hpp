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
//#include "tb2vacutils.hpp"
#include "search/tb2clusters.hpp"

class AllDifferentConstraint : public AbstractNaryConstraint {
    Cost Original_ub; // initial upper bound when creating the constraint
    StoreCost lb; // projected cost to problem lower bound (if it is zero then all deltaCosts must be zero)
    StoreCost assigneddeltas;
    vector<Long> conflictWeights; // used by weighted degree heuristics
    vector<vector<StoreCost>> deltaCosts; // extended costs from unary costs to values in the cost function

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
    {
        if (arity_in > 0) {
            for (int i = 0; i < arity_in; i++) {
                conflictWeights.push_back(0);
                assert((int)scope[i]->getDomainInitSize() == arity_in); //TODO: do not assume identical initial domain size equal to the arity
                deltaCosts.emplace_back(arity_in, MIN_COST);
            }
            propagate();
        } else {
            deconnect();
        }
    }

    virtual ~AllDifferentConstraint() {}

    void read(istream& file) {} // TODO: add a parameter for controlling the level of propagation if necessary

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
        vector<bool> alreadyUsed(arity_, false);
        for (int i = 0; i < arity_; i++) {
            assert(s[i] < arity_);
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

            if (getNonAssigned() <= 2) { // TODO: use NARYPROJECTIONSIZE
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
                //TODO: compute a lower bound and prune forbidden domain values
                Cost hungarian = MIN_COST;
                //TODO: modify unary costs by a set of EPTs using ExtOrProj
                projectLB(hungarian);
            }
        }
    }

    //TODO: checks that the constraint is still satisfiable (called by WCSP::verify in Debug mode at each search node)
    //bool verify() override;

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
            propagate();
            if (revise)
                reviseEACGreedySolution();
        } else {
            assign(index);
        }
    }
    void projectFromZero(int index) override
    {
        // TODO: incremental cost propagation
        bool revise = ToulBar2::FullEAC && (getVar(index)->cannotbe(getVar(index)->getSupport()) || getVar(index)->getCost(getVar(index)->getSupport()) + getCost(index, getVar(index)->getSupport()) > MIN_COST);
        propagate();
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
        os << " arity: " << arity_;
        os << " unassigned: " << getNonAssigned() << "/" << unassigned_ << endl;
    }

    void dump(ostream& os, bool original = true) override
    {
        if (original) {
            os << arity_;
            for (int i = 0; i < arity_; i++)
                os << " " << scope[i]->wcspIndex;
            os << " -1 alldiff" << endl;
        } else {
            os << getNonAssigned();
            for (int i = 0; i < arity_; i++)
                if (scope[i]->unassigned())
                    os << " " << scope[i]->getCurrentVarId();
            os << " -1 alldiff" << endl;
        }
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
            os << "],\"type\":\"alldiff\",\"params\":{}";
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
            os << "],\"type\":\"alldiff\",\"params\":{}";
        }
        os << "},\n";
    }

};
#endif /*TB2ALLDIFFERENT_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

