#ifndef TB2GLOBALWCSP_HPP_
#define TB2GLOBALWCSP_HPP_

#include "tb2abstractconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2wcsp.hpp"

extern void setvalue(int wcspId, int varIndex, Value value, void* solver);
extern void tb2setvalue(int wcspId, int varIndex, Value value, void* solver);
extern void tb2removevalue(int wcspId, int varIndex, Value value, void* solver);
extern void tb2setmin(int wcspId, int varIndex, Value value, void* solver);
extern void tb2setmax(int wcspId, int varIndex, Value value, void* solver);

class WeightedCSPConstraint : public AbstractNaryConstraint {
    Cost lb; // encapsulated slave problem lower bound hard constraint (must be greater or equal to this bound)
    Cost ub; // encapsulated slave problem upper bound hard constraint (must be strictly less than this bound)
    WCSP *problem; // encapsulated slave problem
    WCSP *negproblem; // encapsulated slave problem in negation form
    StoreInt nonassigned; // number of non-assigned variables during search, must be backtrackable!
    vector<int> varIndexes; // copy of scope using integer identifiers inside slave problem (should be equal to [0, 1, 2, ..., arity-1])
    vector<Value> newValues; // used to convert Tuples into variable assignments
    vector<Long> conflictWeights; // used by weighted degree heuristics

public:
    static WCSP* MasterWeightedCSP; // Master problem used by value and variable ordering heuristics
    static map<int, WeightedCSPConstraint *> WeightedCSPConstraints;

    WeightedCSPConstraint(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in, WCSP *problem_in, WCSP *negproblem_in, Cost lb_in, Cost ub_in)
        : AbstractNaryConstraint(wcsp, scope_in, arity_in)
        , lb(lb_in)
        , ub(ub_in)
        , problem(problem_in)
        , negproblem(negproblem_in)
        , nonassigned(arity_in)
    {
        assert(arity_in == (int)problem_in->numberOfVariables());
        assert(arity_in == (int)negproblem_in->numberOfVariables());
        for (int i = 0; i < arity_in; i++) {
            assert(scope_in[i]->getDomainInitSize() == ((EnumeratedVariable *)problem->getVar(i))->getDomainInitSize());
            assert(scope_in[i]->getDomainInitSize() == ((EnumeratedVariable *)negproblem->getVar(i))->getDomainInitSize());
            varIndexes.push_back(i);
            newValues.push_back(scope_in[i]->getInf());
            conflictWeights.push_back(0);
        }
        ToulBar2::setvalue = ::tb2setvalue;
        ToulBar2::removevalue = ::tb2removevalue;
        ToulBar2::setmin = ::tb2setmin;
        ToulBar2::setmax = ::tb2setmax;
        WeightedCSPConstraints[problem->getIndex()] = this;
        WeightedCSPConstraints[negproblem->getIndex()] = this;
        assert(MasterWeightedCSP == NULL || MasterWeightedCSP == wcsp); //FIXME: the slave problem cannot contain a WeightedCSPConstraint inside!
        MasterWeightedCSP = wcsp;
        problem->updateUb(ub);
        negproblem->updateUb(-lb + negproblem->getNegativeLb() + UNIT_COST);
    }

    virtual ~WeightedCSPConstraint() {}

    bool extension() const FINAL { return false; } // this is not a cost function represented by an exhaustive table of costs

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
        if (from == this) {
            if (deconnected() || nonassigned == arity_) {
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

    //FIXME: only valid if all hard constraints in the slave problem are also present in the master
    bool universal() override
    {
        if (problem->getLb() >= lb && -(negproblem->getLb() - negproblem->getNegativeLb()) < ub) {
            return true;
        } else {
            return false;
        }
    }

    Cost eval(const Tuple& s) override
    {
        for (int i=0; i < arity_; i++) {
            newValues[i] = ((EnumeratedVariable *)getVar(i))->toValue(s[i]);
        }
        int depth = Store::getDepth();
        bool unsat = false;
        try {
            Store::store();
            problem->assignLS(varIndexes, newValues); // throw a Contradiction if unsatisfied
            negproblem->assignLS(varIndexes, newValues);
        } catch (const Contradiction&) {
            problem->whenContradiction();
            negproblem->whenContradiction();
            unsat = true;
        }
        Store::restore(depth);
        if (unsat) {
            return MAX_COST;
        } else {
            return MIN_COST;
        }
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

    double computeTightness() override
    {
        double res = 0.; //FIXME: take into account elimBinConstr and elimTernConstr
        for (unsigned int c=0; c < problem->numberOfConstraints(); c++) {
            if (problem->getCtr(c)->connected()) {
                res += problem->getCtr(c)->getTightness();
            }
        }
        return res / problem->numberOfConnectedConstraints();
    }

    Cost getMaxFiniteCost() override
    {
        return MIN_COST;
    }

    void assign(int varIndex) override
    {
        if (connected(varIndex)) {
            deconnect(varIndex);
            nonassigned = nonassigned - 1;
            assert(nonassigned >= 0);

            if (universal()) {
                deconnect();
                return;
            }

            if (nonassigned <= 3) {
                deconnect();
                projectNary();
            }
        }
    }
    void increase(int varIndex) override {}
    void decrease(int varIndex) override {}
    void remove(int varIndex) override {}

    // propagates from scratch
    void propagate() override
    {
        //FIXME: synchronize current domains between master and slave problems at initialization?
        wcsp->revise(this);
        assigns();
        if (connected()) {
            problem->enforceUb();
            problem->propagate();
        }
        if (connected()) {
            negproblem->enforceUb();
            negproblem->propagate();
        }
        assert(problem->getLb() < ub);
        assert(negproblem->getLb() < -lb + negproblem->getNegativeLb() + UNIT_COST);
    }

    void print(ostream& os) override
    {
        os << endl
           << "WeightedWCSPConstraint(";
        int unassigned_ = 0;
        for (int i = 0; i < arity_; i++) {
            if (scope[i]->unassigned())
                unassigned_++;
            os << wcsp->getName(scope[i]->wcspIndex);
            if (i < arity_ - 1)
                os << ",";
        }
        os << ") in [" << lb << "," << ub << "[ ";
        if (ToulBar2::weightedDegree) {
            os << "/" << getConflictWeight();
            for (int i = 0; i < arity_; i++) {
                os << "," << conflictWeights[i];
            }
        }
        os << " arity: " << arity_;
        os << " unassigned: " << (int)nonassigned << "/" << unassigned_ << endl;
        os << *problem << endl << *negproblem << endl;
    }

    //TODO:
//    void dump(ostream& os, bool original = true) override
//    {
//    }

//    void dump_CFN(ostream& os, bool original = true) override
//    {
//    }

    friend void tb2setvalue(int wcspId, int varIndex, Value value, void* solver);
    friend void tb2removevalue(int wcspId, int varIndex, Value value, void* solver);
    friend void tb2setmin(int wcspId, int varIndex, Value value, void* solver);
    friend void tb2setmax(int wcspId, int varIndex, Value value, void* solver);
};

WCSP* WeightedCSPConstraint::MasterWeightedCSP = NULL;
map<int, WeightedCSPConstraint *> WeightedCSPConstraint::WeightedCSPConstraints;

void tb2setvalue(int wcspId, int varIndex, Value value, void* solver)
{
    assert(WeightedCSPConstraint::MasterWeightedCSP);
    Variable *masterVar = NULL;
    if (wcspId != WeightedCSPConstraint::MasterWeightedCSP->getIndex()) { // we came from a slave, wake up the master
        masterVar = WeightedCSPConstraint::WeightedCSPConstraints[wcspId]->getVar(varIndex);
        masterVar->assign(value);
        assert(solver && ((Solver *)solver)->getWCSP() == WeightedCSPConstraint::MasterWeightedCSP);
        setvalue(WeightedCSPConstraint::MasterWeightedCSP->getIndex(), masterVar->wcspIndex, value, solver);
    } else {
        // we came from the master
        masterVar = WeightedCSPConstraint::MasterWeightedCSP->getVar(varIndex);
    }
    if (ToulBar2::verbose >= 2)
        cout << "EVENT: x" << varIndex << "_" << wcspId << " = " << value << endl;
    for (auto gc: WeightedCSPConstraint::WeightedCSPConstraints) if (gc.second->connected()) {
        int varCtrIndex = gc.second->getIndex(masterVar);
        if (varCtrIndex != -1) {
            if (wcspId != gc.second->problem->getIndex()) { // do not reenter inside the same problem as the one we came
                gc.second->problem->getVar(varCtrIndex)->assign(value);
            }
            if (wcspId != gc.second->negproblem->getIndex()) { // do not reenter inside the same problem as the one we came
                gc.second->negproblem->getVar(varCtrIndex)->assign(value);
            }
        }
    }
}

void tb2removevalue(int wcspId, int varIndex, Value value, void* solver)
{
    assert(WeightedCSPConstraint::MasterWeightedCSP);
    Variable *masterVar = NULL;
    if (wcspId != WeightedCSPConstraint::MasterWeightedCSP->getIndex()) { // we came from a slave, wake up the master
        masterVar = WeightedCSPConstraint::WeightedCSPConstraints[wcspId]->getVar(varIndex);
        masterVar->remove(value);
    } else {
        // we came from the master
        masterVar = WeightedCSPConstraint::MasterWeightedCSP->getVar(varIndex);
    }
    if (ToulBar2::verbose >= 2)
        cout << "EVENT: x" << varIndex << "_" << wcspId << " != " << value << endl;
    for (auto gc: WeightedCSPConstraint::WeightedCSPConstraints) if (gc.second->connected()) {
        int varCtrIndex = gc.second->getIndex(masterVar);
        if (varCtrIndex != -1) {
            if (wcspId != gc.second->problem->getIndex()) { // do not reenter inside the same problem as the one we came
                gc.second->problem->getVar(varCtrIndex)->remove(value);
            }
            if (wcspId != gc.second->negproblem->getIndex()) { // do not reenter inside the same problem as the one we came
                gc.second->negproblem->getVar(varCtrIndex)->remove(value);
            }
        }
    }
}

void tb2setmin(int wcspId, int varIndex, Value value, void* solver)
{
    assert(WeightedCSPConstraint::MasterWeightedCSP);
    Variable *masterVar = NULL;
    if (wcspId != WeightedCSPConstraint::MasterWeightedCSP->getIndex()) { // we came from a slave, wake up the master
        masterVar = WeightedCSPConstraint::WeightedCSPConstraints[wcspId]->getVar(varIndex);
        masterVar->increase(value);
    } else {
        // we came from the master
        masterVar = WeightedCSPConstraint::MasterWeightedCSP->getVar(varIndex);
    }
    if (ToulBar2::verbose >= 2)
        cout << "EVENT: x" << varIndex << "_" << wcspId << " >= " << value << endl;
    for (auto gc: WeightedCSPConstraint::WeightedCSPConstraints) if (gc.second->connected()) {
        int varCtrIndex = gc.second->getIndex(masterVar);
        if (varCtrIndex != -1) {
            if (wcspId != gc.second->problem->getIndex()) { // do not reenter inside the same problem as the one we came
                gc.second->problem->getVar(varCtrIndex)->increase(value);
            }
            if (wcspId != gc.second->negproblem->getIndex()) { // do not reenter inside the same problem as the one we came
                gc.second->negproblem->getVar(varCtrIndex)->increase(value);
            }
        }
    }
}

void tb2setmax(int wcspId, int varIndex, Value value, void* solver)
{
    assert(WeightedCSPConstraint::MasterWeightedCSP);
    Variable *masterVar = NULL;
    if (wcspId != WeightedCSPConstraint::MasterWeightedCSP->getIndex()) { // we came from a slave, wake up the master
        masterVar = WeightedCSPConstraint::WeightedCSPConstraints[wcspId]->getVar(varIndex);
        masterVar->decrease(value);
    } else {
        // we came from the master
        masterVar = WeightedCSPConstraint::MasterWeightedCSP->getVar(varIndex);
    }
    if (ToulBar2::verbose >= 2)
        cout << "EVENT: x" << varIndex << "_" << wcspId << " <= " << value << endl;
    for (auto gc: WeightedCSPConstraint::WeightedCSPConstraints) if (gc.second->connected()) {
        int varCtrIndex = gc.second->getIndex(masterVar);
        if (varCtrIndex != -1) {
            if (wcspId != gc.second->problem->getIndex()) { // do not reenter inside the same problem as the one we came
                gc.second->problem->getVar(varCtrIndex)->decrease(value);
            }
            if (wcspId != gc.second->negproblem->getIndex()) { // do not reenter inside the same problem as the one we came
                gc.second->negproblem->getVar(varCtrIndex)->decrease(value);
            }
        }
    }
}
#endif /*TB2GLOBALWCSP_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
