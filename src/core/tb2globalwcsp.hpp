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
    Cost negcost; // sum of negative cost shift from slave problem and its negation form
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
        , negcost(MIN_COST)
        , problem(problem_in)
        , negproblem(negproblem_in)
        , nonassigned(arity_in)
    {
        assert(arity_in == (int)problem_in->numberOfVariables());
        assert(arity_in == (int)negproblem_in->numberOfVariables());
        if (lb >= ub) {
            cerr << "Wrong bounds in WeightedCSPConstraint: " << lb << " < " << ub << endl;
            throw WrongFileFormat();
        }
        for (int i = 0; i < arity_in; i++) {
            assert(!problem || scope_in[i]->getDomainInitSize() == ((EnumeratedVariable *)problem->getVar(i))->getDomainInitSize());
            assert(!negproblem || scope_in[i]->getDomainInitSize() == ((EnumeratedVariable *)negproblem->getVar(i))->getDomainInitSize());
            varIndexes.push_back(i);
            newValues.push_back(scope_in[i]->getInf());
            conflictWeights.push_back(0);
        }
        ToulBar2::setvalue = ::tb2setvalue;
        ToulBar2::removevalue = ::tb2removevalue;
        ToulBar2::setmin = ::tb2setmin;
        ToulBar2::setmax = ::tb2setmax;
        assert(MasterWeightedCSP == NULL || MasterWeightedCSP == wcsp); //FIXME: the slave problem cannot contain a WeightedCSPConstraint inside!
        MasterWeightedCSP = wcsp;
        if (problem) {
            negcost += problem->getNegativeLb();
            WeightedCSPConstraints[problem->getIndex()] = this;
            problem->setSolver(wcsp->getSolver()); // force slave problems to use the same solver as the master
            problem->updateUb(ub);
            problem->enforceUb();
        }
        if (negproblem) {
            negcost += negproblem->getNegativeLb();
            WeightedCSPConstraints[negproblem->getIndex()] = this;
            negproblem->setSolver(wcsp->getSolver());
            negproblem->updateUb(-lb + negcost + UNIT_COST);
            negproblem->enforceUb();
        }
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
        if ((!problem || problem->getLb() >= lb) && (!negproblem || -(negproblem->getLb() - negcost) < ub)) {
            return true;
        } else {
            return false;
        }
    }

    Cost eval(const Tuple& s) override
    {
        assert((int)s.size() == arity_);
        for (int i=0; i < arity_; i++) {
            newValues[i] = ((EnumeratedVariable *)getVar(i))->toValue(s[i]);
        }
        ToulBar2::setvalue = NULL;
        ToulBar2::removevalue = NULL;
        ToulBar2::setmin = NULL;
        ToulBar2::setmax = NULL;
        bool rasps = ToulBar2::RASPS;
        ToulBar2::RASPS = false;
        int depth = Store::getDepth();
        bool unsat = false;
        try {
            Store::store();
            if (problem) problem->assignLS(varIndexes, newValues); // throw a Contradiction if unsatisfied
            if (negproblem) negproblem->assignLS(varIndexes, newValues); // idem
        } catch (const Contradiction&) {
            if (problem) problem->whenContradiction();
            if (negproblem) negproblem->whenContradiction();
            unsat = true;
        }
        Store::restore(depth);
        ToulBar2::setvalue = ::tb2setvalue;
        ToulBar2::removevalue = ::tb2removevalue;
        ToulBar2::setmin = ::tb2setmin;
        ToulBar2::setmax = ::tb2setmax;
        ToulBar2::RASPS = rasps;
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
        if (problem) {
            for (unsigned int c=0; c < problem->numberOfConstraints(); c++) {
                if (problem->getCtr(c)->connected()) {
                    res += problem->getCtr(c)->getTightness();
                }
            }
            return res / problem->numberOfConnectedConstraints();
        } else if (negproblem) {
            for (unsigned int c=0; c < negproblem->numberOfConstraints(); c++) {
                if (negproblem->getCtr(c)->connected()) {
                    res += negproblem->getCtr(c)->getTightness();
                }
            }
            return res / negproblem->numberOfConnectedConstraints();
        } else {
            return 1.;
        }
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

            if (nonassigned <= NARYPROJECTIONSIZE) {
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
        if (problem) problem->enforceUb();
        if (negproblem) negproblem->enforceUb();
        assigns();
        if (connected()) {
            if (problem) problem->propagate();
            if (connected()) {
                if (negproblem) negproblem->propagate();
            }
        }
        assert(!problem || problem->getLb() < ub);
        assert(!negproblem || negproblem->getLb() < -lb + negcost + UNIT_COST);
    }

    void print(ostream& os) override
    {
        os << this << "WeightedWCSPConstraint(";
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
        if (problem) os << *problem << endl;
        if (negproblem) os << *negproblem << endl;
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
    assert(wcspId == WeightedCSPConstraint::MasterWeightedCSP->getIndex() || WeightedCSPConstraint::WeightedCSPConstraints.find(wcspId) != WeightedCSPConstraint::WeightedCSPConstraints.end());
    Variable *masterVar = NULL;
    if (wcspId != WeightedCSPConstraint::MasterWeightedCSP->getIndex()) { // we came from a slave, wake up the master
        masterVar = WeightedCSPConstraint::WeightedCSPConstraints[wcspId]->getVar(varIndex);
        masterVar->assign(value);
    } else {
        // we came from the master
        masterVar = WeightedCSPConstraint::MasterWeightedCSP->getVar(varIndex);
        setvalue(WeightedCSPConstraint::MasterWeightedCSP->getIndex(), masterVar->wcspIndex, value, WeightedCSPConstraint::MasterWeightedCSP->getSolver());
    }
    if (ToulBar2::verbose >= 2)
        cout << "EVENT: x" << varIndex << "_" << wcspId << " = " << value << endl;
    bool rasps = ToulBar2::RASPS;
    ToulBar2::RASPS = false;
    for (auto gc: WeightedCSPConstraint::WeightedCSPConstraints) if (gc.second->connected()) {
        int varCtrIndex = gc.second->getIndex(masterVar);
        if (varCtrIndex != -1) { // only for slave problems which are concerned by this variable
            if (gc.second->problem && wcspId != gc.second->problem->getIndex()) { // do not reenter inside the same problem as the one we came
                gc.second->problem->getVar(varCtrIndex)->assign(value);
            }
            if (gc.second->negproblem && wcspId != gc.second->negproblem->getIndex()) { // do not reenter inside the same problem as the one we came
                gc.second->negproblem->getVar(varCtrIndex)->assign(value);
            }
        }
    }
    ToulBar2::RASPS = rasps;
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
    bool rasps = ToulBar2::RASPS;
    ToulBar2::RASPS = false;
    for (auto gc: WeightedCSPConstraint::WeightedCSPConstraints) if (gc.second->connected()) {
        int varCtrIndex = gc.second->getIndex(masterVar);
        if (varCtrIndex != -1) {
            if (gc.second->problem && wcspId != gc.second->problem->getIndex()) { // do not reenter inside the same problem as the one we came
                gc.second->problem->getVar(varCtrIndex)->remove(value);
            }
            if (gc.second->negproblem && wcspId != gc.second->negproblem->getIndex()) { // do not reenter inside the same problem as the one we came
                gc.second->negproblem->getVar(varCtrIndex)->remove(value);
            }
        }
    }
    ToulBar2::RASPS = rasps;
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
    bool rasps = ToulBar2::RASPS;
    ToulBar2::RASPS = false;
    for (auto gc: WeightedCSPConstraint::WeightedCSPConstraints) if (gc.second->connected()) {
        int varCtrIndex = gc.second->getIndex(masterVar);
        if (varCtrIndex != -1) {
            if (gc.second->problem && wcspId != gc.second->problem->getIndex()) { // do not reenter inside the same problem as the one we came
                gc.second->problem->getVar(varCtrIndex)->increase(value);
            }
            if (gc.second->negproblem && wcspId != gc.second->negproblem->getIndex()) { // do not reenter inside the same problem as the one we came
                gc.second->negproblem->getVar(varCtrIndex)->increase(value);
            }
        }
    }
    ToulBar2::RASPS = rasps;
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
    bool rasps = ToulBar2::RASPS;
    ToulBar2::RASPS = false;
    for (auto gc: WeightedCSPConstraint::WeightedCSPConstraints) if (gc.second->connected()) {
        int varCtrIndex = gc.second->getIndex(masterVar);
        if (varCtrIndex != -1) {
            if (gc.second->problem && wcspId != gc.second->problem->getIndex()) { // do not reenter inside the same problem as the one we came
                gc.second->problem->getVar(varCtrIndex)->decrease(value);
            }
            if (gc.second->negproblem && wcspId != gc.second->negproblem->getIndex()) { // do not reenter inside the same problem as the one we came
                gc.second->negproblem->getVar(varCtrIndex)->decrease(value);
            }
        }
    }
    ToulBar2::RASPS = rasps;
}
#endif /*TB2GLOBALWCSP_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
