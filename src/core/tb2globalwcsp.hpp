#ifndef TB2GLOBALWCSP_HPP_
#define TB2GLOBALWCSP_HPP_

#include "tb2abstractconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2naryconstr.hpp"
#include "tb2vac.hpp"
#include "tb2wcsp.hpp"

extern void setvalue(int wcspId, int varIndex, Value value, void* solver);
extern void tb2setvalue(int wcspId, int varIndex, Value value, void* solver);
extern void tb2removevalue(int wcspId, int varIndex, Value value, void* solver);
extern void tb2setmin(int wcspId, int varIndex, Value value, void* solver);
extern void tb2setmax(int wcspId, int varIndex, Value value, void* solver);

class WeightedCSPConstraint : public AbstractNaryConstraint {
    bool isfinite; // true if any complete assignment of the input problem (or negproblem), before enforcing lb and ub, has a finite cost or if it is forbidden then there is another redundant constraint in the master problem which forbids the same assigment
    bool strongDuality; // if true then it assumes the propagation is complete when all channeling variables in the scope are assigned and the semantic of the constraint enforces that the optimum on the remaining variables is between lb and ub.
    Cost lb; // encapsulated slave problem lower bound hard constraint (must be greater or equal to this bound)
    Cost ub; // encapsulated slave problem upper bound hard constraint (must be strictly less than this bound)
    Cost negCost; // sum of cost shifts from slave problem and from its negative form
    Cost top; // forbidden cost returned if evaluated as unsatisfied
    WCSP* problem; // encapsulated slave problem
    WCSP* negproblem; // encapsulated slave problem in negative form (should be equivalent to -problem)
    WCSP* original_problem; // pointer to the wcsp given as input, necessary to compute the cost of a solution
    WCSP* original_negproblem; // pointer to the wcsp given as input, necessary to compute the cost of a solution
    vector<int> varIndexes; // copy of scope using integer identifiers inside slave problem (should be equal to [0, 1, 2, ..., arity-1])
    vector<Value> newValues; // used to convert Tuples into variable assignments
    vector<Long> conflictWeights; // used by weighted degree heuristics

public:
    static WCSP* MasterWeightedCSP; // Master problem used by value and variable ordering heuristics
    static map<int, WeightedCSPConstraint*> WeightedCSPConstraints;
    static bool _protected_;
    static int preprocessFunctional;
    static int elimDegree;
    static int elimDegree_preprocessing;
    static int elimDegree_;
    static int elimDegree_preprocessing_;
    static int DEE;
    static int DEE_;
    static bool FullEAC;
    static bool RASPS;
    static int useRASPS;
    static void protect(bool master = true) ///< \brief deactivate some preprocessing/propagation features not compatible with our channeling mechanism
    {
        assert(!_protected_);
        if (master) {
            preprocessFunctional = ToulBar2::preprocessFunctional;
            elimDegree = ToulBar2::elimDegree;
            elimDegree_preprocessing = ToulBar2::elimDegree_preprocessing;
            elimDegree_ = ToulBar2::elimDegree_;
            elimDegree_preprocessing_ = ToulBar2::elimDegree_preprocessing_;
            DEE = ToulBar2::DEE;
            DEE_ = ToulBar2::DEE_;
            FullEAC = ToulBar2::FullEAC;
            RASPS = ToulBar2::RASPS;
            useRASPS = ToulBar2::useRASPS;
        }

        _protected_ = true;
        ToulBar2::preprocessFunctional = 0;
        ToulBar2::elimDegree = -1;
        ToulBar2::elimDegree_preprocessing = -1;
        ToulBar2::elimDegree_ = -1;
        ToulBar2::elimDegree_preprocessing_ = -1;
        ToulBar2::DEE = 0;
        ToulBar2::DEE_ = 0;
        ToulBar2::FullEAC = false; // FIXME: remove this restriction?
        ToulBar2::RASPS = false;
        ToulBar2::useRASPS = 0;
    }
    static void unprotect() ///< \brief reactivate preprocessing/propagation features
    {
        if (_protected_) {
            _protected_ = false;
            ToulBar2::preprocessFunctional = preprocessFunctional;
            ToulBar2::elimDegree = elimDegree;
            ToulBar2::elimDegree_preprocessing = elimDegree_preprocessing;
            ToulBar2::elimDegree_ = elimDegree_;
            ToulBar2::elimDegree_preprocessing_ = elimDegree_preprocessing_;
            ToulBar2::DEE = DEE;
            ToulBar2::DEE_ = DEE_;
            ToulBar2::FullEAC = FullEAC;
            ToulBar2::RASPS = RASPS;
            ToulBar2::useRASPS = useRASPS;
        }
    }

    // TODO: add local NARYPROJECTIONSIZE parameter
    WeightedCSPConstraint(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in, WCSP* problem_in, WCSP* negproblem_in, Cost lb_in, Cost ub_in, bool duplicateHard = false, bool strongDuality_ = false)
        : AbstractNaryConstraint(wcsp, scope_in, arity_in)
        , isfinite(true)
        , strongDuality(strongDuality_)
        , lb(lb_in)
        , ub(ub_in)
        , negCost(MIN_COST)
        , top(MAX_COST)
        , problem(problem_in)
        , negproblem(negproblem_in)
        , original_problem(problem_in)
        , original_negproblem(negproblem_in)
    {
        assert(problem);
        assert(!problem || arity_ == (int)problem->numberOfVariables());
        assert(!negproblem || arity_ == (int)negproblem->numberOfVariables());
        if (lb >= ub) {
            cerr << "Wrong bounds in WeightedCSPConstraint: " << lb << " < " << ub << endl;
            throw WrongFileFormat();
        }
        for (int i = 0; i < arity_in; i++) {
            assert(!problem || scope_in[i]->getDomainInitSize() == ((EnumeratedVariable*)problem->getVar(i))->getDomainInitSize());
            assert(!negproblem || scope_in[i]->getDomainInitSize() == ((EnumeratedVariable*)negproblem->getVar(i))->getDomainInitSize());
            varIndexes.push_back(i);
            newValues.push_back(scope_in[i]->getInf());
            conflictWeights.push_back(0);
        }
        assert(ToulBar2::setvalue == NULL || ToulBar2::setvalue == ::tb2setvalue);
        ToulBar2::setvalue = ::tb2setvalue;
        ToulBar2::removevalue = ::tb2removevalue;
        ToulBar2::setmin = ::tb2setmin;
        ToulBar2::setmax = ::tb2setmax;

        if (MasterWeightedCSP != NULL && MasterWeightedCSP != wcsp) {
            WeightedCSPConstraints.clear();
        }
        MasterWeightedCSP = wcsp; // FIXME: the slave problem should not contain a WeightedCSPConstraint inside!
        if (problem) {
            negCost += problem->getNegativeLb();
            WeightedCSPConstraints[problem->getIndex()] = this;
            problem->setSolver(wcsp->getSolver()); // force slave problems to use the same solver as the master
            if (!duplicateHard && !problem->isfinite())
                isfinite = false;

            Cost summaxcost = problem->getLb() + UNIT_COST;
            for (unsigned int i = 0; i < problem->numberOfVariables(); i++) {
                if (problem->enumerated(i)) {
                    Cost maxcost = MIN_COST;
                    EnumeratedVariable* var = (EnumeratedVariable*)problem->getVar(i);
                    for (EnumeratedVariable::iterator iter = var->begin(); iter != var->end(); ++iter) {
                        if (var->getCost(*iter) > maxcost)
                            maxcost = var->getCost(*iter);
                    }
                    summaxcost += maxcost;
                } else {
                    summaxcost += max(problem->getUnaryCost(i, problem->getInf(i)), problem->getUnaryCost(i, problem->getSup(i)));
                }
            }

            if (isfinite && problem->finiteUb() <= ub) {
                WeightedCSPConstraints.erase(problem->getIndex());
                problem = NULL; // no need to check ub anymore
            } else if (problem->numberOfConstraints() == 0 || (duplicateHard && problem->finiteUb() == summaxcost)) { // special case where there are only unary cost functions or there are only hard constraints already duplicated in the master problem
                vector<int> thescope;
                string params = to_string(-ub + problem->getLb() + UNIT_COST);
                for (int i = 0; i < arity_; i++) {
                    thescope.push_back(scope[i]->wcspIndex);
                    vector<pair<Value, Cost>> vc = problem->getEnumDomainAndCost(i);
                    params += to_string(" ") + to_string(vc.size());
                    for (unsigned int j = 0; j < vc.size(); j++) {
                        params += to_string(" ") + to_string(vc[j].first) + to_string(" ") + to_string(-(vc[j].second));
                    }
                }
                WeightedCSPConstraints.erase(problem->getIndex());
                problem = NULL;
                wcsp->postKnapsackConstraint(thescope, params, false, true, false);
            } else {
                problem->updateUb(ub);
                problem->enforceUb();
            }
        }
        if (negproblem) {
            negCost += negproblem->getNegativeLb();
            WeightedCSPConstraints[negproblem->getIndex()] = this;
            negproblem->setSolver(wcsp->getSolver());
            if (!duplicateHard && !negproblem->isfinite())
                isfinite = false;

            Cost summaxcost = negproblem->getLb() + UNIT_COST;
            for (unsigned int i = 0; i < negproblem->numberOfVariables(); i++) {
                if (negproblem->enumerated(i)) {
                    Cost maxcost = MIN_COST;
                    EnumeratedVariable* var = (EnumeratedVariable*)negproblem->getVar(i);
                    for (EnumeratedVariable::iterator iter = var->begin(); iter != var->end(); ++iter) {
                        if (var->getCost(*iter) > maxcost)
                            maxcost = var->getCost(*iter);
                    }
                    summaxcost += maxcost;
                } else {
                    summaxcost += max(negproblem->getUnaryCost(i, negproblem->getInf(i)), negproblem->getUnaryCost(i, negproblem->getSup(i)));
                }
            }

            if (isfinite && (negproblem->finiteUb() <= (-lb + negCost + UNIT_COST))) {
                WeightedCSPConstraints.erase(negproblem->getIndex());
                negproblem = NULL; // no need to check lb anymore
            } else if (negproblem->numberOfConstraints() == 0 || (duplicateHard && negproblem->finiteUb() == summaxcost)) { // special case where there are only unary cost functions or there are only hard constraints already duplicated in the master problem
                vector<int> thescope;
                string params = to_string(lb - negCost + negproblem->getLb());
                for (int i = 0; i < arity_; i++) {
                    thescope.push_back(scope[i]->wcspIndex);
                    vector<pair<Value, Cost>> vc = negproblem->getEnumDomainAndCost(i);
                    params += to_string(" ") + to_string(vc.size());
                    for (unsigned int j = 0; j < vc.size(); j++) {
                        params += to_string(" ") + to_string(vc[j].first) + to_string(" ") + to_string(-(vc[j].second));
                    }
                }
                WeightedCSPConstraints.erase(negproblem->getIndex());
                negproblem = NULL;
                wcsp->postKnapsackConstraint(thescope, params, false, true, false);
            } else {
                negproblem->updateUb(-lb + negCost + UNIT_COST);
                negproblem->enforceUb();
            }
        }
        if (!problem && !negproblem) {
            deconnect();
            if (WeightedCSPConstraints.size() == 0) { // we do not need anymore hook functions
                MasterWeightedCSP = NULL;
                ToulBar2::setvalue = NULL;
                ToulBar2::removevalue = NULL;
                ToulBar2::setmin = NULL;
                ToulBar2::setmax = NULL;
            }
        } else {
            assert(connected());
            if (problem) {
                protect(true);
                try {
                    for (int i = 0; i < arity_in && connected(); i++) {
                        for (unsigned int j = 0; j < scope_in[i]->getDomainInitSize() && connected(); j++) {
                            Value a = scope_in[i]->toValue(j);
                            EnumeratedVariable* vari = ((EnumeratedVariable*)problem->getVar(i));
#ifndef NDEBUG
                            Value b = vari->toValue(j);
                            assert(a == b);
#endif
                            if (scope_in[i]->cannotbe(a) && vari->canbe(a)) {
                                problem->remove(i, a);
                            }
                        }
                    }
                    problem->propagate(true); // preprocessing();
                } catch (const Contradiction&) {
                    deconnect();
                    clearPtrReferences();
                    throw;
                }
                unprotect();
                for (int i = 0; i < arity_in && connected(); i++) {
                    for (unsigned int j = 0; j < scope_in[i]->getDomainInitSize() && connected(); j++) {
                        Value a = scope_in[i]->toValue(j);
                        EnumeratedVariable* vari = ((EnumeratedVariable*)problem->getVar(i));
                        if (scope_in[i]->canbe(a) && vari->cannotbe(a)) {
                            scope_in[i]->remove(a);
                        }
                    }
                }
            }
            if (connected() && negproblem) {
                protect(true);
                try {
                    for (int i = 0; i < arity_in && connected(); i++) {
                        for (unsigned int j = 0; j < scope_in[i]->getDomainInitSize() && connected(); j++) {
                            Value a = scope_in[i]->toValue(j);
                            EnumeratedVariable* vari = ((EnumeratedVariable*)negproblem->getVar(i));
#ifndef NDEBUG
                            Value b = vari->toValue(j);
                            assert(a == b);
#endif
                            if (scope_in[i]->cannotbe(a) && vari->canbe(a)) {
                                negproblem->remove(i, a);
                            }
                        }
                    }
                    negproblem->propagate(true); // preprocessing();
                } catch (const Contradiction&) {
                    deconnect();
                    clearPtrReferences();
                    throw;
                }
                unprotect();
                for (int i = 0; i < arity_in && connected(); i++) {
                    for (unsigned int j = 0; j < scope_in[i]->getDomainInitSize() && connected(); j++) {
                        Value a = scope_in[i]->toValue(j);
                        EnumeratedVariable* vari = ((EnumeratedVariable*)negproblem->getVar(i));
                        if (scope_in[i]->canbe(a) && vari->cannotbe(a)) {
                            scope_in[i]->remove(a);
                        }
                    }
                }
            }
        }
    }

    virtual ~WeightedCSPConstraint()
    {
        clearPtrReferences();
    }

    bool extension() const FINAL { return false; } // this is not a cost function represented by an exhaustive table of costs

    Long getConflictWeight() const override { return Constraint::getConflictWeight(); }
    Long getConflictWeight(int varIndex) const override
    {
        assert(varIndex >= 0);
        assert(varIndex < arity_);
        return conflictWeights[varIndex] + Constraint::getConflictWeight() + ((problem && ToulBar2::weightedDegree > 0) ? problem->getWeightedDegree(varIndex) : 0) + ((negproblem && ToulBar2::weightedDegree > 0) ? negproblem->getWeightedDegree(varIndex) : 0);
    }
    void incConflictWeight(Constraint* from) override
    {
        assert(from != NULL);
        if (from == this) {
            if (getNonAssigned() == arity_ || deconnected()) {
                Constraint::incConflictWeight(1);
            } else {
                for (int i = 0; i < arity_; i++) {
                    if (connected(i)) { // It will favor diversification by branching on new variables
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
        if (problem) {
            problem->resetWeightedDegree();
        }
        if (negproblem) {
            negproblem->resetWeightedDegree();
        }
    }

    Cost getUnaryCost(int wcspIndex, Value value, int sign = 1)
    {
        if (scope_inv.find(wcspIndex) == scope_inv.end())
            return MAX_COST;
        int varIndex = scope_inv[wcspIndex];
        if (problem && negproblem) {
            if (sign >= 0) {
                return problem->getUnaryCost(varIndex, value);
            } else {
                return negproblem->getUnaryCost(varIndex, value);
            }
        } else if (problem) {
            return problem->getUnaryCost(varIndex, value);
        } else if (negproblem) {
            return negproblem->getUnaryCost(varIndex, value);
        } else {
            return MAX_COST;
        }
    }

    Value getSupport(int wcspIndex, int sign = 0, bool mingap = true) ///< \brief returns EAC support value of variable index wcspIndex in problem if sign>0 or negproblem if sign<0 or heuristically chosen between both if sign==0
    {
        if (scope_inv.find(wcspIndex) == scope_inv.end())
            return WRONG_VAL;
        int varIndex = scope_inv[wcspIndex];
        if (problem && negproblem) {
            if (sign > 0 || (sign == 0 && (!mingap != (max(MIN_COST, lb - problem->getLb()) <= max(MIN_COST, -ub + negCost - negproblem->getLb() + UNIT_COST))))) {
                return problem->getSupport(varIndex);
            } else {
                return negproblem->getSupport(varIndex);
            }
        } else if (problem) {
            return problem->getSupport(varIndex);
        } else if (negproblem) {
            return negproblem->getSupport(varIndex);
        } else {
            return WRONG_VAL;
        }
    }

    bool universal(Cost zero = MIN_COST) override
    {
        if (isfinite && problem && negproblem && problem->getLb() >= lb && negproblem->getLb() > -ub + negCost) {
            return true;
        } else {
            return false;
        }
    }

    /// \brief returns true if all remaining (unassigned) variables are only connected to this global constraint
    bool canbeDeconnected() const
    {
        for (int i = 0; i < arity_; i++) {
            if (getVar(i)->unassigned() && getVar(i)->getDegree() > 1) {
                return false;
            }
        }
        return true;
    }

    Cost eval(const Tuple& s) override
    {
        assert((int)s.size() == arity_);
        for (int i = 0; i < arity_; i++) {
            newValues[i] = ((EnumeratedVariable*)getVar(i))->toValue(s[i]);
        }
        ToulBar2::setvalue = NULL;
        ToulBar2::removevalue = NULL;
        ToulBar2::setmin = NULL;
        ToulBar2::setmax = NULL;
        int weightedDegree = ToulBar2::weightedDegree;
        ToulBar2::weightedDegree = 0; // do not update weighted degrees inside slave problems
        protect();
        int depth = Store::getDepth();
        bool unsat = false;
        try {
            Store::store();
            if (problem && original_problem) {
                original_problem->enforceUb();
                assert(original_problem->isactivatePropagate());
                original_problem->assignLS(varIndexes, newValues); // may-be throw a Contradiction if it violates ub
                if (original_problem->getLb() < lb || original_problem->getLb() >= ub) { // checks if the solution violates lb or ub
                    unsat = true;
                }
            } else if (negproblem && original_negproblem) {
                original_negproblem->enforceUb();
                assert(original_negproblem->isactivatePropagate());
                original_negproblem->assignLS(varIndexes, newValues); // may-be throw a Contradiction if it violates lb
                if (original_negproblem->getLb() <= -ub + negCost || original_negproblem->getLb() > -lb + negCost) { // checks if the solution violates ub or lb
                    unsat = true;
                }
            }
        } catch (const Contradiction&) {
            if (original_problem)
                original_problem->whenContradiction();
            if (original_negproblem)
                original_negproblem->whenContradiction();
            unsat = true;
        }
        Store::restore(depth);
        unprotect();
        ToulBar2::weightedDegree = weightedDegree;
        ToulBar2::setvalue = ::tb2setvalue;
        ToulBar2::removevalue = ::tb2removevalue;
        ToulBar2::setmin = ::tb2setmin;
        ToulBar2::setmax = ::tb2setmax;
        if (unsat) {
            return top;
        } else {
            return MIN_COST;
        }
    }

    /*!
     * \brief compute the sum of the cost functions in the constraint from a solution of the wcsp (variables are indexed from the main cfn)
     * \param solution the values of the variables in the solution returned by the solver
     * \return the cost of the solution for the cfn in the constraint expressed as a cost for the problem_in cfn
     */
    Cost computeSolutionCost(vector<Value>& solution)
    {
        assert(solution.size() == wcsp->numberOfVariables());
        Cost cost = MIN_COST;
        for (int i = 0; i < arity_; i++) {
            newValues[i] = solution[((EnumeratedVariable*)getVar(i))->wcspIndex];
        }
        externalevent oldSetValue = ToulBar2::setvalue;
        externalevent oldRemoveValue = ToulBar2::removevalue;
        externalevent oldSetMin = ToulBar2::setmin;
        externalevent oldSetMax = ToulBar2::setmax;
        ToulBar2::setvalue = NULL;
        ToulBar2::removevalue = NULL;
        ToulBar2::setmin = NULL;
        ToulBar2::setmax = NULL;
        int weightedDegree = ToulBar2::weightedDegree;
        ToulBar2::weightedDegree = 0; // do not update weighted degrees inside slave problems
        protect();
        int depth = Store::getDepth();
        try {
            Store::store();
            if (original_problem) {
                assert(original_problem->isactivatePropagate());
                original_problem->assignLS(varIndexes, newValues);
                cost = original_problem->getLb();
            } else if (original_negproblem) {
                assert(original_negproblem->isactivatePropagate());
                original_negproblem->assignLS(varIndexes, newValues);
                cost = negCost - original_negproblem->getLb();
            }
        } catch (const Contradiction&) {
            if (original_problem)
                original_problem->whenContradiction();
            if (original_negproblem)
                original_negproblem->whenContradiction();
            cost = MAX_COST;
        }
        Store::restore(depth);
        unprotect();
        ToulBar2::weightedDegree = weightedDegree;
        ToulBar2::setvalue = oldSetValue;
        ToulBar2::removevalue = oldRemoveValue;
        ToulBar2::setmin = oldSetMin;
        ToulBar2::setmax = oldSetMax;
        return cost;
    }

    void resetTightness() override
    {
        Constraint::resetTightness();
        if (problem) {
            problem->resetTightness();
        }
        if (negproblem) {
            negproblem->resetTightness();
        }
    }

    double computeTightness() override
    {
        double res = 0.;
        if (problem && problem->numberOfConnectedConstraints() > 0) {
            for (unsigned int c = 0; c < problem->numberOfConstraints(); c++) {
                if (problem->getCtr(c)->connected() && !problem->getCtr(c)->isSep()) {
                    res += problem->getCtr(c)->getTightness();
                }
            }
            for (int i = 0; i < problem->getElimBinOrder(); i++) {
                BinaryConstraint *c = (BinaryConstraint *)problem->getElimBinCtr(i);
                if (c->connected() && !c->isSep()) {
                    res += c->getTightness();
                }
            }
            for (int i = 0; i < problem->getElimTernOrder(); i++) {
                TernaryConstraint *c = (TernaryConstraint *)problem->getElimTernCtr(i);
                if (c->connected() && !c->isSep()) {
                    res += c->getTightness();
                }
            }
            return res / problem->numberOfConnectedConstraints();
        } else if (negproblem && negproblem->numberOfConnectedConstraints() > 0) {
            for (unsigned int c = 0; c < negproblem->numberOfConstraints(); c++) {
                if (negproblem->getCtr(c)->connected() && !negproblem->getCtr(c)->isSep()) {
                    res += negproblem->getCtr(c)->getTightness();
                }
            }
            for (int i = 0; i < negproblem->getElimBinOrder(); i++) {
                BinaryConstraint *c = (BinaryConstraint *)negproblem->getElimBinCtr(i);
                if (c->connected() && !c->isSep()) {
                    res += c->getTightness();
                }
            }
            for (int i = 0; i < negproblem->getElimTernOrder(); i++) {
                TernaryConstraint *c = (TernaryConstraint *)negproblem->getElimTernCtr(i);
                if (c->connected() && !c->isSep()) {
                    res += c->getTightness();
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
    void setInfiniteCost(Cost ub) override
    {
        Cost mult_ub = ((wcsp->getUb() < (MAX_COST / MEDIUM_COST)) ? (max(LARGE_COST, wcsp->getUb() * MEDIUM_COST)) : wcsp->getUb());
        if (CUT(top, ub)) {
            top = mult_ub;
        }
    }

    void assign(int varIndex) override
    {
        if ((problem && !problem->isactivatePropagate()) || (negproblem && !negproblem->isactivatePropagate()))
            return; // wait until propagation is done for each subproblem before trying to deconnect or project
        if (connected(varIndex)) {
            deconnect(varIndex);
            assert(getNonAssigned() >= 0);

            if (universal()) {
                deconnect();
                return;
            }

            if (getNonAssigned() <= NARYPROJECTIONSIZE && (getNonAssigned() < 3 || maxInitDomSize <= NARYPROJECTION3MAXDOMSIZE || prodInitDomSize <= NARYPROJECTION3PRODDOMSIZE) && (!strongDuality || getNonAssigned() == 0)) {
                deconnect();
                projectNary();
            } else {
                propagate();
            }
        }
    }
    //    void increase(int varIndex) override {}
    //    void decrease(int varIndex) override {}
    //    void remove(int varIndex) override {}
    bool verify() override
    {
        for (int i = 0; i < arity_; i++) {
            vector<Value> vals = wcsp->getEnumDomain(getVar(i)->wcspIndex);
            if (problem) {
                vector<Value> vals2 = problem->getEnumDomain(i);
                if (vals2 != vals) {
                    cout << "Error WeightedCSPConstraint(" << problem->getIndex() << "): wrong domain values " << vals2 << " in variable " << *getVar(i) << endl;
                    return false;
                }
            }
            if (negproblem) {
                vector<Value> vals2 = negproblem->getEnumDomain(i);
                if (vals2 != vals) {
                    cout << "Error WeightedCSPConstraint(" << negproblem->getIndex() << "): wrong domain values " << vals2 << " in variable " << *getVar(i) << endl;
                    return false;
                }
            }
        }
        return (!problem || problem->verify()) && (!negproblem || negproblem->verify());
    }

    // propagates from scratch
    void propagate() override
    {
        if (ToulBar2::dumpWCSP % 2) // do not propagate if problem is dumped before preprocessing
            return;
        // FIXME: synchronize current domains between master and slave problems at initialization?
        wcsp->revise(this);
        if (problem) {
            problem->enforceUb();
            if (!problem->isactivatePropagate())
                return; // do not propagate during recursive calls of tb2setvalue/tb2removevalue/tb2setmin/tb2setmax
        }
        if (negproblem) {
            negproblem->enforceUb();
            if (!negproblem->isactivatePropagate())
                return; // idem
        }
        assigns();
        if (connected()) {
            try {
                if (problem && problem->isactivatePropagate()) {
                    protect();
                    problem->propagate();
                    unprotect();
                    if (strongDuality && connected() && canbeDeconnected()) {
                        if (problem->getLb() < lb) {
                            assert(!wcsp->vac || wcsp->vac->getThreshold() == UNIT_COST); // FIXME: we should check also variable.myThreshold and constraint.myThreshold are all zero
                            THROWCONTRADICTION;
                        } else {
                            deconnect();
                        }
                    }
                }
                if (connected()) {
                    if (negproblem && negproblem->isactivatePropagate()) {
                        protect();
                        negproblem->propagate();
                        unprotect();
                        if (connected()) {
                            assert(negproblem && negproblem->isactivatePropagate());
                            assert(negproblem->propagated());
                            assert(!problem || problem->isactivatePropagate());
                            if (problem && !problem->propagated())
                                propagate(); // continue recursively if one subproblem is not fully propagated
                        }
                    }
                }
            } catch (const Contradiction&) {
                if (problem)
                    problem->whenContradiction();
                if (negproblem)
                    negproblem->whenContradiction();
                unprotect();
                THROWCONTRADICTION;
            }
        }
        assert(!problem || problem->getLb() < ub);
        assert(!negproblem || negproblem->getLb() < -lb + negCost + UNIT_COST);
    }

    //    void setDACScopeIndex(); //TODO: reorder variables inside problem and negproblem when setDACOrder is called

    bool checkEACGreedySolution(int index = -1, Value supportValue = 0) FINAL
    {
        for (int i = 0; i < arity_; i++) {
            EnumeratedVariable* var = (EnumeratedVariable*)getVar(i);
            evalTuple[i] = var->toIndex((i == index) ? supportValue : var->getSupport());
        }
        return eval(evalTuple) == MIN_COST;
    }

    bool reviseEACGreedySolution(int index = -1, Value supportValue = 0) FINAL
    {
        bool result = checkEACGreedySolution(index, supportValue);
        if (!result) {
            if (index >= 0) {
                getVar(index)->unsetFullEAC();
            } else {
                for (int i = 0; i < arity_; i++) {
                    getVar(i)->unsetFullEAC();
                }
            }
        }
        return result;
    }

    void print(ostream& os) override
    {
        os << this << " WeightedWCSPConstraint(";
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
        os << " isfinite: " << isfinite;
        os << " strongDuality: " << strongDuality;
        os << " arity: " << arity_;
        os << " unassigned: " << getNonAssigned() << "/" << unassigned_ << endl;
        if (problem)
            os << *problem << endl;
        if (negproblem)
            os << *negproblem << endl;
    }

    void dump(ostream& os, bool original = true) override
    {
        if (!problem && !negproblem)
            return;
        if (original) {
            os << arity_;
            for (int i = 0; i < arity_; i++)
                os << " " << scope[i]->wcspIndex;
        } else {
            os << getNonAssigned();
            for (int i = 0; i < arity_; i++)
                if (scope[i]->unassigned())
                    os << " " << scope[i]->getCurrentVarId();
        }
        if (problem) {
            os << " -1 wcsp " << lb << " " << ub << " " << problem->isfinite() << " " << strongDuality << endl;
            problem->dump(os, original);
        } else if (negproblem) {
            os << " -1 wcsp " << -ub + negCost << " " << -lb + negCost << " " << negproblem->isfinite() << " " << strongDuality << endl;
            negproblem->dump(os, original);
        }
    }

    void dump_CFN(ostream& os, bool original = true) override
    {
        if (!problem && !negproblem)
            return;
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
        }
        os << "],\"type\": \"cfnconstraint\",\"params\":{\n\"cfn\":\n";
        if (problem) {
            problem->dump_CFN(os, original);
            os << ",\n\"lb\":" << wcsp->DCost2Decimal(problem->Cost2ADCost(lb)) << ",\"ub\":" << wcsp->DCost2Decimal(problem->Cost2ADCost(ub)) << ",\"duplicatehard\":" << problem->isfinite() << ",\"strongduality\":" << strongDuality << "}},\n";
        } else if (negproblem) {
            negproblem->dump_CFN(os, original);
            os << ",\n\"lb\":" << wcsp->DCost2Decimal(negproblem->Cost2ADCost(-ub + negCost)) << ",\"ub\":" << wcsp->DCost2Decimal(negproblem->Cost2ADCost(-lb + negCost)) << ",\"duplicatehard\":" << negproblem->isfinite() << ",\"strongduality\":" << strongDuality << "}},\n";
        }
    }

    friend void tb2setvalue(int wcspId, int varIndex, Value value, void* solver);
    friend void tb2removevalue(int wcspId, int varIndex, Value value, void* solver);
    friend void tb2setmin(int wcspId, int varIndex, Value value, void* solver);
    friend void tb2setmax(int wcspId, int varIndex, Value value, void* solver);

public:
    void clearPtrReferences()
    {
        for (auto it = WeightedCSPConstraints.begin(); it != WeightedCSPConstraints.end();) {
            if (it->second == this) {
                it = WeightedCSPConstraints.erase(it);
            } else {
                ++it;
            }
        }
    }
};

#endif /*TB2GLOBALWCSP_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
