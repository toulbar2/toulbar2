#ifndef TB2CLQ_HPP_
#define TB2CLQ_HPP_

#include "tb2abstractconstr.hpp"
#include "tb2ternaryconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2wcsp.hpp"

/* Enforce that among the variables in scope, at most rhs_in get any
   of the values appearing in clq_in. This is isomorphic to having
   just two values, 0 (not in the clique) and 1 (in the clique), so
   we use these terms.
*/
class CliqueConstraint : public AbstractNaryConstraint {
public:
    CliqueConstraint(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in,
        vector<vector<int>> clq_in, int rhs_in);
    CliqueConstraint(WCSP* wcsp, EnumeratedVariable** scope_in, int arity_in);
    ~CliqueConstraint();

    void read(istream& file);

    bool extension() const FINAL { return false; }

    void assign(int idx) override;
    void remove(int idx) override;
    void increase(int idx) override;
    void decrease(int idx) override;
    void projectFromZero(int idx) override;

    void propagate() override;

    Cost eval(const Tuple& s) override
    {
        bool iszerotuple = true;
        for (int i = 0; i < arity_; i++) {
            if (inclq[i][s[i]]) {
                iszerotuple = false;
                break;
            }
        }
        if (iszerotuple)
            return all0;
        else
            return MIN_COST;
    }

    vector<Long> conflictWeights; // used by weighted degree heuristics
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
    void resetConflictWeight() override
    {
        conflictWeights.assign(conflictWeights.size(), 0);
        Constraint::resetConflictWeight();
    }
    double computeTightness() override;
    Cost getMaxFiniteCost() override
    {
        if (!CUT(all0, wcsp->getUb())) {
            return all0;
        } else {
            return MIN_COST;
        }
    }
    void setInfiniteCost(Cost ub) override
    {
        Cost mult_ub = ((wcsp->getUb() < (MAX_COST / MEDIUM_COST)) ? (max(LARGE_COST, wcsp->getUb() * MEDIUM_COST)) : wcsp->getUb());
        if (CUT(all0, ub))
            all0 = mult_ub;
    }
    void dump(ostream&, bool) override { cerr << "Warning! Clique constraint cannot be dump." << endl; } // TODO
    void dump_CFN(ostream&, bool) override { cerr << "Warning! Clique constraint cannot be dump." << endl; } // TODO

private:
    // ----------------------------------------------------------------------
    // definition

    // two views of values in the clique: vector<bool> per variable
    // (inclq[var][var->toIndex(val)] == true iff (var,val) is in clique) and array
    // of values in clique per variable
    vector<vector<bool>> inclq;
    vector<vector<int>> clqvals;
    vector<vector<int>> nonclqvals;

    // We require that we use at most rhs values among those in the
    // clique (and of course, even if rhs > 1, no more than one value
    // from each variable)
    int rhs;

    //----------------------------------------------------------------------
    // operation

    // construct current_scope (vector of uninstantiated variables)
    // and current_scope_idx (map from indices of current_scope to
    // indices of scope)
    void get_current_scope(std::vector<EnumeratedVariable*>& s,
        std::vector<int>& si);

    // compute zero (resp., one) cost of a var (given by index into
    // original scope): min Cost among values that do not (resp., do)
    // appear in the clique. var is index into scope
    Cost get_zero_cost(int var);
    Cost get_one_cost(int var);

    void extend_zero_cost(int var, Cost c);
    void project_zero_cost(int var, Cost c);
    void project_one_cost(int var, Cost c);

    // compute binary zero cost of a pair of vars, i.e., the minimum
    // c_{ij}(k,l) where k and l are not in the clique. call this the
    // fat tuple 00 of that constraint. idx and jdx are indices into
    // scope.
    Cost get_binary_zero_cost(int idx, int jdx);

    // extend from the binary constraint (idx,jdx), from the fat tuple
    // 0,0, cost c
    void extend_binary_cost(int idx, int jdx, Cost c);
    BinaryConstraint* project_binary_cost(int idx, int jdx, Cost c);

    // when we have a small number of unassigned vars, project to a
    // binary/unary/nullary constraint.
    void handle_low_arity();

    // gather unary/binary costs to increase the lower bound. We
    // gather the minimum cost of all the 1s separately, as it is
    // useful to do it after gathering binary costs
    void gather_unary_0s();
    void gather_unary_1s();
    void gather_binary();

    void initialize_binary();

    // the remove/projectFromZero etc all come to this.
    void propagate_incremental();

    // amount we have already projected to c_zero
    StoreCost lb;
    // cost of assigning everything to 0
    StoreCost all0;

    // number of variables already assigned to 1
    StoreInt num1;
    // number of variables that remain unassigned
    StoreInt carity;

    // for each variable, a 0-cost 0 value and a zero-cost 1 value
    vector<int> supports0;
    vector<int> supports1;

    // buffer: current_scope (scope excluding instantiated variables),
    // plus indices to original scope
    vector<EnumeratedVariable*> current_scope;
    vector<int> current_scope_idx;

    vector<EnumeratedVariable*> current_scope_asgn;
    vector<int> current_scope_asgn_idx;

    vector<Cost> zerocosts;
    vector<Cost> binary_extra;

    // binary constraints in scope
    vector<vector<BinaryConstraint*>> bc;

    std::ostream& printstate(std::ostream& os);

    int run{ 0 };
    int id{ 0 };
    static int nextid;

public:
    struct state {
        CliqueConstraint* clq;
        std::ostream& print(std::ostream& os) { return clq->printstate(os); }
    };

    void print(ostream& os) override { printstate(os); }
};

inline std::ostream& operator<<(std::ostream& os, CliqueConstraint::state s)
{
    return s.print(os);
}

#endif /* TB2CLQ_HPP_ */

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
