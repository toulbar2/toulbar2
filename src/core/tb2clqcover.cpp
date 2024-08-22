/*----------------------------------------------------------------------
 *
 * Clique constraint
 *
 * Typically redundant constraint. Given a set of a variables, a set
 * of values V(i) for each variable (i) and a constant rhs, ensure
 * that at most rhs variables get a value from the set v(i)
 */

#include "tb2clqcover.hpp"
#include "search/tb2clusters.hpp"

int CliqueConstraint::nextid{ 0 };

CliqueConstraint::CliqueConstraint(WCSP* wcsp, EnumeratedVariable** scope_in,
    int arity_in, vector<vector<int>> clq_in,
    int rhs_in)
    : AbstractNaryConstraint(wcsp, scope_in, arity_in)
    , clqvals(clq_in)
    , rhs(rhs_in)
    , lb(MIN_COST)
    , all0(MIN_COST)
    , num1(0)
    , carity(arity_in)
{
    inclq.resize(arity_);
    for (int i = 0; i != arity_; ++i) {
        conflictWeights.push_back(0);
        inclq[i].resize(scope[i]->getDomainInitSize());
        for (int v : clqvals[i]) {
            assert(static_cast<size_t>(v) < inclq[i].size());
            inclq[i][v] = true;
        }
        auto* x = scope[i];
        if (x->assigned()) {
            if (inclq[i][x->toIndex(x->getValue())]) {
                num1 += 1;
                if (num1 == rhs)
                    deconnect();
                return;
            }
            deconnect(i);
            carity -= 1;
        }
        for (auto val : *x) {
            int vindex = x->toIndex(val);
            if (!inclq[i][vindex])
                nonclqvals[i].push_back(vindex);
        }
    }
    id = nextid++;
    if (rhs != 1) {
        cout << "Unsupported: rhs == " << rhs << "\n";
        throw WrongFileFormat();
    }
}

CliqueConstraint::CliqueConstraint(WCSP* wcsp, EnumeratedVariable** scope_in,
    int arity_in)
    : AbstractNaryConstraint(wcsp, scope_in, arity_in)
    , lb(MIN_COST)
    , all0(MIN_COST)
    , num1(0)
    , carity(arity_in)
{
    id = nextid++;
}

CliqueConstraint::~CliqueConstraint()
{
}

std::ostream& CliqueConstraint::printstate(std::ostream& os)
{
    os << endl
       << this << " clique cut: ";
    os << "all0 = " << all0 << " carity = " << carity << " run = " << run
       << " id = " << id << " connected = " << connected()
       << " depth = " << Store::getDepth() << "\n";
    for (int i = 0; i != arity_; ++i) {
        auto* x = scope[i];
        if (connected(i))
            os << " * ";
        else
            os << "   ";
        os << x->getName();
        x->print(os);
        if (x->assigned() && connected(i))
            os << "*****";
        os << "\n";
    }
    return os;
}

void CliqueConstraint::propagate()
{
    if (ToulBar2::dumpWCSP % 2) // do not propagate if problem is dumped before preprocessing
        return;
    if (bc.empty()) {
        assert(Store::getDepth() == 0);
        initialize_binary();
    }
    propagate_incremental();
}

void CliqueConstraint::propagate_incremental()
{
    ++run;
    static const bool debug{ false };

    if (!connected())
        return;

    if (debug) {
        cout << "--------------------------------------------------\n";
        cout << "Propagator id = " << id << " run = " << run << "\n";
    }

    get_current_scope(current_scope, current_scope_idx);
    if (debug)
        cout << "carity = " << carity << "\n"
             << state{ this } << "\n";

    handle_low_arity();
    if (!connected()) // handle_low_arity may disconnect
        return;
    assert(carity > 0);

    wcsp->revise(this);
    gather_unary_0s();
#ifdef PROPAGATE_CLIQUE_WITH_BINARIES
    TreeDecomposition* td = wcsp->getTreeDec();
    if (!td)
        gather_binary(); // Warning! does not work with tree decomposition-based methods
#endif
    gather_unary_1s();

    if (debug) {
        cout << "After propagate, state\n"
             << state{ this } << "\n";
    }
}

void CliqueConstraint::get_current_scope(std::vector<EnumeratedVariable*>& s,
    std::vector<int>& s_idx)
{
    // recover current scope
    s.clear();
    s_idx.clear();

    num1 = 0;
    for (int i = 0; i < arity_; i++) {
        EnumeratedVariable* x = scope[i];
        if (x->assigned()) {
            if (inclq[i][x->toIndex(x->getValue())])
                num1 += 1;
            deconnect(i);
            continue;
        }
        s.push_back(scope[i]);
        s_idx.push_back(i);
    }
    carity = s.size();
}

void CliqueConstraint::gather_unary_0s()
{
    static const bool debug{ false };

    if (debug)
        cout << "------------------------------\ngather_unary_0s state = \n"
             << state{ this } << "\n";

    zerocosts.clear();
    zerocosts.resize(carity);
    Cost maxc{ MIN_COST }, sumc{ MIN_COST }, secondmax{ MIN_COST };
    for (int i = 0, e = current_scope.size(); i != e; ++i) {
        auto i0 = get_zero_cost(current_scope_idx[i]);
        zerocosts[i] = i0;
        if (maxc < i0) {
            secondmax = maxc;
            maxc = i0;
        } else if (secondmax < i0) {
            secondmax = i0;
        }
        maxc = std::max(maxc, i0);
        sumc += i0;
        if (debug) {
            auto* x = current_scope[i];
            cout << "var " << i << " ";
            x->print(cout);
            cout << " i0 = " << i0 << "\n";
        }
    }

    for (int i = 0, e = current_scope.size(); i != e; ++i) {
        extend_zero_cost(current_scope_idx[i],
            std::min(zerocosts[i], secondmax));
    }
    Cost l0 = sumc - maxc;

    if (debug)
        cout << "sumc = " << sumc << " maxc = " << maxc << " l0 = " << l0
             << " lb = " << wcsp->getLb() << " depth = " << Store::getDepth()
             << " all0 = " << all0 << "\nzerocosts = " << zerocosts
             << " arity = " << carity << "\n";

    Constraint::projectLB(l0);
    all0 += secondmax;
    Cost fixedsumc = sumc - maxc + secondmax;
    for (int i = 0, e = current_scope.size(); i != e; ++i)
        project_one_cost(current_scope_idx[i],
            fixedsumc - std::min(zerocosts[i], secondmax) - l0);
}

void CliqueConstraint::gather_unary_1s()
{
    Cost min1cost{ wcsp->getUb() };
    for (int i = 0, e = current_scope.size(); i != e; ++i)
        min1cost = std::min(get_one_cost(current_scope_idx[i]), min1cost);

    Cost extra{ std::min(min1cost, (Cost)all0) };
    if (extra > MIN_COST) {
        TreeDecomposition* td = wcsp->getTreeDec();
        all0 -= extra;
        assert(all0 >= 0);
        for (int i = 0, e = current_scope.size(); i != e; ++i) {
            auto* x = current_scope[i];
            for (auto vindex : clqvals[current_scope_idx[i]]) {
                Value v = x->toValue(vindex);
                if (x->canbe(v)) {
                    if (td)
                        td->addDelta(cluster, x, v, -extra);
                    x->extend(v, extra);
                }
            }
        }
        Constraint::projectLB(extra);
    }
}

void CliqueConstraint::initialize_binary()
{
    get_current_scope(current_scope, current_scope_idx);
    bc.resize(arity());
    for (auto& bcx : bc)
        bcx.resize(arity());
    using std::begin;
    using std::end;

    std::map<EnumeratedVariable*, int> rmap;
    for (int i = 0, e = current_scope.size(); i != e; ++i)
        rmap[current_scope[i]] = i;

    for (int i = 0, e = current_scope.size(); i != e; ++i) {
        auto& xvar = *current_scope[i];
        for (auto&& cle : *xvar.getConstrs()) {
            auto* cons = cle.constr;
            if (!cons->isBinary())
                continue;
            BinaryConstraint *bincons = (BinaryConstraint *)cons;
            auto& bcons = *bincons;
            auto& yvar = [&]() -> EnumeratedVariable& {
                if (&xvar == bcons.getVar(0))
                    return static_cast<EnumeratedVariable&>(*bcons.getVar(1));
                else
                    return static_cast<EnumeratedVariable&>(*bcons.getVar(0));
            }();
            if (!rmap.count(&yvar))
                continue;
            int j = rmap[&yvar];
            bc[current_scope_idx[i]][current_scope_idx[j]] = bincons;
        }
    }
}

void CliqueConstraint::gather_binary()
{
    if (bc.empty()) {
        assert(Store::getDepth() == 0);
        initialize_binary();
    }

    Cost sum{ MIN_COST };
    Cost maxe{ MIN_COST };
    vector<Cost>& extra = binary_extra;
    extra.clear();
    extra.resize(current_scope.size());
    for (int i = 0, e = current_scope.size(); i != e; ++i) {
        if (!connected(current_scope_idx[i]))
            continue;
        for (int j = i + 1; j != e; ++j) {
            if (!connected(current_scope_idx[j]))
                continue;
            if (!bc[current_scope_idx[i]][current_scope_idx[j]])
                continue;
            auto c00 = get_binary_zero_cost(current_scope_idx[i],
                current_scope_idx[j]);
            extend_binary_cost(current_scope_idx[i], current_scope_idx[j], c00);
            sum += c00;
            extra[i] += c00;
            extra[j] += c00;
            maxe = std::max({ maxe, extra[i], extra[j] });
        }
    }

    Cost l0 = sum - maxe;
    Constraint::projectLB(l0);
    all0 += sum - l0;
    for (int i = 0, e = current_scope.size(); i != e; ++i)
        project_one_cost(current_scope_idx[i], sum - extra[i] - l0);
}

void CliqueConstraint::assign(int idx)
{
    ++run;
    static const bool debug{ false };

    auto* x = scope[idx];

    if (debug)
        cout << "In assign " << idx << "=" << x->getValue() << " run = " << run
             << " id = " << id << "\n";

    if (!connected(idx))
        return;
    deconnect(idx);
    carity -= 1;

    if (inclq[idx][x->toIndex(x->getValue())])
        num1 += 1;

    handle_low_arity();

    if (num1 == rhs) {
        if (debug)
            cout << "disconnecting, state = \n"
                 << state{ this } << "\n";
        deconnect();
    } else {
        if (ToulBar2::FullEAC)
            reviseEACGreedySolution();
    }

    if (debug)
        cout << "After assign of " << idx << " state = \n"
             << state{ this } << "\n";
}

void CliqueConstraint::handle_low_arity()
{
    static const bool debug{ false };

    if (debug)
        cout << "in handle_low_arity state = " << state{ this } << "\n";

    if (carity > 3)
        return;

    deconnect();
    if (num1 != rhs && all0 > MIN_COST) {
        projectNary();
    }

    // // make sure we have the right arity
    // get_current_scope(current_scope_asgn, current_scope_asgn_idx);
    // if (num1 == rhs || all0 == MIN_COST) {
    //     deconnect();
    //     return;
    // }

    // if (carity == 2) {
    //     if (bc.empty()) {
    //         assert(Store::getDepth() == 0);
    //         initialize_binary();
    //     }
    //     auto* cons = project_binary_cost(current_scope_asgn_idx[0],
    //                                      current_scope_asgn_idx[1], all0);
    //     if (!cons->connected())
    //         cout << "OOPS. Binary constraint disconnected " << __FILE__ << ":"
    //              << __LINE__ << "\n";
    //     cons->propagate();
    //     deconnect();
    //     return;
    // }

    // if (carity == 1) {
    //     project_zero_cost(current_scope_asgn_idx[0], all0);
    //     deconnect();
    //     return;
    // }

    // if (carity == 0) {
    //     Constraint::projectLB(all0);
    //     deconnect();
    //     return;
    // }
}

void CliqueConstraint::remove(int idx)
{
    static const bool debug{ false };

    if (!connected(idx))
        return;

    if (debug)
        cout << "In remove " << idx << " run = " << run << "\n";

    // FIXME: be incremental
    propagate_incremental();
}

void CliqueConstraint::increase(int idx)
{
    remove(idx);
}

void CliqueConstraint::decrease(int idx)
{
    remove(idx);
}

void CliqueConstraint::projectFromZero(int idx)
{
    static const bool debug{ false };

    if (!connected(idx))
        return;

    if (debug)
        cout << "In projectFromZero state = " << state{ this } << "\n";

    propagate_incremental();
}

Cost CliqueConstraint::get_zero_cost(int idx) // TODO remember last support and choose between smallest current domain size and clqvalue set
{
    EnumeratedVariable* x = scope[idx];
    return std::accumulate(x->begin(), x->end(), wcsp->getUb(),
        [&](Cost m, Value val) {
            int vindex = x->toIndex(val);
            if (!inclq[idx][vindex])
                return std::min(m, x->getCost(val));
            else
                return m;
        });
}

Cost CliqueConstraint::get_binary_zero_cost(int idx, int jdx) // TODO: avoid iterate if last support still valid and zero cost
{
    EnumeratedVariable* x = scope[idx];
    EnumeratedVariable* y = scope[jdx];
    auto* cons = bc[idx][jdx];
    assert(cons);
    assert(cons->connected());
    Cost c00{ wcsp->getUb() };
    for (auto ivalindex : nonclqvals[idx]) {
        Value ival = x->toValue(ivalindex);
        if (!x->canbe(ival))
            continue;
        for (auto jvalindex : nonclqvals[jdx]) {
            Value jval = y->toValue(jvalindex);
            if (!y->canbe(jval))
                continue;
            c00 = std::min(c00, cons->getCost(x, y, ival, jval));
        }
    }
    return c00;
}

Cost CliqueConstraint::get_one_cost(int idx) // TODO: remember last support
{
    EnumeratedVariable* x = scope[idx];
    return std::accumulate(x->begin(), x->end(), wcsp->getUb(),
        [&](Cost m, Value val) {
            if (inclq[idx][x->toIndex(val)])
                return std::min(m, x->getCost(val));
            else
                return m;
        });
}

void CliqueConstraint::extend_zero_cost(int var, Cost c)
{
    if (c == MIN_COST)
        return;
    TreeDecomposition* td = wcsp->getTreeDec();
    auto* x = scope[var];
    for (auto q = x->begin(), e = x->end(); q != e; ++q) {
        if (!inclq[var][x->toIndex(*q)]) {
            if (td)
                td->addDelta(cluster, x, *q, -c);
            x->extend(*q, c);
        }
    }
}

void CliqueConstraint::project_zero_cost(int var, Cost c)
{
    if (c == MIN_COST)
        return;
    auto* x = scope[var];
    if (x->assigned()) {
        deconnect(var);
        if (!inclq[var][x->toIndex(x->getValue())])
            Constraint::projectLB(c);
        return;
    }
    TreeDecomposition* td = NULL;
    if (!CUT(c + wcsp->getLb(), wcsp->getUb())) {
        td = wcsp->getTreeDec();
    }
    for (auto vindex : nonclqvals[var]) {
        Value v = x->toValue(vindex);
        if (x->canbe(v)) {
            if (td)
                td->addDelta(cluster, x, v, c);
            x->project(v, c, true);
        }
    }
    x->findSupport();
}

void CliqueConstraint::project_one_cost(int var, Cost c)
{
    if (c == MIN_COST)
        return;
    auto* x = scope[var];
    if (x->assigned()) {
        deconnect(var);
        if (inclq[var][x->toIndex(x->getValue())])
            Constraint::projectLB(c);
        return;
    }
    TreeDecomposition* td = NULL;
    if (!CUT(c + wcsp->getLb(), wcsp->getUb())) {
        td = wcsp->getTreeDec();
    }
    for (auto vindex : clqvals[var]) {
        Value v = x->toValue(vindex);
        if (x->canbe(v)) {
            if (td)
                td->addDelta(cluster, x, v, c);
            x->project(v, c, true);
        }
    }
    x->findSupport();
}

void CliqueConstraint::extend_binary_cost(int idx, int jdx, Cost c)
{
    project_binary_cost(idx, jdx, -c);
}

BinaryConstraint* CliqueConstraint::project_binary_cost(int idx, int jdx, Cost c)
{
    EnumeratedVariable* x = scope[idx];
    EnumeratedVariable* y = scope[jdx];
    auto* cons = bc[idx][jdx];
    assert(cons);
    assert(cons->connected());
    for (auto ivalindex : nonclqvals[idx]) {
        Value ival = x->toValue(ivalindex);
        if (!x->canbe(ival))
            continue;
        for (auto jvalindex : nonclqvals[jdx]) {
            Value jval = y->toValue(jvalindex);
            if (!y->canbe(jval))
                continue;
            cons->addcost(x, y, ival, jval, c); // TODO: update deltas for tree decomposition-based methods
        }
    }
    return cons;
}

double CliqueConstraint::computeTightness()
{
    return 1.0 * all0 / std::pow(2, arity_);
}

void CliqueConstraint::read(istream& is)
{
    inclq.resize(arity_);
    clqvals.resize(arity_);
    nonclqvals.resize(arity_);
    is >> rhs;
    for (int i = 0; i != arity_; ++i) {
        conflictWeights.push_back(0);
        int nv{ 0 };
        is >> nv;
        auto* x = scope[i];
        inclq[i].resize(x->getDomainInitSize());
        for (int j = 0; j != nv; ++j) {
            int val{ 0 };
            is >> val;
            assert(static_cast<size_t>(x->toIndex(val)) < inclq[i].size());
            inclq[i][x->toIndex(val)] = true;
            clqvals[i].push_back(x->toIndex(val));
        }
        if (x->assigned()) {
            if (inclq[i][x->toIndex(x->getValue())]) {
                num1 += 1;
                if (num1 == rhs)
                    deconnect();
                return;
            }
            deconnect(i);
            carity -= 1;
        }
        for (auto val : *x) {
            int vindex = x->toIndex(val);
            if (!inclq[i][vindex])
                nonclqvals[i].push_back(vindex);
        }
    }
    if (rhs != 1) {
        cerr << "Unsupported: rhs == " << rhs << " and should be set to one.\n";
        throw WrongFileFormat();
    }
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
