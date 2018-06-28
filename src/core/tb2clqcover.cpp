/*
 * **************** Clique constraint **************************
 *
 */

#include "tb2clqcover.hpp"
#include "search/tb2clusters.hpp"
#include "utils/Bitset.hpp"

//--------------------------------------------------
// output operator for vector, pair
namespace std {

template <typename T> ostream& operator<<(ostream& os, vector<T> const& v)
{
    os << "v(sz=" << v.size() << ")[";
    bool first = true;
    for (auto&& t : v) {
        if (first)
            first = false;
        else
            os << ",";
        os << t;
    }
    os << "]";
    return os;
}

template <typename U, typename T>
ostream& operator<<(ostream& os, pair<U, T> const& p)
{
    return os << "p{" << p.first << "," << p.second << "}";
}

}

//----------------------------------------------------------------------
// a graph structure with optional storing of additional data for each
// vertex

struct vertex_ref {
    int id;
    vertex_ref()
        : id()
    {
    }
    explicit vertex_ref(int i)
        : id(i)
    {
    }
};

bool operator==(vertex_ref u, vertex_ref v) { return u.id == v.id; }
bool operator!=(vertex_ref u, vertex_ref v) { return u.id != v.id; }
bool operator<(vertex_ref u, vertex_ref v) { return u.id < v.id; }

std::ostream& operator<<(std::ostream& os, vertex_ref u)
{
    os << "v" << u.id;
    return os;
}

template<typename VertexData>
struct vertex : public boost::compressed_pair<vertex_ref, VertexData>
{
    using base_type = boost::compressed_pair<vertex_ref, VertexData>;
    using base_type::base_type;

    int ref() const { return this->first(); }
    VertexData& data() { return this->second(); }
    VertexData const& data() const { return this->second(); }
};

/* An undirected graph with both adjacency list
   representations. No data here. */
struct graph_base_mutable {
public:
    void add_edge(vertex_ref v1, vertex_ref v2);

    std::vector<vertex_ref> const& neighbors(vertex_ref v) const { return E[v.id]; }
    std::vector<vertex_ref> const& vertices() const { return V; }


protected:
    std::vector<std::vector<vertex_ref>> E;
    std::vector<vertex_ref> V;

    bool has_vertex(vertex_ref r)
    {
        return 0 <= r.id && r.id < static_cast<int>(vertices().size());
    }
};

inline void graph_base_mutable::add_edge(vertex_ref v1, vertex_ref v2)
{
    assert(has_vertex(v1));
    assert(has_vertex(v2));
    E[v1.id].push_back(v2);
    E[v2.id].push_back(v1);
}

struct graph_base {
    graph_base() = delete;
    graph_base(graph_base_mutable &&g) : gd(g) {
      constructM();
    }

    bool M(vertex_ref u, vertex_ref v) const {
        if (u.id > v.id)
            std::swap(u,v);
        return Mstorage[u.id][v.id-u.id];
    }

    BitSet const& M(vertex_ref u) const {
        return Mstorage[u.id];
    }

    std::vector<vertex_ref> const &neighbors(vertex_ref v) const {
      return gd.neighbors(v);
    }
    std::vector<vertex_ref> const &vertices() const { return gd.vertices(); }

  private:
    std::vector<std::vector<vertex_ref>> E;
    std::vector<vertex_ref> V;

    bool has_vertex(vertex_ref r)
    {
        return 0 <= r.id && r.id < static_cast<int>(gd.vertices().size());
    }

    using boolref = bool_reference;
    boolref M(vertex_ref u, vertex_ref v) {
        return Mstorage[u.id][v.id];
    }
    void constructM() {
        BitSet bs(0, vertices().size()-1, BitSet::empt);
        Mstorage.resize(gd.vertices().size(), bs);
        for(auto v : gd.vertices()) {
            for(auto u : gd.neighbors(v))
                M(u,v)=true;
        }
    }
    std::vector<BitSet> Mstorage;

    graph_base_mutable gd;
};

/* Graph, stores VertexData for each vertex. Can be void and will then
   take no extra space */
template<typename VertexData> struct graph;

template<typename VertexData> struct graph_mutable : public graph_base_mutable {
    friend class graph<VertexData>;
public:
    using vertex_type = vertex<VertexData>;

    vertex_ref get_vertex(VertexData const& vd);
    VertexData const& get_data(vertex_ref u) const;
protected:
    std::vector<vertex_type> vertex_data;
    std::map<VertexData, vertex_ref> rmap;
};

template <typename VD> vertex_ref graph_mutable<VD>::get_vertex(VD const &vd) {
  if (!rmap.count(vd)) {
    auto newid = int{static_cast<int>(vertices().size())};
    auto newref = vertex_ref{newid};
    vertex_data.push_back(vertex<VD>(newref, vd));
    V.push_back(newref);
    E.push_back({});
    rmap[vd] = newref;
    return newref;
  }
  return rmap[vd];
}

template <typename VD>
VD const &graph_mutable<VD>::get_data(vertex_ref u) const {
  assert(u.id >= 0 && u.id < static_cast<int>(vertex_data.size()));
  return vertex_data[u.id].data();
}

template <typename VertexData> struct graph : public graph_base {
public:
    using vertex_type = vertex<VertexData>;

    graph() = delete;
    graph(graph_mutable<VertexData> &&g)
        : graph_base(std::move(g)), vertex_data(std::move(g.vertex_data)),
          rmap(std::move(g.rmap)) {}

    vertex_ref get_vertex(VertexData const& vd) const;
    VertexData const& get_data(vertex_ref u) const;
protected:
    std::vector<vertex_type> vertex_data;
    std::map<VertexData, vertex_ref> rmap;
};

template <typename VD> vertex_ref graph<VD>::get_vertex(VD const &vd) const {
  if (!rmap.count(vd)) {
    throw std::runtime_error("vertex not in graph");
  }
  return rmap[vd];
}

template <typename VD> VD const& graph<VD>::get_data(vertex_ref u) const
{
    assert(u.id >= 0 && u.id < static_cast<int>(vertex_data.size()));
    return vertex_data[u.id].data();
}

//----------------------------------------------------------------------
// Find a degeneracy ordering of the graph: each vertex in the
// ordering has at most d neighbors ahead of it in the
// ordering. Return the ordering and d

inline
std::pair<std::vector<vertex_ref>, int> degeneracy_ordering(graph_base& g)
{
    std::pair<std::vector<vertex_ref>, int> rv;
    auto& order{rv.first};
    auto& d{rv.second};

    std::vector<std::list<vertex_ref>> buckets;
    std::vector<int> degrees(g.vertices().size());
    std::vector<std::list<vertex_ref>::iterator> iterators(g.vertices().size());
    std::vector<bool> ordered(g.vertices().size());
    for (auto v : g.vertices()) {
        auto vd = g.neighbors(v).size();
        if (vd >= buckets.size())
            buckets.resize(vd+1);
        buckets[vd].push_front(v);
        degrees[v.id] = vd;
        iterators[v.id] = buckets[vd].begin();
        ordered[v.id] = false;
    }

    while(true) {
        size_t i{0};
        for (; i != buckets.size(); ++i)
            if (!buckets[i].empty())
                break;
        if (i == buckets.size())
            break;
        d = std::max(d,static_cast<int>(i));
        auto v = buckets[i].back();
        order.push_back(v);
        buckets[i].pop_back();
        ordered[v.id] = true;
        for (auto u : g.neighbors(v)) {
            if (ordered[u.id])
                continue;
            auto &ud = degrees[u.id];
            buckets[ud].erase(iterators[u.id]);
            --ud;
            buckets[ud].push_front(u);
            iterators[u.id] = buckets[ud].begin();
        }
    }

    return rv;
}

//----------------------------------------------------------------------

struct clique_cover
{
    std::vector<BitSet> partitions;
    std::vector<int> partition_of;
    std::vector<BitSet> candidates;
};

clique_cover compute_clique_cover(graph_base const &G,
                                  std::vector<vertex_ref> &V)
{
  clique_cover cc;

  cc.partition_of.resize(G.vertices().size());

  for(auto v : V) {
      int clq = cc.partitions.size();
      for (int i = 0; static_cast<size_t>(i) != cc.partitions.size(); ++i)
          if (cc.candidates[i].contain(v.id)) {
              clq = i;
              break;
          }
      cc.partition_of[v.id] = clq;
      if (static_cast<size_t>(clq) == cc.partitions.size()) {
          // we use e here because otherwise emplace_back will try to get a
          // reference to BitSet::empt and fail to link
          auto e = BitSet::empt;
          cc.partitions.emplace_back(0, G.vertices().size() - 1, e);
          cc.candidates.push_back(G.M(v));
      } else {
          cc.candidates[clq].intersect_with(G.M(v));
      }
      cc.partitions[clq].add(v.id);
  }

  std::cout << "|V| = " << G.vertices().size() << "\n";
  std::cout << "Found " << cc.partitions.size() << " cliques\n";

  for (unsigned i = 0; i != cc.partitions.size(); ++i) {
      auto& p = cc.partitions[i];
      auto& c = cc.candidates[i];
      std::cout << p << "--" << c << "\n";
  }

  return cc;
}

//----------------------------------------------------------------------

struct wcsp_info
{
    int var;
    int val;
};

bool operator==(wcsp_info u, wcsp_info v) {
  return make_pair(u.var, u.val) == make_pair(v.var, v.val);
}
bool operator!=(wcsp_info u, wcsp_info v) {
  return make_pair(u.var, u.val) != make_pair(v.var, v.val);
}
bool operator<(wcsp_info u, wcsp_info v) {
  return make_pair(u.var, u.val) < make_pair(v.var, v.val);
}

std::ostream& operator<<(ostream& os, wcsp_info i)
{
    return os << "X" << i.var << "=" << i.val;
}

void CliqueCoverPropagator::propagate() {
    return; //TODO: why not using this part???
    graph_mutable<wcsp_info> gm;

    std::vector<Value> dom_buffer;

    static int run{0};
    cout << "start clq propagator run = " << ++run << " with "
         << wcsp.numberOfUnassignedVariables() << " unassigned variables "
         << std::endl;

    // construct graph
    for (int var = 0; static_cast<unsigned>(var) != wcsp.numberOfVariables();
         ++var) {
        if (wcsp.assigned(var))
            continue;
        if (!wcsp.enumerated(var))
            continue;
        auto &xvar = *static_cast<EnumeratedVariable *>(wcsp.getVar(var));
        for (auto val : xvar) {
            auto v1 = gm.get_vertex({var, val});
            for (auto val2 : xvar) {
                if (val2 <= val)
                    continue;
                auto v2 = gm.get_vertex({var, val});
                gm.add_edge(v1, v2);
            }
        }
        for (auto &&cle : *xvar.getConstrs()) {
          auto *cons = cle.constr;
          auto *bincons = dynamic_cast<BinaryConstraint *>(cons);
          if (!bincons)
            continue;
          auto &bc = *bincons;
          auto &yvar = [&]() -> EnumeratedVariable & {
              if (&xvar == bc.getVar(0))
                  return static_cast<EnumeratedVariable&>(*bc.getVar(1));
              else
                  return static_cast<EnumeratedVariable&>(*bc.getVar(0));
          }();
          for (auto xval : xvar)
            for (auto yval : yvar)
              if (wcsp.getLb() + xvar.getCost(xval) + yvar.getCost(yval) +
                      bc.getCost(&xvar, &yvar, xval, yval) >=
                  wcsp.getUb()) {
                gm.add_edge(gm.get_vertex({xvar.wcspIndex, xval}),
                            gm.get_vertex({yvar.wcspIndex, yval}));
              }
        }
    }

    graph<wcsp_info> g(std::move(gm));
    auto V = g.vertices(); // copy
    auto cc = compute_clique_cover(g, V);

    BitSet vars(0, wcsp.numberOfVariables()-1, BitSet::empt);
    for (auto& part : cc.partitions) {
        cout << "using clique " << part << "\n";
        if (part.size() < 3)
            continue;
        vars.clear();
        Cost partcost{0};
        for (auto vid : part) {
            vertex_ref v{vid};
            cout << v << " -- " << g.get_data(v) << "\n";
            auto x = g.get_data(v).var;
            if (vars[x])
                continue;
            vars[x] = true;
            auto &xvar = *static_cast<EnumeratedVariable *>(wcsp.getVar(x));
            auto i = std::min_element(
                begin(xvar), end(xvar), [&](int val1, int val2) {
                  return make_pair(!!part[val1], xvar.getCost(val1)) <
                         make_pair(!!part[val2], xvar.getCost(val2));
                });
            assert(i != end(xvar));
            assert(!part[*i]);
            partcost += xvar.getCost(*i);
            {
                // how much could we project to it
                Cost canproject{0};
                vector<int> contrib;
                for( auto &&cle : *xvar.getConstrs()) {
                    auto *cons = cle.constr;
                    auto *bincons = dynamic_cast<BinaryConstraint *>(cons);
                    if (!bincons)
                        continue;
                    auto &bc = *bincons;
                    auto &yvar = [&]() -> EnumeratedVariable & {
                        if (&xvar == bc.getVar(0))
                            return static_cast<EnumeratedVariable&>(*bc.getVar(1));
                        else
                            return static_cast<EnumeratedVariable&>(*bc.getVar(0));
                    }();
                    Cost minbc{wcsp.getUb()};
                    for (auto yval : yvar)
                      minbc =
                          std::min(minbc, bc.getCost(&xvar, &yvar, *i, yval) +
                                              yvar.getCost(yval));
                    if (minbc) {
                        cout << "\t" << minbc << " from ";
                        yvar.print(cout);
                        cout << "\n";
                        contrib.push_back(yvar.wcspIndex);
                    }
                    canproject += minbc;
                }
                cout << "Can project " << canproject << " to " << *i
                     << " from " << contrib << "\n";
            }
            std::cout << v << " is in var " << x << "\n";
            for (auto val : xvar) {
              if (!part[val])
                cout << "out " << val << " cost " << xvar.getCost(val) << "\n";
              else
                cout << "in " << val << "\n";
            }
            cout << "min " << *i << "\n";
            cout << "var " << x << " contributes cost " << xvar.getCost(*i)
                 << "\n";
        }
        cout << "clique cost " << partcost << "\n";
        cout << "--------------------------------------------------\n";
    }
}

/*----------------------------------------------------------------------
 *
 * Clique constraint
 *
 * Typically redundant constraint. Given a set of a variables, a set
 * of values V(i) for each variable (i) and a constant rhs, ensure
 * that at most rhs variables get a value from the set v(i)
 */

int CliqueConstraint::nextid{0};

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
        auto *x = scope[i];
        if (x->assigned()) {
            if (inclq[i][x->getValue()]) {
                num1 += 1;
                if (num1 == rhs)
                    deconnect();
                return;
            }
            deconnect(i);
            carity -= 1;
        }
        for (auto val : *scope[i]) {
            if (!inclq[i][val])
                nonclqvals[i].push_back(val);
        }
    }
    id = nextid++;
    if (rhs != 1) {
        cout << "Unsupported: rhs == " << rhs << "\n";
        exit(1);
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
    os << "all0 = " << all0 << " carity = " << carity << " run = " << run
       << " id = " << id << " connected = " << connected()
       << " depth = " << Store::getDepth() << "\n";
    for (int i = 0; i != arity_; ++i) {
        auto *x = scope[i];
        if (connected(i))
            os << " * ";
        else
            os << "   ";
        os << "var " << i << " ";
        x->print(os);
        if (x->assigned() && connected(i))
            os << "*****";
        os << "\n";
    }
    return  os;
}

void CliqueConstraint::propagate()
{
    if (bc.empty()) {
        assert(Store::getDepth() == 0);
        initialize_binary();
    }
    propagate_incremental();
}

void CliqueConstraint::propagate_incremental()
{
    ++run;
    static const bool debug{false};

    if (!connected())
        return;

    if (debug) {
        cout << "--------------------------------------------------\n";
        cout << "Propagator id = " << id << " run = " << run << "\n";
    }

    get_current_scope(current_scope, current_scope_idx);
    if (debug)
        cout << "carity = " << carity << "\n"
             << state{this} << "\n";

    handle_low_arity();
    if (!connected()) // handle_low_arity may disconnect
        return;
    assert(carity > 0);

    wcsp->revise(this);
    gather_unary_0s();
    TreeDecomposition* td = wcsp->getTreeDec();
    if (!td) gather_binary(); // Warning! does not work with tree decomposition-based methods
    gather_unary_1s();

    if (debug) {
        cout << "After propagate, state\n" << state{this} << "\n";
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
        EnumeratedVariable *x = scope[i];
        if (x->assigned()) {
            if (inclq[i][x->getValue()])
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
    static const bool debug{false};

    if (debug)
        cout << "------------------------------\ngather_unary_0s state = \n"
             << state{this} << "\n";

    zerocosts.clear();
    zerocosts.resize(carity);
    Cost maxc{MIN_COST}, sumc{MIN_COST}, secondmax{MIN_COST};
    for (int i = 0, e = current_scope.size(); i != e; ++i) {
        auto i0{get_zero_cost(current_scope_idx[i])};
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
            auto *x = current_scope[i];
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
    Cost min1cost{wcsp->getUb()};
    for (int i = 0, e = current_scope.size(); i != e; ++i)
        min1cost = std::min(get_one_cost(current_scope_idx[i]), min1cost);

    Cost extra{std::min(min1cost, (Cost)all0)};
    if (extra > MIN_COST) {
        TreeDecomposition* td = wcsp->getTreeDec();
        all0 -= extra;
        assert(all0 >= 0);
        for (int i = 0, e = current_scope.size(); i != e; ++i) {
            auto *x = current_scope[i];
            for (auto v : clqvals[current_scope_idx[i]])
                if (x->canbe(v)) {
                    if(td) td->addDelta(cluster,x,v,-extra);                  
                    x->extend(v, extra);
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
            auto* bincons = dynamic_cast<BinaryConstraint*>(cons);
            if (!bincons)
                continue;
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

    Cost sum{MIN_COST};
    Cost maxe{MIN_COST};
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
            maxe = std::max({maxe, extra[i], extra[j]});
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
    static const bool debug{false};

    auto *x = scope[idx];

    if (debug)
        cout << "In assign " << idx << "=" << x->getValue() << " run = " << run
             << " id = " << id << "\n";

    if (!connected(idx))
        return;
    deconnect(idx);
    carity -= 1;

    if (inclq[idx][x->getValue()])
        num1 += 1;

    handle_low_arity();

    if (num1 == rhs) {
        if (debug)
            cout << "disconnecting, state = \n" << state{this} << "\n";
        deconnect();
    }

    if (debug)
        cout << "After assign of " << idx << " state = \n"
             << state{this} << "\n";
}

void CliqueConstraint::handle_low_arity()
{
    static const bool debug{false};

    if (debug)
        cout << "in handle_low_arity state = " << state{this} << "\n";

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
    static const bool debug{false};

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
    static const bool debug{false};

    if (!connected(idx))
        return;

    if (debug)
        cout << "In projectFromZero state = " << state{this} << "\n";

    propagate_incremental();
}

Cost CliqueConstraint::get_zero_cost(int idx)
{
    EnumeratedVariable *x = scope[idx];
    return std::accumulate(x->begin(), x->end(), wcsp->getUb(),
                           [&](Cost m, int val) {
                               if (!inclq[idx][val])
                                   return std::min(m, x->getCost(val));
                               else
                                   return m;
                           });
}

Cost CliqueConstraint::get_binary_zero_cost(int idx, int jdx)
{
    EnumeratedVariable *x = scope[idx];
    EnumeratedVariable *y = scope[jdx];
    auto *cons = bc[idx][jdx];
    assert(cons);
    assert(cons->connected());
    Cost c00{wcsp->getUb()};
    for (auto ival : nonclqvals[idx]) {
        if (!x->canbe(ival))
            continue;
        for (auto jval : nonclqvals[jdx]) {
            if (!y->canbe(jval))
                continue;
            c00 = std::min(c00, cons->getCost(x,y,ival,jval));
        }
    }
    return c00;
}

Cost CliqueConstraint::get_one_cost(int idx)
{
    EnumeratedVariable *x = scope[idx];
    return std::accumulate(x->begin(), x->end(), wcsp->getUb(),
                           [&](Cost m, int val) {
                               if (inclq[idx][val])
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
    auto *x = scope[var];
    for(auto q = x->begin(), e = x->end(); q != e; ++q) {
        if (!inclq[var][*q]) {
            if(td) td->addDelta(cluster,x,*q,-c);
            x->extend(*q, c);
        }
    }
}

void CliqueConstraint::project_zero_cost(int var, Cost c)
{
    if (c == MIN_COST)
        return;
    auto *x = scope[var];
    if (x->assigned()) {
        deconnect(var);
        if (!inclq[var][x->getValue()])
            Constraint::projectLB(c);
        return;
    }
    TreeDecomposition* td = NULL;
    if (!CUT(c + wcsp->getLb(), wcsp->getUb())) {
        td = wcsp->getTreeDec();
    }
    for (auto v : nonclqvals[var])
        if (x->canbe(v)) {
            if(td) td->addDelta(cluster,x,v,c);
            x->project(v, c, true);
        }
    x->findSupport();
}

void CliqueConstraint::project_one_cost(int var, Cost c)
{
    if (c == MIN_COST)
        return;
    auto *x = scope[var];
    if (x->assigned()) {
        deconnect(var);
        if (inclq[var][x->getValue()])
            Constraint::projectLB(c);
        return;
    }
    TreeDecomposition* td = NULL;
    if (!CUT(c + wcsp->getLb(), wcsp->getUb())) {
        td = wcsp->getTreeDec();
    }
    for (auto v : clqvals[var])
        if (x->canbe(v)) {
            if(td) td->addDelta(cluster,x,v,c);
            x->project(v, c, true);
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
    for (auto ival : nonclqvals[idx]) {
        if (!x->canbe(ival))
            continue;
        for (auto jval : nonclqvals[jdx]) {
            if (!y->canbe(jval))
                continue;
            cons->addcost(x, y, ival, jval, c); //TODO: update deltas for tree decomposition-based methods
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
        int nv{0};
        is >> nv;
        inclq[i].resize(scope[i]->getDomainInitSize());
        for (int j = 0; j != nv; ++j) {
          int val{0};
          is >> val;
          assert(static_cast<size_t>(val) < inclq[i].size());
          inclq[i][val] = true;
          clqvals[i].push_back(val);
        }
        auto* x = scope[i];
        if (x->assigned()) {
            if (inclq[i][x->getValue()]) {
                num1 += 1;
                if (num1 == rhs)
                    deconnect();
                return;
            }
            deconnect(i);
            carity -= 1;
        }
        for (auto val : *x) {
            if (!inclq[i][val])
                nonclqvals[i].push_back(val);
        }
    }
    if (rhs != 1) {
        cerr << "Unsupported: rhs == " << rhs << " and should be set to one.\n";
        exit(1);
    }
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

