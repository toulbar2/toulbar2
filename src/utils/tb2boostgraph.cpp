/*
 * ****** Graph Algorithms using Boost Graph Library ********
 *
 * Variable elimination ordering heuristics developed by Cyril Terrioux <cyril.terrioux@lsis.org>
 */

#include "core/tb2wcsp.hpp"
#include "core/tb2binconstr.hpp"
#include "core/tb2knapsack.hpp"
#include "core/tb2naryconstr.hpp"
#include "core/tb2vacutils.hpp"

#ifdef BOOST
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/graph/minimum_degree_ordering.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>

using namespace boost;

namespace boost {
struct edge_component_t {
    enum { num = 555 };
    typedef edge_property_tag kind;
} edge_component;
} // namespace boost

typedef adjacency_list<setS, vecS, undirectedS> Graph;
typedef adjacency_list<setS, vecS, directedS> DirectedGraph;
typedef adjacency_list<setS, vecS, undirectedS, no_property,
    property<edge_weight_t, int, property<edge_component_t, std::size_t>>>
    IntWeightedGraph;
typedef adjacency_list<setS, vecS, undirectedS, no_property,
    property<edge_weight_t, double, property<edge_component_t, std::size_t>>>
    DoubleWeightedGraph;
typedef adjacency_list<setS, vecS, undirectedS, property<vertex_color_t, default_color_type, property<vertex_degree_t, int>>> ColoredGraph;

template <typename T>
static void addConstraint(Constraint* c, T& g)
{
    int a = c->arity();
    for (int i = 0; i < a; i++) {
        for (int j = i + 1; j < a; j++) {
            Variable* vari = c->getVar(i);
            Variable* varj = c->getVar(j);
            add_edge(vari->wcspIndex, varj->wcspIndex, g);
        }
    }
}

static void addConstraint(Constraint* c, DirectedGraph& g)
{
    int a = c->arity();
    for (int i = 0; i < a; i++) {
        for (int j = i + 1; j < a; j++) {
            Variable* vari = c->getVar(i);
            Variable* varj = c->getVar(j);
            add_edge(vari->wcspIndex, varj->wcspIndex, g);
            add_edge(varj->wcspIndex, vari->wcspIndex, g);
        }
    }
}

static void addConstraint(Constraint* c, IntWeightedGraph& g, int weight = 1)
{
    property_map<IntWeightedGraph, edge_weight_t>::type weights = get(edge_weight, g);
    int a = c->arity();
    for (int i = 0; i < a; i++) {
        for (int j = i + 1; j < a; j++) {
            Variable* vari = c->getVar(i);
            Variable* varj = c->getVar(j);
            weights[add_edge(vari->wcspIndex, varj->wcspIndex, g).first] = weight;
        }
    }
}

static void addConstraint(Constraint* c, DoubleWeightedGraph& g, double maxweight = 1000000)
{
    property_map<DoubleWeightedGraph, edge_weight_t>::type weights = get(edge_weight, g);
    int a = c->arity();
    for (int i = 0; i < a; i++) {
        for (int j = i + 1; j < a; j++) {
            Variable* vari = c->getVar(i);
            Variable* varj = c->getVar(j);
            weights[add_edge(vari->wcspIndex, varj->wcspIndex, g).first] = maxweight - c->getTightness();
        }
    }
}

int WCSP::connectedComponents()
{
    Graph G;
    for (unsigned int i = 0; i < vars.size(); i++)
        add_vertex(G);
    for (unsigned int i = 0; i < constrs.size(); i++)
        if (constrs[i]->connected() && !constrs[i]->universal())
            addConstraint(constrs[i], G);
    for (int i = 0; i < elimBinOrder; i++)
        if (elimBinConstrs[i]->connected())
            addConstraint(elimBinConstrs[i], G);
    for (int i = 0; i < elimTernOrder; i++)
        if (elimTernConstrs[i]->connected())
            addConstraint(elimTernConstrs[i], G);
    vector<int> component(num_vertices(G));
    int num = connected_components(G, &component[0]);
    vector<int> cctruesize(num, 0);
    for (size_t i = 0; i < num_vertices(G); ++i) {
        assert(component[i] >= 0 && component[i] < num);
        if (unassigned(i))
            cctruesize[component[i]]++;
    }
    int res = 0;
    char c = '(';
    for (int i = 0; i < num; ++i) {
        if (cctruesize[i] >= 1) {
            res++;
            cout << c << cctruesize[i];
            c = ' ';
        }
    }
    cout << ")";

    return res;
}

int WCSP::biConnectedComponents()
{
    IntWeightedGraph G;
    for (unsigned int i = 0; i < vars.size(); i++)
        add_vertex(G);
    for (unsigned int i = 0; i < constrs.size(); i++)
        if (constrs[i]->connected())
            addConstraint(constrs[i], G);
    for (int i = 0; i < elimBinOrder; i++)
        if (elimBinConstrs[i]->connected())
            addConstraint(elimBinConstrs[i], G);
    for (int i = 0; i < elimTernOrder; i++)
        if (elimTernConstrs[i]->connected())
            addConstraint(elimTernConstrs[i], G);
    property_map<IntWeightedGraph, edge_component_t>::type component = get(edge_component, G);

    int num = biconnected_components(G, component);

    vector<int> art_points;
    articulation_points(G, back_inserter(art_points));
    cout << "Articulation points: " << art_points.size() << endl;
    if (art_points.size() > 0) {
        for (unsigned int i = 0; i < art_points.size(); i++)
            cout << " " << art_points[i];
        cout << endl;
    }
    return num;
}

int WCSP::diameter()
{
    if (vars.size() >= LARGE_NB_VARS)
        return -1;

    IntWeightedGraph G;
    for (unsigned int i = 0; i < vars.size(); i++)
        add_vertex(G);
    for (unsigned int i = 0; i < constrs.size(); i++)
        if (constrs[i]->connected()) {
            if (constrs[i]->arity() > MAX_ARITY / 10) {
                cerr << "Warning! Cost function arity of " << constrs[i]->arity() << " is too large for diameter computation and will be skipped." << endl;
                continue;
            }
            addConstraint(constrs[i], G);
        }
    for (int i = 0; i < elimBinOrder; i++)
        if (elimBinConstrs[i]->connected())
            addConstraint(elimBinConstrs[i], G);
    for (int i = 0; i < elimTernOrder; i++)
        if (elimTernConstrs[i]->connected())
            addConstraint(elimTernConstrs[i], G);

    typedef int* int_ptr;
    int** D;
    D = new int_ptr[num_vertices(G)];
    for (unsigned int i = 0; i < num_vertices(G); ++i)
        D[i] = new int[num_vertices(G)];
    johnson_all_pairs_shortest_paths(G, D);

    if (ToulBar2::verbose >= 2) {
        cout << "     ";
        for (unsigned int i = 0; i < num_vertices(G); ++i) {
            cout << i << " -> ";
            for (unsigned int j = 0; j < num_vertices(G); ++j) {
                cout << " " << D[i][j];
            }
            cout << endl;
        }
    }

    int maxd = 0;
    double meand = 0;
    for (unsigned int i = 0; i < num_vertices(G); ++i) {
        for (unsigned int j = 0; j < num_vertices(G); ++j) {
            if (D[i][j] > maxd)
                maxd = D[i][j];
            meand += D[i][j];
        }
    }
    meand /= num_vertices(G) * num_vertices(G);
    if (ToulBar2::verbose >= 1) {
        cout << "Mean diameter: " << meand << endl;
    }

    for (unsigned int i = 0; i < num_vertices(G); ++i)
        delete[] D[i];
    delete[] D;

    return maxd;
}

inline bool cmp_vars(Variable* v1, Variable* v2) { return (v1->wcspIndex < v2->wcspIndex); }

/// \brief Minimum Degree Ordering algorithm
/// \warning Output order usually worse than WCSP::minimumDegreeOrdering ???
void WCSP::minimumDegreeOrderingBGL(vector<int>& order_inv)
{
    DirectedGraph G;
    for (unsigned int i = 0; i < vars.size(); i++)
        add_vertex(G);
    for (unsigned int i = 0; i < constrs.size(); i++)
        if (constrs[i]->connected())
            addConstraint(constrs[i], G);
    for (int i = 0; i < elimBinOrder; i++)
        if (elimBinConstrs[i]->connected())
            addConstraint(elimBinConstrs[i], G);
    for (int i = 0; i < elimTernOrder; i++)
        if (elimTernConstrs[i]->connected())
            addConstraint(elimTernConstrs[i], G);

    int n = num_vertices(G);
    int delta = 0;
    typedef vector<int> Vector;
    Vector inverse_perm(n, 0);
    Vector perm(n, 0);

    Vector supernode_sizes(n, 1); // init has to be 1

    property_map<DirectedGraph, vertex_index_t>::type id = get(vertex_index, G);

    Vector degree(n, 0);

    minimum_degree_ordering(G,
        make_iterator_property_map(&degree[0], id, degree[0]),
        &inverse_perm[0],
        &perm[0],
        make_iterator_property_map(&supernode_sizes[0], id, supernode_sizes[0]),
        delta,
        id);

    order_inv = inverse_perm;
    if (ToulBar2::verbose >= 1) {
        cout << "Minimum degree ordering:";
        for (size_t i = 0; i < num_vertices(G); ++i) {
            cout << " " << order_inv[i];
        }
        cout << endl;
    }

    // // \bug reordering of vars array is dubious!!! (invalidates further use of variable indexes)
    // for (size_t i=0; i < num_vertices(G); ++i) {
    //    vars[i]->wcspIndex = num_vertices(G) - perm[i] - 1;
    //  }
    //  stable_sort(vars.begin(), vars.end(), cmp_vars);
    //  for (size_t i=0; i < num_vertices(G); ++i) {
    //    assert(vars[i]->wcspIndex == (int) i);
    //  }

    assert(order_inv.size() == numberOfVariables());
}

void WCSP::spanningTreeOrderingBGL(vector<int>& order_inv)
{
    double alltight = 0;
    double maxt = 0;
    for (unsigned int i = 0; i < constrs.size(); i++)
        if (constrs[i]->connected()) {
            double t = constrs[i]->getTightness();
            alltight += t;
            if (t > maxt)
                maxt = t;
        }
    for (int i = 0; i < elimBinOrder; i++)
        if (elimBinConstrs[i]->connected()) {
            double t = elimBinConstrs[i]->getTightness();
            alltight += t;
            if (t > maxt)
                maxt = t;
        }
    for (int i = 0; i < elimTernOrder; i++)
        if (elimTernConstrs[i]->connected()) {
            double t = elimTernConstrs[i]->getTightness();
            alltight += t;
            if (t > maxt)
                maxt = t;
        }

    DoubleWeightedGraph G;
    for (unsigned int i = 0; i < vars.size(); i++)
        add_vertex(G);
    for (unsigned int i = 0; i < constrs.size(); i++)
        if (constrs[i]->connected())
            addConstraint(constrs[i], G, maxt);
    for (int i = 0; i < elimBinOrder; i++)
        if (elimBinConstrs[i]->connected())
            addConstraint(elimBinConstrs[i], G, maxt);
    for (int i = 0; i < elimTernOrder; i++)
        if (elimTernConstrs[i]->connected())
            addConstraint(elimTernConstrs[i], G, maxt);

    int n = num_vertices(G);

    vector<graph_traits<DoubleWeightedGraph>::vertex_descriptor> p(n);
    prim_minimum_spanning_tree(G, &p[0]);

    double tight = 0;
    bool tightok = true;
    vector<int> roots;
    vector<vector<int>> listofsuccessors(n, vector<int>());
    if (ToulBar2::verbose >= 0)
        cout << "Maximum spanning tree ordering"; // << endl;
    for (size_t i = 0; i != p.size(); ++i) {
        if (p[i] != i) {
            BinaryConstraint* bctr = getVar(i)->getConstr(getVar(p[i]));
            if (bctr) {
                //      cout << "parent[" << i << "] = " << p[i] << " (" << bctr->getTightness() << ")" << endl;
                tight += bctr->getTightness();
            } else {
                tightok = false;
            }
            listofsuccessors[p[i]].push_back(i);
        } else {
            roots.push_back(i);
            //      cout << "parent[" << i << "] = no parent" << endl;
        }
    }
    if (ToulBar2::verbose >= 0) {
        if (tightok)
            cout << " (" << 100.0 * tight / alltight << "%)";
        cout << endl;
    }

    vector<bool> marked(n, false);
    for (int i = roots.size() - 1; i >= 0; i--) {
        visit(roots[i], order_inv, marked, listofsuccessors);
    }
    for (int i = n - 1; i >= 0; i--) {
        if (!marked[i]) {
            visit(i, order_inv, marked, listofsuccessors);
        }
    }

    if (ToulBar2::verbose >= 1) {
        cout << "Maximum spanning tree ordering:";
        for (int i = 0; i < n; i++) {
            cout << " " << order_inv[i];
        }
        cout << endl;
    }

    assert(order_inv.size() == numberOfVariables());
}

void WCSP::DAGOrdering(vector<int>& order_inv)
{
    if (ToulBar2::verbose >= 0)
        cout << "DAG ordering"; // << endl;

    int  n = numberOfVariables();
    set<int> roots; // first it assumes all variables as root
    for (int i = 0; i < n; ++i) {
        roots.insert(roots.end(), i);
    }
    for (int i = 0; i < n; i++) {
        int sz = getListSuccessors()->at(i).size();
        for (int j = 0; j < sz; j++) {
            roots.erase(getListSuccessors()->at(i)[j]);
        }
    }
    vector<bool> marked(n, false);
    for (auto it = roots.rbegin(); it != roots.rend(); ++it) {
        visit(*it, order_inv, marked, listofsuccessors);
    }
    for (int i = n - 1; i >= 0; i--) {
        if (!marked[i]) {
            visit(i, order_inv, marked, listofsuccessors);
        }
    }

    if (ToulBar2::verbose >= 1) {
        cout << ":";
        for (int i = 0; i < n; i++) {
            cout << " " << order_inv[i];
        }
    }
    if (ToulBar2::verbose >= 0) cout << endl;

    assert(order_inv.size() == numberOfVariables());
}

void WCSP::reverseCuthillMcKeeOrderingBGL(vector<int>& order_inv)
{
    ColoredGraph G;
    for (unsigned int i = 0; i < vars.size(); i++)
        add_vertex(G);
    for (unsigned int i = 0; i < constrs.size(); i++)
        if (constrs[i]->connected())
            addConstraint(constrs[i], G);
    for (int i = 0; i < elimBinOrder; i++)
        if (elimBinConstrs[i]->connected())
            addConstraint(elimBinConstrs[i], G);
    for (int i = 0; i < elimTernOrder; i++)
        if (elimTernConstrs[i]->connected())
            addConstraint(elimTernConstrs[i], G);

    int n = num_vertices(G);
    vector<int> inverse_perm(n, 0);

    cuthill_mckee_ordering(G, inverse_perm.rbegin(), get(vertex_color, G), make_degree_map(G));
    order_inv = inverse_perm;
    if (ToulBar2::verbose >= 1) {
        cout << "Reverse Cuthill-McKee ordering:";
        for (size_t i = 0; i < num_vertices(G); ++i) {
            cout << " " << order_inv[i];
        }
        cout << endl;
    }

    assert(order_inv.size() == numberOfVariables());
}

/// \brief Maximum Cardinality Search algorithm (Tarjan & Yannakakis)
/// \note code from Cyril Terrioux
void WCSP::maximumCardinalitySearch(vector<int>& order_inv)
{
    Graph G;
    for (unsigned int i = 0; i < vars.size(); i++)
        add_vertex(G);
    for (unsigned int i = 0; i < constrs.size(); i++)
        if (constrs[i]->connected())
            addConstraint(constrs[i], G);
    for (int i = 0; i < elimBinOrder; i++)
        if (elimBinConstrs[i]->connected())
            addConstraint(elimBinConstrs[i], G);
    for (int i = 0; i < elimTernOrder; i++)
        if (elimTernConstrs[i]->connected())
            addConstraint(elimTernConstrs[i], G);

    int n = num_vertices(G);
    vector<int> inverse_perm(n, 0);

    vector<vector<int>> sets(n, vector<int>(n));
    vector<int> size(n);
    vector<int> card(n);
    vector<int> degree(n);

    Graph::adjacency_iterator neighbourIt, neighbourEnd;

    /* initialize sets, card and size */
    for (int v = 0; v < n; v++) {
        size[v] = 0;

        for (int i = 0; i < n; i++)
            sets[v][i] = 0;

        sets[0][v] = 1;
        card[v] = 0;
        degree[v] = boost::degree(v, G);
    }

    card[0] = n;
    int i = n - 1;
    int j = 0;
    int v = 0;

    while (i >= 0) {
        /* choose a vertex */
        int deg = -1;
        for (int x = 0; x < n; x++)
            if ((sets[j][x] == 1) && (degree[x] > deg)) {
                v = x;
                deg = degree[x];
            }
        sets[j][v] = 0;
        card[j]--;

        /* build the order */
        inverse_perm[i] = v;
        size[v] = -1;

        /* update sets and size */
        boost::tie(neighbourIt, neighbourEnd) = adjacent_vertices(v, G);
        for (; neighbourIt != neighbourEnd; ++neighbourIt) {
            if (size[*neighbourIt] >= 0) {
                sets[size[*neighbourIt]][*neighbourIt] = 0;
                card[size[*neighbourIt]]--;

                size[*neighbourIt]++;

                sets[size[*neighbourIt]][*neighbourIt] = 1;
                card[size[*neighbourIt]]++;
            }
        }

        i--;
        j++;
        while ((j >= 0) && (card[j] == 0))
            j--;
    }

    order_inv = inverse_perm;
    if (ToulBar2::verbose >= 1) {
        cout << "MCS ordering:";
        for (size_t i = 0; i < num_vertices(G); ++i) {
            cout << " " << order_inv[i];
        }
        cout << endl;
    }
    assert(order_inv.size() == numberOfVariables());
}

/// \brief Minimum Fill-In Ordering algorithm
/// \note code from Cyril Terrioux
void WCSP::minimumFillInOrdering(vector<int>& order_inv)
{
    Graph G;
    for (unsigned int i = 0; i < vars.size(); i++)
        add_vertex(G);
    for (unsigned int i = 0; i < constrs.size(); i++)
        if (constrs[i]->connected())
            addConstraint(constrs[i], G);
    for (int i = 0; i < elimBinOrder; i++)
        if (elimBinConstrs[i]->connected())
            addConstraint(elimBinConstrs[i], G);
    for (int i = 0; i < elimTernOrder; i++)
        if (elimTernConstrs[i]->connected())
            addConstraint(elimTernConstrs[i], G);

    int n = num_vertices(G);
    vector<int> order(n, -1);
    order_inv = order;

    vector<int> nb_fillin(n, 0);
    vector<int> degree(n, 0);

    Graph::adjacency_iterator neighbourIt, neighbourEnd;

    for (int v = 0; v < n; v++) {
        degree[v] = boost::degree(v, G);
        /* compute initial number of edges to add for each vertex */
        boost::tie(neighbourIt, neighbourEnd) = adjacent_vertices(v, G);
        for (; neighbourIt != neighbourEnd; ++neighbourIt) {
            Graph::adjacency_iterator neighbourIt2 = neighbourIt;
            for (++neighbourIt2; neighbourIt2 != neighbourEnd; ++neighbourIt2) {
                if (!edge(*neighbourIt, *neighbourIt2, G).second)
                    nb_fillin[v]++;
            }
        }
    }
    for (int i = 0; i < n - 1; i++) {
        /* compute number of fill-in edges to add for each unprocessed vertex */
        /* choose vertex with minimum fill-in */
        int v = 0;
        Long minfill = (Long)n * n;
        int deg = -1;
        for (int x = 0; x < n; x++) {
            if ((order[x] == -1) && ((nb_fillin[x] < minfill) || ((nb_fillin[x] == minfill) && (degree[x] > deg)))) {
                v = x;
                minfill = nb_fillin[x];
                deg = degree[x];
            }
        }
        order[v] = i;
        order_inv[i] = v;

        /* remove vertex v from nb_fillin */
        boost::tie(neighbourIt, neighbourEnd) = adjacent_vertices(v, G);
        for (; neighbourIt != neighbourEnd; ++neighbourIt) {
            if (order[*neighbourIt] == -1) {
                Graph::adjacency_iterator neighbourIt2, neighbourEnd2;
                boost::tie(neighbourIt2, neighbourEnd2) = adjacent_vertices(*neighbourIt, G);
                for (; neighbourIt2 != neighbourEnd2; ++neighbourIt2) {
                    if (order[*neighbourIt2] == -1 && !edge(v, *neighbourIt2, G).second)
                        nb_fillin[*neighbourIt]--;
                }
            }
        }
        /* add fill-in edges to G */
        boost::tie(neighbourIt, neighbourEnd) = adjacent_vertices(v, G);
        for (; neighbourIt != neighbourEnd; ++neighbourIt) {
            if (order[*neighbourIt] == -1) {
                Graph::adjacency_iterator neighbourIt2 = neighbourIt;
                for (++neighbourIt2; neighbourIt2 != neighbourEnd; ++neighbourIt2) {
                    if ((order[*neighbourIt2] == -1) && !edge(*neighbourIt, *neighbourIt2, G).second) {
                        add_edge(*neighbourIt, *neighbourIt2, G);
                        degree[*neighbourIt]++;
                        degree[*neighbourIt2]++;

                        unsigned int x = *neighbourIt;
                        unsigned int y = *neighbourIt2;
                        /* update nb_fillin with missing edges between x and neighbors of y */
                        Graph::adjacency_iterator neighbourItX, neighbourEndX;
                        boost::tie(neighbourItX, neighbourEndX) = adjacent_vertices(x, G);
                        for (; neighbourItX != neighbourEndX; ++neighbourItX) {
                            if ((order[*neighbourItX] == -1) && (*neighbourItX != y)) {
                                if (!edge(y, *neighbourItX, G).second)
                                    nb_fillin[x]++;
                                else
                                    nb_fillin[*neighbourItX]--; /* new added edge between x and y has to be removed from nb_fillin  */
                            }
                        }
                        /* update nb_fillin with missing edges between y and neighbors of x */
                        Graph::adjacency_iterator neighbourItY, neighbourEndY;
                        boost::tie(neighbourItY, neighbourEndY) = adjacent_vertices(y, G);
                        for (; neighbourItY != neighbourEndY; ++neighbourItY) {
                            if ((order[*neighbourItY] == -1) && (*neighbourItY != x) && !edge(x, *neighbourItY, G).second)
                                nb_fillin[y]++;
                        }
                    }
                }
            }
        }
    }

    int v = 0;
    while (order[v] != -1)
        v++;
    order[v] = n - 1;
    order_inv[n - 1] = v;
    if (ToulBar2::verbose >= 1) {
        cout << "Min-fill ordering:";
        for (int j = 0; j < n; ++j) {
            cout << " " << getName(order_inv[j]);
        }
        cout << endl;
    }
    assert(order_inv.size() == numberOfVariables());
}

/// \brief Minimum Degree Ordering algorithm
/// \note code from Cyril Terrioux
void WCSP::minimumDegreeOrdering(vector<int>& order_inv)
{
    Graph G;
    for (unsigned int i = 0; i < vars.size(); i++)
        add_vertex(G);
    for (unsigned int i = 0; i < constrs.size(); i++)
        if (constrs[i]->connected())
            addConstraint(constrs[i], G);
    for (int i = 0; i < elimBinOrder; i++)
        if (elimBinConstrs[i]->connected())
            addConstraint(elimBinConstrs[i], G);
    for (int i = 0; i < elimTernOrder; i++)
        if (elimTernConstrs[i]->connected())
            addConstraint(elimTernConstrs[i], G);
    int n = num_vertices(G);
    vector<int> order(n, -1);
    order_inv = order;
    vector<int> degree(n, 0);

    //    vector<int> preorder;
    //    reverseCuthillMcKeeOrderingBGL(preorder);

    Graph::adjacency_iterator neighbourIt, neighbourEnd;

    for (int v = 0; v < n; v++) {
        degree[v] = boost::degree(v, G);
    }
    for (int i = 0; i < n - 1; i++) {
        /* find vertex with minimum degree */
        int v = 0;
        int deg_min = n + 1;
        for (int x = 0; x < n; x++) {
            //           int x = preorder[xx];
            if ((order[x] == -1) && (degree[x] < deg_min)) {
                v = x;
                deg_min = degree[x];
            }
        }
        order[v] = i;
        order_inv[i] = v;

        boost::tie(neighbourIt, neighbourEnd) = adjacent_vertices(v, G);
        for (; neighbourIt != neighbourEnd; ++neighbourIt) {
            if (order[*neighbourIt] == -1) {
                degree[*neighbourIt]--;
                Graph::adjacency_iterator neighbourIt2 = neighbourIt;
                for (++neighbourIt2; neighbourIt2 != neighbourEnd; ++neighbourIt2) {
                    if ((order[*neighbourIt2] == -1) && !edge(*neighbourIt, *neighbourIt2, G).second) {
                        add_edge(*neighbourIt, *neighbourIt2, G);
                        degree[*neighbourIt]++;
                        degree[*neighbourIt2]++;
                    }
                }
            }
        }
    }
    int v = 0;
    while (order[v] != -1)
        v++;
    order[v] = n - 1;
    order_inv[n - 1] = v;
    if (ToulBar2::verbose >= 1) {
        cout << "Minimum degree ordering:";
        for (int j = 0; j < n; ++j) {
            cout << " " << order_inv[j];
        }
        cout << endl;
    }
    assert(order_inv.size() == numberOfVariables());
}

int cmpValueCost3(const void* p1, const void* p2)
{
    Cost c1 = ((ValueCost*)p1)->cost;
    Cost c2 = ((ValueCost*)p2)->cost;
    Value v1 = ((ValueCost*)p1)->value;
    Value v2 = ((ValueCost*)p2)->value;
    if (c1 < c2)
        return -1;
    else if (c1 > c2)
        return 1;
    else if (v1 < v2)
        return -1;
    else if (v1 > v2)
        return 1;
    else
        return 0;
}

template <typename T>
static vector<vector<pair<int, Value>>> FindClique(vector<int> scope, T& g, vector<int> CorrespVertices)
{
    Graph G;
    for (unsigned int i = 0; i < scope.size(); ++i) {
        add_vertex(G);
        add_vertex(G);
    }
    for (unsigned int i = 0; i < scope.size(); i++) {
        for (unsigned int j = i + 1; j < scope.size(); j++) {
            if (edge(CorrespVertices[scope[i]], CorrespVertices[scope[j]], g).second)
                add_edge(2 * i, 2 * j, G);
            if (edge(CorrespVertices[scope[i]] + 1, CorrespVertices[scope[j]], g).second)
                add_edge(2 * i + 1, 2 * j, G);
            if (edge(CorrespVertices[scope[i]], CorrespVertices[scope[j]] + 1, g).second)
                add_edge(2 * i, 2 * j + 1, G);
            if (edge(CorrespVertices[scope[i]] + 1, CorrespVertices[scope[j]] + 1, g).second)
                add_edge(2 * i + 1, 2 * j + 1, G);
        }
    }
    vector<int> Temp;
    vector<vector<int>> Tempclq;
    vector<int> order;
    for (unsigned int i = 0; i < scope.size(); ++i) {
        order.push_back(i);
    }
    // sort order by degree instead of taking the first variable
    sort(order.begin(), order.end(),
        [&](int x, int y) {
            return (max(degree(2 * x, G), degree(2 * x + 1, G)) > max(degree(2 * y, G), degree(2 * y + 1, G)) || (max(degree(2 * x, G), degree(2 * x + 1, G)) == max(degree(2 * y, G), degree(2 * y + 1, G)) && x < y));
        });
    if (degree(2 * order[0], G) > degree(2 * order[0] + 1, G))
        Tempclq.push_back({ 2 * order[0] });
    else
        Tempclq.push_back({ 2 * order[0] + 1 });
    for (int i = 1; i < (int)order.size(); i++) {
        bool ok = false;
        int j = 0;
        int curr = 2 * order[i] + 1;
        if (degree(2 * order[i], G) > degree(2 * order[i] + 1, G))
            curr = 2 * order[i];
        while (!ok && j < (int)Tempclq.size()) {
            int k = 0;
            ok = true;
            while (ok && k < (int)Tempclq[j].size()) {
                if (!edge(curr, Tempclq[j][k], G).second && !edge(Tempclq[j][k], curr, G).second) // undirected graph should always have both terms equal
                    ok = false;
                k++;
            }
            j++;
        }
        if (ok) {
            assert(j > 0);
            Tempclq[j - 1].push_back(curr);
            // update Tempclq to keep a sorted list of cliques of decreasing size
            if (j - 1 > 0 && Tempclq[j - 1].size() > Tempclq[j - 2].size()) {
                int jnew = j - 2;
                while (jnew > 0 && Tempclq[j - 1].size() > Tempclq[jnew - 1].size()) {
                    jnew--;
                }
                assert(jnew < j - 1);
                swap(Tempclq[jnew], Tempclq[j - 1]);
            }
        } else
            Tempclq.push_back({ curr });
    }
    vector<vector<pair<int, Value>>> clq;
    vector<pair<int, Value>> clq2;
    for (unsigned int i = 0; i < Tempclq.size(); ++i) {
        clq2.clear();
        for (unsigned int l = 0; l < Tempclq[i].size(); ++l) {
            clq2.push_back(make_pair(scope[Tempclq[i][l] / 2], Tempclq[i][l] % 2)); // scope[(int)floor(Tempclq[i][l] / 2.0 + 0.1)]
        }
        if (clq2.size() > 1)
            clq.push_back(clq2);
    }
    return clq;
}

void WCSP::addAMOConstraints()
{
    if (ToulBar2::verbose >= 1)
        cout << "Add AMO constraints to knapsack constraints." << endl;
    double startCpuTime = cpuTime();
    double startRealTime = realTime();
    Graph G;
    vector<int> CorrespVertices{ 0 };
    vector<unsigned int> DomainSizes;
    vector<vector<Value>> Domains;
    vector<map<Value, int>> DomainsInv;
    vector<Cost> ucosts;
    for (int i = 0; i < (int)vars.size(); ++i) {
        EnumeratedVariable* var = (EnumeratedVariable*)vars[i];
        int size = getDomainSize(var->wcspIndex);
        ValueCost sorted[size];
        getEnumDomainAndCost(var->wcspIndex, sorted);
        map<Value, int> DomainInv;
        for (int a = 0; a < size; ++a) {
            add_vertex(G);
            ucosts.push_back(sorted[a].cost);
            DomainInv[sorted[a].value] = a;
        }
        DomainSizes.push_back(size);
        Domains.push_back(getEnumDomain(i));
        DomainsInv.push_back(DomainInv);
        CorrespVertices.push_back(CorrespVertices.back() + size);
    }
    CorrespVertices.pop_back();
    if (ToulBar2::verbose >= 1)
        cout << "Construct a conflict graph of size " << num_vertices(G) << endl;
    int elimDegree_ = ToulBar2::elimDegree_;
    ToulBar2::elimDegree_ = -1; // Warning! constructive disjunction is not compatible with variable elimination!!
    for (unsigned int varIndex = 0; varIndex < vars.size(); varIndex++) {
        EnumeratedVariable* var = (EnumeratedVariable*)vars[varIndex];
        if (var->unassigned()) {
            int size = getDomainSize(var->wcspIndex);
            ValueCost sorted[size];
            getEnumDomainAndCost(var->wcspIndex, sorted);
            qsort(sorted, size, sizeof(ValueCost), cmpValueCost3);
            set<pair<int, Value>> disjunctive;
            bool first = true;
            for (int a = 0; a < size; a++)
                if (canbe(var->wcspIndex, sorted[a].value)) {
                    bool pruned = false;
                    int storedepth = Store::getDepth();
                    try {
                        Store::store();
                        assign(var->wcspIndex, sorted[a].value);
                        propagate();
                        set<pair<int, Value>> removals;
                        for (unsigned int i = 0; i < vars.size(); ++i) {
                            EnumeratedVariable* vari = (EnumeratedVariable*)vars[i];
                            if (i != varIndex && vari->getDomainSize() != DomainSizes[i]) {
                                assert(DomainSizes[i] == Domains[i].size());
                                for (unsigned int j = 0; j < DomainSizes[i]; ++j) {
                                    if (vari->cannotbe(Domains[i][j])) {
                                        add_edge(CorrespVertices[varIndex] + DomainsInv[varIndex][sorted[a].value], CorrespVertices[i] + j, G);
                                        if (first || disjunctive.size() > 0)
                                            removals.insert(make_pair(vari->wcspIndex, Domains[i][j]));
                                    }
                                }
                            }
                        }
                        if (first) {
                            first = false;
                            disjunctive = removals;
                        } else if (disjunctive.size() > 0) {
                            set<pair<int, Value>> intersection;
                            set_intersection(disjunctive.begin(), disjunctive.end(), removals.begin(), removals.end(), inserter(intersection, intersection.end()));
                            disjunctive = intersection;
                        }
                    } catch (const Contradiction&) {
                        whenContradiction();
                        pruned = true;
                    }
                    Store::restore(storedepth);
                    if (pruned) {
                        assert(canbe(var->wcspIndex, sorted[a].value));
                        remove(var->wcspIndex, sorted[a].value); // singleton consistency for free!
                        propagate();
                        if (ToulBar2::verbose >= 0) {
                            cout << ".";
                        }
                    }
                }
            if (disjunctive.size() > 0) { // constructive disjunction for free!
                for (auto elt : disjunctive) {
                    if (canbe(elt.first, elt.second)) {
                        remove(elt.first, elt.second);
                        if (ToulBar2::verbose >= 0) {
                            cout << "+";
                        }
                    }
                }
                propagate();
            }
        }
    }
    ToulBar2::elimDegree_ = elimDegree_;
    for (unsigned int i = 0; i < vars.size(); ++i) {
        EnumeratedVariable* vari = (EnumeratedVariable*)vars[i];
        for (unsigned int j = 0; j < DomainSizes[i]; ++j) {
            if (vari->cannotbe(Domains[i][j]) || vari->getDomainSize() != 2 || vari->toValue(0) != 0 || vari->toValue(vari->getDomainInitSize() - 1) != 1) {
                clear_vertex(CorrespVertices[i] + j, G);
            }
        }
    }
    if (ToulBar2::verbose >= 1) {
        if (ToulBar2::verbose >= 4) {
            for (auto iter = edges(G).first; iter != edges(G).second; ++iter) {
                cout << "Edge: " << *iter << endl;
            }
        }
        cout << "Graph done with " << num_edges(G) << " edges." << endl;
    }
    int count = 0;
    vector<vector<pair<int, Value>>> clq;
    int MaxAMO = 0;
    if (ToulBar2::addAMOConstraints >= 0) {
        bool b = false;
        int total = 0;
        vector<int> scope2;
        unsigned int nbconstrs = constrs.size();
        for (unsigned int i = 0; i < nbconstrs; ++i) {
            if (!constrs[i]->isKnapsack())
                continue;
            else {
                KnapsackConstraint *k = (KnapsackConstraint *)constrs[i];
                if (constrs[i]->arity() > 3 && constrs[i]->connected()) {
                    scope2.clear();
                    clq.clear();
                    for (int j = 0; j < constrs[i]->arity(); j++) {
                        if (constrs[i]->getVar(j)->unassigned() && constrs[i]->getVar(j)->getDomainSize() == 2) { // TODO: find cliques for non-Boolean domains
                            scope2.push_back(constrs[i]->getVar(j)->wcspIndex);
                        }
                    }
                    if (scope2.size() > 3) {
                        clq = FindClique(scope2, G, CorrespVertices);
                        if (!clq.empty()) {
                            b = true;
                            for (unsigned int j = 0; j < clq.size(); ++j) {
                                if ((int)clq[j].size() > MaxAMO)
                                    MaxAMO = (int)clq[j].size();
                                if ((int)clq[j].size() == (int)constrs[i]->arity())
                                    b = false;
                                if ((int)clq[j].size() > 1)
                                    total += 1;
                            }
                            if (b) {
                                count++;
                                k->addAMOConstraints(clq, vars, this);
                            }
                        }
                    }
                }
            }
        }

        if (ToulBar2::verbose >= 0) {
            if (count > 0)
                cout << count << " knapsack constraint" << ((count > 1) ? "s" : "") << " extended with " << total << " AMO constraint" << ((total > 1) ? "s" : "") << " of maximum arity " << MaxAMO;
            else
                cout << "No local AMO constraints added";
            cout << " in " << ((ToulBar2::parallel) ? (realTime() - startRealTime) : (cpuTime() - startCpuTime)) << " seconds." << endl;
        }
    }
    int count2 = 0;
    clq.clear();
    if (ToulBar2::addAMOConstraints != 0) {
        // see Katsirelos et al paper at CP 2017 for a better list of overlapping cliques
        // TODO: try graph coloring heuristics like DSATUR on the complementary graph for a better clique cover
        vector<int> Temp;
        vector<vector<int>> Tempclq;
        vector<int> order;
        for (unsigned int i = 0; i < num_vertices(G); ++i) {
            order.push_back(i);
        }
        // sort order by degree instead of taking the first variable
        assert(ucosts.size() == num_vertices(G));
        sort(order.begin(), order.end(),
            [&](int x, int y) {
                return (degree(x, G) > degree(y, G) || (degree(x, G) == degree(y, G) && ucosts[x] > ucosts[y]) || (degree(x, G) == degree(y, G) && ucosts[x] == ucosts[y] && x < y));
            });
        Tempclq.push_back({ order[0] });
        for (int i = 1; i < (int)num_vertices(G); ++i) {
            bool ok = false;
            int j = 0;
            int curr = order[i];
            while (!ok && j < (int)Tempclq.size()) {
                int k = 0;
                ok = true;
                while (ok && k < (int)Tempclq[j].size()) {
                    if (!edge(curr, Tempclq[j][k], G).second && !edge(Tempclq[j][k], curr, G).second)
                        ok = false;
                    k++;
                }
                j++;
            }
            if (ok) {
                assert(j > 0);
                Tempclq[j - 1].push_back(curr);
                // update Tempclq to keep a sorted list of cliques of decreasing size
                if (j - 1 > 0 && Tempclq[j - 1].size() > Tempclq[j - 2].size()) {
                    int jnew = j - 2;
                    while (jnew > 0 && Tempclq[j - 1].size() > Tempclq[jnew - 1].size()) {
                        jnew--;
                    }
                    assert(jnew < j - 1);
                    swap(Tempclq[jnew], Tempclq[j - 1]);
                }
            } else {
                Tempclq.push_back({ curr });
            }
        }
        if ((int)Tempclq.size() > abs(ToulBar2::addAMOConstraints))
            Tempclq.resize(abs(ToulBar2::addAMOConstraints));
        for (unsigned int i = 0; i < Tempclq.size(); ++i) {
            if (Tempclq[i].size() > 2) {
                clq.push_back({});
                for (unsigned int l = 0; l < Tempclq[i].size(); ++l) {
                    int f = 0;
                    while (CorrespVertices[f] < Tempclq[i][l] && f + 1 < (int)CorrespVertices.size() && Tempclq[i][l] >= CorrespVertices[f + 1]) {
                        f++;
                    }
                    clq.back().push_back(make_pair(vars[f]->wcspIndex, Tempclq[i][l] - CorrespVertices[f]));
                }
                // sort scope of cliques by wcspIndex (lexicographic order)
                sort(clq.back().begin(), clq.back().end(),
                    [&](pair<int, Value> x, pair<int, Value> y) {
                        return (x.first < y.first);
                    });
            }
        }
        AbstractNaryConstraint* cc = NULL;
        vector<EnumeratedVariable*> scopeVars;
        vector<vector<Long>> weights;
        vector<vector<Value>> VarVal, NotVarVal;
        vector<int> CorrAMO;
        MaxAMO = 0;
        for (int i = 0; i < (int)clq.size(); ++i) {
            int arity_ = (int)clq[i].size();
            if (arity_ > 2) {
                count2++;
                if (arity_ > MaxAMO)
                    MaxAMO = arity_;
                weights.resize(arity_, { 0, 1 });
                CorrAMO.resize(arity_, 0);
                VarVal.resize(arity_);
                NotVarVal.resize(arity_);
                scopeVars.resize(arity_);
                for (int j = 0; j < arity_; j++) {
                    scopeVars[j] = (EnumeratedVariable*)vars[clq[i][j].first];
                    int size = scopeVars[j]->getDomainSize();
                    Value* VV = new Value[size];
                    scopeVars[j]->getDomain(VV);
                    NotVarVal[j] = {};
                    for (int l = 0; l < size; l++) {
                        if (l != clq[i][j].second)
                            NotVarVal[j].push_back(VV[l]);
                    }
                    VarVal[j] = { clq[i][j].second, NotVarVal[j].front() };
                }
                cc = new KnapsackConstraint(this, scopeVars.data(), arity_, arity_ - 1, weights, arity_, VarVal,
                    NotVarVal, {}, weights, CorrAMO, CorrAMO, arity_);
                BinaryConstraint* bctr;
                TernaryConstraint* tctr = new TernaryConstraint(this);
                elimTernConstrs.push_back(tctr);
                for (int j = 0; j < 3; j++) {
                    if (!ToulBar2::vac && !ToulBar2::vac_prev)
                        bctr = new BinaryConstraint(this);
                    else
                        bctr = new VACBinaryConstraint(this);
                    elimBinConstrs.push_back(bctr);
                }
                cc->propagate();
            }
        }
        if (ToulBar2::verbose >= 0) {
            if (count2 > 0)
                cout << count2 << " AMO constraint" << ((count2 > 1) ? "s" : "") << " of maximum arity " << MaxAMO << " directly added to the problem";
            else
                cout << "No global AMO constraints added";
            cout << " in " << ((ToulBar2::parallel) ? (realTime() - startRealTime) : (cpuTime() - startCpuTime)) << " seconds." << endl;
        }
    }
    if (count > 0 || count2 > 0) {
        propagate();
    }
}

#endif

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
