/*
 * ****** Graph Algorithms using Boost Graph Library ********
 *
 * Variable elimination ordering heuristics developed by Cyril Terrioux <cyril.terrioux@lsis.org>
 */

#include "core/tb2wcsp.hpp"
#include "core/tb2binconstr.hpp"
#include "core/tb2knapsack.hpp"
#include "core/tb2naryconstr.hpp"

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
}

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
        if (constrs[i]->connected())
            addConstraint(constrs[i], G);
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
static vector<vector<pair<int, int>>> FindClique(vector<int> scope, T& g)
{
    Graph G;
    for (unsigned int i = 0; i < scope.size(); ++i) {
        add_vertex(G);
        add_vertex(G);
    }
    for (unsigned int i = 0; i < scope.size(); i++) {
        for (unsigned int j = i + 1; j < scope.size(); j++) {
            if (edge(2 * scope[i], 2 * scope[j], g).second)
                add_edge(2 * i, 2 * j, G);
            if (edge(2 * scope[i] + 1, 2 * scope[j], g).second)
                add_edge(2 * i + 1, 2 * j, G);
            if (edge(2 * scope[i], 2 * scope[j] + 1, g).second)
                add_edge(2 * i, 2 * j + 1, G);
            if (edge(2 * scope[i] + 1, 2 * scope[j] + 1, g).second)
                add_edge(2 * i + 1, 2 * j + 1, G);
        }
    }
    vector<int> Temp;
    vector<vector<int>> Tempclq;
    vector<int> order;
    for (unsigned int i = 0; i < scope.size(); ++i) {
        order.push_back(i);
    }
    if (G.m_vertices[2 * order[0]].m_out_edges.size() > G.m_vertices[2 * order[0] + 1].m_out_edges.size())
        Temp.push_back(0);
    else
        Temp.push_back(1);
    Tempclq.push_back(Temp);
    bool ok;
    unsigned j, k, curr;
    for (int i = 1; i < (int)order.size(); i++) {
        ok = false;
        j = 0;
        curr = 2 * order[i] + 1;
        if (G.m_vertices[2 * order[i]].m_out_edges.size() > G.m_vertices[2 * order[i] + 1].m_out_edges.size())
            curr = 2 * order[order[i]];
        while (!ok && j < Tempclq.size()) {
            k = 0;
            ok = true;
            while (ok && k < Tempclq[j].size()) {
                if (!edge(curr, Tempclq[j][k], G).second && !edge(Tempclq[j][k], curr, G).second)
                    ok = false;
                k++;
            }
            j++;
        }
        if (ok) {
            Tempclq[j - 1].push_back(curr);
        } else {
            Temp.clear();
            Temp.push_back(curr);
            Tempclq.push_back(Temp);
        }
    }
    vector<vector<pair<int, int>>> clq;
    vector<pair<int, int>> clq2;
    for (unsigned int i = 0; i < Tempclq.size(); ++i) {
        clq2.clear();
        for (unsigned int l = 0; l < Tempclq[i].size(); ++l) {
            clq2.push_back(pair(scope[(int)floor(Tempclq[i][l] / 2.0 + 0.1)], Tempclq[i][l] % 2));
        }
        if (clq2.size() > 1)
            clq.push_back(clq2);
    }
    return clq;
}

void WCSP::addAMOConstraints()
{
    if (ToulBar2::verbose >= 1)
        cout << "Add AMO constraints to knapsack contraints." << endl;
    double startCpuTime = cpuTime();
    double startRealTime = realTime();

    vector<int> Var;
    Graph G;
    int count = 0;
    if (ToulBar2::verbose >= 1)
        cout << "Construct Graph of size " << vars.size() << endl;
    for (int i = 0; i < (int)vars.size(); ++i) {
        if (vars[i]->unassigned()) {
            Var.push_back(i);
            assert(i == vars[i]->wcspIndex);
        }
        add_vertex(G);
        add_vertex(G);
    }
    for (unsigned int varIndex = 0; varIndex < Var.size(); varIndex++) {
        int size = getDomainSize(Var[varIndex]);
        ValueCost sorted[size];
        getEnumDomainAndCost(Var[varIndex], sorted);
        qsort(sorted, size, sizeof(ValueCost), cmpValueCost3);
        for (int a = 0; a < size; a++) {
            int storedepth = Store::getDepth();
            try {
                Store::store();
                assign(Var[varIndex], sorted[a].value);
            } catch (const Contradiction&) {
            }
            for (unsigned int i = 0; i < Var.size(); ++i) {
                if (vars[Var[i]]->assigned() && i != varIndex) {
                    add_edge(2 * Var[varIndex] + sorted[a].value, 2 * Var[i] + 1 - vars[Var[i]]->getValue(), G);
                }
            }
            Store::restore(storedepth);
        }
    }
    if (ToulBar2::verbose >= 1)
        cout << "Graph done." << endl;
    count = 0;
    vector<int> scope;
    vector<int> scope2;
    vector<vector<pair<int, int>>> clq;
    int MaxAMO = 0;
    int total = 0;
    unsigned int nbconstrs = constrs.size();
    for (unsigned int i = 0; i < nbconstrs; ++i) {
        auto* k = dynamic_cast<KnapsackConstraint*>(constrs[i]);
        if (!k)
            continue;
        else {
            if (constrs[i]->arity() > 3 && constrs[i]->connected()) {
                scope.clear();
                scope2.clear();
                clq.clear();
                //scope=k->GetOrder();
                for (int j = 0; j < constrs[i]->arity(); j++) {
                    if (constrs[i]->getVar(j)->unassigned()) {
                        scope2.push_back(constrs[i]->getVar(j)->wcspIndex);
                    }
                }
                if (scope2.size() > 3) {
                    clq = FindClique(scope2, G);
                    if (clq.size() > 0) {
                        for (unsigned int j = 0; j < clq.size(); ++j) {
                            if ((int)clq[j].size() > MaxAMO)
                                MaxAMO = clq[j].size();
                            total += clq[j].size();
                        }
                        count++;
                        k->addAMOConstraints(clq, vars, this);
                    }
                }
            }
        }
    }
    if (ToulBar2::verbose >= 0) {
        if (count>0) cout << count << " AMO constraint" << ((count>1)?"s":"") << " added with max size " << MaxAMO;
        else cout << "No AMO constraint added";
        cout << " in " << ((ToulBar2::parallel) ? (realTime() - startRealTime) : (cpuTime() - startCpuTime)) << " seconds." << endl;
    }
}

#endif

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
