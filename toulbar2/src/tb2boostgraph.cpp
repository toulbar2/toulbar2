
#ifdef BOOST
#include <boost/config.hpp>
#include <boost/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/graph/minimum_degree_ordering.hpp>

using namespace boost;

namespace boost
{
  struct edge_component_t
  {
    enum
    { num = 555 };
    typedef edge_property_tag kind;
  }
  edge_component;
}

typedef adjacency_list< vecS, vecS, undirectedS, no_property, 
                        property< edge_weight_t, int, property < edge_component_t, std::size_t > > > Graph;
typedef graph_traits < Graph >::vertex_descriptor Vertex;
typedef graph_traits < Graph >::edge_descriptor Edge;
#endif

#include "tb2wcsp.hpp"

#ifdef BOOST
static void addConstraint(Constraint *c, Graph& g)
{
    property_map<Graph, edge_weight_t>::type weight = get(edge_weight, g);
    int a = c->arity();
    for(int i=0;i<a;i++) {
        for(int j=i+1;j<a;j++) {
            Variable* vari = c->getVar(i);
            Variable* varj = c->getVar(j);
            weight[add_edge( vari->wcspIndex, varj->wcspIndex, g).first] = 1;
        }
    }
}
#endif

#ifdef BOOST
int WCSP::connectedComponents()
{
    Graph G;
    for (unsigned int i=0; i<vars.size(); i++) add_vertex(G);
    for (unsigned int i=0; i<constrs.size(); i++) if (constrs[i]->connected() && !constrs[i]->universal()) addConstraint(constrs[i], G);
	for (int i=0; i<elimBinOrder; i++) if (elimBinConstrs[i]->connected()) addConstraint(elimBinConstrs[i], G);
	for (int i=0; i<elimTernOrder; i++) if (elimTernConstrs[i]->connected()) addConstraint(elimTernConstrs[i], G);
    vector<int> component(num_vertices(G));
    int num = connected_components(G, &component[0]);
    vector<int> cctruesize(num, 0);
    for (size_t i = 0; i < num_vertices(G); ++i) {
        assert(component[i]>=0 && component[i]<num);
        if (unassigned(i)) cctruesize[component[i]]++;
    }
    int res = 0;
    char c = '(';
    for (int  i = 0; i < num; ++i) {
        if (cctruesize[i] >= 1) {
            res++;
            cout << c << cctruesize[i];
            c = ' ';
        }
    }
    cout << ")";
    
    return res;
}
#endif

#ifdef BOOST 
int WCSP::biConnectedComponents()
{
    Graph G;
    for (unsigned int i=0; i<vars.size(); i++) add_vertex(G);
    for (unsigned int i=0; i<constrs.size(); i++) if (constrs[i]->connected()) addConstraint(constrs[i], G);
	for (int i=0; i<elimBinOrder; i++) if (elimBinConstrs[i]->connected()) addConstraint(elimBinConstrs[i], G);
	for (int i=0; i<elimTernOrder; i++) if (elimTernConstrs[i]->connected()) addConstraint(elimTernConstrs[i], G);
    property_map < Graph, edge_component_t >::type component = get(edge_component, G);

    int num = biconnected_components(G, component);

    vector<Vertex> art_points;
    articulation_points(G, back_inserter(art_points));
    cout << "Articulation points: " << art_points.size() <<  endl;
    for (unsigned int i=0; i<art_points.size(); i++) cout << " " << art_points[i];
	cout << endl;
    return num;
}
#endif


#ifdef BOOST     
int WCSP::diameter()
{
  Graph G;
  for (unsigned int i=0; i<vars.size(); i++) add_vertex(G);
  for (unsigned int i=0; i<constrs.size(); i++) if (constrs[i]->connected()) addConstraint(constrs[i], G);
  for (int i=0; i<elimBinOrder; i++) if (elimBinConstrs[i]->connected()) addConstraint(elimBinConstrs[i], G);
  for (int i=0; i<elimTernOrder; i++) if (elimTernConstrs[i]->connected()) addConstraint(elimTernConstrs[i], G);

  typedef int *int_ptr;
  int **D;
  D = new int_ptr[num_vertices(G)];
  for (unsigned int i = 0; i < num_vertices(G); ++i) D[i] = new int[num_vertices(G)];
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
        if (D[i][j] > maxd) maxd = D[i][j];
        meand += D[i][j];
    }
  }
  meand /= num_vertices(G)*num_vertices(G);
  if (ToulBar2::verbose >= 1) {
    cout << "Mean diameter: " << meand << endl;
  }

  for (unsigned int i = 0; i < num_vertices(G); ++i) delete[] D[i];
  delete[] D;

  return maxd;
}
#endif

inline bool cmp_vars(Variable *v1, Variable *v2) { return (v1->wcspIndex < v2->wcspIndex); }


#ifdef BOOST   
/// \bug reordering of vars array is dubious!!! (invalidates further use of variable indexes)    
void WCSP::minimumDegreeOrdering()
{
  Graph G;
  for (unsigned int i=0; i<vars.size(); i++) add_vertex(G);
  for (unsigned int i=0; i<constrs.size(); i++) if (constrs[i]->connected()) addConstraint(constrs[i], G);
  for (int i=0; i<elimBinOrder; i++) if (elimBinConstrs[i]->connected()) addConstraint(elimBinConstrs[i], G);
  for (int i=0; i<elimTernOrder; i++) if (elimTernConstrs[i]->connected()) addConstraint(elimTernConstrs[i], G);
  
  int delta = 0;
  int n = num_vertices(G);
  typedef vector<int> Vector;
  Vector inverse_perm(n, 0);
  Vector perm(n, 0);

  Vector supernode_sizes(n, 1); // init has to be 1

  property_map<Graph, vertex_index_t>::type id = get(vertex_index, G);

  Vector degree(n, 0);

  minimum_degree_ordering(G,
     make_iterator_property_map(&degree[0], id, degree[0]),
     &inverse_perm[0],
     &perm[0],
     make_iterator_property_map(&supernode_sizes[0], id, supernode_sizes[0]), 
     delta,
     id);
     
  cout << "Minimum Degree Order:";
  for (size_t i=0; i < num_vertices(G); ++i) {
    cout << " " << inverse_perm[i];
  }
  cout << endl;
  
  for (size_t i=0; i < num_vertices(G); ++i) {
    vars[i]->wcspIndex = num_vertices(G) - perm[i] - 1;
  }
  stable_sort(vars.begin(), vars.end(), cmp_vars);
  for (size_t i=0; i < num_vertices(G); ++i) {
    assert(vars[i]->wcspIndex == (int) i);
  }
  // update DAC ordering
  setDACOrder(inverse_perm);
}
#endif

