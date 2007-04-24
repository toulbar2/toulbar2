#include <boost/config.hpp>
#include <boost/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
//#include <boost/graph/dijkstra_shortest_paths.hpp>

using namespace boost;

#include "tb2wcsp.hpp"

typedef adjacency_list< vecS, vecS, undirectedS, no_property, 
                        property< edge_weight_t, int > > Graph;
typedef graph_traits< Graph >::vertex_descriptor vertex_ptr_t;
typedef graph_traits< Graph >::edge_descriptor   edge_ptr_t;

typedef int *int_ptr;
    
static Graph *WcspGraph;                             // link to Boost Graph internal graph structure

void WCSP::boostGraphConnection()
{
   if (::WcspGraph != NULL) delete ::WcspGraph;
   
   int V = numberOfVariables();
   int C = numberOfConstraints();
   ::WcspGraph = new Graph;
   // add nodes dynamically instead of using "const int n; Graph G(n);"
   for (int i=0; i < V; i++) {
        add_vertex(*::WcspGraph);
   }
   // fill edges with original constraints only
   property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, *::WcspGraph);
   for (int ictr=0; ictr < C; ictr++) {
        Constraint* c = constrs[ictr];

        int a = c->arity();
        for(int i=0;i<a;i++) {
            for(int j=i+1;j<a;j++) {
                Variable* vari = c->getVar(i);
                Variable* varj = c->getVar(j);
                edge_ptr_t e; 
                bool inserted;
                tie(e, inserted) = add_edge( vari->wcspIndex, varj->wcspIndex, *::WcspGraph);
                weightmap[e] = 1;
            }
        }
   }

//  output graph into file   
//   ofstream f("problem.dot");
//   write_graphviz(f, *::WcspGraph);
}

int WCSP::diameter()
{
   boostGraphConnection();
   
   int V = numberOfVariables();

//  vector<int> D(numberOfVariables());
//  vertex_ptr_t s = *(vertices(*::WcspGraph).first);
//  dijkstra_shortest_paths(*::WcspGraph, s, distance_map(&D[0]));

  int **D;
  D = new int_ptr[V];
  for (int i = 0; i < V; ++i) D[i] = new int[V];
  johnson_all_pairs_shortest_paths(*::WcspGraph, D);

  if (ToulBar2::verbose >= 2) {
      cout << "     ";
      for (int i = 0; i < V; ++i) {
        cout << i << " -> ";
        for (int j = 0; j < V; ++j) {
            cout << " " << D[i][j];
        }
        cout << endl;
      }
  }
  
  int maxd = 0;
  double meand = 0;
  for (int i = 0; i < V; ++i) {
    for (int j = 0; j < V; ++j) {
        if (D[i][j] > maxd) maxd = D[i][j];
        meand += D[i][j];
    }
  }
  meand /= V*V;
  if (ToulBar2::verbose >= 1) {
    cout << "Mean diameter: " << meand << endl;
  }

  for (int i = 0; i < V; ++i) delete[] D[i];
  delete[] D;

  return maxd;
}
	
	
