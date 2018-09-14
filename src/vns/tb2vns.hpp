/*
 * \file tb2vns.hpp
 * \brief various neighborhood structures for VNS-like algorithms
 *
 *  Created on: 3 mars 2015
 *      Author: Mathieu Fontaine, Abdelkader Ouali
 *      PhD Student: LITIO, University of Oran ; GREYC, University of Caen.
 */

#ifndef TB2VNS_HPP_
#define TB2VNS_HPP_

#include "tb2localsearch.hpp"
#include <set>
#include <queue>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <iostream>
#include <sstream>

#include <algorithm>
#include <boost/tokenizer.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/foreach.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/connected_components.hpp>

//TODO: is it need for other versions???
#define BOOSTGRAPH134
using namespace boost;

/**
 * Basic structure
 */

struct separator {
    string name;
    set<int> vars;
    float size;
};

struct cluster {
    string name;
    set<int> vars;
    set<Constraint*> consts;
    map<int, float> absorptions;
    int degree;
    Cost lastCost;
    float absorption;
    bool mark;
};

typedef property<vertex_index_t, int> variable_vertex;
typedef adjacency_list<vecS, vecS, undirectedS, variable_vertex, no_property, graph_name_t> TGraph;
typedef adjacency_list<vecS, vecS, undirectedS, cluster, separator> TCDGraph;
typedef graph_traits<TCDGraph>::vertex_descriptor TDCluster;
typedef graph_traits<TCDGraph>::edge_descriptor Cluster_edge;
typedef set<int> zone;

#ifdef BOOSTGRAPH134
namespace boost {
inline bool operator<(const Cluster_edge& __x, const Cluster_edge&__y)
{
    return __x.get_property() < __y.get_property();
}
}
#endif

class NeighborhoodStructure {
protected:
    WeightedCSP* wcsp;
    LocalSearch* l;
public:
    //initialization
    virtual void init(WeightedCSP* wcsp_, LocalSearch* l_)=0;
    virtual const zone getNeighborhood(size_t neighborhood_size) =0;
    virtual const zone getNeighborhood(size_t neighborhood_size, zone z) =0;
    virtual ~NeighborhoodStructure() {}
    virtual const bool incrementK() { return true; }
};

// for vns/lds-cp
class RandomNeighborhoodChoice: public NeighborhoodStructure {

protected:
	//
	zone last_zone;
	// Toulbar2::AA_vector index of the last first element chosen for the last neighborhood
	// Used for a diversification step (if < 0 this value is not considered) 
	unsigned int last_first_selection;
	// Choose an index according to probabilities
	unsigned int getIndexWithProba(vector<double> proba);
	// Select a point randomly and uniformly 
	unsigned int randFirstChoice(zone z);
	// Select a point according to its =1/distance from the molecule's centroid
	unsigned int probaFirstChoice(zone z);
	// Select a point according to its =distance from the molecule's centroid
	unsigned int diversifyFirstChoice(zone z);
	// Select the first point for the algorithms which build the neighborhoods
	unsigned int getFirstElt(size_t neighborhood_size, zone z); 
	// Build a neighborhood selecting point randomly and uniformly in the zone
	const zone directShuffle(size_t neighborhood_size, zone z);
	// Build a neighborhood selecting points in the zone according to their distances from the first selected element
	const zone probaDistance(size_t neighborhood_size, zone z);
	// Build a neighborhood with the nearest points of the first selected element
	const zone kNearest(size_t neighborhood_size, zone z);
	// Display data of a neighborhood and its first selected element
	void printNeighborhood(zone z, unsigned int first_selected_elt);

public:
    virtual void init(WeightedCSP* wcsp_, LocalSearch* l_);
    virtual const zone getNeighborhood(size_t neighborhood_size);
    virtual const zone getNeighborhood(size_t neighborhood_size, zone z);

};

class ClustersNeighborhoodStructure: public NeighborhoodStructure {
protected:
    vector<int> clusters;
    TCDGraph m_graph;
    vector<int> file;
    NeighborhoodStructure* insideHeuristic;
    uint counter;
    uint precK;
public:
    void load_decomposition();
    const int getSize()
    {
        return clusters.size();
    }
    TCDGraph getGraph()
    {
        return m_graph;
    }
};

// for dgvns
class RandomClusterChoice: public ClustersNeighborhoodStructure {
public:
    virtual void init(WeightedCSP* wcsp_, LocalSearch* l_);
    virtual const zone getNeighborhood(size_t neighborhood_size);
    virtual const zone getNeighborhood(size_t neighborhood_size, zone z);
    virtual const bool incrementK();

};

// for rpdgvns
class ParallelRandomClusterChoice: public ClustersNeighborhoodStructure {
public:
    virtual void init(WeightedCSP* wcsp_, LocalSearch* l_);
    virtual const zone getNeighborhood(size_t neighborhood_size);
    virtual const zone getNeighborhood(size_t neighborhood_size, zone z);
    // Master / Slave
    virtual const zone SlaveGetNeighborhood(uint CurrentCluster, size_t neighborhood_size);
    virtual const zone SlaveGetNeighborhood(uint CurrentCluster, uint number, size_t NeighborhoodSize);
    virtual const bool incrementK(); // for master process
    virtual vector<int> getClustersIndex();
    virtual uint getClustersSize(uint c, uint number);
};

#endif /* TB2VNS_HPP_ */

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
