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
inline bool operator<(const Cluster_edge& __x, const Cluster_edge& __y)
{
    return __x.get_property() < __y.get_property();
}
}
#endif

extern void cluster_graph_absorption(TCDGraph& input, TCDGraph& output);
extern void print_decomposition(ostream& os, TCDGraph& cg);

class NeighborhoodStructure {
protected:
    WeightedCSP* wcsp;
    LocalSearch* l;

public:
    //initialization
    virtual void init(WeightedCSP* wcsp_, LocalSearch* l_) = 0;
    virtual const zone getNeighborhood(size_t neighborhood_size) = 0;
    virtual const zone getNeighborhood(size_t neighborhood_size, zone z) const = 0;
    virtual ~NeighborhoodStructure() {}
    virtual const bool incrementK() { return true; }
};

// for vns/lds-cp
class RandomNeighborhoodChoice : public NeighborhoodStructure {
public:
    virtual void init(WeightedCSP* wcsp_, LocalSearch* l_);
    virtual const zone getNeighborhood(size_t neighborhood_size);
    virtual const zone getNeighborhood(size_t neighborhood_size, zone z) const;
};

class ClustersNeighborhoodStructure : public NeighborhoodStructure {
protected:
    vector<int> clusters;
    TCDGraph m_graph;
    vector<int> file;
    NeighborhoodStructure* insideHeuristic;
    uint counter;
    uint precK;
    uint maxClusterSize;
    uint minClusterSize;

public:
    void load_decomposition();
    void printClusters(ostream& os) { print_decomposition(os, m_graph); }
    int getSize() const
    {
        return clusters.size();
    }
    int getMaxClusterSize() const
    {
        return maxClusterSize;
    }
    int getMinClusterSize() const
    {
        return minClusterSize;
    }
    double getMeanClusterSize() const;
    uint getMedianClusterSize() const;
};

// for dgvns
class RandomClusterChoice : public ClustersNeighborhoodStructure {
public:
    virtual void init(WeightedCSP* wcsp_, LocalSearch* l_);
    virtual const zone getNeighborhood(size_t neighborhood_size);
    virtual const zone getNeighborhood(size_t neighborhood_size, zone z) const;
    virtual const bool incrementK();
};

// for rpdgvns
class ParallelRandomClusterChoice : public ClustersNeighborhoodStructure {
public:
    virtual void init(WeightedCSP* wcsp_, LocalSearch* l_);
    virtual const zone getNeighborhood(size_t neighborhood_size);
    virtual const zone getNeighborhood(size_t neighborhood_size, zone z) const;
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
