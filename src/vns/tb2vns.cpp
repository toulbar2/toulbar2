/*
 * tb2vns.cpp
 *
 *  Created on: 3 mars 2015
 *      Authors: Mathieu Fontaine, Abdelkader Ouali
 *      Phd. Student : LITIO, University of Oran. GREYC, University of Caen.
 */

#ifdef BOOST
#include "tb2vns.hpp"
#include "core/tb2wcsp.hpp"
#include "search/tb2clusters.hpp"
#include <random>

using namespace boost;

// loading decomposition for only a given file
void ClustersNeighborhoodStructure::load_decomposition()
{
    if (ToulBar2::clusterFile != "") {
        bool cover = false;
        if (ToulBar2::clusterFile.find(".cov") != string::npos)
            cover = true;
        std::fstream file(ToulBar2::clusterFile.c_str());
        // file >> clustersNumber;
        set<int> nbvars;
        set<int> nbunvars;
        while (!file.eof()) {
            string cluster;
            getline(file, cluster);
            if (!cluster.empty()) {
                set<int> tmp;
                stringstream ss(cluster);
                if (cover) {
                    int clusterid, parentid;
                    ss >> clusterid;
                    ss >> parentid;
                }
                while (ss.good()) {
                    unsigned int var;
                    ss >> var;
                    nbvars.insert(var);
                    if (var >= wcsp->numberOfVariables()) {
                        cerr << "Error: cluster decomposition contains bad variable index!" << endl;
                        throw WrongFileFormat();
                    }
                    if (wcsp->unassigned(var)) {
                        tmp.insert(var);
                        nbunvars.insert(var);
                    }
                };
                if (tmp.size() > 0) {
                    TDCluster c = add_vertex(m_graph);
                    m_graph[c].vars = tmp;
                }
            }
        }
        file.close();
        if (nbvars.size() < wcsp->numberOfVariables()) {
            cout << "Warning: cluster decomposition has missing variables! (" << nbvars.size() << "!=" << wcsp->numberOfVariables() << ")" << endl;
        }
        assert(nbunvars.size() == wcsp->numberOfUnassignedVariables());
        TCDGraph::vertex_iterator v, vend, v2;
        int num = 0;
        for (tie(v, vend) = vertices(m_graph); v != vend; ++v) {
            num++;
            string name;
            vector<int> cl;
            ostringstream ss(name);
            ss << *v << ":(";
            bool first = true;
            for (set<int>::iterator i = m_graph[*v].vars.begin();
                 i != m_graph[*v].vars.end(); ++i) {
                if (not first)
                    ss << ",";
                ss << *i;
                first = false;
            }
            ss << ")";
            m_graph[*v].name = ss.str();
            for (v2 = v + 1; v2 != vend; ++v2) {
                set<int> separator;
                set<int> v_vars = m_graph[*v].vars;
                set<int> v2_vars = m_graph[*v2].vars;
                set_intersection(v_vars.begin(), v_vars.end(), v2_vars.begin(),
                    v2_vars.end(), inserter(separator, separator.begin()));

                if (separator.size() > 0) {
                    Cluster_edge sep;
                    tie(sep, tuples::ignore) = add_edge(*v, *v2, m_graph);
                    m_graph[sep].vars = separator;
                    m_graph[sep].size = (float)1 / separator.size();
                }
            }
        }
    } else if (ToulBar2::varOrder) {
        set<int> nbunvars;
        TreeDecomposition* td = new TreeDecomposition((WCSP*)wcsp);
        double time = cpuTime();
        td->buildFromOrder();
        int nc = td->getNbOfClusters();
        for (int i = 0; i < nc; i++) {
            Cluster* ct = td->getCluster(i);
            if (ct->getSep())
                ct->getSep()->deconnect(); // deconnect separator constraints
            set<int> unassignedvars;
            set<int> cvars = ct->getVars();
            for (TVars::iterator iter = cvars.begin(); iter != cvars.end(); ++iter)
                if (wcsp->unassigned(*iter)) {
                    unassignedvars.insert(*iter);
                    nbunvars.insert(*iter);
                }
            if (unassignedvars.size() > 0) {
                TDCluster c = add_vertex(m_graph);
                m_graph[c].name = to_string(i);
                m_graph[c].vars = unassignedvars;
            }
        }
        assert(nbunvars.size() == wcsp->numberOfUnassignedVariables());
        TCDGraph::vertex_iterator v, vend, v2;
        for (tie(v, vend) = vertices(m_graph); v != vend; ++v) {
            for (v2 = v + 1; v2 != vend; ++v2) {
                set<int> separator;
                td->intersection(m_graph[*v].vars, m_graph[*v2].vars, separator);
                if (separator.size() > 0) {
                    Cluster_edge sep;
                    tie(sep, tuples::ignore) = add_edge(*v, *v2, m_graph);
                    m_graph[sep].vars = separator;
                    m_graph[sep].size = (float)1 / separator.size();
                }
            }
        }
        if (ToulBar2::verbose >= 0)
            cout << "Tree decomposition time: " << cpuTime() - time << " seconds." << endl;
    } else {
        // make only one big cluster with all the problem variables
        TDCluster c = add_vertex(m_graph);
        m_graph[c].name = to_string(0);
        m_graph[c].vars = l->getUnassignedVars();
    }
    if (std::abs(ToulBar2::boostingBTD) > 0. && std::abs(ToulBar2::boostingBTD) < 1.) {
        TCDGraph abs_graph; // graph absorption
        cluster_graph_absorption(m_graph, abs_graph);
        m_graph = abs_graph;
    }
}

double ClustersNeighborhoodStructure::getMeanClusterSize() const
{
    uint nbc = getSize();
    if (nbc == 0)
        return 0;
    uint totalsize = 0;
    for (uint i = 0; i < nbc; i++)
        totalsize += m_graph[i].vars.size();
    return (double)totalsize / nbc;
}

uint ClustersNeighborhoodStructure::getMedianClusterSize() const
{
    uint nbc = getSize();
    if (nbc == 0)
        return 0;
    int csize[nbc];
    for (uint i = 0; i < nbc; i++)
        csize[i] = m_graph[i].vars.size();
    return stochastic_selection<int>(csize, 0, nbc - 1, nbc / 2);
}

void RandomNeighborhoodChoice::init(WeightedCSP* wcsp_, LocalSearch* l_)
{
    this->l = l_;
    wcsp = wcsp_;
}
const zone RandomNeighborhoodChoice::getNeighborhood(size_t neighborhood_size)
{
    zone neighborhood;
    vector<int> z(l->unassignedVars->getSize());
    unsigned int j = 0;

    for (BTList<Value>::iterator iter = l->unassignedVars->begin(); iter != l->unassignedVars->end(); ++iter) {
        z[j] = *iter;
        ++j;
    }
    shuffle(z.begin(), z.end(), myrandom_generator);
    assert(neighborhood_size <= z.size());
    neighborhood.insert(z.begin(), z.begin() + neighborhood_size);
    return neighborhood;
}
const zone RandomNeighborhoodChoice::getNeighborhood(size_t neighborhood_size, zone z) const
{
    zone neighborhood;
    vector<int> zv(z.begin(), z.end());

    shuffle(zv.begin(), zv.end(), myrandom_generator);
    assert(neighborhood_size <= zv.size());
    neighborhood.insert(zv.begin(), zv.begin() + neighborhood_size);
    return neighborhood;
}

void RandomClusterChoice::init(WeightedCSP* wcsp_, LocalSearch* l_)
{
    this->l = l_;
    this->wcsp = wcsp_;
    maxClusterSize = 0;
    minClusterSize = wcsp_->numberOfVariables();
    load_decomposition();
    TCDGraph::vertex_iterator v, vend;
    tie(v, vend) = vertices(m_graph);
    for (; v != vend; ++v) {
        TCDGraph::out_edge_iterator e, eend;
        float sum = 0.0;
        tie(e, eend) = out_edges(*v, (m_graph));
        for (; e != eend; ++e) {
            if (m_graph[*v].vars.size() > 0) {
                m_graph[*v].absorptions[target(*e, m_graph)] = m_graph[*e].vars.size()
                    / m_graph[*v].vars.size();
                sum += m_graph[*e].vars.size() / m_graph[*v].vars.size();
            }
        }
        m_graph[*v].absorption = sum / degree(*v, m_graph);
        clusters.push_back(*v);
        if (m_graph[*v].vars.size() > maxClusterSize)
            maxClusterSize = m_graph[*v].vars.size();
        if (m_graph[*v].vars.size() < minClusterSize)
            minClusterSize = m_graph[*v].vars.size();
        set<Constraint*> csts = m_graph[*v].consts;
        m_graph[*v].lastCost = 0;
        for (set<Constraint*>::iterator it = csts.begin(); it != csts.end();
             ++it) {
            m_graph[*v].lastCost = m_graph[*v].lastCost
                + (*it)->getConflictWeight();
        }
    }
    if (clusters.size() >= 1 && m_graph[clusters[clusters.size() - 1]].vars.empty()) {
        clusters.pop_back();
    }
    file = clusters;
    shuffle(file.begin(), file.end(), myrandom_generator);
    insideHeuristic = new RandomNeighborhoodChoice();
    insideHeuristic->init(wcsp, l);
}

const zone RandomClusterChoice::getNeighborhood(size_t neighborhood_size)
{
    assert(neighborhood_size <= wcsp->numberOfUnassignedVariables());
    set<int> selclusters;
    if (file.size() == 0) {
        file = clusters;
        shuffle(file.begin(), file.end(), myrandom_generator);
    }
    assert(file.size() > 0);
    int c = file.back();
    file.pop_back();
    selclusters.insert(c);
    if (ToulBar2::verbose >= 1)
        cout << "Select cluster " << c << endl;
    zone z = m_graph[c].vars;
    // if variables are missing
    if (z.size() < neighborhood_size) {
        queue<int> fifo;
        TCDGraph::adjacency_iterator v, vend;
        int currclu = c;
        int i = 0;
        do {
            tie(v, vend) = adjacent_vertices(currclu, m_graph);
            if (v != vend) {
                vector<int> neighbors(v, vend);
                // merge neighbors of current cluster
                shuffle(neighbors.begin(), neighbors.end(), myrandom_generator);
                // add them to list
                for (vector<int>::iterator it = neighbors.begin();
                     it != neighbors.end(); ++it) {
                    if (selclusters.count(*it) == 0) {
                        fifo.push(*it);
                        selclusters.insert(*it);
                    }
                }
            }
            if (fifo.size() == 0) {
                assert(selclusters.size() < clusters.size());
                while (selclusters.count(i) > 0) {
                    i = (i + 1) % clusters.size(); // TODO: use file instead in order to shuffle clusters randomly
                }
                fifo.push(clusters[i]);
                selclusters.insert(clusters[i]);
            }
            currclu = fifo.front();
            fifo.pop();
            z.insert(m_graph[currclu].vars.begin(),
                m_graph[currclu].vars.end());
        } while (z.size() < neighborhood_size);
    }
    return insideHeuristic->getNeighborhood(neighborhood_size, z);
}

const zone RandomClusterChoice::getNeighborhood(size_t neighborhood_size, zone z) const
{
    assert("not implemented!!!");
    return zone();
}

bool RandomClusterChoice::incrementK()
{
    if (file.size() == 0) {
        file = clusters;
        shuffle(file.begin(), file.end(), myrandom_generator);
        return true;
    }

    return true;
}

void ParallelRandomClusterChoice::init(WeightedCSP* wcsp_, LocalSearch* l_)
{
    this->l = l_;
    this->wcsp = wcsp_;
    maxClusterSize = 0;
    minClusterSize = wcsp_->numberOfVariables();
    load_decomposition();
    TCDGraph::vertex_iterator v, vend;
    tie(v, vend) = vertices(m_graph);
    for (; v != vend; ++v) {
        TCDGraph::out_edge_iterator e, eend;
        float sum = 0.0;
        tie(e, eend) = out_edges(*v, (m_graph));
        for (; e != eend; ++e) {
            if (m_graph[*v].vars.size() > 0) {
                m_graph[*v].absorptions[target(*e, m_graph)] = m_graph[*e].vars.size()
                    / m_graph[*v].vars.size();
                sum += m_graph[*e].vars.size() / m_graph[*v].vars.size();
            }
        }
        m_graph[*v].absorption = sum / degree(*v, m_graph);
        clusters.push_back(*v);
        if (m_graph[*v].vars.size() > maxClusterSize)
            maxClusterSize = m_graph[*v].vars.size();
        if (m_graph[*v].vars.size() < minClusterSize)
            minClusterSize = m_graph[*v].vars.size();
        set<Constraint*> csts = m_graph[*v].consts;
        m_graph[*v].lastCost = 0;
        for (set<Constraint*>::iterator it = csts.begin(); it != csts.end();
             ++it) {
            m_graph[*v].lastCost = m_graph[*v].lastCost
                + (*it)->getConflictWeight();
        }
    }
    if (clusters.size() >= 1 && m_graph[clusters[clusters.size() - 1]].vars.empty()) {
        clusters.pop_back();
    }
    file = clusters;
    shuffle(file.begin(), file.end(), myrandom_generator);
    insideHeuristic = new RandomNeighborhoodChoice();
    insideHeuristic->init(wcsp, l);
}

const zone ParallelRandomClusterChoice::getNeighborhood(size_t neighborhood_size)
{
    set<int> selclusters;
    int c = file.back();
    file.pop_back();
    selclusters.insert(c);
    zone z = m_graph[c].vars;
    // if variables are missing
    if (z.size() < neighborhood_size) {
        queue<int> fifo;
        TCDGraph::adjacency_iterator v, vend;
        int currclu = c;
        int i = 0;
        do {
            tie(v, vend) = adjacent_vertices(currclu, m_graph);
            if (v != vend) {
                vector<int> neighbors(v, vend);
                // merge neighbors of current cluster
                shuffle(neighbors.begin(), neighbors.end(), myrandom_generator);
                // add them to list
                for (vector<int>::iterator it = neighbors.begin();
                     it != neighbors.end(); ++it) {
                    if (selclusters.count(*it) == 0) {
                        fifo.push(*it);
                        selclusters.insert(*it);
                    }
                }
            }
            if (fifo.size() == 0) {
                assert(selclusters.size() < clusters.size());
                while (selclusters.count(i) > 0) {
                    i = (i + 1) % clusters.size(); // TODO: use file instead in order to shuffle clusters randomly
                }
                fifo.push(clusters[i]);
                selclusters.insert(clusters[i]);
            }
            currclu = fifo.front();
            fifo.pop();
            z.insert(m_graph[currclu].vars.begin(),
                m_graph[currclu].vars.end());
        } while (z.size() < neighborhood_size);
    }
    return insideHeuristic->getNeighborhood(neighborhood_size, z);
}

const zone ParallelRandomClusterChoice::getNeighborhood(size_t neighborhood_size, zone z) const
{
    assert("Not implemented!!!");
    return zone();
}

const zone ParallelRandomClusterChoice::SlaveGetNeighborhood(unsigned int CurrentCluster, size_t neighborhood_size)
{
    set<int> selclusters;
    selclusters.insert(CurrentCluster);
    zone z = m_graph[CurrentCluster].vars;
    // if variables are missing
    if (z.size() < neighborhood_size) {
        queue<int> fifo;
        TCDGraph::adjacency_iterator v, vend;
        int currclu = CurrentCluster;
        int i = 0;
        do {
            tie(v, vend) = adjacent_vertices(currclu, m_graph);
            if (v != vend) {
                vector<int> neighbors(v, vend);
                // merge neighbors of current cluster
                shuffle(neighbors.begin(), neighbors.end(), myrandom_generator);
                // add them to list
                for (vector<int>::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
                    if (selclusters.count(*it) == 0) {
                        fifo.push(*it);
                        selclusters.insert(*it);
                    }
                }
            }
            if (fifo.size() == 0) {
                assert(selclusters.size() < clusters.size());
                while (selclusters.count(i) > 0) {
                    i = (i + 1) % clusters.size(); // TODO: use file instead in order to shuffle clusters randomly
                }
                fifo.push(clusters[i]);
                selclusters.insert(clusters[i]);
            }
            currclu = fifo.front();
            fifo.pop();
            z.insert(m_graph[currclu].vars.begin(), m_graph[currclu].vars.end());
        } while (z.size() < neighborhood_size);
    }
    return insideHeuristic->getNeighborhood(neighborhood_size, z);
}

const zone ParallelRandomClusterChoice::SlaveGetNeighborhood(unsigned int CurrentCluster, uint number, size_t NeighborhoodSize)
{
    set<int> selclusters;
    selclusters.insert(CurrentCluster);
    zone z = m_graph[CurrentCluster].vars;
    // if variables are missing
    if (z.size() < NeighborhoodSize && number > 0) {
        queue<int> fifo;
        TCDGraph::adjacency_iterator v, vend;
        int currclu = CurrentCluster;
        uint numclu = 0;
        int i = 0;
        do {
            tie(v, vend) = adjacent_vertices(currclu, m_graph);
            if (v != vend) {
                vector<int> neighbors(v, vend);
                // merge neighbors of current cluster
                shuffle(neighbors.begin(), neighbors.end(), myrandom_generator);
                // add them to list
                for (vector<int>::iterator it = neighbors.begin();
                     it != neighbors.end(); ++it) {
                    if (selclusters.count(*it) == 0) {
                        fifo.push(*it);
                        selclusters.insert(*it);
                    }
                }
            }
            if (fifo.size()) {
                assert(selclusters.size() < clusters.size());
                while (selclusters.count(i) > 0) {
                    i = (i + 1) % clusters.size(); // TODO: use file instead in order to shuffle clusters randomly
                }
                fifo.push(clusters[i]);
                selclusters.insert(clusters[i]);
            }
            currclu = fifo.front();
            fifo.pop();
            z.insert(m_graph[currclu].vars.begin(),
                m_graph[currclu].vars.end());
            numclu++;
        } while (z.size() < NeighborhoodSize && numclu <= number);
    }

    if (z.size() < NeighborhoodSize) {
        NeighborhoodSize = z.size();
    }
    return insideHeuristic->getNeighborhood(NeighborhoodSize, z);
}

bool ParallelRandomClusterChoice::incrementK()
{
    if (file.size() == 0) {
        file = clusters;
        shuffle(file.begin(), file.end(), myrandom_generator);
        return true;
    }

    return true;
}

vector<int> ParallelRandomClusterChoice::getClustersIndex()
{
    return clusters;
}

uint ParallelRandomClusterChoice::getClustersSize(uint c, uint number)
{
    zone z = m_graph[c].vars;
    uint numclu = 0;
    if (number > 0) {
        set<int> selclusters;
        queue<int> fifo;
        TCDGraph::adjacency_iterator v, vend;
        int currclu = c;
        int i = 0;
        do {
            tie(v, vend) = adjacent_vertices(currclu, m_graph);
            if (v != vend) {
                vector<int> neighbors(v, vend);
                // merge neighbors of current cluster
                shuffle(neighbors.begin(), neighbors.end(), myrandom_generator);
                // add them to list
                for (vector<int>::iterator it = neighbors.begin();
                     it != neighbors.end(); ++it) {
                    if (selclusters.count(*it) == 0) {
                        fifo.push(*it);
                        selclusters.insert(*it);
                    }
                }
            }
            if (fifo.size() == 0) {
                assert(selclusters.size() < clusters.size());
                while (selclusters.count(i) > 0) {
                    i = (i + 1) % clusters.size(); // TODO: use file instead in order to shuffle clusters randomly
                }
                fifo.push(clusters[i]);
                selclusters.insert(clusters[i]);
            }
            currclu = fifo.front();
            fifo.pop();
            z.insert(m_graph[currclu].vars.begin(),
                m_graph[currclu].vars.end());
            numclu++;
        } while (numclu <= number);
    }
    return z.size();
}

#endif

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
