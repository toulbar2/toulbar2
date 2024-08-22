/*
 * tb2vnsutils.cpp
 *
 *  Created on: 12 December 2016
 *      Author: Abdelkader Ouali
 *      PhD Student: LITIO, University of Oran ; GREYC, University of Caen.
 */

#ifdef BOOST
#include "tb2vnsutils.hpp"
#include "core/tb2wcsp.hpp"

using namespace boost;

void fusionCluster(TCDGraph::vertex_descriptor v, TCDGraph::vertex_descriptor p, TCDGraph& cg)
{
    set<int> varsfusion;
    set_union(cg[v].vars.begin(), cg[v].vars.end(), cg[p].vars.begin(),
        cg[p].vars.end(), inserter(varsfusion, varsfusion.begin()));
    cg[p].vars = varsfusion;
    cg[v].mark = true;
}

void treeClusterFusion(TCDGraph::vertex_descriptor p, TCDGraph::vertex_descriptor v, TCDGraph& cg)
{
    TCDGraph::vertex_iterator vend;
    tie(tuples::ignore, vend) = vertices(cg);
    if (degree(v, cg) > 1 || p == *vend) // on est pas sur une feuille
    {
        TCDGraph::adjacency_iterator a, aend;

        tie(a, aend) = adjacent_vertices(v, cg);
        for (; a != aend; ++a) {
            if (*a != p) {

                treeClusterFusion(v, *a, cg);
            }
        }
    }
    if (p == *vend)
        return;
    set<int> varsinter;
    set_intersection(cg[v].vars.begin(), cg[v].vars.end(), cg[p].vars.begin(),
        cg[p].vars.end(), inserter(varsinter, varsinter.begin()));
    //    int size_res = cg[v].vars.size() + cg[p].vars.size() - varsinter.size();

    // si absorber à 90% merge whatever
    if ((float)varsinter.size() / (float)cg[v].vars.size() >= std::abs(ToulBar2::boostingBTD)
        || (float)varsinter.size() / (float)cg[p].vars.size() >= std::abs(ToulBar2::boostingBTD)) {
        fusionCluster(v, p, cg);
        //    } else {
        //        // si absorber à 70% merge whatever merged size <= 100
        //        if ((float)varsinter.size() / (float)cg[v].vars.size() >= 0.7
        //            || (float)varsinter.size() / (float)cg[p].vars.size() >= 0.7) {
        //            if (size_res <= 100)
        //                fusionCluster(v, p, cg);
        //        } else {
        //            // si absorber à 50% merge whatever merged size <= 50
        //            if ((float)varsinter.size() / (float)cg[v].vars.size() >= 0.5
        //                || (float)varsinter.size() / (float)cg[p].vars.size() >= 0.5) {
        //                if (size_res <= 50)
        //                    fusionCluster(v, p, cg);
        //            }
        //        }
    }
}

bool TreeDecRefinement::solve(bool first)
{
    load_decomposition();
    cluster_graph_absorption(m_graph, abs_graph);
    print_decomposition(ToulBar2::vnsOutput, abs_graph);
    return true;
}

void print_decomposition(ostream& os, TCDGraph& m_graph)
{
    TCDGraph::vertex_iterator v, vend, v2;
    unsigned int i = 0;
    for (tie(v, vend) = vertices(m_graph); v != vend; ++v) {
        zone z = m_graph[i].vars;
        for (zone::iterator it = z.begin(); it != z.end(); ++it)
            os << (int)*it << " ";
        os << endl;
        i++;
    }
}

void TreeDecRefinement::print_dec_satistics()
{
    ToulBar2::vnsOutput << "Number of clusters = " << nbClusters << endl;
    ToulBar2::vnsOutput
        << "Number of separator variables = " << nbSeparators << " => "
        << double(nbSeparators) / double(wcsp->numberOfVariables()) << endl;
    int imax = 0, imin = 0;
    int max_size = clusters[0].size;
    int min_size = clusters[0].size;
    vector<int> data_sizes;
    vector<int> data_propers;
    vector<int> data_separators;
    for (int i = 1; i < nbClusters; i++) {
        if (clusters[i].size > max_size) {
            imax = i;
            max_size = clusters[i].size;
        }
        if (clusters[i].size < min_size) {
            imin = i;
            min_size = clusters[i].size;
        }
        data_sizes.push_back(clusters[i].size);
        data_propers.push_back(clusters[i].sizeProper);
        data_separators.push_back(clusters[i].sizeSeparator);
    }
    double ecartype_size = ecart_type(data_sizes);
    double ecartype_propers = ecart_type(data_propers);
    double ecartype_separators = ecart_type(data_separators);
    ToulBar2::vnsOutput << "Cluster with maximum size:" << endl;
    ToulBar2::vnsOutput << "Number of variables: " << clusters[imax].size
                        << endl;
    ToulBar2::vnsOutput << "Number of proper variables: " << clusters[imax].sizeProper
                        << endl;
    ToulBar2::vnsOutput << "Number of separators: "
                        << clusters[imax].sizeSeparator << endl;

    ToulBar2::vnsOutput << "Cluster with minimum size:" << endl;
    ToulBar2::vnsOutput << "Number of variables: " << clusters[imin].size
                        << endl;
    ToulBar2::vnsOutput << "Number of proper variables: " << clusters[imin].sizeProper
                        << endl;
    ToulBar2::vnsOutput << "Number of separators: "
                        << clusters[imin].sizeSeparator << endl;

    ToulBar2::vnsOutput << "Clusters deviation (standard deviation):" << endl;
    ToulBar2::vnsOutput << "Deviation on variables size: " << ecartype_size
                        << endl;
    ToulBar2::vnsOutput << "Deviation on proper variables size: " << ecartype_propers
                        << endl;
    ToulBar2::vnsOutput << "Deviation on separators size: "
                        << ecartype_separators << endl;

    ToulBar2::vnsOutput << endl;
}

double TreeDecRefinement::ecart_type(vector<int>& data)
{
    double moy = 0.0; // mean
    double ecart_type = 0.0; // standard deviation
    vector<double> diff; // array of differences between data and mean
    diff.reserve(data.size());
    diff.resize(data.size());

    for (uint i = 0; i < data.size(); i++)
        moy += data[i];
    moy /= data.size();

    // compute standard deviation
    for (uint i = 0; i < data.size(); i++) {
        diff[i] = pow(data[i] - moy, 2);
        ecart_type += diff[i];
    }
    ecart_type = sqrt(ecart_type / data.size());
    return ecart_type;
}

void TreeDecRefinement::load_decomposition()
{
    if (ToulBar2::clusterFile != "") {
        bool cover = false;
        if (ToulBar2::clusterFile.find(".cov") != string::npos)
            cover = true;
        std::fstream file(ToulBar2::clusterFile.c_str());
        set<int> nbvars;
        while (!file.eof()) {
            string cluster;
            getline(file, cluster);
            if (!cluster.empty()) {
                TDCluster c = add_vertex(m_graph);
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
                    tmp.insert(var);
                    nbvars.insert(var);
                    if (var >= wcsp->numberOfVariables()) {
                        cerr << "Error: cluster decomposition contains bad variable index!" << endl;
                        throw WrongFileFormat();
                    }
                };
                m_graph[c].vars = tmp;
            }
        }
        file.close();
        if (nbvars.size() < wcsp->numberOfVariables()) {
            cout << "Warning: cluster decomposition has missing variables! (" << nbvars.size() << "!=" << wcsp->numberOfVariables() << ")" << endl;
        }
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
    } else {
        cerr << "No cluster decomposition file!" << endl;
        throw WrongFileFormat();
    }
}

void cluster_graph_absorption(TCDGraph& m_graph, TCDGraph& abs_graph)
{
    TGraph::vertex_iterator vi, viend;
    TCDGraph::vertex_iterator v, vend, v2;
    TCDGraph::edge_iterator e, eend;

    for (tie(v, vend) = vertices(m_graph); v != vend; ++v) {
        set<int> cl;
        for (set<int>::iterator i = m_graph[*v].vars.begin();
             i != m_graph[*v].vars.end(); ++i) {
            cl.insert(*i);
        }
        TDCluster c = add_vertex(abs_graph);
        abs_graph[c].vars = cl;
        abs_graph[c].mark = false;
    }
    for (tie(v, vend) = vertices(abs_graph); v != vend; ++v) {
        string name;
        vector<int> cl;
        ostringstream ss(name);
        ss << *v << ":(";
        bool first = true;
        for (set<int>::iterator i = abs_graph[*v].vars.begin();
             i != abs_graph[*v].vars.end(); ++i) {
            if (not first)
                ss << ",";
            ss << *i;
            first = false;
        }
        ss << ")";
        abs_graph[*v].name = ss.str();
        for (v2 = v + 1; v2 != vend; ++v2) {
            set<int> separator;
            set<int> v_vars = abs_graph[*v].vars;
            set<int> v2_vars = abs_graph[*v2].vars;
            set_intersection(v_vars.begin(), v_vars.end(), v2_vars.begin(),
                v2_vars.end(), inserter(separator, separator.begin()));

            if (separator.size() > 0) {
                Cluster_edge sep;
                tie(sep, tuples::ignore) = add_edge(*v, *v2, abs_graph);
                abs_graph[sep].vars = separator;
                abs_graph[sep].size = (float)1 / separator.size();
            }
        }
    }

    // creation du spanning tree pour la fusion
    tie(e, eend) = edges(abs_graph);
    set<Cluster_edge> edges;
    for (; e != eend; ++e)
        edges.insert(*e);
    set<Cluster_edge> spanning_tree;
    set<Cluster_edge> removable_edges;
    kruskal_minimum_spanning_tree(abs_graph,
        inserter(spanning_tree, spanning_tree.begin()),
        weight_map(get(&separator::size, abs_graph)));
    set_difference(edges.begin(), edges.end(), spanning_tree.begin(),
        spanning_tree.end(),
        inserter(removable_edges, removable_edges.begin()));
    remove_edge_if(IsRemovable(removable_edges), abs_graph);
    std::vector<int> component(num_vertices(abs_graph));
    int num = connected_components(abs_graph, &component[0]);
    vector<int> max(num, 0);
    vector<int> maxv(num, 0);
    tie(v, vend) = vertices(abs_graph);
    for (; v != vend; ++v)
        if (max[component[*v]] < (int)abs_graph[*v].vars.size()) {
            max[component[*v]] = abs_graph[*v].vars.size();
            maxv[component[*v]] = *v;
        }
    for (int i = 0; i < num; ++i)
        treeClusterFusion(*vend, maxv[i], abs_graph);

    remove_edge_if(TrueRem(), abs_graph);
    tie(v, vend) = vertices(abs_graph);

    while (v != vend) {
        if (abs_graph[*v].mark) {
            remove_vertex(*v, abs_graph);
            tie(v, vend) = vertices(abs_graph);
        } else {
            ++v;
        }
    }
    num = 0;
    for (tie(v, vend) = vertices(abs_graph); v != vend; ++v) {
        num++;
        string name;
        vector<int> cl;
        ostringstream ss(name);
        ss << *v << ":(";
        bool first = true;
        for (set<int>::iterator i = abs_graph[*v].vars.begin();
             i != abs_graph[*v].vars.end(); ++i) {
            if (not first)
                ss << ",";
            ss << *i;
            first = false;
        }
        ss << ")";
        abs_graph[*v].name = ss.str();
        for (v2 = v + 1; v2 != vend; ++v2) {
            set<int> separator;
            set<int> v_vars = abs_graph[*v].vars;
            set<int> v2_vars = abs_graph[*v2].vars;
            set_intersection(v_vars.begin(), v_vars.end(), v2_vars.begin(),
                v2_vars.end(), inserter(separator, separator.begin()));

            if (separator.size() > 0) {
                Cluster_edge sep;
                tie(sep, tuples::ignore) = add_edge(*v, *v2, abs_graph);
                abs_graph[sep].vars = separator;
                abs_graph[sep].size = (float)1 / separator.size();
            }
        }
    }
}

#endif

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
