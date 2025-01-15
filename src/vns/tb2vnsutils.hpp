/*
 * \file tb2vnsutils.hpp
 * \brief absorption of clusters for VNS-like algorithms
 *
 *  Created on: 12 December 2016
 *      Author: Abdelkader Ouali
 *      PhD Student: LITIO, University of Oran ; GREYC, University of Caen.
 */

#ifndef TB2VNSUTILS_HPP_
#define TB2VNSUTILS_HPP_
#ifdef BOOST

#include "tb2vns.hpp"
#include "tb2localsearch.hpp"
#include "tb2cpdgvns.hpp"

void fusionCluster(TCDGraph::vertex_descriptor v, TCDGraph::vertex_descriptor p, TCDGraph& cg);
void treeClusterFusion(TCDGraph::vertex_descriptor p, TCDGraph::vertex_descriptor v, TCDGraph& cg);
void cluster_graph_absorption(TCDGraph& input, TCDGraph& output);
void print_decomposition(ostream& os, TCDGraph& cg);

class TrueRem {
public:
    bool operator()(Cluster_edge e)
    {
        return true;
    }
};

class IsRemovable {
public:
    set<Cluster_edge> removable;
    IsRemovable(set<Cluster_edge>& s)
    {
        removable = s;
    }
    bool operator()(Cluster_edge e)
    {
        return removable.count(e) > 0;
    }
};

struct StructCluster {
    // Variables in the cluster
    int size;
    int* variablesIn;
    // Proper Variables
    int sizeProper;
    int* variablesProper;
    // Separators
    int sizeSeparator;
    int* variablesSeparator;
    // Shortcut to test the belonging
    int* belong;

    float BestTime;

    bool isIn(int variable)
    {
        return belong[variable];
    }
    bool isSeparator(int variable)
    {
        return (belong[variable] == 2);
    }
    bool isProper(int variable)
    {
        return (belong[variable] == 1);
    }

    void display()
    {
        cout << size << "\t(s=" << sizeSeparator << "/p=" << sizeProper
             << ") \t all{";
        for (int i = 0; i < size; i++) {
            if (i)
                cout << ",";
            if (belong[variablesIn[i]] == 2)
                cout << "\033[31m\033[1m -";
            cout << variablesIn[i];
            if (belong[variablesIn[i]] == 2)
                cout << "- \033[0m";
        }
        cout << "} ";
        // cout << "sep{";
        // for (int i = 0 ; i < sizeSeparator ; i++) {
        // if (i) cout << ",";
        // cout << variablesSeparator[i];
        // }
        // cout << "}";
        // cout << "sep{";
        // for (int i = 0 ; i < sizeProper ; i++) {
        // if (i) cout << ",";
        // cout << variablesProper[i];
        // }
        // cout << "}";
        cout << endl;
    }
    void display2()
    {
        ToulBar2::vnsOutput << size << "\t(s=" << sizeSeparator << "/p="
                            << sizeProper << ") \t all{";
        for (int i = 0; i < size; i++) {
            if (i)
                ToulBar2::vnsOutput << ",";
            if (belong[variablesIn[i]] == 2)
                ToulBar2::vnsOutput << " -";
            ToulBar2::vnsOutput << variablesIn[i];
            if (belong[variablesIn[i]] == 2)
                ToulBar2::vnsOutput << "- ";
        }
        ToulBar2::vnsOutput << "} ";
        ToulBar2::vnsOutput << endl;
    }
};

class TreeDecRefinement : public LocalSearch {
protected:
    // Decomposition
    TCDGraph m_graph;
    TCDGraph abs_graph; // graph absorption
    int nbClusters;
    int nbSeparators;
    vector<int> separators;
    vector<StructCluster> clusters;

public:
    TreeDecRefinement(Cost initUpperBound)
        : LocalSearch(initUpperBound)
        , nbClusters(0)
        , nbSeparators(0)
    {
    }
    ~TreeDecRefinement() {}

    bool solve(bool first = true) override;

    // decomposition tools
    void load_decomposition();

    //
    void print_dec_satistics();
    double ecart_type(vector<int>& data);
};

#endif
#endif /* TB2VNSUTILS_HPP_ */

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
