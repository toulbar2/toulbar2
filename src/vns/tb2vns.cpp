/*
 * tb2vns.cpp
 *
 *  Created on: 3 mars 2015
 *      Authors: Mathieu Fontaine, Abdelkader Ouali
 *      Phd. Student : LITIO, University of Oran. GREYC, University of Caen.
 */

#include "tb2config.hpp"

#ifdef BOOST
#include "tb2vns.hpp"
#include "core/tb2wcsp.hpp"
#include "search/tb2clusters.hpp"
#include <random>
#include <iomanip>
#include <fstream>
#include <algorithm>

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
        for (tie(v, vend) = vertices(m_graph); v != vend; ++v) {
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



// RandomNeighborhoodChoice
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


//NaturelNeighborhoodChoice
void NaturelNeighborhoodChoice::init(WeightedCSP* wcsp_, LocalSearch* l_)
{
    this->l = l_;
    wcsp = wcsp_;
    counter = 0; // on repart toujours du début de la liste au démarrage
}
const zone NaturelNeighborhoodChoice::getNeighborhood(size_t neighborhood_size)
{
    zone neighborhood;
    vector<int> z(l->unassignedVars->getSize());
    unsigned int j = 0;


    //BTList, elle contient les variables non affectées dans l'ordre de leur index.
    for (BTList<Value>::iterator iter = l->unassignedVars->begin(); iter != l->unassignedVars->end(); ++iter) {
        z[j] = *iter;
        ++j;
    }
    
    
    if (counter + neighborhood_size > z.size())
        counter = 0; // on a parcouru toutes les variables, on repart au début
    assert(neighborhood_size <= z.size());
    
    //cout << "counter: " << counter << " neighborhood_size: " << neighborhood_size << " z.size(): " << z.size() << endl;
    neighborhood.insert(z.begin() + counter, z.begin() + counter + neighborhood_size);
    counter += neighborhood_size;  // on avance d'un bloc pour le prochain appel
    return neighborhood;
}
const zone NaturelNeighborhoodChoice::getNeighborhood(size_t neighborhood_size, zone z) const
{
    zone neighborhood;
    vector<int> zv(z.begin(), z.end());

    

    assert(neighborhood_size <= zv.size());
    neighborhood.insert(zv.begin(), zv.begin() + neighborhood_size);
    return neighborhood;
}
// ProteinNeighborhoodChoice

void ProteinNeighborhoodChoice::printBilan()
{
    if (ToulBar2::showvns >= 1 || ToulBar2::verbose >= 1) {
        cout << "[Vns geode] Summary: clusters visited for optimum=" << clustersVisitedAtBest 
             << "/" << clusters.size() 
             << " | number improvement=" << nbImprovements 
             << " | time on best cluster=" << std::fixed << std::setprecision(3) << timeOnBestCluster << "s" << endl;
    }

}

void ProteinNeighborhoodChoice::init(WeightedCSP* wcsp_, LocalSearch* l_)
{
    this->wcsp = wcsp_;
    this->l = l_;

    if (clustersBuilt) {
        needsKReset = false;
        cycleComplete = false;
        clustersVisitedAtBest = nbVisitedClusters + 1;
        nbVisitedClusters = 0;
        timeOnBestCluster = cpuTime() - clusterEntryTime;
        clusterEntryTime = cpuTime();
        nbImprovements++;
        if (ToulBar2::vnsAdaptive) {
            lastAggregatedCluster = currentClusterIdx;
            currentZoneSize = (int)clusters[currentClusterIdx].size();
        }

        if (ToulBar2::showvns >= 1 || ToulBar2::verbose >= 1) {
            cout << "[Vns geode] Improvement on cluster_idx=" << currentClusterIdx 
                 << " (var=" << clusterRootWcspIdx[currentClusterIdx] << ")"
                 << " | k=kmin | number improvement=" << nbImprovements << endl;
        }
        return;
    }

    if (ToulBar2::vnsGeode <= 0) {
        cerr << "Error: with -vns=13 (PROTEIN mode), you must provide -geode=<radius> with radius > 0." << endl;
        throw BadConfiguration();
    }

    buildClusters(ToulBar2::vnsGeode);

    if (clusters.empty()) {
        cerr << "Error: no unassigned variable found, cannot build protein clusters." << endl;
        throw BadConfiguration();
    }

    size_t minSize = clusters[0].size();
    size_t maxSize = 0;
    size_t totalSize = 0;
    for (const auto& c : clusters) {
        if (c.size() < minSize) minSize = c.size();
        if (c.size() > maxSize) maxSize = c.size();
        totalSize += c.size();
    }

    if (ToulBar2::showvns >= 1 || ToulBar2::verbose >= 1) {
        cout << "[Vns geode] " << clusters.size() << " clusters | min_size=" << minSize
             << " mean_size=" << std::fixed << std::setprecision(1) << (totalSize / (double)clusters.size())
             << " max_size=" << maxSize << " | radius=" << ToulBar2::vnsGeode 
             << " | reverse=" << (ToulBar2::vnsReverseOrder ? "true" : "false") << endl;
        cout << "[Vns geode] Start: cluster=0 (var=" << clusterRootWcspIdx[0] << ")" << endl;
    }

    currentClusterIdx = 0;
    needsKReset = false;
    cycleComplete = false;
    nbVisitedClusters = 0;
    timeOnBestCluster = 0.0;
    clusterEntryTime = cpuTime();
    if (ToulBar2::vnsAdaptive) {
        currentZoneSize = (int)clusters[0].size();
    }

    if (ToulBar2::showvns >= 2 || ToulBar2::verbose >= 1) {
        cout << "[Vns geode] Mapping cluster_idx -> variable WCSP root" << endl;
        for (size_t i = 0; i < clusters.size(); i++) {
            cout << "[Vns geode] cluster_idx=" << i 
                 << " -> var=" << clusterRootWcspIdx[i] 
                 << " (size ball=" << clusters[i].size() << ") | content: ";
            for (size_t j = 0; j < clusters[i].size(); j++) {
                cout << clusters[i][j] << (j < clusters[i].size() - 1 ? ";" : "");
            }
            cout << endl;
        }
    }
    
    clustersBuilt = true;
}

const zone ProteinNeighborhoodChoice::getNeighborhood(size_t neighborhood_size)
{
    assert(neighborhood_size <= wcsp->numberOfUnassignedVariables());
    assert(currentClusterIdx >= 0 && currentClusterIdx < (int)clusters.size());

    set<int> z;
    int idx = currentClusterIdx;
    int nbClusters = (int)clusters.size();

    int aggregated = 0;
    int savedLastAggregated = lastAggregatedCluster;
    lastAggregatedCluster = currentClusterIdx;
    if (ToulBar2::vnsAdaptive) {
        int targetClusters = (savedLastAggregated - currentClusterIdx + nbClusters) % nbClusters + 1;
    while (aggregated < targetClusters) {
        for (int v : clusters[idx]) z.insert(v);
        aggregated++;
        idx = (idx + 1) % nbClusters;
    }
    lastAggregatedCluster = (currentClusterIdx + targetClusters - 1) % nbClusters;
} else {
    while (z.size() < neighborhood_size && aggregated < nbClusters) {
        for (int v : clusters[idx]) {
            z.insert(v);
        }
        lastAggregatedCluster = idx;
        aggregated++;
        idx = (idx + 1) % nbClusters;
    }
}

    if (aggregated >= nbClusters) {
        needsKReset = true;
    }

    if (ToulBar2::showvns >= 2 || ToulBar2::verbose >= 1) {
        cout << "[Vns geode] cluster_idx=" << currentClusterIdx
             << " (var=" << clusterRootWcspIdx[currentClusterIdx] << ")"
             << " | k=" << neighborhood_size
             << " | cluster size=" << z.size()
             << " | cluster aggregate=";
        for (int c = 0; c < aggregated; c++) {
            int aggIdx = (currentClusterIdx + c) % (int)clusters.size();
            cout << aggIdx << (c < aggregated - 1 ? ";" : "");
        }
        cout << endl;
    }


    lastRepairTime = cpuTime();
    zone neighborhood;
    if (ToulBar2::vnsAdaptive) {
        neighborhood = z;
        ToulBar2::vnsKcur = (int)z.size();
    } else {
        size_t actual_size = min(neighborhood_size, z.size());
        neighborhood.insert(z.begin(), next(z.begin(), actual_size));
    }
    return neighborhood;
}

const zone ProteinNeighborhoodChoice::getNeighborhood(size_t neighborhood_size, zone z) const
{
    assert("not implemented!!!");
    return zone();
}

bool ProteinNeighborhoodChoice::incrementK()
{
    double timeSpent = cpuTime() - lastRepairTime;

    if (ToulBar2::vnsTLimit > 0 && timeSpent >= ToulBar2::vnsTLimit) {

        int nextCluster;
        if (lastAggregatedCluster == currentClusterIdx) {
            // pas d'agrégation externe : on saute au cluster suivant
            nextCluster = (currentClusterIdx + 1) % (int)clusters.size();
        } else {
            nextCluster = lastAggregatedCluster - 1;
            if (nextCluster < 0) {
                nextCluster += (int)clusters.size();
            }
            // garde-fou : si l'avant-dernier == cluster racine, on prend lastAggregated
            if (nextCluster == currentClusterIdx) {
                nextCluster = lastAggregatedCluster;
            }
        }
        int jumpSize = nextCluster - currentClusterIdx;
        if (jumpSize <= 0) jumpSize += (int)clusters.size();

        currentClusterIdx = nextCluster;
        nbVisitedClusters += jumpSize;
        clusterEntryTime = cpuTime();
        needsKReset = false;
        
        if (ToulBar2::vnsAdaptive) {
            lastAggregatedCluster = currentClusterIdx;
            currentZoneSize = (int)clusters[currentClusterIdx].size();
        }


        if (nbVisitedClusters >= (int)clusters.size()) {
            cycleComplete = true;
            nbVisitedClusters = 0;
        }
        if (ToulBar2::showvns >= 1) {
            cout << "[Vns geode] tLimit reached (" << std::fixed << std::setprecision(1)
                 << timeSpent << "s > " << ToulBar2::vnsTLimit
                 << "s). Jumping to cluster_idx=" << currentClusterIdx;
            if (ToulBar2::vnsAdaptive) {
                cout << " (k adapted to new root size)" << endl;
            } else {
                cout << " (k remains unchanged)" << endl;
            }
        }
        return false;

    }

    if (ToulBar2::vnsAdaptive) {
        // agrège le cluster suivant et met à jour currentZoneSize
        int nextIdx = (lastAggregatedCluster + 1) % (int)clusters.size();
        set<int> newZone;
        // reconstruire la zone courante
        int idx = currentClusterIdx;
        int aggregated = 0;
        while (aggregated <= (lastAggregatedCluster - currentClusterIdx + (int)clusters.size()) % (int)clusters.size()) {
            for (int v : clusters[idx]) newZone.insert(v);
            idx = (idx + 1) % (int)clusters.size();
            aggregated++;
        }
        // ajouter le prochain cluster
        for (int v : clusters[nextIdx]) newZone.insert(v);
        currentZoneSize = (int)newZone.size();
        lastAggregatedCluster = nextIdx;
        needsKReset = true;

        if (nextIdx == currentClusterIdx) {
            if (ToulBar2::showvns >= 1) {
                cout << "[Vns geode] Full topological cycle completed. Increasing LDS." << endl;
            }
            ToulBar2::vnsKcur = ToulBar2::vnsKmax + 1;
            needsKReset = false;
            return true;
        } else if ((int)newZone.size() >= ToulBar2::vnsKmax) {
            if (ToulBar2::showvns >= 1) {
                cout << "[Vns geode] Zone size reached kmax. Increasing LDS." << endl;
            }
            ToulBar2::vnsKcur = ToulBar2::vnsKmax + 1;
            needsKReset = false;
            return true;
        }
        if (ToulBar2::showvns >= 1) {
            cout << "[Vns geode] adaptive: zone size=" << currentZoneSize
                 << " | last_aggregated=" << lastAggregatedCluster << endl;
        }
    }
    

    else if (needsKReset) {
        if (!ToulBar2::vnsAdaptive) {
            currentClusterIdx = (currentClusterIdx + 1) % (int)clusters.size();
            nbVisitedClusters++;
            clusterEntryTime = cpuTime();

            if (nbVisitedClusters >= (int)clusters.size()) {
                cycleComplete = true;
                nbVisitedClusters = 0;
            }
        }
        if (ToulBar2::showvns >= 2 || ToulBar2::verbose >= 1) {
            cout << "[Vns geode] Next cluster: cluster_idx=" << currentClusterIdx
                 << " (var=" << clusterRootWcspIdx[currentClusterIdx] << ")" << endl;
        } else if (ToulBar2::showvns == 1) {
            cout << "[Vns geode] -> cluster_idx=" << currentClusterIdx << endl;
        }
    }
    return true;
}


bool ProteinNeighborhoodChoice::shouldResetK()
{

    if (ToulBar2::vnsAdaptive) {
        needsKReset = false;
        return false;
    }
    if (ToulBar2::vnsTLimit > 0) {
        cycleComplete = false;
        bool reset = needsKReset;
        needsKReset = false;
        return reset;
    }
    if (cycleComplete) {
        cycleComplete = false;
        needsKReset = false;
        if (ToulBar2::showvns >= 1) {
            cout << "[Vns geode] Full cycle completed. Increasing LDS." << endl;
        }
        return false;
    }
    bool reset = needsKReset;
    needsKReset = false;
    return reset;
}




void ProteinNeighborhoodChoice::buildClusters(int radius)
{
    // BFS depuis chaque variable non-affectée, jusqu'à profondeur "radius".
    // CHOIX B (vecteur compact) :
    // - clusters[i] = boule géodésique du i-ème cluster non vide
    // - clusterRootWcspIdx[i] = indice WCSP de la variable racine
    clusters.clear();
    clusterRootWcspIdx.clear();

    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
        if (!wcsp->unassigned(i))
            continue;  // on saute les variables affectées

        set<int> ball;
        queue<pair<int, int>> bfsQueue;
        ball.insert(i);
        bfsQueue.push(make_pair((int)i, 0));
        set<int> neighbors;
        while (!bfsQueue.empty()) {
            int currVar = bfsQueue.front().first;
            int depth = bfsQueue.front().second;
            bfsQueue.pop();
            if (depth >= radius)
                continue;

            getDirectNeighbors(currVar, neighbors);
            for (int n : neighbors) {

                if (wcsp->unassigned(n) && ball.find(n) == ball.end()) {
                    ball.insert(n);
                    bfsQueue.push(make_pair(n, depth + 1));
                }
            }
        }

        // Ajouter au vecteur compact (pas de trous)
        clusters.push_back(vector<int>(ball.begin(), ball.end()));
        clusterRootWcspIdx.push_back((int)i);
    }

    if (ToulBar2::vnsOrderFile != "") {
        ifstream file(ToulBar2::vnsOrderFile.c_str());
        if (!file) {
            cerr << "Error: Cannot open TSP file: " << ToulBar2::vnsOrderFile << endl;
            exit(EXIT_FAILURE);
        }
        vector<vector<int>> tspClusters;
        vector<int> tspClusterRootWcspIdx;
        int varIdx;
        while (file >> varIdx) {
            auto it = std::find(clusterRootWcspIdx.begin(), clusterRootWcspIdx.end(), varIdx);
            if (it != clusterRootWcspIdx.end()) {
                int originalIndex = std::distance(clusterRootWcspIdx.begin(), it);
                tspClusters.push_back(clusters[originalIndex]);
                tspClusterRootWcspIdx.push_back(clusterRootWcspIdx[originalIndex]);
            }
        }
        clusters = tspClusters;
        clusterRootWcspIdx = tspClusterRootWcspIdx;
        if (ToulBar2::showvns >= 1) {
            cout << "[Vns geode] TSP order applied from file: " << ToulBar2::vnsOrderFile << endl;
        }
    }

    // reverse .
    // Si l'option -reverse est activée, inverser l'ordre des clusters.
    // Cela permet de tester l'impact du point de départ sur la convergence.
    if (ToulBar2::vnsReverseOrder) {
        std::reverse(clusters.begin(), clusters.end());
        std::reverse(clusterRootWcspIdx.begin(), clusterRootWcspIdx.end());
        if (ToulBar2::verbose >= 0) {
            cout << "[Vns geode] Reverse order activated: starting from last cluster (root_var="
                 << clusterRootWcspIdx[0] << ")" << endl;
        }
    }
}


void ProteinNeighborhoodChoice::getDirectNeighbors(int varIdx, set<int>& neighbors) const
{
   // Retourne l'ensemble des variables connectées à varIdx par une fonction de coût
    // binaire (= voisins directs dans le graphe primal du WCSP).
    // Pour les SCP/CPD, cela correspond aux résidus dont la distance physique
    // est inférieure au cutoff utilisé pour générer le CFN.
    neighbors.clear();
    EnumeratedVariable* var = (EnumeratedVariable*)((WCSP*)wcsp)->getVar(varIdx);
    ConstraintList* ctrlist = var->getConstrs();
    for (ConstraintList::iterator it = ctrlist->begin(); it != ctrlist->end(); ++it) {
        Constraint* ctr = (*it).constr;
        if (ctr->arity() != 2)
            continue;
        BinaryConstraint* bctr = (BinaryConstraint*)ctr;
        TSCOPE scope;
        bctr->getScope(scope);
        for (auto elt : scope) {
            if (elt.first != varIdx) {
                neighbors.insert(elt.first);
            }
        }
    }

}



// GraphNeighborhoodChoice

void GraphNeighborhoodChoice::init(WeightedCSP* wcsp_, LocalSearch* l_)
{
    this->wcsp = wcsp_;
    this->l = l_;

    // SOFT-RESET : si les clusters sont déjà construits (cas où init() est rappelée
    // depuis tb2dgvns.cpp ligne 210 après une amélioration), on évite de tout
    // reconstruire ; on réinitialise juste l'ordre de parcours.
    if (clustersBuilt) {
        file = clusters;
        reverse(file.begin(), file.end());  // back() = indice le plus petit
        return;
    }

    maxClusterSize = 0;
    minClusterSize = wcsp_->numberOfVariables();

    //  info : Choix de la source de la décomposition 
    if (ToulBar2::clusterFile != "") {
        // si l'utilisateur fournit un fichier, on l'utilise.
        if (ToulBar2::vnsGeode > 0 && ToulBar2::verbose >= 0) {
            cout << "Warning: -geode parameter ignored because a decomposition file is provided ("
                 << ToulBar2::clusterFile << ")." << endl;
        }
        load_decomposition();  // fichier fourni par l'utilisateur, méthode héritée de ClustersNeighborhoodStructure .
    } else {
        // Mode interne : construction des boules géodésiques
        if (ToulBar2::vnsGeode <= 0) {
            cerr << "Error: with -vns=11 (GRAPH mode), you must provide either"
                 << " --decfile=<file> or -geode=<radius> with radius > 0." << endl;
            throw BadConfiguration();
        }
        // info : geodesic decomposition.
        // ma classe construit bien sa propre décomposition géodésique .
        if (ToulBar2::verbose >= 0) {
            cout << "Building geodesic decomposition with radius "
                 << ToulBar2::vnsGeode << " (binary constraints only)." << endl;
        }
        buildGeodesicClusters(ToulBar2::vnsGeode); // je declenche ma decomposition qui sera stocké dans m_graph .
    }

    
    // En mode géodésique, on préserve volontairement l'intégrité des boules :
    // l'absorption fusionnerait des clusters et casserait la sémantique du rayon.

    //  Calcul des stats sur les clusters et préparation de file ----
    TCDGraph::vertex_iterator v, vend;
    tie(v, vend) = vertices(m_graph);
    for (; v != vend; ++v) {
        unsigned int sz = m_graph[*v].vars.size();
        if (sz > maxClusterSize)
            maxClusterSize = sz;
        if (sz < minClusterSize)
            minClusterSize = sz;
        clusters.push_back(*v);
    }

    // Ordre déterministe : indice croissant.
    // file est dépilée par back() => on trie en ordre décroissant.
    file = clusters;
    reverse(file.begin(), file.end());
    // info c'est ici que je prépare l'ordre de parcours des clusters pour le moteur VNS, en commençant par le premier qui est entré (celui avec l'indice le plus petit) .   
    // Heuristique interne : Naturel (sélection par indice croissant dans la zone)
    insideHeuristic = new NaturelNeighborhoodChoice();
    insideHeuristic->init(wcsp_, l_);

    clustersBuilt = true;
}

void GraphNeighborhoodChoice::buildGeodesicClusters(int radius)
{
    // BFS depuis chaque variable non-affectée, jusqu'à profondeur vnsGeode.
    // on s'est limité aux contraintes binaires, il faudra remplacer ctr->isBinary()par ctr->arity() si on veut les n-aires .
    
    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
        if (!wcsp->unassigned(i))
            continue;

        set<int> ball;
        queue<pair<int, int>> bfsQueue;  // <variable_index, depth>
        ball.insert(i);
        bfsQueue.push(make_pair(i, 0));

        while (!bfsQueue.empty()) {
            int currVarIndex = bfsQueue.front().first;
            int currentDepth = bfsQueue.front().second;
            bfsQueue.pop();

            if (currentDepth >= radius)
                continue;

            Variable* var = ((WCSP*)wcsp)->getVar(currVarIndex);
            for (ConstraintList::iterator it = var->getConstrs()->begin();
                 it != var->getConstrs()->end(); ++it) {
                Constraint* ctr = (*it).constr;
                if (!ctr->isBinary())
                    continue;
                BinaryConstraint* bctr = (BinaryConstraint*)ctr;
                int neighborIndex = bctr->getVarDiffFrom(var)->wcspIndex;
                if (wcsp->unassigned(neighborIndex)
                    && ball.find(neighborIndex) == ball.end()) {
                    ball.insert(neighborIndex);
                    bfsQueue.push(make_pair(neighborIndex, currentDepth + 1));
                }
            }
        }

        // Création du cluster (Point 2 : on garde même les singletons)
        TDCluster c = add_vertex(m_graph);
        m_graph[c].name = to_string(i);
        m_graph[c].vars = ball;
    }

    // Construction des arêtes par intersection (cf. load_decomposition lignes 81-94)
    TCDGraph::vertex_iterator v, vend, v2;
    for (tie(v, vend) = vertices(m_graph); v != vend; ++v) {
        for (v2 = v + 1; v2 != vend; ++v2) {
            set<int> separator;
            set_intersection(m_graph[*v].vars.begin(), m_graph[*v].vars.end(),
                             m_graph[*v2].vars.begin(), m_graph[*v2].vars.end(),
                             inserter(separator, separator.begin()));
            if (separator.size() > 0) {
                Cluster_edge sep;
                tie(sep, tuples::ignore) = add_edge(*v, *v2, m_graph);
                m_graph[sep].vars = separator;
                m_graph[sep].size = (float)1 / separator.size();
            }
        }
    }
}

const zone GraphNeighborhoodChoice::getNeighborhood(size_t neighborhood_size)
{
    // Identique à RandomClusterChoice::getNeighborhood, sauf que :
    // - L'ordre de parcours est déterministe (pas de shuffle).
    // - Les clusters sont des boules géodésiques.
    assert(neighborhood_size <= wcsp->numberOfUnassignedVariables());
    set<int> selclusters;
    if (file.size() == 0) {
        file = clusters;
        reverse(file.begin(), file.end());  // on inverse l'ordre pour commencer par le premier qui est entré.
    }
    assert(file.size() > 0);
    int c = file.back();
    file.pop_back();
    selclusters.insert(c);
    if (ToulBar2::verbose >= 1)
        cout << "Select geodesic cluster " << c << endl;
    zone z = m_graph[c].vars;
    // info : Select geodesic cluster .
    // on explore les clusters dans l'ordre déterministe

    // Mécanisme d'union (k-jump) : si la boule est trop petite, on agrège
    // les clusters voisins jusqu'à atteindre neighborhood_size.
    if (z.size() < neighborhood_size) {
        queue<int> fifo;
        TCDGraph::adjacency_iterator av, avend;
        int currclu = c;
        int i = 0;
        do {
            tie(av, avend) = adjacent_vertices(currclu, m_graph); // fonction boost pour retourner les voisins direct d'un sommet dans le graphe. 
            if (av != avend) {
                vector<int> neighbors(av, avend);
                // info :Tri par proximité topologique (taille du séparateur décroissante).
                // Tie-breaker par indice croissant : OBLIGATOIRE pour la reproductibilité,
                // car sort n'est pas stable et les clusters géodésiques ont
                // fréquemment des séparateurs de taille égale.
                sort(neighbors.begin(), neighbors.end(), [&](int a, int b) {
                    Cluster_edge ea, eb; // identifiant d'aretes .
                    bool found_a, found_b;
                    tie(ea, found_a) = edge(currclu, a, m_graph); // existante d'arete entre currclu et a dans m_graph ?
                    tie(eb, found_b) = edge(currclu, b, m_graph);
                    int taille_a = found_a ? m_graph[ea].vars.size() : 0;
                    int taille_b = found_b ? m_graph[eb].vars.size() : 0;
                    if (taille_a == taille_b) return a < b;  // tie-breaker déterministe
                    return taille_a > taille_b;              // proximité topologique décroissante
                });
                for (vector<int>::iterator it = neighbors.begin();
                     it != neighbors.end(); ++it) {
                    if (selclusters.count(*it) == 0) {
                        fifo.push(*it);
                        selclusters.insert(*it);
                    }
                }
            }
            if (fifo.size() == 0) {
                // Le voisinage de currclu est épuisé : sauter sur un cluster
                // non encore sélectionné.
                if (selclusters.size() >= clusters.size())
                    break;  // sécurité : tous les clusters ont été agrégés
                while (selclusters.count(i) > 0) {
                    i = (i + 1) % clusters.size();
                }
                fifo.push(clusters[i]);
                selclusters.insert(clusters[i]);
            }
            currclu = fifo.front();
            fifo.pop();
            // le z est un set,donc on aura pas de doublons.
            z.insert(m_graph[currclu].vars.begin(), m_graph[currclu].vars.end());
        } while (z.size() < neighborhood_size);
    }
    return insideHeuristic->getNeighborhood(neighborhood_size, z);
}

const zone GraphNeighborhoodChoice::getNeighborhood(size_t neighborhood_size, zone z) const
{
    assert("not implemented!!!");
    return zone();
}

bool GraphNeighborhoodChoice::incrementK()
{
    // Identique à RandomClusterChoice::incrementK : si la file est vide,
    // on signale au moteur VNS que k peut être incrémenté.
    if (file.size() == 0) {
        file = clusters;
        reverse(file.begin(), file.end());
        return true;
    }
    return false;
}
// KgeodeNeighborhoodChoice .
void KgeodeNeighborhoodChoice::init(WeightedCSP* wcsp_, LocalSearch* l_)
{
    this->wcsp = wcsp_;
    this->l = l_;
    cerr << "[KgeodeNeighborhoodChoice::init] Not implemented yet (étape 1)" << endl;
    throw BadConfiguration();
}

const zone KgeodeNeighborhoodChoice::getNeighborhood(size_t neighborhood_size)
{
    cerr << "[KgeodeNeighborhoodChoice::getNeighborhood] Not implemented yet (étape 2)" << endl;
    throw BadConfiguration();
}

const zone KgeodeNeighborhoodChoice::getNeighborhood(size_t neighborhood_size, zone z) const
{
    assert("not implemented!!!");
    return zone();
}

bool KgeodeNeighborhoodChoice::incrementK()
{
    return false;
}

void KgeodeNeighborhoodChoice::buildGeodesicClusters(int radius)
{
    cerr << "[KgeodeNeighborhoodChoice::buildGeodesicClusters] Not implemented yet (étape 1)" << endl;
}

int KgeodeNeighborhoodChoice::findClosestUnmergedCluster()
{
    return -1;
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
