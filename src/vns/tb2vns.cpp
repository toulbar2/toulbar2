/*
 * tb2vns.cpp
 *
 *  Created on: 3 mars 2015
 *      Authors: Mathieu Fontaine, Abdelkader Ouali
 *      Phd. Student : LITIO, University of Oran. GREYC, University of Caen.
 */

#include "tb2vns.hpp"
#include "tb2wcsp.hpp"
#include "tb2pdb.hpp" 
#include "mpi.h"

// loading decomposition for only a given file
void ClustersNeighborhoodStructure::load_decomposition()
{
    if (ToulBar2::clusterFile != "") {
        bool cover = false;
        if (ToulBar2::clusterFile.find(".cov") != string::npos) cover =true;
        fstream file(ToulBar2::clusterFile.c_str());
        //file >> clustersNumber;
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
                do {
                    unsigned int var;
                    ss >> var;
                    nbvars.insert(var);
                    if (var >= wcsp->numberOfVariables()) {
                        cerr << "Error: cluster decomposition contains bad variable index!" << endl;
                        exit(EXIT_FAILURE);
                    }
                    tmp.insert(var);
                } while (ss.good());
                m_graph[c].vars = tmp;
            }
        }
        file.close();
        if (nbvars.size() != wcsp->numberOfVariables()) {
            cerr << "Error: cluster decomposition has missing variables! (" << nbvars.size() << "!=" << wcsp->numberOfVariables() << ")" << endl;
            exit(EXIT_FAILURE);
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
                    m_graph[sep].size = (float) 1 / separator.size();
                }
            }
        }
    } else {
        cerr << "Cluster decomposition file is missing!" << endl;
        exit(EXIT_FAILURE);
    }
}

void RandomNeighborhoodChoice::init(WeightedCSP* wcsp_, LocalSearch* l_)
{
    this->l = l_;
    wcsp = wcsp_;
	//
	last_first_selection = 0;
	// check if the euclidean distances are availabled
	if(!ToulBar2::check_edrsn) {
		ToulBar2::random_selection_neighborhood = 0;
		ToulBar2::random_selection_first_element = 0;
		cout << "Warning : check edrsn fail -> use default algorithm f0 n0." << endl;
	}
	// centre/periphery for the fisrt selection doesn't work
	// in the case of a new neighborhood build on the previous one
	if(ToulBar2::incr_k_algo == 2 && ToulBar2::random_selection_first_element == 1) {
		ToulBar2::random_selection_first_element = 0;
		cout << "Warning : root algo = 1 is incompatible with kpp algo = 2 -> use default root algo" << endl;
	}
	//
	int myrank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if(ToulBar2::seed < 1) {
		srand(time(0) + getpid());
	} else {
		srand(ToulBar2::seed + myrank);
	}	
}



const zone RandomNeighborhoodChoice::getNeighborhood(size_t neighborhood_size)
{
	const zone neighborhood = getNeighborhood(neighborhood_size, zone ());
	return neighborhood;
}
const zone RandomNeighborhoodChoice::getNeighborhood(size_t neighborhood_size, zone z)
{
	
	// The returned neighborhood		
	zone neighborhood;
		
	// Select an algorithm to build the neighborhood
	switch (ToulBar2::random_selection_neighborhood)
	{
        case 2:
			// k nearest neighbors algorithm
			neighborhood = kNearest(neighborhood_size, z);
			break;
            
        case 1:
			// (1/d) euclidean distance probabilites algorithm
			neighborhood = probaDistance(neighborhood_size, z);
			break;
            
        case 0:
			// shuffle uniformly and randomly (default)
			neighborhood = directShuffle(neighborhood_size, z);
			break;
            
        default:
			// shuffle uniformly and randomly (default)
			neighborhood = directShuffle(neighborhood_size, z);
	}
	
	// Update last zone
	last_zone = neighborhood;
	
	return neighborhood;
}

// Select an element (index in the zone z or the direct variable id) as root for the neighborhood
unsigned int RandomNeighborhoodChoice::getFirstElt(
	size_t neighborhood_size,
	zone z) {

	// the returned first element
	unsigned int first_element = 0;
	
	// Select an algorithm to choose the first element
	switch (ToulBar2::random_selection_first_element)
	{
		case 1:
			// diversification and intensification
			if (neighborhood_size == (unsigned int) ToulBar2::vnsKmin) {
				// intensification : get the first elt near the molecule's centroid
				first_element = probaFirstChoice(z);
			} else {
				// diversification : get the first elt away from the last first point 
				first_element = diversifyFirstChoice(z);
			}
			break;
		
		case 0:
			// full random
			first_element = randFirstChoice(z);
			break;
		
		default:
			// full random
			first_element = randFirstChoice(z);
	}

	return first_element;
}

// Choose an element in z uniformly and randomly 
// return the chosen position in z not the corresponding value (return i but not z[i])
unsigned int RandomNeighborhoodChoice::randFirstChoice(zone z) {
	// Check the input zone
	const unsigned int size_zone = z.size();
	unsigned int select_nb_var = ToulBar2::nbvar; // if z is empty all the variables is considered 
	if(size_zone > 0) {
		select_nb_var = size_zone;
	}
	// Select the first point uniformly and randomly between all points (z or full)
	const unsigned int first_elt = rand() % select_nb_var;
	
	return first_elt;
}

// Choose an element in z uniformly and randomly according to distances between elements in z and the molecule's centroid
// return the chosen position in z not the corresponding value (return i but not z[i])
unsigned int RandomNeighborhoodChoice::probaFirstChoice(zone z) {
	// Copy the input zone in a vector
	vector<int> zv(z.begin(), z.end());
	// Check the input zone
	const unsigned int size_zone = z.size();
	unsigned int select_nb_var = ToulBar2::nbvar; // if z is empty all the variables is considered 
	if(size_zone > 0) {
		select_nb_var = size_zone;
	}
	// 1. Compute probabilities
	// Vector to store probabilities 1/d
	vector<double> proba(select_nb_var, 0.0);
	// 1. Compute probabilites for all poitns (= 1/d)
	for (unsigned int i = 0; i < select_nb_var; ++i) {
		unsigned int true_i = i;
		// Get the corresponding index in ToulBar2::AA_vector
		if(size_zone > 0) {
			true_i = zv[i];
		}
		// the distances from the molecule's centroid are in the diagonal of the distance matrix
		proba[i] = 1.0 / (0.001 + (* ToulBar2::distance_matrix)[true_i][true_i]);	
	}	
		
	// 2. Select an elt according to the probabilities
	const unsigned int first_elt = getIndexWithProba(proba);
	
	return first_elt;
}

// Choose an element in z uniformly and randomly according to distances between elements in z and the molecule's centroid
// return the chosen position in z not the corresponding value (return i but not z[i])
unsigned int RandomNeighborhoodChoice::diversifyFirstChoice(zone z) {
	// Copy the input zone in a vector
	vector<int> zv(z.begin(), z.end());
	// Check the input zone
	const unsigned int size_zone = z.size();
	unsigned int select_nb_var = ToulBar2::nbvar; // if z is empty all the variables is considered 
	if(size_zone > 0) {
		select_nb_var = size_zone;
	}
	// 1. Compute probabilities
	// Vector to store probabilities = d (diversification)
	vector<double> proba(select_nb_var, 0.0);
	// 1. Compute probabilites for all poitns (=d)
	for (unsigned int i = 0; i < select_nb_var; ++i) {
		double dist = 0.0;
		unsigned int true_i = i;
		// Get the corresponding index in ToulBar2::AA_vector
		if(size_zone > 0) {
			true_i = zv[i];
		}
		if (true_i < last_first_selection) {
			dist = (* ToulBar2::distance_matrix)[last_first_selection][true_i];
		} else {
			dist = (* ToulBar2::distance_matrix)[true_i][last_first_selection];
		}
		
		// the distances from the molecule's centroid are in the diagonal of the distance matrix
		proba[i] = dist;	
	}	
		
	// 2. Select an elt according to the probabilities
	const unsigned int first_elt = getIndexWithProba(proba);
	
	return first_elt;
}

// Choose an index in [0, proba.size()-1] according to probabilities in the input vector
unsigned int RandomNeighborhoodChoice::getIndexWithProba(vector<double> proba) {

	// the returned index
	unsigned int selected_elt = 0;
	// 1. Compute the sum of probabilities
	double sum_proba = 0.0;
	for (unsigned int i = 0; i < proba.size(); ++i) {
		sum_proba += proba[i];
	}
	// Select an elt according to a random and uniform value in [0, sum_proba] 
	// A point is selected if random value is in [cur_sum_proba - proba[i], cur_sum_proba]
	const double random_value = sum_proba * (rand() / (double) RAND_MAX); 
	double current_sum_proba = 0.0;
	for (unsigned int i = 0; i < proba.size(); ++i) {
		current_sum_proba += proba[i];
		if (random_value < current_sum_proba) {
			selected_elt = i;
			break;
		}
	}

	return selected_elt;
}

// Select randomly and uniformly points in the input zone z
const zone RandomNeighborhoodChoice::directShuffle(size_t neighborhood_size, zone z) {
	
	// The returned neighborhood		
	zone neighborhood;
	// Use the last neighborhoods to build a new one if the algorithm using memory is selected
	if ((ToulBar2::incr_k_algo == 2)
			&& (neighborhood_size > (unsigned int) ToulBar2::vnsKmin)) {
		neighborhood = last_zone;
	}
	
	//
	const unsigned int size_zone = z.size();
	// Copy the input zone in a vector
	vector<int> zv(z.begin(), z.end());
	if (size_zone == 0) { // if z is empty all the variables is considered 
		zv.resize(ToulBar2::nbvar);
		for (unsigned int i = 0; i < wcsp->numberOfVariables(); ++i) {
			zv[i] = i;
		}
	}
	
	// Get the first point for the neighborhood if empty 
	unsigned int true_first_elt = 0;
	if(neighborhood.size() == 0) {
		// Get the first point for the neighborhood
		const unsigned int first_elt = getFirstElt(neighborhood_size, z);		
		// Get the corresponding index in ToulBar2::AA_vector
		true_first_elt = first_elt;
		if (size_zone > 0) {
			true_first_elt = zv[first_elt];
		}
		//
		neighborhood.insert(true_first_elt);
	} else {
		true_first_elt = last_first_selection;
	}
	
	// shuffle variables
	random_shuffle(zv.begin(), zv.end());
	
	// add variables in the set until the requested size is reached
	for (unsigned int i = 0; i < zv.size(); ++i) {
		if(neighborhood_size > neighborhood.size()) {
			neighborhood.insert(zv[i]);
		} else {
			break;
		}
	}

	// print data about this neighborhood
	printNeighborhood(neighborhood, true_first_elt);
	// update variable
	last_first_selection = true_first_elt;
	
	return neighborhood;	
}

// Select the first point uniformly and randomly in z, 
// other points are chosen considering their distances from this first selected point
const zone RandomNeighborhoodChoice::probaDistance(size_t neighborhood_size, zone z) {
	
	// The returned neighborhood		
	zone neighborhood;
	// Use the last neighborhoods to build a new one if the algorithm using memory is selected
	if ((ToulBar2::incr_k_algo == 2)
			&& (neighborhood_size > (unsigned int) ToulBar2::vnsKmin)) {
		neighborhood = last_zone;
	}
            
	// Copy the input zone in a vector
	vector<int> zv(z.begin(), z.end());
	// Check the input zone
	const unsigned int size_zone = z.size();
	unsigned int select_nb_var = ToulBar2::nbvar; // if z is empty all the variables is considered 
	if(size_zone > 0) {
		select_nb_var = size_zone;
	}
		
	// Get the first point for the neighborhood if empty 
	unsigned int true_first_elt = 0;
	if(neighborhood.size() == 0) {
		// Get the first point for the neighborhood
		const unsigned int first_elt = getFirstElt(neighborhood_size, z);		
		// Get the corresponding index in ToulBar2::AA_vector
		true_first_elt = first_elt;
		if (size_zone > 0) {
			true_first_elt = zv[first_elt];
		}
		//
		neighborhood.insert(true_first_elt);
	} else {
		true_first_elt = last_first_selection;
	}
	
	// Selected other points according to 1/d		 
	if (neighborhood_size > neighborhood.size()) {
		// Vector to identify the selected points (true : not selected, false : selected)
		vector<bool> not_selected(select_nb_var, true);
		// check which variables are already selected
		for(unsigned int i=0; i<select_nb_var; ++i) {
			unsigned int true_i = i;
			if(size_zone > 0) {
				true_i = zv[i];
			}
			if (neighborhood.find(true_i) !=  neighborhood.end()){
				not_selected[i]=false;
			}
		}
		
		// Vector to store probabilities 1/d
		vector<double> proba(select_nb_var, 0.0);
		// 1. Compute probabilites for all poitns (= 1/d), =0 if already selected
		for (unsigned int i = 0; i < select_nb_var; ++i) {
			if (not_selected[i]) {
				double dist = 0.0;
				unsigned int true_i = i;
				// Get the corresponding index in ToulBar2::AA_vector
				if(size_zone > 0) {
					true_i = zv[i];
				}
				if (true_i < true_first_elt) {
					dist = (* ToulBar2::distance_matrix)[true_first_elt][true_i];
				} else {
					dist = (* ToulBar2::distance_matrix)[true_i][true_first_elt];
				}
				proba[i] = 1.0 / (0.001 + dist);
			}
		}
       
		// number of points to add in the neighborhood	
        const unsigned int delta_elt = neighborhood_size-neighborhood.size();
		// 2. Select k-1 elements 
		for (unsigned int i = 0; i < delta_elt; ++i) {
			// Choose an index according to probabilities				
			const unsigned int proba_index = getIndexWithProba(proba);
			// update data	
			not_selected[proba_index] = false;
			proba[proba_index] = 0.0;
			// Get the corresponding index in ToulBar2::AA_vector
			if(size_zone > 0) {
				neighborhood.insert(zv[proba_index]);
			} else {
				neighborhood.insert(proba_index);
			}
		}
	}
	// print data about this neighborhood
	printNeighborhood(neighborhood, true_first_elt);
	// update variable
	last_first_selection = true_first_elt;
			
	return	neighborhood;
}

// Select the first point uniformly and randomly in z, 
// other points are the k-1 nearest neighbors of the first selected point
const zone RandomNeighborhoodChoice::kNearest(size_t neighborhood_size, zone z) {
	// to do ?: select randomly and uniformly last neighbors if the distances are equals ?
	
	// The returned neighborhood		
	zone neighborhood;
	// Use the last neighborhoods to build a new one if the algorithm using memory is selected
	if ((ToulBar2::incr_k_algo == 2)
			&& (neighborhood_size > (unsigned int) ToulBar2::vnsKmin)) {
		neighborhood = last_zone;
	}
	
	// Check the input zone
	const unsigned int size_zone = z.size();
	// Copy the input zone in a vector
	vector<int> zv(z.begin(), z.end());
	
	// Get the first point for the neighborhood if empty 
	unsigned int true_first_elt = 0;
	if(neighborhood.size() == 0) {
		// Get the first point for the neighborhood			
		const unsigned int first_elt = getFirstElt(neighborhood_size, z);
		// Get the corresponding index in ToulBar2::AA_vector
		true_first_elt = first_elt;
		if (size_zone > 0) {
			true_first_elt = zv[first_elt];
		}
		//
		neighborhood.insert(true_first_elt);
	} else {
		true_first_elt = last_first_selection;
	}
		
	// Select the k-1 nearest neighbors of the first element
	if (neighborhood_size > neighborhood.size()) {
		vector<tb2_point *> & neighbors = (* ToulBar2::AA_vector)[true_first_elt]->neighbors;
		for(unsigned int i = 0; i < neighbors.size(); ++i) {
			if(neighborhood_size > neighborhood.size()) {
				const unsigned int new_elt = neighbors[i]->var_id;
				neighborhood.insert(new_elt);
			} else {
				break;
			}
		}
	}
	// print data about this neighborhood
	printNeighborhood(neighborhood, true_first_elt);
	
	// update variable
	last_first_selection = true_first_elt;
	
	return neighborhood;
}

// Print information about the selected neighborhood z
void RandomNeighborhoodChoice::printNeighborhood(zone z, unsigned int first_selected_elt) {
		
	// level 0
	if (ToulBar2::vns_verbosity_level > 0) {
			
		const unsigned int neighborhood_size = z.size();
		if (neighborhood_size > 0) {
			int myrank = 0;
			MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
			cout << "Thread = " << myrank 
				<< " ; t(s) = " << time(0)                                 
				<< " ; algo(KFN) = " 
				<< ToulBar2::incr_k_algo
				<< ToulBar2::random_selection_first_element
				<< ToulBar2::random_selection_neighborhood
				<< " ; k = " << neighborhood_size;
		
			// level 1
			if (ToulBar2::vns_verbosity_level > 1) {
		
				if (ToulBar2::check_edrsn) {
					// Get the amino acid from the input zone
					vector<tb2_amino_acid *> * AA_neighborhood = new vector<tb2_amino_acid *>;
					for (zone::iterator it=z.begin(); it!=z.end(); ++it) {
						const int index_elt = (* it);
						AA_neighborhood->push_back((* ToulBar2::AA_vector)[index_elt]);	
					} 
					// Compute neighborhood's centroid (NC)
					tb2_point * centroid_neighborhood = compute_centroid(AA_neighborhood);
					// Get the first selected amino acid (F) 
					tb2_amino_acid * first_AA = ((* ToulBar2::AA_vector)[first_selected_elt]);
		
					//
					const double distFNC =  first_AA->distance(centroid_neighborhood);
					//
					const double distFC =  first_AA->distance(ToulBar2::AA_centroid);
					// Compute distance between the molecule's centroid (C) and neighborhood's centroid (NC)
					const double distNCC = ToulBar2::AA_centroid->distance(centroid_neighborhood);
		
					cout << " ; 1erAAid = id:" << first_AA->AA_id
						<< " x:" << first_AA->x 
						<< " y:" << first_AA->y 
						<< " z:" << first_AA->z
						<< " distToNC:" << distFNC
						<< " distToC:" << distFC
						<< " ; newCentro = x:" << centroid_neighborhood->x 
						<< " y:" << centroid_neighborhood->y 
						<< " z:" << centroid_neighborhood->z
						<< " distToC:" << distNCC;
			
					// level 2
					if (ToulBar2::vns_verbosity_level > 2) {
						// 0=average_dist_to_centro_neighb, 1=var_dist_to_centro_neighb,
						// 2=stdev_dist_to_centro_neighb, 3=coef_variation, 4=dispersion, 5=distorsion
						vector<double> analyze = analyze_points(AA_neighborhood, centroid_neighborhood);
		
						cout << " ; aver = " << analyze[0]
							<< " ; var = " << analyze[1]
							<< " ; stdev = " << analyze[2]
							<< " ; coefVar = " << analyze[3]
							<< " ; disper = " << analyze[4]
							<< " ; distor = " << analyze[5]
							<< " ; neighbAAid =";
						for (vector<tb2_amino_acid *>::iterator it=AA_neighborhood->begin(); it!=AA_neighborhood->end(); ++it) {
							cout << " " << (* it)->AA_id;
						}
					}
			
					// clean vns_point_neighborhood
					AA_neighborhood->clear();
					delete AA_neighborhood;
					AA_neighborhood = NULL;
					// clean centroid_neighborhood
					delete centroid_neighborhood;
					centroid_neighborhood = NULL;
				}
			}
			cout << endl; 
		} 
	}
}

void RandomClusterChoice::init(WeightedCSP* wcsp_, LocalSearch* l_)
{

    this->l = l_;
    this->wcsp = wcsp_;
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
        set<Constraint*> csts = m_graph[*v].consts;
        m_graph[*v].lastCost = 0;
        for (set<Constraint*>::iterator it = csts.begin(); it != csts.end();
                ++it) {
            m_graph[*v].lastCost = m_graph[*v].lastCost
                                   + (*it)->getConflictWeight();
        }
    }
    if (clusters.size()>=1 && m_graph[clusters[clusters.size() - 1]].vars.empty()) {
        clusters.pop_back();
    }
    file = clusters;
    random_shuffle(file.begin(), file.end());
    insideHeuristic = new RandomNeighborhoodChoice();
    precK = -1;
    insideHeuristic->init(wcsp, l);
}

const zone RandomClusterChoice::getNeighborhood(size_t neighborhood_size)
{
    precK = neighborhood_size;
    set<int> selclusters;
    if (file.size() == 0) {
        file = clusters;
        random_shuffle(file.begin(), file.end());
    }
    int c = file.back();
    file.pop_back();
    selclusters.insert(c);
    if (ToulBar2::verbose >= 1) cout << "Select cluster " << c << endl;
    zone z = m_graph[c].vars;
    //if variables are missing
    assert(neighborhood_size <= wcsp->numberOfVariables());
    if (z.size() < neighborhood_size) {
        queue<int> fifo;
        TCDGraph::adjacency_iterator v, vend;
        int currclu = c;
        int i = 0;
        do {
            tie(v, vend) = adjacent_vertices(currclu, m_graph);
            if (v != vend) {
                vector<int> neighbors(v, vend);
                //merge neighbors of current cluster
                random_shuffle(neighbors.begin(), neighbors.end());
                //add them to list
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
                    i = (i + 1) % clusters.size(); //TODO: use file instead in order to shuffle clusters randomly
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

const zone RandomClusterChoice::getNeighborhood(size_t neighborhood_size, zone z)
{
    assert("not implemented!!!");
    return zone();
}

const bool RandomClusterChoice::incrementK()
{
    if (file.size() == 0) {
        file = clusters;
        random_shuffle(file.begin(), file.end());
        return true;
    }

    return true;
}

void ParallelRandomClusterChoice::init(WeightedCSP* wcsp_, LocalSearch* l_)
{
    this->l = l_;
    this->wcsp = wcsp_;
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
        set<Constraint*> csts = m_graph[*v].consts;
        m_graph[*v].lastCost = 0;
        for (set<Constraint*>::iterator it = csts.begin(); it != csts.end();
                ++it) {
            m_graph[*v].lastCost = m_graph[*v].lastCost
                                   + (*it)->getConflictWeight();
        }
    }
    if (clusters.size()>=1 && m_graph[clusters[clusters.size() - 1]].vars.empty()) {
        clusters.pop_back();
    }
    file = clusters;
    random_shuffle(file.begin(), file.end());
    insideHeuristic = new RandomNeighborhoodChoice();
    precK = -1;
    insideHeuristic->init(wcsp, l);
}

const zone ParallelRandomClusterChoice::getNeighborhood(size_t neighborhood_size)
{
    precK = neighborhood_size;
    set<int> selclusters;
    int c = file.back();
    file.pop_back();
    selclusters.insert(c);
    zone z = m_graph[c].vars;
    //if variables are missing
    if (z.size() < neighborhood_size) {
        queue<int> fifo;
        TCDGraph::adjacency_iterator v, vend;
        int currclu = c;
        int i = 0;
        do {
            tie(v, vend) = adjacent_vertices(currclu, m_graph);
            if (v != vend) {
                vector<int> neighbors(v, vend);
                //merge neighbors of current cluster
                random_shuffle(neighbors.begin(), neighbors.end());
                //add them to list
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
                    i = (i + 1) % clusters.size(); //TODO: use file instead in order to shuffle clusters randomly
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

const zone ParallelRandomClusterChoice::getNeighborhood(size_t neighborhood_size, zone z)
{
    assert("Not implemented!!!");
    return zone();
}

const zone ParallelRandomClusterChoice::SlaveGetNeighborhood(unsigned int CurrentCluster, size_t neighborhood_size)
{
    precK = neighborhood_size;
    set<int> selclusters;
    selclusters.insert(CurrentCluster);
    zone z = m_graph[CurrentCluster].vars;
    //if variables are missing
    if (z.size() < neighborhood_size) {
        queue<int> fifo;
        TCDGraph::adjacency_iterator v, vend;
        int currclu = CurrentCluster;
        int i = 0;
        do {
            tie(v, vend) = adjacent_vertices(currclu, m_graph);
            if (v != vend) {
                vector<int> neighbors(v, vend);
                //merge neighbors of current cluster
                random_shuffle(neighbors.begin(), neighbors.end());
                //add them to list
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
                    i = (i + 1) % clusters.size(); //TODO: use file instead in order to shuffle clusters randomly
                }
                fifo.push(clusters[i]);
                selclusters.insert(clusters[i]);
            }
            currclu = fifo.front();
            fifo.pop();
            z.insert(m_graph[currclu].vars.begin(),m_graph[currclu].vars.end());
        } while (z.size() < neighborhood_size);

    }
    return insideHeuristic->getNeighborhood(neighborhood_size, z);
}

const zone ParallelRandomClusterChoice::SlaveGetNeighborhood(unsigned int CurrentCluster, uint number, size_t NeighborhoodSize)
{
    precK = NeighborhoodSize;
    set<int> selclusters;
    selclusters.insert(CurrentCluster);
    zone z = m_graph[CurrentCluster].vars;
    //if variables are missing
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
                //merge neighbors of current cluster
                random_shuffle(neighbors.begin(), neighbors.end());
                //add them to list
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
                    i = (i + 1) % clusters.size(); //TODO: use file instead in order to shuffle clusters randomly
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

const bool ParallelRandomClusterChoice::incrementK()
{
    if (file.size() == 0) {
        file = clusters;
        random_shuffle(file.begin(), file.end());
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
                //merge neighbors of current cluster
                random_shuffle(neighbors.begin(), neighbors.end());
                //add them to list
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
                    i = (i + 1) % clusters.size(); //TODO: use file instead in order to shuffle clusters randomly
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

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
