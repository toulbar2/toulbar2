/*
 * tb2pdb.cpp
 *
 *  Created on: 8 november 2016
 *      Author: Charpentier Antoine
 */
 
#include "tb2pdb.hpp"

#include <sstream>  
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>

// Constructor tb2_point
tb2_point::tb2_point (
	string _label,	
	unsigned int _var_id,
	double _x,			
	double _y,			
	double _z) 
{
	label = _label;	
	var_id = _var_id;				
	x = _x;		
	y = _y;			
	z = _z;
}

// Compute the distance between this tb2_point and another tb2_point
double tb2_point::distance(
	tb2_point * pt) 
{
	double dist = 0.0;
	if (pt != NULL) {
		dist = sqrt(
			(x - pt->x)*(x - pt->x) 
			+ (y - pt->y)*(y - pt->y)
			+ (z - pt->z)*(z - pt->z)
			);
	} else {
		cout << "ERROR : in tb2pdb.cpp/tb2_point::distance the input point to compute distance is NULL" << endl;
	}	
	return dist;
}

// Compute the distance between this tb2_point and origin (0,0,0)
double tb2_point::distance() 
{
	return sqrt(x*x + y*y + z*z);
}

// Print informations on terminal
void tb2_point::print() 
{
	cout << "label: " << label << " / "
		<< "var_id: " << var_id << " / "
		<< "x: " << x << " / " 
		<< "y: " << y << " / " 
		<< "z: " << z << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// fusion sort
///////////////////////////////////////////////////////////////////////////////////////////

// PRIVATE
vector<tb2_point *> tb2_point::fusionSort(
	vector<tb2_point *> T, 
	const vector<vector<double> > * distM) 
{
	if (T.size() <= 1) {
		return T;
	} else {
		const unsigned int middle = T.size()/2;
		vector<tb2_point *> A;
		for (unsigned int i = 0; i < middle; ++i) {
			A.push_back(T[i]);
		}
		vector<tb2_point *> B;
		for (unsigned int i = middle; i < T.size(); ++i) {
			B.push_back(T[i]);
		}
		return fusion(fusionSort(A, distM), fusionSort(B, distM), distM);
	}
}

// PRIVATE
vector<tb2_point *> tb2_point::fusion(
	vector<tb2_point *> A, 
	vector<tb2_point *> B, 
	const vector<vector<double> > * distM) 
{
	if (A.size()==0) {
		return B;
	}
	if (B.size()==0) {
		return A;
	}
	
	double dA = 0.0;
	const unsigned int indexA = A[0]->var_id;
	if (indexA < var_id) {
		dA = (* distM)[var_id][indexA];	
	} else {
		dA = (* distM)[indexA][var_id];
	}
		
	double dB = 0.0;
	const unsigned int indexB = B[0]->var_id;
	if (indexB < var_id) {
		dB = (* distM)[var_id][indexB];	
	} else {
		dB = (* distM)[indexB][var_id];
	}
	
	if (dA <= dB) {
		vector<tb2_point *> supA;
		for (unsigned int i = 1; i < A.size(); ++i) {
			supA.push_back(A[i]);
		}
		vector<tb2_point *> F = fusion(supA, B, distM);
		F.insert(F.begin(), A[0]);
		return F;
	} else {
		vector<tb2_point *> supB;
		for (unsigned int i = 1; i < B.size(); ++i) {
			supB.push_back(B[i]);
		}
		vector<tb2_point *> F = fusion(A, supB, distM);
		F.insert(F.begin(), B[0]);
		return F;
	}
}

// Sort the neighbors of a vns point according to the distance matrix
void tb2_point::sortNeighbors(
	const vector<vector<double> > * distM) 
{
	if (distM != NULL) {
		// sort neighbors according to distances
		neighbors = fusionSort(neighbors, distM);
	} else {
		cout << "ERROR : in tb2pdb.cpp/tb2_point::sortNeighbors the input distance matrix to sort neighbors is NULL" << endl;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
// END fusion sort
///////////////////////////////////////////////////////////////////////////////////////////

// Constructor tb2_amino_acid
tb2_amino_acid::tb2_amino_acid (
	string _label,
	unsigned int _var_id,
	string _AA_name,			
	unsigned int _AA_id,	
	string _atom_name,
	unsigned int _atom_id,						
	double _x,			
	double _y,			
	double _z)
	: tb2_point(_label, _var_id, _x, _y, _z)
{
	AA_name = _AA_name;		
	AA_id = _AA_id;	
	atom_name = _atom_name;
	atom_id = _atom_id;				
}
	
// Print informations on terminal
void tb2_amino_acid::print() {
	cout << "label: " << label << " / " 
		<< "var_id: " << var_id << " / " 
		<< "AA_name: " << AA_name << " / "
		<< "AA_id: " << AA_id << " / "
		<< "atom_name: " << atom_name << " / "
		<< "atom_id: " << atom_id << " / "
		<< "x: " << x << " / "
		<< "y: " << y << " / "
		<< "z: " << z << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Read pdb file and store data corresponding to a specific atom name
///////////////////////////////////////////////////////////////////////////////////////////

vector<tb2_amino_acid *> * extract_pdb(
	string path_to_pdb,
	string selected_atom,
	string forbiden_AA)
{
	// Parse the forbiden amino acid string (delimiter = ':')
	vector<string> forbiden_AA_vector;
    istringstream iss(forbiden_AA);
    string token;
    while (std::getline(iss, token, ':'))
    {
        forbiden_AA_vector.push_back(token);
    }
	
	// vector returned by the function
	vector<tb2_amino_acid *> * AA_vector = NULL; 
	
	// get the pdb name
	string pdb_name = path_to_pdb.c_str();
	
	// open the file (read)
	ifstream pdb_file(path_to_pdb.c_str(), ifstream::in); 
	
	// check if the file is opened (read)
	if(pdb_file) {
		
		AA_vector = new vector<tb2_amino_acid *>;
		// vns id for str_amino_acid (read order)
		unsigned int var_id = 0;
		
		// read every lines in pdb file
		string line = ""; 
			
		while(getline(pdb_file, line)) {
			
			// check if the line is an atom
			string first_word = line.substr(0,4); // ATOM
			if (first_word == "ATOM") {
				
				// get the atom_id
				unsigned int atom_id = atoi((line.substr(6,5)).c_str());
				
				// get the atom_name
				string atom_name = line.substr(13,2);
				
				// check if the atom name is the selected atom
				if (atom_name == selected_atom) {
					
					// get the amino acid name
					string AA_name = line.substr(17,3);
				
					// skip forbiden AA
					bool is_forbiden_AA = false;
					for(unsigned int i = 0; i < forbiden_AA_vector.size(); ++i) {
						if (AA_name.compare(forbiden_AA_vector[i]) == 0) {
							is_forbiden_AA = true;
							break;
						}
					}
				
					if (!is_forbiden_AA) {
				
						// get the AA id
						unsigned int AA_id = atoi((line.substr(22,4)).c_str());
				
						// get the coord x
						double x = atof((line.substr(30,8)).c_str());
				
						// get the coord y
						double y = atof((line.substr(38,8)).c_str());
				
						// get the coord z
						double z = atof((line.substr(46,8)).c_str());
					
						tb2_amino_acid * my_amino_acid = new tb2_amino_acid(
							pdb_name,
							var_id,
							AA_name,
							AA_id,
							atom_name,
							atom_id,
							x,
							y,
							z);
					
						AA_vector->push_back(my_amino_acid);	
						//my_amino_acid->print();
						++var_id;
					}
				}
			}
		}
    } else {
		cout << "ERROR : in tb2pdb.cpp : Couldn't open pdb file " << path_to_pdb << endl;
	}
	pdb_file.close();
	
	return AA_vector;
} 

///////////////////////////////////////////////////////////////////////////////////////////
// Read xyz file and store data in vns points
///////////////////////////////////////////////////////////////////////////////////////////

vector<tb2_point *> * extract_xyz(
	string path_to_xyz)
{
	// vector returned by the function
	vector<tb2_point *> * tb2_points_vector = NULL; 
	
	// open the file (read)
	ifstream xyz_file(path_to_xyz.c_str(), ifstream::in); 
	
	// check if the file is opened (read)
	if(xyz_file) {
		
		tb2_points_vector = new vector<tb2_point *>;
		// vns id for str_amino_acid (read order)
		unsigned int var_id = 0;
		
		// read every lines in pdb file
		string line = ""; 
		while(getline(xyz_file, line)) {
			
			// store the line in a buffer
			istringstream buffer(line);
			
			// check if the line is an atom
			string label = "";
			buffer >> label;
				
			// get the coord x
			double x = 0.0;
			buffer >> x;
				
			// get the coord y
			double y = 0.0;
			buffer >> y;
				
			// get the coord z
			double z = 0.0;
			buffer >> z;
					
			tb2_point * my_point = new tb2_point(
				label,
				var_id,
				x,
				y,
				z);
					
			tb2_points_vector->push_back(my_point);	
			//my_point->print();
			++var_id;
			
		}
    } else {
		cout << "ERROR : in tb2pdb.cpp : Couldn't open xyz file " << path_to_xyz << endl;
	}
	xyz_file.close();
	
	return tb2_points_vector;
} 

///////////////////////////////////////////////////////////////////////////////////////////
// Compute the euclidean distance matrix according to a tb2_point vector
// If a centroid is given the diagonal is the distance between the point and the centroid
///////////////////////////////////////////////////////////////////////////////////////////

vector<vector<double> > * compute_distances(
	const vector<tb2_point *> * points,
	tb2_point * centroid = NULL)
{
	// the returned matrix
	vector<vector<double> > * distances_matrix = NULL;
	
	if (points != NULL) {
		distances_matrix = new vector<vector<double> >;
		// get the size of the input vector
		const unsigned int elt_size = points->size();
		distances_matrix->resize(elt_size);
		// Compute every distances (triangular matrix)
		for (unsigned int i = 0; i < elt_size; ++i) {
			(* distances_matrix)[i].resize(i+1);
			const double pi_x = (* points)[i]->x;
			const double pi_y = (* points)[i]->y;
			const double pi_z = (* points)[i]->z;
			for (unsigned int j = 0; j <= i; ++j) {
				if (i==j) {
					// the diagonal is the distance to the input centroid
					if (centroid == NULL) {
						(* distances_matrix)[i][i] = 0.0;
					} else {
						(* distances_matrix)[i][i] = centroid->distance((* points)[i]);
					}
				} else {
					const double delta_x = (pi_x - (* points)[j]->x);
					const double delta_y = (pi_y - (* points)[j]->y);
					const double delta_z = (pi_z - (* points)[j]->z);	
					(* distances_matrix)[i][j] = sqrt((delta_x)*(delta_x) 
					+ (delta_y)*(delta_y) + (delta_z)*(delta_z));
					// insert neighbors for tb2_points
					(* points)[i]->neighbors.push_back((* points)[j]);
					(* points)[j]->neighbors.push_back((* points)[i]);
				}
				//cout << (* distances_matrix)[i][j] << " ";
			}
			//cout << endl;
		}
		// sort neighbors for all points
		for (unsigned int i = 0; i < elt_size; ++i) {
			(* points)[i]->sortNeighbors(distances_matrix);
		}
	} else {
		cout << "ERROR : in tb2pdb.cpp/compute_distances the input vector to compute distances is NULL" << endl;
	}
	return distances_matrix;
}

// Same function but for tb2_amino_acid (solve direct conversion mother/daughter structures with vectors)
vector<vector<double> > * compute_distances(
	const vector<tb2_amino_acid *> * points,
	tb2_point * centroid = NULL) 
{
	// the returned matrix
	vector<vector<double> > * distances_matrix = NULL;
	if (points != NULL) {
		vector<tb2_point *> * cp_points = new vector<tb2_point *>(points->begin(), points->end());
		distances_matrix = compute_distances(cp_points, centroid);	
		// clean
		cp_points->clear();
		delete cp_points;
		cp_points = NULL;
	} else {
		cout << "ERROR : in tb2pdb.cpp/compute_distances the input vector to compute distances is NULL" << endl;
	}
	return distances_matrix;
}
	
///////////////////////////////////////////////////////////////////////////////////////////
// Compute the centroid of a tb2_point vector
///////////////////////////////////////////////////////////////////////////////////////////

tb2_point * compute_centroid(
	const vector<tb2_point *> * points)
{
	// the returned centroid
	tb2_point * centro = NULL; 
	if (points != NULL) {
		centro = new tb2_point("centroid", 0, 0.0, 0.0, 0.0);
		// get the size of the input vector
		const unsigned int elt_size = points->size();
		if (elt_size > 0) { 
			for (unsigned int i = 0; i < elt_size; ++i) {
				centro->x += (* points)[i]->x;
				centro->y += (* points)[i]->y;
				centro->z += (* points)[i]->z;
			}
			centro->x /= (double) elt_size;
			centro->y /= (double) elt_size;
			centro->z /= (double) elt_size;
		}
	} else {
		cout << "ERROR : in tb2pdb.cpp/compute_centroid the input vector to compute centroid is NULL" << endl;
	}
	return centro;
}

// Same function but for tb2_amino_acid (solve direct conversion mother/daughter structures with vectors)
tb2_point * compute_centroid(
	const vector<tb2_amino_acid *> * points) 
{
	// the returned centroid
	tb2_point * centro = NULL; 
	if (points != NULL) {
		vector<tb2_point *> * cp_points = new vector<tb2_point *>(points->begin(), points->end());
		centro = compute_centroid(cp_points);
		// clean
		cp_points->clear();
		delete cp_points;
		cp_points = NULL;
	} else {
		cout << "ERROR : in tb2pdb.cpp/compute_centroid the input vector to compute centroid is NULL" << endl;
	}
	return centro;
}

// Analyse neighborhood
// return a vector with : 0=average_dist_to_centro_neighb, 1=var_dist_to_centro_neighb,
// 2=stdev_dist_to_centro_neighb, 3=coef_variation, 4=dispersion, 5=distorsion
vector<double> analyze_points(
	vector<tb2_point *> * neighborhood, 
	tb2_point * neighborhood_centroid) 
{
	// the returned vector
	vector<double> analyze;
	
	if ((neighborhood != NULL) && (neighborhood_centroid != NULL)) {
		//
		const unsigned int neighborhood_size = neighborhood->size();
	
		if(neighborhood_size>0) {

			//
			analyze.resize(6);
			
			// Store distances between the neighborhood's center and points in neighborhood
			vector<double> dist_to_centro_neighb;
			// Compute average distances 
			double average_dist_to_centro_neighb = 0.0;
			for (vector<tb2_point *>::iterator it=neighborhood->begin(); it!=neighborhood->end(); ++it) {
				const double dist = (* it)->distance(neighborhood_centroid);
				average_dist_to_centro_neighb += dist;
				dist_to_centro_neighb.push_back(dist);
			}
			average_dist_to_centro_neighb /= (double) neighborhood_size;
			analyze[0]=average_dist_to_centro_neighb;
	
			// Compute variance
			double var_dist_to_centro_neighb = 0.0;
			for (unsigned int i = 0; i < dist_to_centro_neighb.size(); ++i) {
				const double delta_dist = dist_to_centro_neighb[i] - average_dist_to_centro_neighb;
				var_dist_to_centro_neighb += delta_dist * delta_dist;
			}
			var_dist_to_centro_neighb /= (double) neighborhood_size;
			analyze[1] = var_dist_to_centro_neighb;
	
			// Compute standard deviation
			const double stdev_dist_to_centro_neighb = sqrt(var_dist_to_centro_neighb);
			analyze[2] = stdev_dist_to_centro_neighb;
	
			// Compute the coefficient of variation
			const double coef_variation = stdev_dist_to_centro_neighb / average_dist_to_centro_neighb;
			analyze[3] = coef_variation;
	
			// Compute the dispersion : dmax - dmin
			const double dispersion = *std::max_element(dist_to_centro_neighb.begin(), dist_to_centro_neighb.end()) 
										- *std::min_element(dist_to_centro_neighb.begin(), dist_to_centro_neighb.end());
			analyze[4] = dispersion;								
										
			// Compute the distorsion 							
			double distorsion = 0.0;
			for (unsigned int i = 0; i < dist_to_centro_neighb.size(); ++i) {
				distorsion += dist_to_centro_neighb[i] * dist_to_centro_neighb[i];
			}
			analyze[5] = distorsion;
		}
	} else {
		cout << "ERROR : in tb2pdb.cpp/analyze_points the input vector or centroid to compute analyze is NULL" << endl;
	}
	
	return analyze;
}

// Same function but for tb2_amino_acid (solve direct conversion mother/daughter structures with vectors)
vector<double> analyze_points(
	vector<tb2_amino_acid *> * neighborhood, 
	tb2_point * neighborhood_centroid) 
{
	// the returned centroid
	vector<double> analyze; 
	if (neighborhood != NULL) {
		vector<tb2_point *> * cp_neighborhood = new vector<tb2_point *>(neighborhood->begin(), neighborhood->end());
		analyze = analyze_points(cp_neighborhood, neighborhood_centroid);
		// clean
		cp_neighborhood->clear();
		delete cp_neighborhood;
		cp_neighborhood = NULL;
	} else {
		cout << "ERROR : in tb2pdb.cpp/analyze_points the input vector to compute analyze is NULL" << endl;
	}
	return analyze;
}
