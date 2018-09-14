/*
 * tb2pdb.hpp
 *
 *  Created on: 8 november 2016
 *      Author: Charpentier Antoine
 */
 
 
#ifndef SRC_TB2PDB_HPP_
#define SRC_TB2PDB_HPP_
 
#include <string>
#include <vector>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////
// Structure used to store pdb and vns neighborhood informations (ids)
///////////////////////////////////////////////////////////////////////////////////////////

// tb2_point : label / var_id / x / y / z
struct tb2_point {
	// point name
	string label;	
	// variable index [0,n-1]
	unsigned int var_id;
	// coordinates of the point
	double x;			
	double y;			
	double z;
	//
	vector<tb2_point *> neighbors; 
	
	// Constructor 
	tb2_point (
		string _label,	
		unsigned int _var_id,
		double _x,			
		double _y,			
		double _z);

    // Compute the distance between this tb2_point and another tb2_point
	double distance(
		tb2_point * pt);
	// Compute the distance between this tb2_point and origin (0,0,0)
	double distance();

	//
	void sortNeighbors(
		const vector<vector<double> > * distM);
	
	// Print informations on terminal
	void print();
	
	virtual ~tb2_point() {};
	
private:
	//
	vector<tb2_point *> fusionSort(
		vector<tb2_point *> T, 
		const vector<vector<double> > * distM); 
	//
	vector<tb2_point *> fusion(
		vector<tb2_point *> A, 
		vector<tb2_point *> B, 
		const vector<vector<double> > * distM); 
};
typedef struct tb2_point tb2_point;

// tb2_amino_acid : as pdb data
struct tb2_amino_acid : tb2_point {
	// Amino acid data (PDB informations)
	string AA_name;			
	unsigned int AA_id;	
	// Atom data (PDB informations)
	string atom_name;
	unsigned int atom_id;		

	// Constructor 
	tb2_amino_acid (
		string _label,
		unsigned int _var_id,
		string _AA_name,			
		unsigned int _AA_id,	
		string _atom_name,
		unsigned int _atom_id,						
		double _x,			
		double _y,			
		double _z);
	
	// Print informations on terminal
	void print();
	
	virtual ~tb2_amino_acid() {};
};
typedef struct tb2_amino_acid tb2_amino_acid;

// Read xyz file and store data in vns points
vector<tb2_point *> * extract_xyz(
	string path_to_xyz);

// Read pdb file and store data corresponding to a specific atom name
vector<tb2_amino_acid *> * extract_pdb(
	string path_to_pdb,
	string selected_atom,
	string forbiden_AA);
	
// Compute the euclidean distance matrix according to a tb2_point vector
// If a centroid is given the diagonal is the distance between the point and the centroid
vector<vector<double> > * compute_distances(
	const vector<tb2_point *> * points,
	tb2_point * centroid);
	
vector<vector<double> > * compute_distances(
	const vector<tb2_amino_acid *> * points,
	tb2_point * centroid);

// Compute the centroid of a tb2_point vector
tb2_point * compute_centroid(
	const vector<tb2_point *> * points);
	
tb2_point * compute_centroid(
	const vector<tb2_amino_acid *> * points);
	
// Analyse neighborhood
// return a vector with : 0=average_dist_to_centro_neighb, 1=var_dist_to_centro_neighb,
// 2=stdev_dist_to_centro_neighb, 3=coef_variation, 4=dispersion, 5=distorsion
vector<double> analyze_points(
	vector<tb2_point *> * neighborhood, 
	tb2_point * neighborhood_centroid);
	
vector<double> analyze_points(
	vector<tb2_amino_acid *> * neighborhood, 
	tb2_point * neighborhood_centroid);
	
#endif /* SRC_TB2PDB_HPP_ */
