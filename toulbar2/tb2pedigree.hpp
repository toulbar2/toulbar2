/** \file tb2pedigree.hpp
 *  \brief Pedigree data structure
 * 
 */
 
#ifndef TB2PEDIGREE_HPP_
#define TB2PEDIGREE_HPP_

#include "tb2wcsp.hpp"

#include <map>

class Genotype {
public:
  int allele1;
  int allele2;
};

typedef enum {MALE=1, FEMALE=2} Sex;

class Individual {
public:
  int individual;
  int varindex;
  int father;
  int mother;
  int sex;
  Genotype genotype;
  bool typed; // true if one of its descendant children (or itself) is typed
  
  Individual(int ind);

  void print(ostream& os);
};

class Pedigree {
  int locus;    // same locus for all the genotypes
  vector<Individual> pedigree;  // list of individuals
  vector<int> genotypes;    // list of genotyped individuals id.
  vector<Genotype> genoconvert; // convert domain value to genotype
  map<int, int> individuals;    // sorted list of pair<individual id, pedigree id>
  map<int, int> alleles;    // sorted list of pair<allele number, encoding consecutive number>  
  int nbtyped;  // number of individuals with a genotyped descendant
  
  void typeAscendants(int individual);

public:
  Pedigree() : locus(-1), nbtyped(0) {alleles[0] = 0;}
  
  void read(const char *fileName, WCSP *wcsp);
  void save(const char *fileName, WCSP *wcsp, bool corrected);

  void printCorrection(WCSP *wcsp);
  
  void printGenotype(ostream& os, Value value);
};

#endif /*TB2PEDIGREE_HPP_*/
