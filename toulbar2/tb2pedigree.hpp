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
};

class Pedigree {
  vector<Individual> pedigree;
  vector<int> genotypes;
  vector<Genotype> genoconvert;
  map<int, int> individuals;
  int nbtyped;
  
  void typeAscendants(int individual);

public:
  Pedigree() : nbtyped(0) {}
  
  void readPedigree(const char *fileName, WCSP *wcsp);

  void printCorrection(WCSP *wcsp);
};

#endif /*TB2PEDIGREE_HPP_*/
