/** \file tb2pedigree.hpp
 *  \brief Pedigree data structure
 * 
 */
 
#ifndef TB2PEDIGREE_HPP_
#define TB2PEDIGREE_HPP_

#include "tb2wcsp.hpp"

class Genotype {
public:
  int allele1;
  int allele2;
};

class Individual {
public:
  int individual;
  int varindex;
  Genotype genotype;
};

class Pedigree {
  vector<Individual> genotypes;
  vector<Genotype> myconvert;

public:
  void readPedigree(const char *fileName, WCSP *wcsp);

  void printCorrection(WCSP *wcsp);
};

#endif /*TB2PEDIGREE_HPP_*/
