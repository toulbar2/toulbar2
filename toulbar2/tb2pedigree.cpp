/*
 * **************** Read pedigree in pre format files **************************
 * 
 */

#include <map>

#include "tb2system.hpp"
#include "tb2wcsp.hpp"
#include "tb2enumvar.hpp"

typedef struct {
    EnumeratedVariable *var;
    vector<Cost> costs;
} TemporaryUnaryConstraint;


void WCSP::read_pedigree(const char *fileName)
{
  int nbindividuals = 0;
  int nballeles = 0;
  int nbtypings = 0;
  int nbfounders = 0;
  map<int, int> individuals;
  map<int, int> individualsReverse;
  map<int, bool> notfounders;
  vector<TemporaryUnaryConstraint> unaryconstrs;

  // open the file
  ifstream file(fileName);
  if (!file) {
    cerr << "Could not open file " << fileName << endl;
    exit(EXIT_FAILURE);
  }
  
  while (file) {
    int locus = -1;
    int individual = 0;
    int father = 0;
    int mother = 0;
    int sex = -1;
    int allele1 = 0;
    int allele2 = 0;

    file >> locus;
    if (!file) break;

    file >> individual;
    if (!file) {
      cerr << "Wrong data after individual " << individual << endl;
      abort();
    }
    assert(individual != 0);
    if (individuals.count(individual)==0) {
      individuals[individual] = nbindividuals;
      individualsReverse[nbindividuals] = individual;
      nbindividuals++;
    }
    file >> father;
    if (!file) {
      cerr << "Wrong data after individual " << individual << endl;
      abort();
    }
    if (father > 0 && individuals.count(father)==0) {
      individuals[father] = nbindividuals;
      individualsReverse[nbindividuals] = father;
      nbindividuals++;
    }
    file >> mother;
    if (!file) {
      cerr << "Wrong data after individual " << individual << endl;
      abort();
    }
    if (mother > 0 && individuals.count(mother)==0) {
      individuals[mother] = nbindividuals;
      individualsReverse[nbindividuals] = mother;
      nbindividuals++;
    }
    if (father > 0 || mother > 0) notfounders[individual] = true;
    file >> sex;
    if (!file) {
      cerr << "Wrong data after individual " << individual << endl;
      abort();
    }
    file >> allele1;
    if (!file) {
      cerr << "Wrong data after individual " << individual << endl;
      abort();
    }
    nballeles = max(nballeles,allele1);
    file >> allele2;
    if (!file) {
      cerr << "Wrong data after individual " << individual << endl;
      abort();
    }
    nballeles = max(nballeles,allele2);
    if (allele1>0 || allele2>0) nbtypings++;
  }

  assert(vars.empty());
  assert(constrs.empty());
  updateUb(nbtypings+1);
  NCBucketSize = cost2log2(getUb()) + 1;
  
  /* create variables */
  for (int i=0; i<nbindividuals; i++) {
    char varname_[128];
    sprintf(varname_, "%d", individualsReverse[i]);
    string varname(varname_);
    makeEnumeratedVariable(varname, 0, nballeles*(nballeles+1)/2 - 1);
    if (notfounders.count(individualsReverse[i])==0) nbfounders++;
  }

  /* create ternary Mendelian hard constraint table */
  vector<Cost> costs3;
  for (int i=1; i<=nballeles; i++) { /* i = first allele of child */
    for (int j=i; j<=nballeles; j++) { /* j = second allele of child */
      for (int k=1; k<=nballeles; k++) { /* k = first allele of father */
	for (int l=k; l<=nballeles; l++) { /* l = second allele of father */
	  for (int m=1; m<=nballeles; m++) { /* m = first allele of mother */
	    for (int n=m; n<=nballeles; n++) { /* n = second allele of mother */
	      costs3.push_back(((i==k && j==m) || (i==k && j==n) || (i==l && j==m) || (i==l && j==n) || (i==m && j==k) || (i==m && j==l) || (i==n && j==k) || (i==n && j==l))?0:nbtypings+1);
	    }
	  }
	}
      }
    }
  }

  /* create binary Mendelian hard constraint table */
  vector<Cost> costs2;
  for (int i=1; i<=nballeles; i++) { /* i = first allele of child */
    for (int j=i; j<=nballeles; j++) { /* j = second allele of child */
      for (int k=1; k<=nballeles; k++) { /* k = first allele of father or mother */
	for (int l=k; l<=nballeles; l++) { /* l = second allele of father or mother */
	  costs2.push_back((i==k || i==l || j==k || j==l)?0:nbtypings+1);
	}
      }
    }
  }

  // re-read the file
  file.clear();
  file.seekg(0, ios::beg);

  /* create constraint network */
  while (file) {
    int locus = -1;
    int individual = 0;
    int father = 0;
    int mother = 0;
    int sex = -1;
    int allele1 = 0;
    int allele2 = 0;

    file >> locus;
    if (!file) break;
    file >> individual;
    file >> father;
    file >> mother;
    file >> sex;
    file >> allele1;
    file >> allele2;

    /* add unary costs (soft constraint) if genotyping is given */
    if (allele1 > 0 || allele2 > 0) {
      EnumeratedVariable *var = (EnumeratedVariable *) vars[individuals[individual]];
      TemporaryUnaryConstraint unaryconstr;
      unaryconstr.var = var;
      for (int i=1; i<=nballeles; i++) { /* i = first allele of child */
	for (int j=i; j<=nballeles; j++) { /* j = second allele of child */
	  if ((allele1>0 && allele2>0 && ((i==allele1 && j==allele2) || (i==allele2 && j==allele1)))
		|| ((allele1==0 || allele2==0) && (i==allele1 || i==allele2 || j==allele1 || j==allele2))) {
	    unaryconstr.costs.push_back(0);
	  } else {
	    unaryconstr.costs.push_back(1);
	  }
	}
      }
      unaryconstrs.push_back(unaryconstr);
      var->queueNC();
    }

    /* add ternary or binary Mendelian hard constraint */
    if (father > 0 || mother > 0) {
      if (father > 0 && mother > 0) {
	postTernaryConstraint(individuals[individual],individuals[father],individuals[mother],costs3);
      } else if (father > 0) {
	postBinaryConstraint(individuals[individual],individuals[father],costs2);
      } else {
	postBinaryConstraint(individuals[individual],individuals[mother],costs2);
      }
    }
  }
  sortConstraints();

  // apply basic initial propagation AFTER complete network loading
  for (unsigned int u=0; u<unaryconstrs.size(); u++) {
    for (unsigned int a = 0; a < unaryconstrs[u].var->getDomainInitSize(); a++) {
      if (unaryconstrs[u].costs[a] > 0) unaryconstrs[u].var->project(a, unaryconstrs[u].costs[a]);
    }
    unaryconstrs[u].var->findSupport();
  }
  
  if (ToulBar2::verbose >= 0) {
    cout << "Read pedigree with " << nbindividuals << " individuals, " << nbfounders << " founders, " << nballeles << " alleles and " << nbtypings << " typings." << endl;
  }
}
