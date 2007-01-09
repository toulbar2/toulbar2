/*
 * **************** Read pedigree in pre format files **************************
 * 
 */

#include "tb2system.hpp"
#include "toulbar2.hpp"
#include "tb2enumvar.hpp"
#include "tb2pedigree.hpp"


void Pedigree::iniProb( WCSP* wcsp ) { 
   Cost TopProb = 0;

   ToulBar2::NormFactor = (-log(10.)/log1p( - pow(10., -(double)ToulBar2::resolution)));
   if (ToulBar2::NormFactor > (pow(2,(float)INTEGERBITS)-1)/-log10(pow(10, -(float)ToulBar2::resolution))) {
	  fprintf(stderr,"This resolution cannot be ensured on the data type used to represent costs.\n");
	  fprintf(stderr,"Please use toolbarl. Aborting.\n");
	  abort();
   }
	
   int nballeles = alleles.size()-1;
   bool overflow = false;

   TopProb = 0;
   
   int ngenotyped = genotypes.size();
   while(ngenotyped) {
   		if(TopProb < MAX_COST) TopProb += -log10(ToulBar2::errorg / (TProb)(nballeles-1)) * ToulBar2::NormFactor;
   		else overflow = true; 
   		ngenotyped--;
   }
  
   for (vector<Individual>::iterator iter = pedigree.begin(); iter != pedigree.end(); ++iter) {
  	Individual& individual = *iter;
	if(individual.typed) {
		if(individual.mother && individual.father) {
			if(TopProb < MAX_COST) TopProb += -log10(0.25) * ToulBar2::NormFactor;
			else overflow = true;
		}
		else if(individual.mother || individual.father) {
			if(TopProb < MAX_COST) TopProb += -log10(0.50) * ToulBar2::NormFactor;
			else overflow = true;
		}
		else
		{
			TProb minp = 1;
			switch(ToulBar2::foundersprob_class)
			{
				case 0:  if(TopProb < MAX_COST) TopProb += -log10(1./nballeles) * ToulBar2::NormFactor;
						 else overflow = true;			
						 break;
						
				case 1:  for (map<int,int>::iterator iter = freqalleles.begin(); iter != freqalleles.end(); ++iter) {
							TProb p = (TProb) (iter->second * iter->second)/ (TProb)(genotypes.size() * genotypes.size() * 4);
							if(p < minp) minp = p;
						 }
						
						 if(TopProb < MAX_COST) TopProb += -log10(minp) * ToulBar2::NormFactor;
						 else overflow = true;			
						 break;
						 
				default: for (vector<TProb>::iterator iter = foundersprob.begin(); iter != foundersprob.end(); ++iter) {
				        	if(*iter < minp) minp = *iter;
				         }
						 if(TopProb < MAX_COST) TopProb += -log10(minp) * ToulBar2::NormFactor;
						 else overflow = true;			
			}
		}
	}
  }	
  wcsp->updateUb(TopProb);
  if(overflow) cout << "Possible risk of overflow: product of min probabilities < size of used datatype" << endl;
}







typedef struct {
    EnumeratedVariable *var;
    vector<Cost> costs;
} TemporaryUnaryConstraint;

Individual::Individual(int ind)
{
    individual = ind;
    varindex = -1;
    father = 0;
    mother = 0;
    sex = MALE;
    genotype.allele1 = 0;
    genotype.allele2 = 0;
    typed = false;
}

void Pedigree::typeAscendants(int individual)
{
    if (individual > 0) {
        assert(individuals.count(individual)!=0);
        int index = individuals[individual];
        if (!pedigree[index].typed) {
            pedigree[index].typed = true;
            nbtyped++;
            typeAscendants(pedigree[index].father);
            typeAscendants(pedigree[index].mother);
        }
    }
}



void Pedigree::read(const char *fileName, WCSP *wcsp)
{
	readPedigree(fileName, wcsp);
	buildWCSP(fileName, wcsp);
}

void Pedigree::read_bayesian(const char *fileName, WCSP *wcsp )
{
	readPedigree(fileName, wcsp);
	buildWCSP_bayesian(fileName, wcsp);
}



// warning! locus information is not used: assume that only one locus is defined in the pedigree file
void Pedigree::readPedigree(const char *fileName, WCSP *wcsp)
{
  int nbindividuals = 0;
  int nballeles = 0;
  int nbtypings = 0;
  map<int, int> allelesInv;
  
  // open the file
  ifstream file(fileName);
  if (!file) {
    cerr << "Could not open file " << fileName << endl;
    exit(EXIT_FAILURE);
  }

  while (file) {
    int cur_locus = -1;
    int individual = 0;

    file >> cur_locus;
    if (!file) break;
    if (locus == -1) locus = cur_locus;
    if (locus != cur_locus) {
        cerr << "Pedigree datafile contains more than one locus!" << endl;
        exit(EXIT_FAILURE);
    } 

    file >> individual;
    if (!file) {
      cerr << "Wrong data after individual " << individual << endl;
      abort();
    }
    assert(individual != 0);
    if (individuals.count(individual)==0) {
      individuals[individual] = nbindividuals;
      Individual geno(individual);
      pedigree.push_back(geno);
      nbindividuals++;
    }

    file >> pedigree[individuals[individual]].father;
    if (!file) {
      cerr << "Wrong data after individual " << individual << endl;
      abort();
    }
    if (pedigree[individuals[individual]].father > 0 && individuals.count(pedigree[individuals[individual]].father)==0) {
      individuals[pedigree[individuals[individual]].father] = nbindividuals;
      Individual geno(pedigree[individuals[individual]].father);
      geno.sex = MALE;
      pedigree.push_back(geno);
      nbindividuals++;
    }

    file >> pedigree[individuals[individual]].mother;
    if (!file) {
      cerr << "Wrong data after individual " << individual << endl;
      abort();
    }
    if (pedigree[individuals[individual]].mother > 0 && individuals.count(pedigree[individuals[individual]].mother)==0) {
      individuals[pedigree[individuals[individual]].mother] = nbindividuals;
      Individual geno(pedigree[individuals[individual]].mother);
      geno.sex = FEMALE;
      pedigree.push_back(geno);
      nbindividuals++;
    }

    file >> pedigree[individuals[individual]].sex;
    if (!file) {
      cerr << "Wrong data after individual " << individual << endl;
      abort();
    }

    file >> pedigree[individuals[individual]].genotype.allele1;    
    if (!file) { cerr << "Wrong data after individual " << individual << endl; abort(); }

    file >> pedigree[individuals[individual]].genotype.allele2;
    if (!file) { cerr << "Wrong data after individual " << individual << endl; abort(); }
    
    int allele1 = pedigree[individuals[individual]].genotype.allele1;
    int allele2 = pedigree[individuals[individual]].genotype.allele2;
    
    if (alleles.count(allele1)==0) {
        nballeles++;
        alleles[allele1] = nballeles;
        freqalleles[ allele1 ] = 0;
    }
    else { freqalleles[ allele1 ]++; }
   

    if (alleles.count(allele2)==0) {
        nballeles++;
        alleles[allele2] = nballeles;
        freqalleles[ allele2 ] = 0;
    }
    else { freqalleles[ allele2 ]++; }

    
    if (allele1>0 || allele2>0) {
      nbtypings++;
      genotypes.push_back(individual);
    }
  }

  /* re-encoding of alleles */
  int nb = 0;
  if (ToulBar2::verbose >= 2) cout << "Alleles encoding:" << endl;
  for (map<int,int>::iterator iter = alleles.begin(); iter != alleles.end(); ++iter) {
    if ((*iter).first == 0) continue;
    nb++;
    if (ToulBar2::verbose >= 2) cout << (*iter).first << ": " << nb << endl;
    alleles[(*iter).first] = nb;
    allelesInv[nb] = (*iter).first;
  }
  assert(nballeles == nb);
  
  if (ToulBar2::verbose) cout << "Genotype encoding:" << endl;
  for (int i=1; i<=nballeles; i++) { /* i = first allele of child */
    for (int j=i; j<=nballeles; j++) { /* j = second allele of child */
      Genotype geno;
      geno.allele1 = allelesInv[i];
      geno.allele2 = allelesInv[j];
      genoconvert.push_back(geno);
      if (ToulBar2::verbose) {
        cout << genoconvert.size()-1 << ": ";
        printGenotype(cout, genoconvert.size()-1);
        cout << endl;
      }
    }
  }
    
  for (unsigned int i=0; i<genotypes.size(); i++) {
    typeAscendants(genotypes[i]);
  }
  cout << nbtyped << " individuals found with a genotyped descendant.." << endl;
  
  assert(wcsp->numberOfVariables() == 0);
  assert(wcsp->numberOfConstraints() == 0);
}  



  
void Pedigree::buildWCSP(const char *fileName, WCSP *wcsp)
{
  ifstream file(fileName);
  
  vector<TemporaryUnaryConstraint> unaryconstrs;
 
  int nbindividuals = individuals.size();
  int nballeles = alleles.size() - 1;       // a demander pourquoi le -1
  int nbfounders = 0;
  int nbtypings = genotypes.size();

  wcsp->updateUb(nbtypings+1);

  /* create variables */
  int nbvar = 0;
  for (int i=0; i<nbindividuals; i++) {
    if (pedigree[i].father == 0 && pedigree[i].mother == 0) nbfounders++;
    if (pedigree[i].typed) {
        string varname;
        varname = to_string(pedigree[i].individual);
        wcsp->makeEnumeratedVariable(varname, 0, nballeles*(nballeles+1)/2 - 1);
        pedigree[i].varindex = nbvar;
        nbvar++;
    }
  }

  /* create ternary Mendelian hard constraint table */
  vector<Cost> costs3;
  for (int k=1; k<=nballeles; k++) { /* k = first allele of father */
	for (int l=k; l<=nballeles; l++) { /* l = second allele of father */
	  for (int m=1; m<=nballeles; m++) { /* m = first allele of mother */
	    for (int n=m; n<=nballeles; n++) { /* n = second allele of mother */
            for (int i=1; i<=nballeles; i++) { /* i = first allele of child */
                for (int j=i; j<=nballeles; j++) { /* j = second allele of child */
	               costs3.push_back(((i==k && j==m) || (i==k && j==n) || (i==l && j==m) || (i==l && j==n) || (i==m && j==k) || (i==m && j==l) || (i==n && j==k) || (i==n && j==l))?0:wcsp->getUb());
                }
            }
        }
      }
    }
  }

  /* create binary Mendelian hard constraint table */
  vector<Cost> costs2;
  for (int k=1; k<=nballeles; k++) { /* k = first allele of father or mother */
     for (int l=k; l<=nballeles; l++) { /* l = second allele of father or mother */
        for (int i=1; i<=nballeles; i++) { /* i = first allele of child */
            for (int j=i; j<=nballeles; j++) { /* j = second allele of child */
	           costs2.push_back((i==k || i==l || j==k || j==l)?0:wcsp->getUb());
	        }
        }
     }
  }

  /* create constraint network */
  while (file) {
    int cur_locus = -1;
    int individual = 0;
    int father = 0;
    int mother = 0;
    int sex = -1;
    int allele1 = 0;
    int allele2 = 0;

    file >> cur_locus;
    if (!file) break;
    file >> individual;
    file >> father;
    file >> mother;
    file >> sex;
    file >> allele1;
    allele1 = alleles[allele1];
    file >> allele2;
    allele2 = alleles[allele2];
    if (!pedigree[individuals[individual]].typed) continue;
    
    /* add unary costs (soft constraint) if genotyping is given */
    if (allele1 > 0 || allele2 > 0) {
      EnumeratedVariable *var = (EnumeratedVariable *) wcsp->getVar(pedigree[individuals[individual]].varindex);
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
           assert(pedigree[individuals[pedigree[individuals[individual]].father]].typed);
           assert(pedigree[individuals[pedigree[individuals[individual]].mother]].typed);
	       wcsp->postTernaryConstraint(pedigree[individuals[pedigree[individuals[individual]].father]].varindex,pedigree[individuals[pedigree[individuals[individual]].mother]].varindex,pedigree[individuals[individual]].varindex,costs3);
      } else if (father > 0) {
	       wcsp->postBinaryConstraint(pedigree[individuals[pedigree[individuals[individual]].father]].varindex,pedigree[individuals[individual]].varindex,costs2);
      } else {
	       wcsp->postBinaryConstraint(pedigree[individuals[pedigree[individuals[individual]].mother]].varindex,pedigree[individuals[individual]].varindex,costs2);
      }
    }
  }
  wcsp->sortConstraints();

  // apply basic initial propagation AFTER complete network loading
  for (unsigned int u=0; u<unaryconstrs.size(); u++) {
    for (unsigned int a = 0; a < unaryconstrs[u].var->getDomainInitSize(); a++) {
      if (unaryconstrs[u].costs[a] > 0) unaryconstrs[u].var->project(a, unaryconstrs[u].costs[a]);
    }
    unaryconstrs[u].var->findSupport();
  }
  
  if (ToulBar2::verbose >= 0) {
    cout << "Read pedigree with " << nbindividuals << " individuals, " << nbfounders << " founders, " << nballeles << " alleles and " << nbtypings << " genotypings." << endl;
  }
}



void Pedigree::buildWCSP_bayesian( const char *fileName, WCSP *wcsp )
{
  ifstream file(fileName);

  vector<TemporaryUnaryConstraint> unaryconstrs;
  int nbindividuals = individuals.size();
  int nballeles = alleles.size() - 1;
  int domsize = nballeles * ( nballeles + 1 ) / 2;
  int nbfounders = 0;
  int nbtypings = genotypes.size();

  iniProb( wcsp ); 

  /* create variables */
  int nbvar = 0;
  for (int i=0; i<nbindividuals; i++) {
    if (pedigree[i].father == 0 && pedigree[i].mother == 0) nbfounders++;
    if (pedigree[i].typed) {
        string varname;
        varname = to_string(pedigree[i].individual);
        wcsp->makeEnumeratedVariable(varname, 0, nballeles*(nballeles+1)/2 - 1);
        pedigree[i].varindex = nbvar;
        nbvar++;
    }
  }

  /* create ternary Mendelian hard constraint table */
  vector<Cost> costs3;
  for (int k=1; k<=nballeles; k++) { /* k = first allele of father */
	for (int l=k; l<=nballeles; l++) { /* l = second allele of father */
	  for (int m=1; m<=nballeles; m++) { /* m = first allele of mother */
	    for (int n=m; n<=nballeles; n++) { /* n = second allele of mother */
            for (int i=1; i<=nballeles; i++) { /* i = first allele of child */
                for (int j=i; j<=nballeles; j++) { /* j = second allele of child */
				   TProb p = 0.;
				   if((i==k && j==m) || (i==m && j==k)) p += 0.25; 
				   if((i==k && j==n) || (i==n && j==k)) p += 0.25;
				   if((i==l && j==m) || (i==m && j==l)) p += 0.25;
				   if((i==l && j==n) || (i==n && j==l)) p += 0.25;
	               costs3.push_back( wcsp->Prob2Cost(p) );
                }
            }
        }
      }
    }
  }

  map<int, int> allelesInv;

  int i,j;
  switch(ToulBar2::foundersprob_class) {
  	case 0: foundersprob.clear();
		  	for (i=1; i<=nballeles; i++) { /* i = first allele of child  */
		      for (j=i; j<=nballeles; j++) { /* j = second allele of child */
		         foundersprob.push_back( 1./(TProb)domsize );
		      }
		  	}
		    break;
		     
     case 1: foundersprob.clear();
             for (map<int,int>::iterator iter = alleles.begin(); iter != alleles.end(); ++iter) {
		    	allelesInv[iter->second] = iter->first;
		  	 }
		   	 for (i=1; i<=nballeles; i++) { /* i = first allele of child  */
		       for (j=i; j<=nballeles; j++) { /* j = second allele of child */
				  foundersprob.push_back( (TProb)freqalleles[ allelesInv[i] ]*(TProb)freqalleles[ allelesInv[j] ] / (TProb)(nbtypings * nbtypings * 4) );		  
		       }
		   	 }    
			 break;
			 
	 default:;
  }


  /* create binary Mendelian hard constraint table */
  vector<Cost> costs2;
  for (int k=1; k<=nballeles; k++) { /* k = first allele of father or mother */
     for (int l=k; l<=nballeles; l++) { /* l = second allele of father or mother */
        for (int i=1; i<=nballeles; i++) { /* i = first allele of child */
            for (int j=i; j<=nballeles; j++) { /* j = second allele of child */
			   TProb p = 0; 	
			   if(i==k || i==l || j==k || j==l) {
				   if(k == l)   p += 1. / (TProb)nballeles;
				   else 		p += 1. / ((TProb)nballeles + (TProb)nballeles - 1); 
			   }
               costs2.push_back( wcsp->Prob2Cost(p) );
	        }
        }
     }
  }


  /* create constraint network */
  while (file) {
    int cur_locus = -1;
    int individual = 0;
    int father = 0;
    int mother = 0;
    int sex = -1;
    int allele1 = 0;
    int allele2 = 0;

    file >> cur_locus;
    if (!file) break;
    file >> individual;
    file >> father;
    file >> mother;
    file >> sex;
    file >> allele1;
    allele1 = alleles[allele1];
    file >> allele2;
    allele2 = alleles[allele2];
    if (!pedigree[individuals[individual]].typed) continue;
    
    EnumeratedVariable *var = (EnumeratedVariable *) wcsp->getVar(pedigree[individuals[individual]].varindex);
    
    /* add unary costs (soft constraint) if genotyping is given */
    if (allele1 > 0 || allele2 > 0) {
	    TemporaryUnaryConstraint unaryconstr;
	    unaryconstr.var = var;
	    for (int i=1; i<=nballeles; i++) { /* i = first allele of child  */
		  for (int j=i; j<=nballeles; j++) { /* j = second allele of child */
		     bool typed = allele1>0 && allele2>0;
		     bool halftyped = allele1>0 || allele2>0;

		     bool theone = (i==allele1 && j==allele2) || (i==allele2 && j==allele1); 
		     bool posible = (i==allele1 || i==allele2 || j==allele1 || j==allele2);
		     
		     TProb p = 0;
		     if (typed) {
		     	if(theone) p = 1. - ToulBar2::errorg;
		        else	   p = ToulBar2::errorg / (TProb)(domsize-1);
		  	 }
		  	 else if(halftyped) {
		  	 	if(posible) p = (1. - ToulBar2::errorg) / (TProb)nballeles;
		  	 	else p = ToulBar2::errorg / (TProb)(domsize - nballeles);
  	  	 	 }	
 	  	 	 unaryconstr.costs.push_back( wcsp->Prob2Cost(p) );	
		  }
	    }
	    unaryconstrs.push_back(unaryconstr);
	    var->queueNC();
    }

    int id_father = individuals[pedigree[individuals[individual]].father];
    int id_mother = individuals[pedigree[individuals[individual]].mother];
 
    /* add ternary or binary Mendelian hard constraint */
    if (father > 0 || mother > 0) {
      if (father > 0 && mother > 0) {
           assert(pedigree[id_father].typed);
           assert(pedigree[id_mother].typed);
	       wcsp->postTernaryConstraint(pedigree[id_father].varindex,pedigree[id_mother].varindex,pedigree[individuals[individual]].varindex,costs3);
      } else if (father > 0) {
	       wcsp->postBinaryConstraint(pedigree[id_father].varindex,pedigree[individuals[individual]].varindex,costs2);
      } else {
	       wcsp->postBinaryConstraint(pedigree[id_mother].varindex,pedigree[individuals[individual]].varindex,costs2);
      }
    }
    else {
      TemporaryUnaryConstraint unaryconstr;
      unaryconstr.var = var;
      for (vector<TProb>::iterator iter = foundersprob.begin(); iter != foundersprob.end(); ++iter) {
		unaryconstr.costs.push_back( wcsp->Prob2Cost( *iter) );
      }
      unaryconstrs.push_back(unaryconstr);
      var->queueNC();
    }
  }
  wcsp->sortConstraints();

  // apply basic initial propagation AFTER complete network loading
  for (unsigned int u=0; u<unaryconstrs.size(); u++) {
    for (unsigned int a = 0; a < unaryconstrs[u].var->getDomainInitSize(); a++) {
      if (unaryconstrs[u].costs[a] > 0) unaryconstrs[u].var->project(a, unaryconstrs[u].costs[a]);
    }
    unaryconstrs[u].var->findSupport();
  }
  
  if (ToulBar2::verbose >= 0) {
    int nbtypings = genotypes.size();
    cout << "Read pedigree with " << nbindividuals << " individuals, " << nbfounders << " founders, " << nballeles << " alleles and " << nbtypings << " genotypings." << endl;
    cout << "Bayesian translated WCSP (errorg=" <<  ToulBar2::errorg << ", resolution=" << ToulBar2::resolution << ", ub:" << wcsp->getUb() << ")" << endl;
   
  }  
}




void Pedigree::printCorrection(WCSP *wcsp)
{
  cout << "Correction:";
  for (unsigned int i=0; i<genotypes.size(); i++) {
    int sol = wcsp->getValue(pedigree[individuals[genotypes[i]]].varindex);
    int a1 = genoconvert[sol].allele1;
    int a2 = genoconvert[sol].allele2;
    int allele1 = pedigree[individuals[genotypes[i]]].genotype.allele1;
    int allele2 = pedigree[individuals[genotypes[i]]].genotype.allele2;
    if (!((allele1>0 && allele2>0 && ((a1==allele1 && a2==allele2) || (a1==allele2 && a2==allele1)))
		|| ((allele1==0 || allele2==0) && (a1==allele1 || a1==allele2 || a2==allele1 || a2==allele2)))) {
      cout << " " << genotypes[i];
    }
  }
  cout << endl;
}

void Pedigree::printGenotype(ostream& os, Value value)
{
    os << genoconvert[value].allele1 << "/" << genoconvert[value].allele2;
}

void Individual::print(ostream& os)
{
    os << individual << " " << father << " " << mother << " " << sex << " " << genotype.allele1 << " " << genotype.allele2 << endl; 
}

void Pedigree::save(const char *fileName, WCSP *wcsp, bool corrected)
{
    // open the file
    ofstream file(fileName);
    if (!file) {
      cerr << "Could not open file " << fileName << endl;
      exit(EXIT_FAILURE);
    }

    for (map<int,int>::iterator iter = individuals.begin(); iter != individuals.end(); ++iter) {
        file << locus << " ";
        if (corrected && pedigree[(*iter).second].varindex >= 0 && pedigree[(*iter).second].varindex < (int) wcsp->numberOfVariables() && wcsp->assigned(pedigree[(*iter).second].varindex)) {
            int sol = wcsp->getValue(pedigree[(*iter).second].varindex);
            int a1 = genoconvert[sol].allele1;
            int a2 = genoconvert[sol].allele2;
            int allele1 = pedigree[(*iter).second].genotype.allele1;
            int allele2 = pedigree[(*iter).second].genotype.allele2;
            if (!((allele1>0 && allele2>0 && ((a1==allele1 && a2==allele2) || (a1==allele2 && a2==allele1)))
                || ((allele1==0 || allele2==0) && (a1==allele1 || a1==allele2 || a2==allele1 || a2==allele2)))) {
              pedigree[(*iter).second].genotype.allele1 = 0;
              pedigree[(*iter).second].genotype.allele2 = 0;
            }
            pedigree[(*iter).second].print(file);
            pedigree[(*iter).second].genotype.allele1 = allele1;
            pedigree[(*iter).second].genotype.allele2 = allele2;
        } else {
            pedigree[(*iter).second].print(file);
        }
    }
}
