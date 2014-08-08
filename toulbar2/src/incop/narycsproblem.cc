#include <cerrno>
#include <stdio.h>
#include <list>
#include <vector>
#include <string>
#include <set>
#include <algorithm>


using namespace std;
#include <iostream>
#include <fstream>

#include "incop.h"
#include "incoputil.h"
#include "csproblem.h"
#include "narycsproblem.h"
#include "autotuning2.h"

extern ofstream* ofile;  // le fichier de sortie

extern Stat_GWW * Statistiques; 







NaryCSProblem::NaryCSProblem (int nbvar, int nbconst) : CSProblem (nbvar,nbconst) {;}

NaryConstraint::NaryConstraint ( int arit) { arity=arit;}

NaryVariable::NaryVariable () {;}

/** code optimisé pour configuration semi-incrementale IncrCSPConfiguration*/
/*
int NaryCSProblem::move_evaluation 
                    (Configuration* configuration,Move* move)
{int var_changee = ((CSPMove*)move)-> variable;
 int val_changee = ((CSPMove*)move)-> value;
 return(configuration->valuation
       +compute_conflict(configuration,var_changee,val_changee)
               -((IncrCSPConfiguration*)configuration)->tabconflicts[var_changee]);
}
*/




// optimisation pour IncrCSPConfiguration 
/*
int NaryCSProblem::config_evaluation(Configuration* configuration)
{ 
  configuration->init_conflicts();
  int value=0;
  for (int i=0; i< (int) naryconstraints->size();i++)
    {int nbconf = (*naryconstraints)[i]->constraint_value (configuration);
      value+= nbconf;
     for (int j=0 ; j< (*naryconstraints)[i]->arity ; j++)
       ((IncrCSPConfiguration*)configuration)->tabconflicts[(*naryconstraints)[i]->constrainedvariables[j]] +=nbconf;
    }
  return value;
}
*/


Long NaryCSProblem::config_evaluation(Configuration* configuration)
{ 
  configuration->init_conflicts();
  Long value=0;
  for (int i=0; i< (int) naryconstraints->size();i++)
    {Long nbconf = (*naryconstraints)[i]->constraint_value (configuration);
      value+= nbconf;
    }
  for (int i=0; i< nbvar ; i++)
     for (int j=0; j< variable_domainsize(i);j++)
       configuration->incr_conflicts(i,j,j,compute_conflict(configuration,i,j));
  return value;
}




Long NaryConstraint::constraint_value(Configuration* configuration)
{ int index=0;
 for (int  i=0; i<arity;i++)
   index+= configuration->config[constrainedvariables[i]] * multiplyers[i];
 return tuplevalues[index];
}


void NaryCSProblem::incr_update_conflicts (IncrCSPConfiguration* configuration, Move* move)
 { int var = ((CSPMove*)move)-> variable;
   int value = ((CSPMove*)move)-> value;
   int aval = configuration->config[var];
   Long actvalue,nctvalue;
   NaryVariable* varobjct= (*naryvariables)[var];
   NaryConstraint* ct;
   for (int i=0; i< (int) (varobjct->constraints).size() ; i++)
     {ct=(varobjct->constraints)[i];
      actvalue= ct->constraint_value (configuration);
      configuration->config[var]=value;
      nctvalue= ct->constraint_value (configuration);
      configuration->config[var]=aval;
      for (int j=0 ;  j< ct->arity; j++)
       {configuration->tabconflicts[ct->constrainedvariables[j]] += nctvalue - actvalue;
       }
     }
 }


void NaryCSProblem::fullincr_update_conflicts (FullincrCSPConfiguration* configuration, Move* move)
 { int var = ((CSPMove*)move)-> variable;
   int value = ((CSPMove*)move)-> value;
   int aval = configuration->config[var];
   Long actvalue,nctvalue;
   int var1, aval1;
   NaryVariable* varobjct= (*naryvariables)[var];
   NaryConstraint* ct;
   for (int i=0; i< (int) (varobjct->constraints).size() ; i++)
     {ct=(varobjct->constraints)[i];
      for (int j=0 ;  j< ct->arity; j++)
       { 
         var1= ct->constrainedvariables[j];
         if (var1 != var)
	   {aval1= configuration->config[var1];
	   for (int k=0; k< variable_domainsize(var1); k++)
	     {
	       configuration->config[var1]=k;
	       actvalue= ct->constraint_value (configuration);
	       configuration->config[var]=value;
	       nctvalue= ct->constraint_value (configuration);
	       configuration->config[var]=aval;
	       //	       configuration->incr_conflicts(var1,k,k,nctvalue - actvalue);
	       configuration->tabconflicts[var1][k]+= nctvalue - actvalue;
	     }
	   configuration->config[var1]=aval1;
	   }
       }
     }
 }
       

int NaryConstraint::compute_index(int* values, vector<int>* tabdomaines)
{ int index=0;
 for (int  i=0; i<arity;i++)
   index+= compute_indexpart( i, values[i], tabdomaines);
 return index;
}


int NaryConstraint::compute_indexpart (int i, int vali, vector<int>* tabdomaines)
{ int factor=1;
 for (int j=i+1; j< arity; j++)
   {factor = factor * tabdomaines[constrainedvariables[j]].size();}
 return vali*factor;
}


/** nombre de n-uplets d'une contrainte */
/* number of tuples of a constraint */
int NaryConstraint::nbtuples( vector<int>* tabdomaines)
{int nbtuples=1;
for (int j=0; j< arity; j++)
  {nbtuples = nbtuples * tabdomaines[constrainedvariables[j]].size();}
return nbtuples;
}

void NaryConstraint::compute_indexmultiplyers(vector<int>* tabdomaines)
{ 
 for (int i=0; i< arity ; i++)
   multiplyers.push_back(compute_indexmultiplyer(i ,tabdomaines));
}

int NaryConstraint::compute_indexmultiplyer(int i, vector<int>* tabdomaines)     
{ int factor=1;
 for (int j=i+1; j< arity; j++)
   {factor = factor * tabdomaines[constrainedvariables[j]].size();}
 return factor;
}

   




/** calcul du nombre de conflits d'une affectation - appele par l'évaluation d'un mouvement (cas incr)*/

Long NaryCSProblem::compute_conflict (Configuration* configuration, int var , int val)
{Long value=0;
int aval=configuration->config[var];
configuration->config[var]=val;
for (int i=0; i< (int) (*naryvariables)[var]->constraints.size() ; i++)
  value+= (*naryvariables)[var]->constraints[i]->constraint_value (configuration);
configuration->config[var]=aval;
return value;
}


/** utilisation des configurations "semi-incrementales"IncrCSPConfiguration - les conflits des valeurs courantes des variables
    sont stockés dans le tableau tabconflicts 
    ou tout-incrémentales  FullincrCSPConfiguration  : les conflits de toutes les valeurs avec la configuration courante
sont maintenus dans tabconflicts */
Configuration* NaryCSProblem::create_configuration()
    {
      return (new FullincrCSPConfiguration(nbvar,domainsize));           
      // return (new IncrCSPConfiguration(nbvar,domainsize));    

}

NaryCSProblem* weighted_narycsp_creation (int nbvar, int nbconst, int maxdomsize, 
 vector<NaryVariable*>* vv,vector<NaryConstraint*>* vct     )
{NaryCSProblem*  p1 =new  NaryCSProblem (nbvar,nbconst);
 p1->domainsize=maxdomsize;
 p1->naryconstraints  = vct;
 p1->naryvariables  = vv;
 return p1;
}


/** lecture du debut du fichier : le probleme et les variables */
void  wcspdomaines_file_read (ifstream & file, int nbvar, vector<int>* tabdomaines)
{
  //  cout << "nbvar=" << nbvar;
 int size=0;
 for (int i=0; i<nbvar; i++)
   {file >> size;
	 //	 cout << " " << size;
   for (int j=0 ; j<size ; j++)
     tabdomaines[i].push_back(j);
   }
 // cout << endl;
}

/** lecture des contraintes */
void  wcspdata_constraint_read (ifstream & file, int nbconst, vector<NaryVariable*>* vv, vector<NaryConstraint*>* vct, 
				vector <int>* connexions, vector<int> * tabdomaines)
{ 

  for (int i =0; i< nbconst; i++)
  {
    int arity=0; int numvar=0; Long defaultcost=0;
    file >> arity;

    NaryConstraint* ct = new NaryConstraint(arity);
    vct->push_back(ct);
    for (int j=0 ; j< arity ; j++)
      {file >> numvar;
      ct->constrainedvariables.push_back(numvar);
      (*vv)[numvar]->constraints.push_back(ct);
      }
    ct->compute_indexmultiplyers(tabdomaines);
	
    file >> defaultcost;
    int nbtuples = ct->nbtuples(tabdomaines);
    for (int t=0;t<nbtuples;t++)
      {ct->tuplevalues.push_back(defaultcost);}
    int nbcosttuples=0;
    Long tuplevalue=0;
    file >> nbcosttuples;
	//	cout << "Constraint " << i << " " << arity << " " << defaultcost << " " << nbcosttuples << " " << nbtuples << endl;
    int* values = new int[arity];
    for (int t =0; t< nbcosttuples; t++)
      {for (int j=0 ; j< arity; j++) file >> values[j];
      file >> tuplevalue;
      ct->tuplevalues[ct->compute_index(values,tabdomaines)]= tuplevalue;
      }

    delete [] values;
  }

}

//SdG: toulbar2 system call should be replace by narycsp("/dev/stdout", "problem.wcsp", 1, 3, "idwa", 100000, "cv", "v", 0, 200, 1, 0, 0)
//int narycsp(char *filename, char *problem, int graine1, int nbessais, string & *method,  int tuningmode)
int narycsp(int argc, char **argv,  int tuningmode)
{


  NaryCSProblem* problem ;          // pointeur sur le probleme 

  // les divers arguments lus dans la ligne de commande
  int nbvar,nbconst, domsize;
  Long lbound;
  int taille,nbessais;
  int graine1;
  int narg = 2;  // compteur des arguments

  // le nom du fichier de sortie : pour les tests : version logiciel + concaténation des arguments
  char filename [1000];
  if ((string)argv[1] == "arg")
    ofile_name(filename, argc, argv);
  else sprintf(filename,"%s",argv[1]);


  ofstream ofile1 (filename);
  ofile = & ofile1;
  ifstream file (argv[2]); // le fichier de données 

  arguments_borneinf(argv,narg,lbound);
  // lecture des paramètres de l'algo et création de l'objet algo
  IncompleteAlgorithm* algo = algo_creation (argv, narg, taille, graine1, nbessais);

  // allocation de l'objet pour les stats
  Statistiques=new Stat_GWW (1, nbessais);
  

  // argument pour la trace
  arguments_tracemode(argv,narg);
  // pour la recuperation du signal 10
  sigaction();

 // argument de temps maximum 
  double maxtime;
  if (tuningmode) arguments_tempscpu (argv,narg,maxtime);


  // Declaration des variables contenant les structures de données des problemes


  
  string pbname;
  Long upperbound;
  file >> pbname; // nom du  probleme
  file >> nbvar;
  file >> domsize;
  file >> nbconst;
  file >> upperbound;

  vector<int> tabdomaines[nbvar] ; // les différents types de domaines 
  wcspdomaines_file_read (file,nbvar, tabdomaines);

  int domaines[nbvar];  // 1 domaine par variable
  for(int i=0;i<nbvar;i++)
    {domaines[i]=i;}

  // Initialisation des structures de données des problèmes


  vector<NaryConstraint*> constraints[nbconst];
  vector<NaryVariable*> variables[nbvar];
  vector<int> connexions [nbvar];

  for (int i=0;i<nbvar;i++)
    {NaryVariable* nv = new NaryVariable();
    variables->push_back(nv);}

  
  wcspdata_constraint_read (file, nbconst, variables, constraints, connexions, tabdomaines);
  int pbnumber=0;
  Statistiques->init_pb(pbnumber);
  problem = weighted_narycsp_creation (nbvar,nbconst,domsize,variables,constraints);

  problem->lower_bound=lbound;
  // mise en place des domaines 
  problem->set_domains_connections(domaines,tabdomaines,connexions);

    
  // creation de la population et initialisation 
  // La population : tableau de configurations
  Configuration* population[taille];
 
  problem->init_population(population,taille);
 
  problem->allocate_moves();

  if (tuningmode)
    autosolving((LSAlgorithm*)algo,population,problem,0,graine1,nbessais,maxtime,1000000);
  else
    {
    *ofile << " pb  " << pbname;
    // boucle sur les essais 
    for(int nessai = 0;nessai< nbessais ; nessai++)
      executer_essai (problem,algo,population,taille,graine1,nessai);

    // ecriture statistiques 
    Statistiques->current_try++; 
    ecriture_stat_probleme();
    }
  delete problem;


      
  cout << "INCOP time : " << Statistiques->total_execution_time << endl;
  return 0;
  
  
}
