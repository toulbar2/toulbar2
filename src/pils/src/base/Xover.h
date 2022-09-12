/*
  xover.h

*/

//#ifndef _Xover_h
//#define _Xover_h

#include <string>
#include <iostream>
#include <streambuf>
#include <stdlib.h>
#include <set>


#include <base/solution.h>
#include <base/costFunction.h>
#include <base/crossover.h>

using namespace std;

class Xover : public Crossover {
public : 

  Xover(CostFunction & _eval) : eval(_eval) {

    nComponents = 1;

    component_id.resize(eval.n_variables);

    components.resize(eval.n_variables + 1);

  } 

  bool operator()(Solution & x1,Solution & x2, Solution & child) {

    child.resize(x1.size());

    compute_connected_components(x1, x2);

    Cost g1, g2;

    bool n1 = false;
    bool n2 = false;

    if (ToulBar2::verbose >= 1) {
        cout << nComponents << " ";
    }
    // first common component : copy x1 into child
    Cost child_fit = 0;
    child_fit = partial_func(0, x1) ;
    for(unsigned i : components[0])
      child[i] = x1[i];

    // Others true components
    for(unsigned c = 1; c <= nComponents; c++) {
      g1 = partial_func(c, x1);
      g2 = partial_func(c, x2);

      if (g1 >= g2) {
        // x2 is better for that component
        n2 = true;
        child_fit += g2;

        for(unsigned i : components[c])
          child[i] = x2[i];
      } else {
        // x1 is better for that component
        n1 = true;
        child_fit += g1;				

        for(unsigned i : components[c])
          child[i] = x1[i];
      }
    }

    child.fitness(child_fit);
    
    if (ToulBar2::verbose >= 1) {
        cout << (n1 && n2) << endl;
    }
    return (n1 && n2); // true when child is different from x1 and x2
  }


  void print_state() {
    std::cout << endl << "ids:" ;
    for(unsigned i : component_id) {
      std::cout << " "  << i ;
    }
    std::cout << std::endl;

    std::cout << "nb components : " << nComponents << std::endl;

    for(unsigned c = 0; c <= nComponents; c++) {
      std::cout << c << ":" ;
      for(unsigned i : components[c])
        {
          std::cout << " " << i << " links:" << endl ;
          for (unsigned j : eval.links[i])
            cout << j << " ";
          cout << endl;
        }
                                
      std::cout << std::endl;
    }

    std::cout << std::endl;
  }

  void afficheComponents() {

    for(unsigned i=0; i <= nComponents; i++){
      cout << endl << "comp" << i << endl;
      for(unsigned j : components[i]){
        cout << components[i][j] << " ";
      }
    }
    cout << endl;

  }


  //number of connected components
  unsigned nComponents;

  //number of the connected components of each variables
  std::vector<int> component_id;

  //list of variables of each connected component
  std::vector < std::vector<unsigned> > components;

  //costfunction to evaluate the variables
  CostFunction & eval;

protected:

  void compute_connected_components (Solution & x1, Solution & x2) {
    nComponents = 1;

    std::fill(component_id.begin(), component_id.end(), -1);

    for(unsigned i = 0; i < eval.n_variables; i++) // avoid to copy the references on x1 and x2 in set_component function
      if (x1[i] == x2[i])
        component_id[i] = 0; // in the share variable component

    //for all variables
    for (unsigned i = 0; i < eval.n_variables; ++i)
      {
        if (component_id[i] == -1) {
          nComponents++;
          set_component(i);

        }
      }

    for(unsigned c = 0; c <= nComponents; c++)
      components[c].resize( 0 );

    for(unsigned i = 0; i < eval.n_variables; i++) {
        assert(component_id[i] >= 0);
        assert(component_id[i] <= (int)nComponents );
        components[ component_id[i] ].push_back(i);
    }
  }

  void set_component(unsigned i) {
    if (component_id[i] == -1) {
      component_id[i] = nComponents;

      for(unsigned j : eval.links[i]) // all neighbors in the variable links matrix
        set_component(j);
      for(unsigned j : eval.backlinks[i]) // check out backlinks too
        set_component(j);
    }
  }
	
  Cost partial_func(unsigned c, Solution & x){
    Cost res= 0;
                
    vector<bool> inthere;
    inthere.resize(x.size());
    fill(inthere.begin(), inthere.end(), false);
    for(unsigned k : components[c])
      {
        inthere[k]=true;
      }
    for(unsigned k : components[c]) {
      res+= eval.energy[k][x[k]];  // add the energy of each AA
      for(unsigned i = 0; i < eval.links[k].size(); i++) {
        unsigned l = eval.links[k][i];
        if (inthere[l])
          res += eval.energy2[k][l][x[k]][x[l]];
      }
    }                    
                
    if (c!=0)
      {
        for(unsigned k : components[c]) {
          for(unsigned i = 0; i < eval.links[k].size(); i++) {
            unsigned l = eval.links[k][i];
            if (!inthere[l]) {
                assert(component_id[l] == 0);
                res += eval.energy2[k][l][x[k]][x[l]];
            } else {
                assert(component_id[k] == component_id[l]);
            }
          }
          for(unsigned i = 0; i < eval.backlinks[k].size(); i++) {
            unsigned l = eval.backlinks[k][i];
            if (!inthere[l]) {
                assert(component_id[l] == 0);
                res += eval.energy2[l][k][x[l]][x[k]];
            } else {
                assert(component_id[k] == component_id[l]);
            }
          }
        }
      }
    return res;
  }
};
