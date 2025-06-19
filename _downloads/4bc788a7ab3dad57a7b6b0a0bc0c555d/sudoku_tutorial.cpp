#include <iostream>
#include <toulbar2lib.hpp>

using namespace std;

// print a solution as a grid
void printSolution(const std::vector<Value>& solution) {

   cout << "-------------------------" << endl;

   size_t var_index = 0;

   for(size_t row = 0; row < 9; row ++) {
      for(size_t col = 0; col < 9; col ++) {
            
         if(col % 3 == 0) {
            cout << "|";
         }
         cout << " " << to_string(solution[var_index]+1);
         if(col % 3 == 2) {
            cout << " ";
         }
         var_index += 1;
      
      }

      cout << "|" << endl;
      if(row % 3 == 2) {
         cout << "-------------------------" << endl;
      }

   }
}

int main() {

    // grid data
   std::vector<std::vector<Value> > input_grid {   {5, 3, 0, 0, 7, 0, 0, 0, 0},
                                                   {6, 0, 0, 1, 9, 5, 0, 0, 0},
                                                   {0, 9, 8, 0, 0, 0, 0, 6, 0},
                                                   {8, 0, 0, 0, 6, 0, 0, 0, 3},
                                                   {4, 0, 0, 8, 0, 3, 0, 0, 1},
                                                   {7, 0, 0, 0, 2, 0, 0, 0, 6},
                                                   {0, 6, 0, 0, 0, 0, 2, 8, 0},
                                                   {0, 0, 0, 4, 1, 9, 0, 0, 5},
                                                   {0, 0, 0, 0, 8, 0, 0, 7, 9} };

   // initialisation
   tb2init();
   initCosts();

   Cost top = 1;

   // creation of the solver object
   WeightedCSPSolver* solver = WeightedCSPSolver::makeWeightedCSPSolver(top);

   // access to the wcsp object created by the solver
   WeightedCSP* wcsp = solver->getWCSP();

   //-----------------------------------
   // problem variables
   //-----------------------------------

   // variable creation
   for(size_t row = 0; row < 9; row ++) {
      for(size_t col = 0; col < 9; col ++) {
         wcsp->makeEnumeratedVariable("Cell_" + to_string(row) + "," + to_string(col), 0, 8);
      }
   }

   // input values initialisation
   size_t var_ind = 0;
   for(size_t row = 0; row < 9; row ++) {
      for(size_t col = 0; col < 9; col ++) {
         
         if(input_grid[row][col] != 0) {
            wcsp->assign(var_ind, input_grid[row][col]-1);
         }

         var_ind ++;
      }
   }

   //-----------------------------------
   // constraints definition
   //-----------------------------------

   // all different constraint parameters
   std::string semantics = "hard";
   std::string prop = "knapsack";

   // add one "all different" constraint for each row
   for(int row = 0; row < 9; row ++) {
      std::vector<int> row_scope;
      for(int col = 0; col < 9; col ++) {
         row_scope.emplace_back(row*9+col);
      }
      wcsp->postWAllDiff(row_scope, semantics, prop, top);
   }

   // add one "all different" constraint for each column
   for(int col = 0; col < 9; col ++) {
      std::vector<int> col_scope;
      for(int row = 0; row < 9; row ++) {
         col_scope.emplace_back(row*9+col);
      }
      wcsp->postWAllDiff(col_scope, semantics, prop, top);
   }

   // add one "all different" constraint for each 3x3 sub square
   for(int sub_ind1 = 0; sub_ind1 < 3; sub_ind1 ++) {
      for(int sub_ind2 = 0; sub_ind2 < 3; sub_ind2 ++) {
         std::vector<int> sub_scope;
         for(int row_ind = 0; row_ind < 3; row_ind ++) { // iterate inside the 3x3 sub grid
            for(int col_ind = 0; col_ind < 3; col_ind ++) {
               sub_scope.emplace_back((sub_ind1*3+row_ind)*9+sub_ind2*3+col_ind);
            }
         }
         wcsp->postWAllDiff(sub_scope, semantics, prop, top);
      }
   }

   //-----------------------------------
   // solution
   //-----------------------------------

   // close the model definition
   wcsp->sortConstraints();

   // solve the problem
   bool hasSolution = solver->solve();

   if(hasSolution) {
      std::vector<Value> solution = solver->getSolution();
      printSolution(solution);
   }

   // clean the solver object
   delete solver;

   return 0;
}
