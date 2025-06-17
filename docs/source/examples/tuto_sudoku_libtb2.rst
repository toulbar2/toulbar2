.. _tuto_sudoku_libtb2:

===========================================
|cpp_logo| Sudoku puzzle with :ref:`libtb2 <ref_cpp>` in C++
===========================================

.. |cpp_logo| image::  https://raw.githubusercontent.com/Benio101/cpp-logo/master/cpp_logo.png
   :width: 40

.. include:: menu_backto.rst

The Sudoku is a widely known puzzle game that consists in filling a grid with numbers.
The typical grid has 9x9 cells in total and each cell must contain an integer between 1 and 9.

Additional constraints need to be enforced to solve the puzzle.
In each line, the same number cannot appear twice.
The same type of constraint occurs for each column.
Finally, each sub-square of size 3x3 located every 3 cells must also not contain duplicates.
The image below shows an example of a Sudoku grid given with its initial values.
The goal is to deduce the other values while verifying the different constraints.

.. only:: html

   .. figure:: ../../../web/IMAGES/sudoku_cst.svg
      :height: 240px
      :align: center
      :alt: A Sudoku grid with its constraints.

      A Sudoku grid partially filled where some of its constraints have been highlighted. Grid source: `Wikipedia <https://en.wikipedia.org/wiki/Sudoku/>`_

.. only:: latex

   .. figure:: ../../../web/IMAGES/sudoku_cst.pdf
      :height: 240px
      :align: center
      :alt: A Sudoku grid with its constraints.

      A Sudoku grid partially filled where some of its constraints have been highlighted. Grid source: `Wikipedia <https://en.wikipedia.org/wiki/Sudoku/>`_



Getting started
=================

Before starting, make sure the ToulBar2 C++ library binaries are installed in your system (:code:`libtb2.so`, see installation section from :ref:`sources <_README_5>` or :ref:`binaries <install-binaries>` for more instructions).

We first create a `WeightedCSPSolver <WCSPSolverClass_>`_ object.
This object is in charge of executing the algorithm to solve the Sudoku grid and internally creates a `WeightedCSP <WCSPClass_>`_ object to store the problem.

.. _WCSPSolverClass: ../ref/ref_cpp.html#weightedcspsolver-class
.. _WCSPClass: ../ref/ref_cpp.html#weightedcsp-class


:code:`WeightedCSP` objects are used in ToulBar2 to define the optimization or decision problems.
The problem is expressed as a set of discrete variables that are connected to each other through cost functions (or constraints).

As ToulBar2 is an optimization framework, an optional upper bound can be provided to the solver object in order to exclude any solution whose value exceeds this bound.
In the case of a Sudoku puzzle, since the problem does not contain a numerical objective, an upper bound of 1 can be chosen.

In order to compile the following code, assuming we are in the main ToulBar2 source repository, the same compilation flags as used to compile :code:`libtb2.so` must be used:
:code:`g++ -DBOOST -DLONGDOUBLE_PROB -DLONGLONG_COST -I./src -o sudoku sudoku_tutorial.cpp libtb2.so

.. highlight:: c++
   
.. code-block:: c++

   #include <iostream>
   #include <toulbar2lib.hpp>

   using namespace std;

   int main() {

      // initialization
      tb2init();
      initCosts();

      Cost top = 1;
      
      // creation of the solver object
      WeightedCSPSolver* solver = WeightedCSPSolver::makeWeightedCSPSolver(top);

      // access to the wcsp object created by the solver
      WeightedCSP* wcsp = solver->getWCSP();

      delete solver;

      return 0;
   }

Representing the grid in ToulBar2
==================================

To represent our problem in pytoulbar2, it is necessary to define discrete decision variables.
The variables will represent the various choices that can be made to build a solution to the problem.
In the Sudoku puzzle, decision variables are typically the different cells of the grid.
Their values would be the possible integers they can be assigned to, from 1 to 9.
We use the makeEnumeratedVariable function of the wcsp object to make the variables.

.. code-block:: c++

   // variable creation
   for(size_t row = 0; row < 9; row ++) {
      for(size_t col = 0; col < 9; col ++) {
         wcsp->makeEnumeratedVariable("Cell_" + to_string(row) + "," + to_string(col), 0, 8);
      }
   }

Solving first the grid
========================

It is already possible to solve the puzzle with ToulBar2, as the cfn object contains the variables of the problem.
The `WeightedCSPSolver::solve <solveFunc_>`_ function is used to run the solving algorithm.

.. _solveFunc: ../ref/ref_cpp.html#_CPPv4N17WeightedCSPSolver5solveEb

.. code-block:: c++

   // close the model definition
   wcsp->sortConstraints();

   // solve the problem
   bool hasSolution = solver->solve();


.. warning::
   It is important to systematically call the function `WeightedCSP::sortConstraints <sortCstFunc_>`_ before solving the problem to close the model definition.

.. _sortCstFunc: ../ref/ref_cpp.html#_CPPv4N11WeightedCSP15sortConstraintsEv

Problems may sometimes be especially hard to solve.
In such cases, a timeout can be set on the solver to stop the search before the optimality proof is complete, or before a solution is found at all.
When doing so, ToulBar2 throws an exception that must be caught to properly clean the problem data structures:

.. code-block:: c++

   try {
      bool hasSolution = solver->solve();
   } catch(const exception& ex) {
      cout << "no solution found: " << ex.what() << endl;
   }

When ToulBar2 returns a solution, the solution can be accessed as a std::vector of Value, specifying the value that is assigned to each variable of the problem (in the same order they were defined):

.. code-block:: c++

   if(hasSolution) {
        std::vector<Value> solution = solver->getSolution();
        cout << "the first value is " << solution[0] << endl;
    }

The above code would display :code:`0`, meaning the chosen value of the first variable (upper left cell) is 1 (first value of the list).
We then define a function to print the solution as a grid :

.. code-block:: c++

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

This function helps to visualize the variables' values as a real Sudoku grid, as follows :

.. code-block:: c++

   // print the first solution
   if(hasSolution) {
        std::vector<Value> solution = solver->getSolution();
        printSolution(solution)
    }

Which helps to display the first solution:

.. code-block:: text

   -------------------------
   | 1 1 1 | 1 1 1 | 1 1 1 |
   | 1 1 1 | 1 1 1 | 1 1 1 |
   | 1 1 1 | 1 1 1 | 1 1 1 |
   -------------------------
   | 1 1 1 | 1 1 1 | 1 1 1 |
   | 1 1 1 | 1 1 1 | 1 1 1 |
   | 1 1 1 | 1 1 1 | 1 1 1 |
   -------------------------
   | 1 1 1 | 1 1 1 | 1 1 1 |
   | 1 1 1 | 1 1 1 | 1 1 1 |
   | 1 1 1 | 1 1 1 | 1 1 1 |
   -------------------------


Assignment to the input values
===============================

The next step consists in initializing the variables that correspond to cells for which the value is known.
We will use the values in the grid example above, defined as a two-dimensional vector (where 0 means the value is unspecified):

.. code-block:: c++

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

Variables can be assigned with the `WeightedCSP::assign <assignFunc_>`_ function.
The variable and its value can be specified as integer indexes or as strings.

.. _assignFunc: ../ref/ref_cpp.html#_CPPv4N11WeightedCSP6assignEi5Value


.. code-block:: c++

   // input values initialization
   size_t var_ind = 0;
   for(size_t row = 0; row < 9; row ++) {
      for(size_t col = 0; col < 9; col ++) {
            
         if(input_grid[row][col] != 0) {
            wcsp->assign(var_ind, input_grid[row][col]-1);
         }

         var_ind ++;

      }
   }

.. warning::
   Although we already solved the problem once, a :code:`WeightedCSP` object cannot execute its :code:`WeightedCSPSolver::solve` function twice in a row. The object must be recreated or the function must be called only once.

The solution returned by the algorithm this time looks like this:

.. code-block:: text

   -------------------------
   | 5 3 1 | 1 7 1 | 1 1 1 |
   | 6 1 1 | 1 9 5 | 1 1 1 |
   | 1 9 8 | 1 1 1 | 1 6 1 |
   -------------------------
   | 8 1 1 | 1 6 1 | 1 1 3 |
   | 4 1 1 | 8 1 3 | 1 1 1 |
   | 7 1 1 | 1 2 1 | 1 1 6 |
   -------------------------
   | 1 6 1 | 1 1 1 | 2 8 1 |
   | 1 1 1 | 4 1 9 | 1 1 5 |
   | 1 1 1 | 1 8 1 | 1 7 9 |
   -------------------------

The initial values are now correctly specified in the variables.

Adding constraints and solving the grid
=========================================

The missing part to be able to generate a solution is the constraints.
Starting with the row constraints, we must ensure that none of the variables in the same row will be assigned to the same values.
This constraint is usually called *all different* and can be added with the `WeightedCSPSolver::postWAllDiff <AllDiffFunc_>`_ function.
The function takes as arguments a list of indices of the variables that must differ (the scope of the constraint) as well as two parameters specifying how the constraint is encoded, which we do not further describe in this tutorial.
We start by adding a constraint for the first row:

.. _AllDiffFunc: ../ref/ref_cpp.html#_CPPv4N11WeightedCSP12postWAllDiffE6vectorIiERK6stringRK6string4Cost


.. code-block:: c++
   
   std::vector<int> scope = {0,1,2,3,4,5,6,7,8};
   std::string semantics = "hard";
   std::string prop = "knapsack";

   wcsp->postWAllDiff(scope, semantics, prop, top);

Which generates the following first row in the solution :

.. code-block:: text

   -------------------------
   | 5 3 1 | 2 7 4 | 6 8 9 |

Constraints for each row can be added by varying the column index for each row :

.. code-block:: c++

   // add one "all different" constraint for each row
   for(int row = 0; row < 9; row ++) {
      std::vector<int> row_scope;
      for(int col = 0; col < 9; col ++) {
         row_scope.emplace_back(row*9+col);
      }
      wcsp->postWAllDiff(row_scope, semantics, prop, top);
   }

Constraints for each column are obtained similarly:

.. code-block:: c++

   // add one "all different" constraint for each column
   for(int col = 0; col < 9; col ++) {
      std::vector<int> col_scope;
      for(int row = 0; row < 9; row ++) {
         col_scope.emplace_back(row*9+col);
      }
      wcsp->postWAllDiff(col_scope, semantics, prop, top);
   }

At this point, the solution is not correct yet since sub-grids of size 3x3 may contain duplicates, such as the values :code:`9` and :code:`3` in the example below:

.. code-block:: text

   ---------
   | 9 7 6 |
   | 1 9 5 |
   | 5 4 2 |
   ---------

A set of 9 additional :code:`allDifferent` constraints can be defined to finalize our Sudoku model definition:

.. code-block:: c++

   // add one "all different" constraint for each 3x3 sub-grid
   for(int sub_ind1 = 0; sub_ind1 < 3; sub_ind1 ++) {
      for(int sub_ind2 = 0; sub_ind2 < 3; sub_ind2 ++) {
         std::vector<int> sub_scope;
         for(int row_ind = 0; row_ind < 3; row_ind ++) { // iterate inside the 3x3 sub-grid
            for(int col_ind = 0; col_ind < 3; col_ind ++) {
               sub_scope.emplace_back((sub_ind1*3+row_ind)*9+sub_ind2*3+col_ind);
            }
         }
         wcsp->postWAllDiff(sub_scope, semantics, prop, top);
      }
   }

These last constraints allow to finally obtain a consistent solution to the Sudoku puzzle :

.. code-block:: text

   -------------------------
   | 5 3 4 | 6 7 8 | 9 1 2 |
   | 6 7 2 | 1 9 5 | 3 4 8 |
   | 1 9 8 | 3 4 2 | 5 6 7 |
   -------------------------
   | 8 5 9 | 7 6 1 | 4 2 3 |
   | 4 2 6 | 8 5 3 | 7 9 1 |
   | 7 1 3 | 9 2 4 | 8 5 6 |
   -------------------------
   | 9 6 1 | 5 3 7 | 2 8 4 |
   | 2 8 7 | 4 1 9 | 6 3 5 |
   | 3 4 5 | 2 8 6 | 1 7 9 |
   -------------------------

Conclusion
=================

This short introduction shows how to represent a problem in ToulBar2 via its C++ library :ref:`libtb2 <ref_cpp>` by defining the problem variables and constraints and how to obtain a solution to this problem.
Below is the complete C++ source code from this tutorial. 

:download:`sudoku_tutorial.cpp<../../../web/TUTORIALS/sudoku_tutorial.cpp>`

.. literalinclude:: ../../../web/TUTORIALS/sudoku_tutorial.cpp
