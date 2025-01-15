#include <iostream>
#include <vector>

#include "core/tb2wcsp.hpp"

using namespace std;

// an alias for storing the variable costs
// first dim is the grid rows and second is the columns
typedef std::vector<std::vector<std::vector<Cost>>> LatinCostArray;

/*!
    \brief generate random costs for each variable (cell)
 */
void initLatinCosts(size_t N, LatinCostArray& costs) {

    // N*N*N values, costs for each cell
    costs.resize(N);
    for(auto& col: costs) {
        col.resize(N);
        for(auto& cell: col) {
            cell.resize(N);
            for(size_t val_ind = 0; val_ind < N; val_ind += 1) {
                cell[val_ind] = (rand()%N)+1;
            }
        }
    }
}

/*!
    \brief print the costs for each unary variabl (cell)
 */
void printCosts(LatinCostArray& costs) {
    
    for(size_t row_ind = 0; row_ind < costs.size(); row_ind ++) {
        for(size_t col_ind = 0; col_ind < costs[row_ind].size(); col_ind ++) {
            cout << "cell " << row_ind << "_" << col_ind;
            cout << " : ";
            for(auto& cost: costs[row_ind][col_ind]) {
                cout << cost << ", ";
            }
            cout << endl;
        }
    }
}

/*!
    \brief fill in a WCSP object with a latin square problem
 */
void buildWCSP(WeightedCSP& wcsp, LatinCostArray& costs, size_t N, Cost top) {

    // variables
    for(size_t row = 0; row < N; row ++) {
        for(unsigned int col = 0; col < N; col ++) {
            wcsp.makeEnumeratedVariable("Cell_" + to_string(row) + "," + to_string(col), 0, N-1);
        }
    }

    cout << "number of variables: " << wcsp.numberOfVariables() << endl;

    /* costs for all different constraints (top on diagonal) */
    vector<Cost> alldiff_costs;
    for(unsigned int i = 0; i < N; i ++) {
        for(unsigned int j = 0; j < N; j ++) {
            if(i == j) {
                alldiff_costs.push_back(top);
            } else {
                alldiff_costs.push_back(0);
            }
        }
    }

    /* all different constraints */
    for(unsigned int index = 0; index < N; index ++) {
        for(unsigned int var_ind1 = 0; var_ind1 < N; var_ind1 ++) {
            for(unsigned int var_ind2 = var_ind1+1; var_ind2 < N; var_ind2 ++) {
                /* row constraints */
                wcsp.postBinaryConstraint(N*index+var_ind1, N*index+var_ind2, alldiff_costs);
                /* col constraints */
                wcsp.postBinaryConstraint(index+var_ind1*N, index+var_ind2*N, alldiff_costs);
            }
        }
    }
  
    /* unary costs */
    size_t var_ind = 0;
    for(size_t row = 0; row < N; row ++) {
        for(size_t col = 0; col < N; col ++) {
            wcsp.postUnaryConstraint(var_ind, costs[row][col]);
            var_ind += 1;
        }
    }

}

int main() {

    srand(123456789);

    size_t N = 5;
    Cost top = N*N*N + 1;

    // N*N*N values, costs for each cell
    LatinCostArray objective_costs;
    
    // init the costs for each cell
    initLatinCosts(N, objective_costs);

    cout << "Randomly genereated costs : " << endl;
    printCosts(objective_costs);
    cout << endl;

    tb2init();

    ToulBar2::verbose = 0;

    WeightedCSPSolver* solver = WeightedCSPSolver::makeWeightedCSPSolver(top);

    // fill in the WeightedCSP object
    WeightedCSP* wcsp = solver->getWCSP();

    buildWCSP(*wcsp, objective_costs, N, top);

    bool result = solver->solve();

    if(result) {

        Cost bestCost = solver->getSolutionValue();
        Cost bestLowerBound = solver->getDDualBound();

        if(!ToulBar2::limited) {
            cout << "Optimal solution found with cost " << bestCost << endl;
        } else {
            cout << "Best solution found with cost " << bestCost << " and best lower bound of " << bestLowerBound << endl;
        }

        // retrieve the solution
        std::vector<Value> solution = solver->getSolution();


        cout << endl << "Best solution : " << endl;
        for(size_t var_ind = 0; var_ind < solution.size(); var_ind ++) {
            cout << solution[var_ind] << " ; ";
            if((var_ind+1) % N == 0) {
                cout << endl;
            }
        }

    } else {
        cout << "No solution has been found !" << endl;
    }

    delete solver;     

    return 0;
}