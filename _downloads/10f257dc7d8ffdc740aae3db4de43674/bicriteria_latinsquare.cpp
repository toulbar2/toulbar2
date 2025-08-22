#include <iostream>
#include <vector>

#include "core/tb2wcsp.hpp"
#include "mcriteria/multicfn.hpp"
#include "mcriteria/bicriteria.hpp"

using namespace std;

// an alias for storing the variable costs
// first dim is the grid rows and second is the columns
typedef std::vector<std::vector<std::vector<Cost>>> LatinCostArray;


// generate random costs for each variable (cell)
// param N grid size
// param costs the matrix costs
void createCostMatrix(size_t N, LatinCostArray& costs) {

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



// print the costs for each unary variabl (cell)
// param costs the cost matrix
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

// fill in a WCSP object with a latin square problem
// param wcsp the wcsp object to fill
// param LatinCostArray the cost matrix
// param N grid size
// top the top value, problem upper bound (the objective is always lower than top)
void buildLatinSquare(WeightedCSP& wcsp, LatinCostArray& costs, size_t N, Cost top) {

    // variables
    for(size_t row = 0; row < N; row ++) {
        for(unsigned int col = 0; col < N; col ++) {
            wcsp.makeEnumeratedVariable("Cell_" + to_string(row) + "," + to_string(col), 0, N-1);
        }
    }

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

// print a solution as a grid
// param N the size of the grid
// param solution the multicfn solution (dict)
// param point the objective costs (objective space point)
void printSolution(size_t N, MultiCFN::Solution& solution, Bicriteria::Point& point) {
    for(size_t row = 0; row < N; row ++) {
        for(size_t col = 0; col < N; col ++) {
            string var_name = "Cell_" + to_string(row) + "," + to_string(col);
            cout << solution[var_name].substr(1) << " ";
        }
        cout << endl;
    }
    cout << "obj_1 = " << point.first << " ; obj2 = " << point.second << endl;

}

// main function
int main() {

    srand(123456789);

    size_t N = 4;
    Cost top = N*N*N + 1;

    // two cost matrice
    LatinCostArray costs_obj1, costs_obj2;
    
    // init the objective with random costs
    createCostMatrix(N, costs_obj1);
    createCostMatrix(N, costs_obj2);

    // cout << "Randomly genereated costs : " << endl;
    // printCosts(costs_obj1);
    // cout << endl << endl;
    // printCosts(costs_obj2);

    tb2init();
    initCosts();

    // create the two wcsp objects
    WeightedCSP* wcsp1 = WeightedCSP::makeWeightedCSP(top);
    WeightedCSP* wcsp2 = WeightedCSP::makeWeightedCSP(top);

    // initialize the objects as a latin square problem objectives with two different objectves
    buildLatinSquare(*wcsp1, costs_obj1, N, top);
    buildLatinSquare(*wcsp2, costs_obj2, N, top);

    // creation of the multicfn
    MultiCFN mcfn;
    mcfn.push_back(dynamic_cast<WCSP*>(wcsp1));
    mcfn.push_back(dynamic_cast<WCSP*>(wcsp2));

    // computation iof the supported points of the biobjective problem
    Bicriteria::computeSupportedPoints(&mcfn, std::make_pair(Bicriteria::OptimDir::Optim_Min, Bicriteria::OptimDir::Optim_Min));

    // access to the computed solutions and their objective values
    std::vector<MultiCFN::Solution> solutions = Bicriteria::getSolutions();
    std::vector<Bicriteria::Point> points = Bicriteria::getPoints();

    // print all solutions computed
    cout << "Resulting solutions: " << endl;
    for(size_t sol_index = 0; sol_index < solutions.size(); sol_index ++) {
        printSolution(N, solutions[sol_index], points[sol_index]);
        cout << endl;
    }

    // delete the wcsp objects
    delete wcsp1;
    delete wcsp2;

    return 0;
}