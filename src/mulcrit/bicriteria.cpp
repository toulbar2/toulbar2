
#include "toulbar2lib.hpp"

#include "bicriteria.hpp"

#include <algorithm>

using namespace std;

using namespace mulcrit;

//--------------------------------------------------------------------------------------------
bool Bicriteria::notEqual(Bicriteria::Point p1, Bicriteria::Point p2) {
  return fabs(p1.first-p2.first) >= MultiWCSP::epsilon || fabs(p1.second-p2.second) >= MultiWCSP::epsilon;
}

//--------------------------------------------------------------------------------------------
bool Bicriteria::equal(Point p1, Point p2) {
  return fabs(p1.first-p2.first) <= MultiWCSP::epsilon && fabs(p1.second-p2.second) <= MultiWCSP::epsilon;
}

//--------------------------------------------------------------------------------------------
Bicriteria::Point Bicriteria::solveScalarization(MultiWCSP* multiwcsp, pair<Double,Double> weights, Solution* solution) {

  cout << "current weights: " << weights.first << ", " << weights.second << endl;

  multiwcsp->setWeight(0, weights.first);
  multiwcsp->setWeight(1, weights.second);

  tb2init();

  // WeightedCSPSolver* solver = WeightedCSPSolver::makeWeightedCSPSolver(MAX_COST);

  // WCSP* pb = dynamic_cast<WCSP*>(solver->getWCSP());
  // combiner.exportToWCSP(pb);

  /* debug */
  // ofstream file("latin_combined.cfn");
  // pb->dump_CFN(file);
  // file.close();

  WeightedCSP* wcsp = multiwcsp->makeWeightedCSP();

  cout << "n variables in the final wcsp: " << wcsp->numberOfVariables() << ", " << wcsp->numberOfConstraints() << endl;

  WeightedCSPSolver* solver = WeightedCSPSolver::makeWeightedCSPSolver(MAX_COST, wcsp);
  solver->getWCSP()->sortConstraints();

  vector<Double> sol_values;

  if(solver->solve()) {
    cout << "solution found" << endl;
    // combiner.getSolution(solver, &sol_values, solution);
    sol_values = multiwcsp->getSolutionValues();

    // solver->getMultiWCSPSolution(sol_values);
  } else {
    cout << "no solution !" << endl;
  }

  delete solver;
  delete wcsp;

  return make_pair(sol_values[0], sol_values[1]); 

}

//--------------------------------------------------------------------------------------------
void Bicriteria::computeSupportedPoints(mulcrit::MultiWCSP* multiwcsp, pair<Bicriteria::OptimDir, Bicriteria::OptimDir> optim_dir, std::vector<Point>* supported_points, vector<Solution>* solutions) {

  tb2init();

  ToulBar2::verbose = -1;

  ToulBar2::cfn = true;

  vector<Point> sol_points;
  
  Point point1, point2;

  Double lambda1 = 1., lambda2 = 1.;
  if(optim_dir.first == Optim_Max) {
    lambda1 = -1.;
  }
  if(optim_dir.second == Optim_Max) {
    lambda2 = -1.;
  }

  Solution sol1, sol2;

  // cout << endl << endl;
  // cout << "optimizing 1 separately: " << endl;
  // // point1 = solve_scalarization(multiwcsp, make_pair(lambda1,-0.01), &sol1);
  // point1 = solve_scalarization(multiwcsp, make_pair(lambda1,0), &sol1);
  // cout << "Optimal point for 1: " << point1.first << ";" << point1.second << endl;

  // cout << endl << endl;
  // cout << "optimizing 2 separately: " << endl;
  // // point2 = solve_scalarization(multiwcsp, make_pair(-0.01,lambda2), &sol2);
  // point2 = solve_scalarization(multiwcsp, make_pair(0,lambda2), &sol2);
  // cout << "Optimal point for 2: " << point2.first << ";" << point2.second << endl;
  // cout << endl << endl;

  if(log10(Bicriteria::delta) > ToulBar2::decimalPoint) {
    cerr << "Error: delta constant (" << Bicriteria::delta << ") is incompatible with decimalPoint (" << ToulBar2::decimalPoint << ")" << endl; 
  }

  cout << "optimizing 1 separately: " << endl;
  if(optim_dir.second == Optim_Min) {
    point1 = solveScalarization(multiwcsp, make_pair(lambda1,-Bicriteria::delta), &sol1);
  } else {
    point1 = solveScalarization(multiwcsp, make_pair(lambda1,Bicriteria::delta), &sol1);
  }
  cout << "Optimal point for 1: " << point1.first << ";" << point1.second << endl;

  cout << "optimizing 2 separately: " << endl;
  if(optim_dir.first == Optim_Min) {
    point2 = solveScalarization(multiwcsp, make_pair(-Bicriteria::delta,lambda2), &sol2);    
  } else {
    point2 = solveScalarization(multiwcsp, make_pair(Bicriteria::delta,lambda2), &sol2);    
  }
  cout << "Optimal point for 2: " << point2.first << ";" << point2.second << endl;
  cout << endl << endl;

  stack<pair<Point, Point>> pending;

  sol_points.push_back(point1);
  if(solutions != nullptr) {
    solutions->push_back(sol1);
  }
  if(notEqual(point1, point2)) {
    pending.push(make_pair(point1, point2));
    sol_points.push_back(point2);
    if(solutions != nullptr) {
      solutions->push_back(sol2);
    }
  }

  unsigned int iter = 0;

  mulcrit::Solution new_sol;

  while(!pending.empty()) {

    iter ++;

    pair<Point,Point> top = pending.top();
    pending.pop();

    // compute the new weights

    if(optim_dir.first == optim_dir.second) {
      lambda1 = -top.second.second + top.first.second;
      lambda2 = -top.first.first + top.second.first;
    } else {
      lambda1 = top.second.second - top.first.second;
      lambda2 = top.first.first - top.second.first;
    }


    cout << "weights : " << lambda1 << ", " << lambda2 << endl;

    auto new_point = solveScalarization(multiwcsp, make_pair(lambda1,lambda2), &new_sol);


    cout << "new point: " << new_point.first << ", " << new_point.second << endl;
    cout << "from " << top.first.first << "," << top.first.second << " and " << top.second.first << ", " << top.second.second << endl;  

    // cut the search if the point was alrady encountered
    auto it = std::find_if(sol_points.begin(), sol_points.end(), [new_point](Point& point) { return equal(point, new_point); });
    if(it != sol_points.end()) {
      continue;
    }

    if(notEqual(new_point, top.first) && notEqual(new_point, top.second)) {

      sol_points.push_back(new_point);
      if(solutions != nullptr) {
        solutions->push_back(new_sol);
      }

      // add two new scalarizations
      if(fabs(top.first.first - new_point.first) >= MultiWCSP::epsilon && fabs(top.first.second - new_point.second) >= MultiWCSP::epsilon) {
        pending.push(make_pair(top.first, new_point));
      }
      if(fabs(top.second.first - new_point.first) >= MultiWCSP::epsilon && fabs(top.second.second - new_point.second) >= MultiWCSP::epsilon) {
        pending.push(make_pair(new_point, top.second));
      }

    }

    // cout << endl << endl << "debug" << endl << endl << endl;

    // cout << "points found so far: ";
    // for(auto& point: sol_points) {
    //   cout << "(" << point.first << "," << point.second << "), ";
    // }
    // cout << endl;

  }

  cout << "Supported points: " << endl;
  for(unsigned int point_ind = 0; point_ind < sol_points.size(); point_ind ++) {
    cout << sol_points[point_ind].first << ", " << sol_points[point_ind].second << endl;
  }

  if(supported_points != nullptr) {
    *supported_points = sol_points;
  

    // the points are sorted from left to right
    vector<unsigned int> sol_indexes(sol_points.size());
    for(unsigned int ind = 0; ind < sol_indexes.size(); ind ++) {
      sol_indexes[ind] = ind;
    }

    if(optim_dir.first == optim_dir.second) {
      sort(sol_indexes.begin(), sol_indexes.end(), [supported_points](unsigned int& ind1, unsigned int& ind2) { return ((*supported_points)[ind1].first < (*supported_points)[ind2].first) || ((*supported_points)[ind1].first == (*supported_points)[ind2].first && (*supported_points)[ind1].second > (*supported_points)[ind2].second); } );
    } else {
      sort(sol_indexes.begin(), sol_indexes.end(), [supported_points](unsigned int& ind1, unsigned int& ind2) { return ((*supported_points)[ind1].first < (*supported_points)[ind2].first) || ((*supported_points)[ind1].first == (*supported_points)[ind2].first && (*supported_points)[ind1].second < (*supported_points)[ind2].second); } );
    }

    vector<Point> temp_points = *supported_points;
    for(unsigned int ind = 0; ind < supported_points->size(); ind ++) {
      (*supported_points)[ind] = temp_points[sol_indexes[ind]];
    }

    if(solutions != nullptr) {
      vector<mulcrit::Solution> temp_sol = *solutions;
      for(unsigned int ind = 0; ind < solutions->size(); ind ++) {
        (*solutions)[ind] = temp_sol[sol_indexes[ind]];
      }
    }

  }


}