
#include "toulbar2lib.hpp"

#include "bicriteria.hpp"

#include <algorithm>

using namespace std;

using namespace mulcrit;

//--------------------------------------------------------------------------------------------
vector<Bicriteria::Point> Bicriteria::_points = vector<Point>();
vector<mulcrit::Solution> Bicriteria::_solutions = vector<mulcrit::Solution>();

//--------------------------------------------------------------------------------------------
bool Bicriteria::notEqual(Bicriteria::Point p1, Bicriteria::Point p2) {
  return fabs(p1.first-p2.first) >= MultiWCSP::epsilon || fabs(p1.second-p2.second) >= MultiWCSP::epsilon;
}

//--------------------------------------------------------------------------------------------
bool Bicriteria::equal(Point p1, Point p2) {
  return fabs(p1.first-p2.first) <= MultiWCSP::epsilon && fabs(p1.second-p2.second) <= MultiWCSP::epsilon;
}

//--------------------------------------------------------------------------------------------
bool Bicriteria::solveScalarization(MultiWCSP* multiwcsp, pair<Double,Double> weights, Solution* solution, Bicriteria::Point* point) {

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

  bool result = solver->solve();

  if(result) {
    // cout << "solution found" << endl;
    // combiner.getSolution(solver, &sol_values, solution);
    sol_values = multiwcsp->getSolutionValues();

    if(solution != nullptr) {
      *solution = multiwcsp->getSolution();
    }

    if(point != nullptr) {
      *point = make_pair(sol_values[0], sol_values[1]);
    }

    // solver->getMultiWCSPSolution(sol_values);
  } else {
    // cout << "no solution !" << endl;
  }

  delete solver;
  delete wcsp; 

  return result;

}

//--------------------------------------------------------------------------------------------
void Bicriteria::computeSupportedPoints(mulcrit::MultiWCSP* multiwcsp, pair<Bicriteria::OptimDir, Bicriteria::OptimDir> optim_dir) {

  _solutions.clear();
  _points.clear();

  tb2init();

  ToulBar2::verbose = -1;

  ToulBar2::cfn = true;

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

  // delta is compared to an internal decimalPoint because ToulBar2::decimalPoint is only initialized when created the final wcsp (makeMultiWCSP)
  if(log10(Bicriteria::delta) > multiwcsp->getDecimalPoint()) {
    cerr << "Error: delta constant (" << Bicriteria::delta << ") is incompatible with decimalPoint (" << ToulBar2::decimalPoint << ")" << endl; 
  }

  bool result1, result2;

  // cout << "optimizing 1 separately: " << endl;
  if(optim_dir.second == Optim_Min) {
    result1 = solveScalarization(multiwcsp, make_pair(lambda1,-Bicriteria::delta), &sol1, &point1);
  } else {
    result1 = solveScalarization(multiwcsp, make_pair(lambda1,Bicriteria::delta), &sol1, &point1);
  }
  // cout << "Optimal point for 1: " << point1.first << ";" << point1.second << endl;

  // cout << "optimizing 2 separately: " << endl;
  if(optim_dir.first == Optim_Min) {
    result2 = solveScalarization(multiwcsp, make_pair(-Bicriteria::delta,lambda2), &sol2, &point2);    
  } else {
    result2 = solveScalarization(multiwcsp, make_pair(Bicriteria::delta,lambda2), &sol2, &point2);    
  }
  // cout << "Optimal point for 2: " << point2.first << ";" << point2.second << endl;
  // cout << endl << endl;

  stack<pair<Point, Point>> pending;

  if(result1) { // make sure there is a solution
    _points.push_back(point1);
    _solutions.push_back(sol1);
  }
  
  if(result2 && notEqual(point1, point2)) {
    pending.push(make_pair(point1, point2));
    _points.push_back(point2);
    _solutions.push_back(sol2);
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


    // cout << "weights : " << lambda1 << ", " << lambda2 << endl;
    Point new_point;
    bool result = solveScalarization(multiwcsp, make_pair(lambda1,lambda2), &new_sol, &new_point);

    // jump to the next weights if there is no solution
    if(!result) {
      continue;
    }

    // cout << "new point: " << new_point.first << ", " << new_point.second << endl;
    // cout << "from " << top.first.first << "," << top.first.second << " and " << top.second.first << ", " << top.second.second << endl;  

    // cut the search if the point was alrady encountered
    auto it = std::find_if(_points.begin(), _points.end(), [new_point](Point& point) { return equal(point, new_point); });
    if(it != _points.end()) {
      continue;
    }

    if(notEqual(new_point, top.first) && notEqual(new_point, top.second)) {

      _points.push_back(new_point);
      _solutions.push_back(new_sol);

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

  // cout << "Supported points: " << endl;
  // for(unsigned int point_ind = 0; point_ind < sol_points.size(); point_ind ++) {
  //   cout << sol_points[point_ind].first << ", " << sol_points[point_ind].second << endl;
  // }

  // the points are sorted from left to right
  vector<unsigned int> sol_indexes(_points.size());
  for(unsigned int ind = 0; ind < sol_indexes.size(); ind ++) {
    sol_indexes[ind] = ind;
  }

  if(optim_dir.first == optim_dir.second) {
    sort(sol_indexes.begin(), sol_indexes.end(), [](unsigned int& ind1, unsigned int& ind2) { return (_points[ind1].first < _points[ind2].first) || (_points[ind1].first == _points[ind2].first && _points[ind1].second > _points[ind2].second); } );
  } else {
    sort(sol_indexes.begin(), sol_indexes.end(), [](unsigned int& ind1, unsigned int& ind2) { return (_points[ind1].first < _points[ind2].first) || (_points[ind1].first == _points[ind2].first && _points[ind1].second < _points[ind2].second); } );
  }

  vector<Point> temp_points = _points;
  for(unsigned int ind = 0; ind < _points.size(); ind ++) {
    _points[ind] = temp_points[sol_indexes[ind]];
  }

  vector<mulcrit::Solution> temp_sol = _solutions;
  for(unsigned int ind = 0; ind < _solutions.size(); ind ++) {
    _solutions[ind] = temp_sol[sol_indexes[ind]];
  }



}

//--------------------------------------------------------------------------------------------
std::vector<mulcrit::Solution> Bicriteria::getSolutions() {
  return _solutions;
}

//--------------------------------------------------------------------------------------------
std::vector<Bicriteria::Point> Bicriteria::getPoints() {
  return _points;
}