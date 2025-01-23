
#include "toulbar2lib.hpp"

#include "bicriteria.hpp"

#include <algorithm>

#include <iomanip>

using namespace std;

//--------------------------------------------------------------------------------------------
vector<Bicriteria::Weights> Bicriteria::_weights = vector<Weights>();
vector<Bicriteria::Point> Bicriteria::_points = vector<Point>();
vector<MultiCFN::Solution> Bicriteria::_solutions = vector<MultiCFN::Solution>();
vector<Double> Bicriteria::_lower_bounds = vector<Double>();

unsigned int Bicriteria::_first_cfn_index = 0;
unsigned int Bicriteria::_second_cfn_index = 0;

// parameters default values
int Bicriteria::_showSolutions = 0;
int Bicriteria::_vac = 0;
int Bicriteria::_seed = 1;
int Bicriteria::_verbose = -1;
int Bicriteria::_sol_timeout = 0;
int Bicriteria::_global_timeout = 0;
int Bicriteria::_nsol_max = 0;
double Bicriteria::_startGlobalCpuTime = 0.;

//--------------------------------------------------------------------------------------------
void Bicriteria::sortSolutions(pair<OptimDir, OptimDir> optim_dir)
{

    // the points are sorted from left to right
    vector<unsigned int> sol_indexes(_points.size());
    for (unsigned int ind = 0; ind < sol_indexes.size(); ind++) {
        sol_indexes[ind] = ind;
    }

    if (optim_dir.first == optim_dir.second) {
        sort(sol_indexes.begin(), sol_indexes.end(), [](unsigned int& ind1, unsigned int& ind2) { return (_points[ind1].first < _points[ind2].first) || (_points[ind1].first == _points[ind2].first && _points[ind1].second > _points[ind2].second); });
    } else {
        sort(sol_indexes.begin(), sol_indexes.end(), [](unsigned int& ind1, unsigned int& ind2) { return (_points[ind1].first < _points[ind2].first) || (_points[ind1].first == _points[ind2].first && _points[ind1].second < _points[ind2].second); });
    }

    vector<Weights> temp_weights = _weights;
    for (unsigned int ind = 0; ind < _weights.size(); ind++) {
        _weights[ind] = temp_weights[sol_indexes[ind]];
    }

    vector<Point> temp_points = _points;
    for (unsigned int ind = 0; ind < _points.size(); ind++) {
        _points[ind] = temp_points[sol_indexes[ind]];
    }

    vector<Double> temp_lb = _lower_bounds;
    for (unsigned int ind = 0; ind < _lower_bounds.size(); ind++) {
        _lower_bounds[ind] = temp_lb[sol_indexes[ind]];
    }

    vector<MultiCFN::Solution> temp_sol = _solutions;
    for (unsigned int ind = 0; ind < _solutions.size(); ind++) {
        _solutions[ind] = temp_sol[sol_indexes[ind]];
    }
}

//--------------------------------------------------------------------------------------------
bool Bicriteria::notEqual(Bicriteria::Point p1, Bicriteria::Point p2)
{
    return fabs(p1.first - p2.first) >= ToulBar2::epsilon || fabs(p1.second - p2.second) >= ToulBar2::epsilon;
}

//--------------------------------------------------------------------------------------------
bool Bicriteria::equal(Point p1, Point p2)
{
    return fabs(p1.first - p2.first) <= ToulBar2::epsilon && fabs(p1.second - p2.second) <= ToulBar2::epsilon;
}

//--------------------------------------------------------------------------------------------
bool Bicriteria::dominates(Point p1, Point p2, pair<OptimDir, OptimDir> optim_dir)
{

    unsigned int obj_1 = 0;
    if ((optim_dir.first == Optim_Min && p1.first < p2.first) || (optim_dir.first == Optim_Max && p1.first > p2.first)) {
        obj_1 = 2;
    } else if (fabs(p1.first - p2.first) <= ToulBar2::epsilon) {
        obj_1 = 1;
    }

    unsigned int obj_2 = 0;
    if ((optim_dir.second == Optim_Min && p1.second < p2.second) || (optim_dir.second == Optim_Max && p1.second > p2.second)) {
        obj_2 = 2;
    } else if (fabs(p1.second - p2.second) <= ToulBar2::epsilon) {
        obj_2 = 1;
    }

    return obj_1 + obj_2 >= 4;
}

//--------------------------------------------------------------------------------------------
bool Bicriteria::solveScalarization(MultiCFN* multicfn, pair<Double, Double> weights, MultiCFN::Solution* solution, Bicriteria::Point* point, Double* lb)
{

    // cout << "current weights: " << std::setprecision(10) << weights.first << ", " << weights.second << endl;

    multicfn->setWeight(_first_cfn_index, weights.first);
    multicfn->setWeight(_second_cfn_index, weights.second);

    // /!\ this initialization restores the state of the solver but also disrupt the user options (such as verbosity level)
    tb2init();
    ToulBar2::vac = _vac;
    ToulBar2::seed = _seed;
    ToulBar2::verbose = _verbose;
    ToulBar2::showSolutions = _showSolutions;

    // WeightedCSPSolver* solver = WeightedCSPSolver::makeWeightedCSPSolver(MAX_COST);

    // WCSP* pb = dynamic_cast<WCSP*>(solver->getWCSP());
    // combiner.exportToWCSP(pb, {}, {}, {});

    /* debug */
    // ofstream file("latin_combined.cfn");
    // pb->dump_CFN(file);
    // file.close();

    WeightedCSP* wcsp = multicfn->makeWeightedCSP();

    // cout << "n variables in the final wcsp: " << wcsp->numberOfVariables() << ", " << wcsp->numberOfConstraints() << endl;

    WeightedCSPSolver* solver = WeightedCSPSolver::makeWeightedCSPSolver(MAX_COST, wcsp);
    solver->getWCSP()->sortConstraints();

    vector<Double> sol_values;

    // start the local or global timer if required
    if(Bicriteria::_sol_timeout > 0 || Bicriteria::_global_timeout > 0) {
        ToulBar2::startCpuTime = cpuTime();
        if (Bicriteria::_sol_timeout > 0 && Bicriteria::_global_timeout > 0) {
            timer(min(Bicriteria::_sol_timeout, max(1, Bicriteria::_global_timeout - static_cast<int>(ToulBar2::startCpuTime - _startGlobalCpuTime))));
        } else if (Bicriteria::_global_timeout > 0) {
            timer(max(1, Bicriteria::_global_timeout - static_cast<int>(ToulBar2::startCpuTime - _startGlobalCpuTime)));
        } else {
            timer(Bicriteria::_sol_timeout);
        }
    }

    bool result = solver->solve();

    // stop the timer
    if(Bicriteria::_sol_timeout > 0 || Bicriteria::_global_timeout > 0) {
        timerStop();
    }

    *lb = solver->getDDualBound();

    if (result) {
        // cout << "solution found" << endl;
        // combiner.getSolution(solver, &sol_values, solution);
        sol_values = multicfn->getSolutionValues();

        if (solution != nullptr) {
            *solution = multicfn->getSolution();
        }

        if (point != nullptr) {
            *point = make_pair(sol_values[0], sol_values[1]);
        }

        // solver->getMultiFNSolution(sol_values);
    } else {
        // cout << "no solution !" << endl;
    }

    delete solver;
    delete wcsp;

    return result;
}

//--------------------------------------------------------------------------------------------
void Bicriteria::computeAdditionalSolutions(MultiCFN* multicfn, pair<Bicriteria::OptimDir, Bicriteria::OptimDir> optim_dir, unsigned int solIndex, unsigned int nbLimit, Double pct)
{

    tb2init();
    ToulBar2::vac = _vac;
    ToulBar2::seed = _seed;
    ToulBar2::verbose = _verbose;
    ToulBar2::showSolutions = _showSolutions;
    ToulBar2::allSolutions = nbLimit;

    Weights weights = make_pair(fabs(_points[solIndex + 1].second - _points[solIndex].second), fabs(_points[solIndex + 1].first - _points[solIndex].first));

    if (optim_dir.first == Optim_Max && weights.first > 0) {
        weights.first = -weights.first;
    }

    if (optim_dir.second == Optim_Max && weights.second > 0) {
        weights.second = -weights.second;
    }

    // cout << "selected weights: " << weights.first << ", " << weights.second << endl;

    // compute the "corner" to determine a new ub
    Point corner;
    if (optim_dir.first == Optim_Min) {
        corner.first = max(_points[solIndex].first, _points[solIndex + 1].first);
    } else {
        corner.first = min(_points[solIndex].first, _points[solIndex + 1].first);
    }
    if (optim_dir.second == Optim_Min) {
        corner.second = max(_points[solIndex].second, _points[solIndex + 1].second);
    } else {
        corner.second = min(_points[solIndex].second, _points[solIndex + 1].second);
    }

    Double new_ub = corner.first * weights.first + corner.second * weights.second;

    Double new_lb = _points[solIndex].first * weights.first + _points[solIndex].second * weights.second;

    // lower the upper bound if the percentage is not one
    if (fabs(1. - pct) > ToulBar2::epsilon) {
        new_ub -= (new_ub - new_lb) * (1. - pct);
    }

    multicfn->setWeight(0, weights.first);
    multicfn->setWeight(1, weights.second);

    WeightedCSP* wcsp = multicfn->makeWeightedCSP();
    wcsp->setUb(wcsp->DoubletoCost(new_ub));

    WeightedCSPSolver* solver = WeightedCSPSolver::makeWeightedCSPSolver(MAX_COST, wcsp);
    solver->getWCSP()->sortConstraints();

    bool result = solver->solve();

    if (result) {

        vector<pair<Double, vector<Value>>> tb2_sol = solver->getSolutions();

        // cout << "number of solutions: " << tb2_sol.size() << endl;

        /* compute the solutions and points */
        vector<MultiCFN::Solution> sol;
        for (unsigned int ind = 0; ind < tb2_sol.size(); ind++) {
            sol.push_back(multicfn->convertToSolution(tb2_sol[ind].second));
        }

        vector<Point> points;
        for (unsigned int ind = 0; ind < sol.size(); ind++) {
            vector<Double> values = multicfn->computeSolutionValues(sol[ind]);
            points.push_back(make_pair(values[0], values[1]));
        }

        // build a list of indexes on the solutions
        vector<unsigned int> sol_indexes_temp(sol.size());
        for (unsigned int ind = 0; ind < sol.size(); ind++) {
            sol_indexes_temp[ind] = ind;
        }

        // filter the points by indexes
        vector<unsigned int> sol_indexes;
        for (unsigned int ind : sol_indexes_temp) {

            bool add = true;

            // check if the point correspond to one the the supported points
            if (Bicriteria::equal(points[ind], _points[solIndex])) {
                add = false;
            } else if (Bicriteria::equal(points[ind], _points[solIndex + 1])) {
                add = false;
            }

            // check if the point is dominated by one of the point given as index
            if (add && Bicriteria::dominates(points[solIndex], points[ind], optim_dir)) {
                add = false;
            } else if (Bicriteria::dominates(points[solIndex + 1], points[ind], optim_dir)) {
                add = false;
            }

            // make sure the point dominates the corner
            if (add && !Bicriteria::dominates(points[ind], corner, optim_dir)) {
                add = false;
            }

            // look for dominating points
            if (add && find_if(sol_indexes_temp.begin(), sol_indexes_temp.end(), [ind, points, optim_dir](unsigned int ind2) { return ind != ind2 && dominates(points[ind2], points[ind], optim_dir); }) != sol_indexes_temp.end()) {
                add = false;
            }

            if (add) {
                sol_indexes.push_back(ind);
            }
        }

        // non filtered solutions are added to the list
        vector<MultiCFN::Solution> temp_sol;
        vector<Point> temp_points;

        for (unsigned int index : sol_indexes) {
            temp_sol.push_back(sol[index]);
            temp_points.push_back(points[index]);
        }

        sol = temp_sol;
        points = temp_points;

        // cout << "n solutions: " << sol.size() << endl;

        /* solutions are added to the list, then sorted */
        for (unsigned int ind = 0; ind < sol.size(); ind++) {
            _solutions.push_back(sol[ind]);
            _points.push_back(points[ind]);
            _weights.push_back(make_pair(0., 0.));
        }

        sortSolutions(optim_dir);

    } else {
    }

    delete solver;
    delete wcsp;
}

//--------------------------------------------------------------------------------------------
void Bicriteria::computeNonSupported(MultiCFN* multicfn, pair<Bicriteria::OptimDir, Bicriteria::OptimDir> optim_dir, unsigned int nbLimit)
{

    // for all nondominated triangles, compute the list of triangles that will be solved together by solution enumeration

    // all triangle corners
    vector<Point> corners(_solutions.size() - 1);
    vector<bool> isTriangle(_solutions.size() - 1, true); // false if the triangle is flat
    for (unsigned int ind = 0; ind < corners.size(); ind++) {

        if (fabs(_points[ind].first - _points[ind + 1].first) <= ToulBar2::epsilon) {
            isTriangle[ind] = false;
            continue;
        } else if (fabs(_points[ind].second - _points[ind + 1].second) <= ToulBar2::epsilon) {
            isTriangle[ind] = false;
            continue;
        }

        if (optim_dir.first == Optim_Min) {
            corners[ind].first = max(_points[ind].first, _points[ind + 1].first);
        } else {
            corners[ind].first = min(_points[ind].first, _points[ind + 1].first);
        }
        if (optim_dir.second == Optim_Min) {
            corners[ind].second = max(_points[ind].second, _points[ind + 1].second);
        } else {
            corners[ind].second = min(_points[ind].second, _points[ind + 1].second);
        }
    }

    // list of the points of each triangle (excluding the corner)
    vector<pair<Point, Point>> trianglePoints(_solutions.size() - 1);
    for (unsigned int ind = 0; ind < corners.size(); ind++) {
        if (!isTriangle[ind]) {
            continue;
        }
        trianglePoints[ind] = make_pair(_points[ind], _points[ind + 1]);
    }

    // compute the triangle weights
    vector<Weights> triangleWeights(corners.size());
    for (unsigned int ind = 0; ind < triangleWeights.size(); ind++) {

        if (!isTriangle[ind]) {
            continue;
        }

        Weights weights = make_pair(fabs(_points[ind + 1].second - _points[ind].second), fabs(_points[ind + 1].first - _points[ind].first));

        if (optim_dir.first == Optim_Max && weights.first > 0) {
            weights.first = -weights.first;
        }
        if (optim_dir.second == Optim_Max && weights.second > 0) {
            weights.second = -weights.second;
        }

        triangleWeights[ind] = weights;
    }

    // determine a list of "dominated" triangles, i.e. triangles that will be solved by solving other ones
    vector<vector<int>> subTriangles(triangleWeights.size(), vector<int>());

    for (unsigned int ind = 0; ind < corners.size(); ind++) {

        if (!isTriangle[ind]) {
            continue;
        }

        Double cornerVal = triangleWeights[ind].first * corners[ind].first + triangleWeights[ind].second * corners[ind].second;

        for (unsigned int ind2 = 0; ind2 < corners.size(); ind2++) {
            if (ind2 == ind || !isTriangle[ind2]) {
                continue;
            }
            Double compVal = triangleWeights[ind].first * corners[ind2].first + triangleWeights[ind].second * corners[ind2].second;

            if (compVal < cornerVal) {
                subTriangles[ind].push_back(ind2);
            }
        }
    }

    // for(unsigned int ind = 0; ind < corners.size(); ind ++) {
    //   if(!isTriangle[ind]) {
    //     continue;
    //   }
    //   cout << "triangle dom" << ind << ": ";
    //   for(auto elt: subTriangles[ind]) {
    //     cout << elt << ", ";
    //   }
    //   cout << endl;
    // }

    // sort the triangles by the number of dominated ones
    vector<unsigned int> triangleInd;
    for (unsigned int indTri = 0; indTri < isTriangle.size(); indTri++) {
        if (isTriangle[indTri]) {
            triangleInd.push_back(indTri);
        }
    }

    sort(triangleInd.begin(), triangleInd.end(), [subTriangles](unsigned int ind1, unsigned int ind2) { return subTriangles[ind1].size() > subTriangles[ind2].size(); });

    // cout << "triangles order: " << endl;
    // for(unsigned int ind: triangleInd) {
    //   cout << ind << ", ";
    // }
    // cout << endl;

    // list of triangles processed
    vector<bool> processed(corners.size(), false);

    // list of all non supported solutions and points
    vector<MultiCFN::Solution> all_nonsup_solutions;
    vector<Point> all_nonsup_points;

    // computation of the nonsupported nondominated points
    for (unsigned int indTri : triangleInd) {

        if (processed[indTri]) {
            continue;
        }

        // compute and filter all the points
        tb2init();
        ToulBar2::vac = _vac;
        ToulBar2::seed = _seed;
        ToulBar2::verbose = _verbose;
        ToulBar2::showSolutions = _showSolutions;
        ToulBar2::allSolutions = nbLimit;

        multicfn->setWeight(0, triangleWeights[indTri].first);
        multicfn->setWeight(1, triangleWeights[indTri].second);

        WeightedCSP* wcsp = multicfn->makeWeightedCSP();
        wcsp->setUb(wcsp->DoubletoCost(corners[indTri].first * triangleWeights[indTri].first + corners[indTri].second * triangleWeights[indTri].second));

        WeightedCSPSolver* solver = WeightedCSPSolver::makeWeightedCSPSolver(MAX_COST, wcsp);
        solver->getWCSP()->sortConstraints();

        // solve the problem
        bool result = solver->solve();

        if (result) {

            vector<pair<Double, vector<Value>>> tb2_sol = solver->getSolutions();

            if (tb2_sol.size() == nbLimit && ToulBar2::verbose >= 0) {
                cerr << "Warning! The limit on maximum number of solutions has been reached, some solutions might have been missed!" << endl;
            }

            /* compute the solutions and points */
            vector<MultiCFN::Solution> sol;
            for (unsigned int ind = 0; ind < tb2_sol.size(); ind++) {
                sol.push_back(multicfn->convertToSolution(tb2_sol[ind].second));
            }

            vector<Point> points;
            for (unsigned int indSol = 0; indSol < sol.size(); indSol++) {
                vector<Double> values = multicfn->computeSolutionValues(sol[indSol]);
                points.push_back(make_pair(values[0], values[1]));
            }

            // build a list of indexes on the solutions
            vector<unsigned int> sol_indexes_temp(sol.size());
            for (unsigned int ind = 0; ind < sol.size(); ind++) {
                sol_indexes_temp[ind] = ind;
            }

            // filter the points by indexes
            vector<unsigned int> sol_indexes;
            for (unsigned int ind_sol : sol_indexes_temp) {

                // check if the point corresponds to one the the supported points
                if (Bicriteria::equal(points[ind_sol], trianglePoints[indTri].first)) {
                    continue;
                } else if (Bicriteria::equal(points[ind_sol], trianglePoints[indTri].second)) {
                    continue;
                } else {
                    for (auto& indTri2 : subTriangles[indTri]) {
                        if (Bicriteria::equal(points[ind_sol], trianglePoints[indTri2].first)) {
                            continue;
                        } else if (Bicriteria::equal(points[ind_sol], trianglePoints[indTri2].second)) {
                            continue;
                        }
                    }
                }

                // check if the point is dominated by one of the points from the triangles
                if (Bicriteria::dominates(trianglePoints[indTri].first, points[ind_sol], optim_dir)) {
                    continue;
                } else if (Bicriteria::dominates(trianglePoints[indTri].second, points[ind_sol], optim_dir)) {
                    continue;
                }
                bool dom = false;
                for (auto& indTri2 : subTriangles[indTri]) {
                    if (Bicriteria::dominates(trianglePoints[indTri2].first, points[ind_sol], optim_dir)) {
                        dom = true;
                        break;
                    } else if (Bicriteria::dominates(trianglePoints[indTri2].second, points[ind_sol], optim_dir)) {
                        dom = true;
                        break;
                    }
                }
                if (dom) {
                    continue;
                }

                // make sure the point dominates at least one corner
                dom = false;
                if (Bicriteria::dominates(points[ind_sol], corners[indTri], optim_dir)) {
                    dom = true;
                } else {
                    for (auto& indTri2 : subTriangles[indTri]) {
                        if (Bicriteria::dominates(points[ind_sol], corners[indTri2], optim_dir)) {
                            dom = true;
                            break;
                        }
                    }
                }
                if (!dom) {
                    continue;
                }

                // look for dominating points
                if (find_if(sol_indexes_temp.begin(), sol_indexes_temp.end(), [ind_sol, points, optim_dir](unsigned int ind_sol2) { return ind_sol != ind_sol2 && dominates(points[ind_sol2], points[ind_sol], optim_dir); }) != sol_indexes_temp.end()) {
                    continue;
                }

                // they made it here
                sol_indexes.push_back(ind_sol);
            }

            vector<MultiCFN::Solution> temp_sol;
            vector<Point> temp_points;

            for (unsigned int index : sol_indexes) {
                temp_sol.push_back(sol[index]);
                temp_points.push_back(points[index]);
            }

            all_nonsup_points.insert(all_nonsup_points.end(), temp_points.begin(), temp_points.end());
            all_nonsup_solutions.insert(all_nonsup_solutions.end(), temp_sol.begin(), temp_sol.end());
        }

        // mark the triangle and the sub triangles as processed
        processed[indTri] = true;
        for (auto& indTri2 : subTriangles[indTri]) {
            processed[indTri2] = true;
        }
    }

    // cout << "n solutions: " << sol.size() << endl;

    // filter the points by domination one more time
    vector<unsigned int> nondom_indexes;
    for (unsigned int sol_ind = 0; sol_ind < all_nonsup_points.size(); sol_ind++) {
        bool dom = false;
        for (unsigned int sol_ind2 = 0; sol_ind2 < all_nonsup_points.size(); sol_ind2++) {
            if (sol_ind == sol_ind2) {
                continue;
            }
            if (dominates(all_nonsup_points[sol_ind2], all_nonsup_points[sol_ind], optim_dir)) {
                dom = true;
                break;
            }
        }
        if (!dom) {
            nondom_indexes.push_back(sol_ind);
        }
    }

    vector<MultiCFN::Solution> all_nonsup_solutions_temp = all_nonsup_solutions;
    vector<Point> all_nonsup_points_temp = all_nonsup_points;

    all_nonsup_points.clear();
    all_nonsup_solutions.clear();

    for (unsigned int sol_ind : nondom_indexes) {
        all_nonsup_points.push_back(all_nonsup_points_temp[sol_ind]);
        all_nonsup_solutions.push_back(all_nonsup_solutions_temp[sol_ind]);
    }

    /* solutions are added to the list, then sorted */
    for (unsigned int ind = 0; ind < all_nonsup_points.size(); ind++) {
        _solutions.push_back(all_nonsup_solutions[ind]);
        _points.push_back(all_nonsup_points[ind]);
        _weights.push_back(make_pair(0., 0.));
    }

    sortSolutions(optim_dir);
}

//--------------------------------------------------------------------------------------------
void Bicriteria::computeSupportedPoints(MultiCFN* multicfn, pair<Bicriteria::OptimDir, Bicriteria::OptimDir> optim_dir, Double delta)
{
    Bicriteria::computeSupportedPoints(multicfn, 0, 1, optim_dir, delta);
}

//--------------------------------------------------------------------------------------------
void Bicriteria::computeSupportedPoints(MultiCFN* multicfn, unsigned int first_cfn_index, unsigned int second_cfn_index, pair<Bicriteria::OptimDir, Bicriteria::OptimDir> optim_dir, Double delta)
{
    // init the global time
    _startGlobalCpuTime = cpuTime();

    _first_cfn_index = first_cfn_index;
    _second_cfn_index = second_cfn_index;

    _solutions.clear();
    _points.clear();
    _weights.clear();
    _lower_bounds.clear();

    //ToulBar2::cfn = true;

    Point point1, point2;

    Double lambda1 = 1., lambda2 = 1.;
    if (optim_dir.first == Optim_Max) {
        lambda1 = -1.;
    }
    if (optim_dir.second == Optim_Max) {
        lambda2 = -1.;
    }

    MultiCFN::Solution sol1, sol2;

    // cout << endl << endl;
    // cout << "optimizing 1 separately: " << endl;
    // // point1 = solve_scalarization(multicfn, make_pair(lambda1,-0.01), &sol1);
    // point1 = solve_scalarization(multicfn, make_pair(lambda1,0), &sol1);
    // cout << "Optimal point for 1: " << point1.first << ";" << point1.second << endl;

    // cout << endl << endl;
    // cout << "optimizing 2 separately: " << endl;
    // // point2 = solve_scalarization(multicfn, make_pair(-0.01,lambda2), &sol2);
    // point2 = solve_scalarization(multicfn, make_pair(0,lambda2), &sol2);
    // cout << "Optimal point for 2: " << point2.first << ";" << point2.second << endl;
    // cout << endl << endl;

    // delta is compared to an internal decimalPoint because ToulBar2::decimalPoint is only initialized when created the final wcsp (makeMultiCFN)
    if (log10(delta) > multicfn->getDecimalPoint() && ToulBar2::verbose >= 0) {
        cerr << "Error: delta constant (" << delta << ") is incompatible with decimalPoint (" << ToulBar2::decimalPoint << ")" << endl;
    }

    bool result1, result2;
    Weights weights1, weights2;
    Double lb1, lb2;

    // cout << "optimizing 1 separately: " << endl;
    if (optim_dir.second == Optim_Max) {
        weights1 = make_pair(lambda1, -delta);
    } else {
        weights1 = make_pair(lambda1, delta);
    }

    result1 = solveScalarization(multicfn, weights1, &sol1, &point1, &lb1);

    // cout << "Optimal point for 1: " << point1.first << ";" << point1.second << endl;

    // cout << "optimizing 2 separately: " << endl;
    if (optim_dir.first == Optim_Max) {
        weights2 = make_pair(-delta, lambda2);
    } else {
        weights2 = make_pair(delta, lambda2);
    }
    result2 = solveScalarization(multicfn, weights2, &sol2, &point2, &lb2);

    // cout << "Optimal point for 2: " << point2.first << ";" << point2.second << endl;
    // cout << endl << endl;

    stack<pair<Point, Point>> pending;

    if (result1) { // make sure there is a solution
        _weights.push_back(weights1);
        _points.push_back(point1);
        _solutions.push_back(sol1);
        _lower_bounds.push_back(lb1);
    }

    if (result2 && notEqual(point1, point2)) {
        pending.push(make_pair(point1, point2));
        _weights.push_back(weights2);
        _points.push_back(point2);
        _solutions.push_back(sol2);
        _lower_bounds.push_back(lb2);
    }

    unsigned int iter = 0;

    MultiCFN::Solution new_sol;

    while (!pending.empty()) {

        iter++;

        pair<Point, Point> top = pending.top();
        pending.pop();

        // compute the new weights

        if (optim_dir.first == optim_dir.second) {
            lambda1 = -top.second.second + top.first.second;
            lambda2 = -top.first.first + top.second.first;
        } else {
            lambda1 = top.second.second - top.first.second;
            lambda2 = top.first.first - top.second.first;
        }

        // normalize the weights
        Double length = sqrt(lambda1 * lambda1 + lambda2 * lambda2);
        lambda1 *= 2. / length;
        lambda2 *= 2. / length;

        if (optim_dir.first == optim_dir.second && lambda1*lambda2 < 0) {
            continue;
        }
        if (optim_dir.first != optim_dir.second && lambda1*lambda2 > 0) {
            continue;
        }

        // cout << "weights : " << lambda1 << ", " << lambda2 << endl;
        Point new_point;
        Double new_lb;

        Weights new_weights = make_pair(lambda1, lambda2);
        bool result = solveScalarization(multicfn, new_weights, &new_sol, &new_point, &new_lb);

        // jump to the next weights if there is no solution
        if (!result) {
            continue;
        }

        // check the signs of the weights and skip if not consistent
        // may happen when soltuions are not optimal

        // cout << "new point: " << new_point.first << ", " << new_point.second << endl;
        // cout << "from " << top.first.first << "," << top.first.second << " and " << top.second.first << ", " << top.second.second << endl;

        // cut the search if the point was alrady encountered
        auto it = std::find_if(_points.begin(), _points.end(), [new_point](Point& point) { return equal(point, new_point); });
        if (it != _points.end()) {
            continue;
        }

        if (notEqual(new_point, top.first) && notEqual(new_point, top.second)) {

            _weights.push_back(new_weights);
            _points.push_back(new_point);
            _solutions.push_back(new_sol);
            _lower_bounds.push_back(new_lb);

            // check if the scalarized cost is better before adding new weights
            Double cost_new = new_weights.first*new_point.first + new_weights.second*new_point.second;
            Double cost_prev = new_weights.first*top.first.first+new_weights.second*top.first.second;

            if(cost_new < cost_prev) {

                // add two new scalarizations
                if (fabs(top.first.first - new_point.first) >= ToulBar2::epsilon && fabs(top.first.second - new_point.second) >= ToulBar2::epsilon) {
                    pending.push(make_pair(top.first, new_point));
                }
                if (fabs(top.second.first - new_point.first) >= ToulBar2::epsilon && fabs(top.second.second - new_point.second) >= ToulBar2::epsilon) {
                    pending.push(make_pair(new_point, top.second));
                }

            }
        }

        // cout << endl << endl << "debug" << endl << endl << endl;

        // cout << "points found so far: ";
        // for(auto& point: sol_points) {
        //   cout << "(" << point.first << "," << point.second << "), ";
        // }
        // cout << endl;
    
        // leave the loop if max number of solution reached
        // or if global timeout reached
        if((Bicriteria::_global_timeout > 0 && cpuTime() - _startGlobalCpuTime > Bicriteria::_global_timeout) || (Bicriteria::_nsol_max > 0 && static_cast<int>(_solutions.size()) >= Bicriteria::_nsol_max) ) {
            break;
        }
    
    }

    // cout << "Supported points: " << endl;
    // for(unsigned int point_ind = 0; point_ind < sol_points.size(); point_ind ++) {
    //   cout << sol_points[point_ind].first << ", " << sol_points[point_ind].second << endl;
    // }

    // the points are sorted from left to right
    sortSolutions(optim_dir);
}

//--------------------------------------------------------------------------------------------
std::vector<MultiCFN::Solution> Bicriteria::getSolutions()
{
    return _solutions;
}

//--------------------------------------------------------------------------------------------
std::vector<Bicriteria::Point> Bicriteria::getPoints()
{
    return _points;
}

//--------------------------------------------------------------------------------------------
std::vector<Bicriteria::Weights> Bicriteria::getWeights()
{
    return _weights;
}

//--------------------------------------------------------------------------------------------
std::vector<Double> Bicriteria::getLowerBounds()
{
    return _lower_bounds;
}
