#define _USE_MATH_DEFINES

#include <vector>
//#include <random>
#include <algorithm>
#include <cmath>
#include <float.h>
#include <eigen3/Eigen/Dense>

#include "objectiveMatrix.hpp"
#include <core/tb2wcsp.hpp>
#include <search/tb2solver.hpp>
#include <core/tb2binconstr.hpp>

using namespace std;

// typedef for dense vectors and matrices using Double entries
typedef Eigen::Matrix<Double, Eigen::Dynamic, 1> DnVec;
typedef Eigen::Matrix<Double, Eigen::Dynamic, Eigen::Dynamic> DnMat;

/////////////////////////  LR-BCD code  /////////////////////////
//                                                             //
//      Compute primal lower bound using LR-BCD                //
//      Compute feasible integer solution using GW heuristic   //
//      and local 1-opt search                                 //
//                                                             //
/////////////////////////////////////////////////////////////////

// build rank for the low rank factorization
// In :
//      - WeightedCSP* wcsp
// Out :
//      - size_t k : rank
size_t getRank(WeightedCSP* wcsp)
{
    assert(wcsp);

    int n = wcsp->getDomainSizeSum() + 1;

    // add the N exactly one constraints
    n += wcsp->numberOfUnassignedVariables();
    size_t k = ceil(sqrt(2 * n));

    assert(k > 0);
    return k;
}

// draw a random vector from the standard normal distribution N(0,1)
// and normalize it
// In :
//      - DnVec &v
//      - size_t k
//      - default_random_engine &generator
//      - normal_distribution<Double> &distribution
void randUnit(DnVec& v, size_t k) //, default_random_engine& generator, normal_distribution<Double>& distribution)
{
    assert(k > 0);

    for (size_t i = 0; i < k; i++) {
        v(i) = mynrandl(); // distribution(generator);
    }

    v = v / v.norm();
}

// Initialize low rank matrix
// In :
//      - WeightedCSP* wcsp
//      - size_t k: rank of the low rank matrix
// Out :
//      - DnMat V: low rank matrix of size n * k
DnMat mixingInit(WeightedCSP* wcsp, size_t k)
{
    assert(wcsp);
    assert(k > 0);

    size_t n = wcsp->getDomainSizeSum() + 1;
    DnMat V = DnMat::Zero(k, n);

    // default_random_engine generator;
    // generator.seed(time(0));
    //  normal_distribution<Double> distribution(0, 1); // normal distribution
    // uniform_real_distribution<Double> distribution(-1.0, 1.0); // uniform distribution

    for (size_t j = 0; j < n; j++) {
        for (size_t i = 0; i < k; i++) {
            V(i, j) = myurandl(); // distribution(generator);
        }
        V.col(j) = V.col(j) / (V.col(j).norm());
    }

    return V;
}

// Bias introduced by the change of variable {0,1} -> {-1,1}
// In :
//      - WeightedCSP* wcsp
// Out :
//      - Double bias
Double bias(WeightedCSP* wcsp)
{
    assert(wcsp);

    Double unaryBias = 0;
    Double binaryBias = 0;
    Double lb = wcsp->getLb();
    Double acc = 0;

    // sum of unary costs
    for (size_t i = 0; i < wcsp->numberOfVariables(); i++) {
        if (wcsp->unassigned(i)) {
            vector<pair<Value, Cost>> domcosts = wcsp->getEnumDomainAndCost(i); // get variable values and corresponding unary cost

            for (size_t k = 0; k < domcosts.size(); k++) {
                acc += (Double)domcosts[k].second;
            }

            unaryBias += 0.5 * acc;
            acc = 0;
        }
    }

    // sum of binary costs
    for (size_t k = 0; k < wcsp->numberOfConstraints(); k++) {
        auto* ctr = dynamic_cast<BinaryConstraint*>(((WCSP*)wcsp)->getCtr(k));
        if (!ctr)
            continue;
        else if (ctr->connected() && !ctr->isSep() && ctr->isBinary()) {
            vector<Value> domi = wcsp->getEnumDomain(ctr->getVar(0)->wcspIndex);
            vector<Value> domj = wcsp->getEnumDomain(ctr->getVar(1)->wcspIndex);

            for (size_t i = 0; i < domi.size(); i++) {
                for (unsigned j = 0; j < domj.size(); j++) {
                    acc += (Double)ctr->getCost(domi[i], domj[j]);
                }
            }
            binaryBias += 0.25 * acc;
            acc = 0;
        }
    }

    for (int i = 0; i < ((WCSP*)wcsp)->getElimBinOrder(); i++) {
        BinaryConstraint* ctr = (BinaryConstraint*)((WCSP*)wcsp)->getElimBinCtr(i);
        if (ctr->connected() && !ctr->isSep()) {
            vector<Value> domi = wcsp->getEnumDomain(ctr->getVar(0)->wcspIndex);
            vector<Value> domj = wcsp->getEnumDomain(ctr->getVar(1)->wcspIndex);

            for (size_t i = 0; i < domi.size(); i++) {
                for (unsigned j = 0; j < domj.size(); j++) {
                    acc += (Double)ctr->getCost(domi[i], domj[j]);
                }
            }
            binaryBias += 0.25 * acc;
            acc = 0;
        }
    }

    return binaryBias + unaryBias + lb;
}

// Rounding with GW heuristic
// In :
//      - const DnMat &V
//      - const vector<size_t> &domains
//      - size_t n
//      - size_t k
//      - int i : choice of the seed for random vector generation
DnVec rounding(const DnMat& V, const vector<size_t>& domains, size_t n, size_t k, int i)
{
    assert(n > 0);
    assert(k > 0);

    // default_random_engine generator;
    // generator.seed(i + 1);
    // normal_distribution<Double> distribution(0, 1);
    mysrand(abs(ToulBar2::seed) + i);

    DnVec r = DnVec::Zero(k);
    DnVec intSol = DnVec::Constant(n, -1.0);
    intSol(n - 1) = 1; // last value for homogenization must be equal to 1

    randUnit(r, k); //, generator, distribution);

    for (size_t i = 0; i < domains.size() - 1; i++) {
        size_t first = domains[i];
        size_t last = domains[i + 1];

        int ind = first;
        Double max = fabs(r.dot(V.col(first)));
        Double previous_max = max;

        for (size_t j = first + 1; j < last; j++) {
            // compute new trial point
            max = fabs(r.dot(V.col(j)));

            if (max > previous_max) {
                previous_max = max;
                ind = j;
            }
        }

        intSol(ind) = 1;
    }

    return intSol;
}

// return vector of assignments corresponding to the Max-Cut solution
// In :
//      - DnVec &intSol
//      - vector<size_t> &domains
// Out :
//      - vector<int> wcspSol
vector<int> mc2tb2(DnVec& intSol, vector<size_t>& domains)
{

    vector<int> wcspSol;

    for (size_t i = 0; i < domains.size() - 1; i++) {
        int value = 0;
        for (size_t j = domains[i]; j < domains[i + 1]; j++) {
            if (intSol(j) == -1) {
                value++;
            } else {
                break;
            }
        }
        wcspSol.push_back(value);
    }

    return wcspSol;
}

// Evaluate objective function
// In :
//      - const DnMat &C: objective matrix
//      - const DnMat &V: low rank matrix
// Out :
//      - Double eval : function value
Double evalFun(const DnMat& C, const DnMat& V)
{

    DnMat A = C * V.transpose();
    A = V * A;
    Double eval = A.trace();

    return eval;
}

// Evaluate vector objective function
// In :
//      - const DnMat &C: objective matrix
//      - const DnVec &V: vector
// Out :
//      - Double eval : function value
Double evalFun(const DnMat& C, const DnVec& V)
{

    DnVec A = C * V;
    Double eval = V.dot(A);

    return eval;
}

// Build the objective vectors for the block-optimization
// sub-problems
// In :
//    - size_t first : first block index
//    - size_t last : last block index
//    - const DnMat& C : cost matrix
//    - const DnMat& V : SDP variable
//    - size_t k : rank of V
// Out :
//    - DnMat G : (k,last-first)-matrix of objective vectors
DnMat objectiveVec(size_t first, size_t last, const DnMat& C, const DnMat& V, size_t k)
{
    assert(0 <= first && first <= last);
    assert(k > 0);

    size_t d = last - first;
    DnMat G = DnMat::Zero(k, d);
    DnVec g = DnVec::Zero(k);

    for (size_t i = first; i < last; i++) {
        g = V * C.col(i);
        G.col(i - first) = g;
    }

    return G;
}

// Compute the angles between the objective vectors
// and V_{d+1}
// In :
//    - const DnMat& G
//    - const DnVec& u
//    - size_t d
// Out
//    - vector<Double> angles
vector<Double> angles(const DnMat& G, const DnVec& u, size_t d)
{

    vector<Double> angles(d, 0);

    for (size_t i = 0; i < d; i++) {
        DnVec g = G.col(i);
        Double a = acos(g.dot(u) / g.norm());
        angles[i] = a;
    }

#ifndef NDEBUG
    // g_i's and V_{d+1} should not be colinear
    for (size_t i = 0; i < d; i++) {
        assert(angles[i] != 0 && angles[i] != M_PI);
    }
#endif

    return angles;
}

// Compute the solution of the relaxation
// In :
//      - const DnMat &G
//      - const DnVec &u
//      - size_t d
//      - size_t k
// Out :
//      - DnVec res
DnVec relaxedSol(const DnMat& G, const DnVec& u, size_t d, size_t k)
{
    assert(d > 0);
    assert(k > 0);

    size_t n = k * d;
    DnVec g = DnVec::Zero(n);
    DnVec w = DnVec::Zero(n);

    for (size_t i = 0; i < n; i++) {
        size_t q = i / k;
        size_t r = i % k;

        g(i) = G(r, q);
        w(i) = u(r);
    }

    Double r = 4 - (4. / d);

    Double gamma = g.squaredNorm() - (1. / d) * pow(w.dot(g), 2);
    Double beta = -sqrt(r / gamma);
    Double alpha = -(1. / d) * (beta * w.dot(g) + d - 2);
    DnVec res = alpha * w + beta * g;

    return res;
}

// Compute a starting point for the dual optimization
// In :
//      - const DnMat &G
//      - const DnMat &u
//      - const DnVec &rSol
//      - size_t k
// Out :
//      - Double dualInit
Double dualInit(const DnMat& G, const DnVec& u, const DnVec& rSol, size_t k)
{
    assert(k > 0);

    DnVec x = rSol.head(k);
    DnVec g = G.col(0);

    Double nrm_g = g.norm();
    Double nrm_x = x.norm();
    Double angle_g_x = acos(g.dot(x) / (nrm_g * nrm_x));
    Double angle_u_x = acos(u.dot(x) / nrm_x);
    Double dualInit = -nrm_g * sin(angle_g_x) / sin(angle_u_x);

    return dualInit;
}

// Compute the solution angles to the primal problem
// In :
//      - const DnMat &G
//      - const vector<Double> &angles
//      - Double dualOpti
//      - size_t d
// Out :
//      - vector<Double> solAngle
vector<Double> solAngle(const DnMat& G, const vector<Double>& angles, Double dualOpti, size_t d)
{
    assert(d > 0);

    vector<Double> solAngle(d, 0);

    for (size_t i = 0; i < d; i++) {
        DnVec g = G.col(i);
        Double nrm_g = g.norm();
        Double theta = angles[i];

        solAngle[i] = fmod(M_PI - atan(nrm_g * sin(theta) / (nrm_g * cos(theta) + dualOpti)), M_PI);
    }

    return solAngle;
}

// Compute the solution vectors to the primal problem
// In :
//      - const DnMat &G
//      - const vector<Double> &solAngle
//      - const vector<Double> &angles
//      - const DnVec &u
//      - vector<Double> &Vb
//      - size_t d
//      - size_t k
// Out :
//      - vector<Double> res
DnMat vecSol(const DnMat& G, const vector<Double>& solAngle,
    const vector<Double>& angles, const DnVec& u, size_t d, size_t k)
{
    assert(d > 0);
    assert(k > 0);

    DnMat res(k, d);
    for (size_t i = 0; i < d; i++) {
        DnVec g = G.col(i);
        Double nrm_g = g.norm();
        Double vdot = u.dot(g);
        Double x = solAngle[i];
        Double theta = angles[i];
        Double c_x = cos(x);
        Double c_tx = cos(x + theta);
        Double denom = g.squaredNorm() - pow(vdot, 2);

        Double b = (nrm_g * c_tx - c_x * vdot) / denom;
        Double a = c_x - vdot * b;

        assert(b < 0);

        res.col(i) = a * u + b * g;
    }

    return res;
}

// Compute the objective difference after the update of a block
// In :
//      - const DnMat &G
//      - const DnMat &V
//      - const DnMat &Vb
//      - size_t first
//      - size_t last
// Out :
//      - Double currDelta
Double currentDelta(const DnMat& G, const DnMat& V, const DnMat& Vb, size_t first,
    size_t last)
{
    assert(0 <= first && first <= last);

    Double currDelta = 0;
    for (size_t i = first; i < last; i++) {
        currDelta += 2 * (V.col(i) - Vb.col(i - first)).dot(G.col(i - first));
    }

    return currDelta;
}

// update the rows of V corresponding to the block indexed by first, last
// In :
//      - const DnMat &V
//      - const DnMat &Vb
//      - size_t first
//      - size_t last
void solUpdate(DnMat& V, const DnMat& Vb, size_t first, size_t last)
{
    assert(0 <= first && first <= last);

    for (size_t i = first; i < last; i++) {
        V.col(i) = Vb.col(i - first);
    }
}

// compute jacobian of the dual fonction
// In :
//      - const DnMat &G
//      - const DnVec &u
//      - Double x
//      - size_t d
// Out :
//      - Double res : value of the jacobian at x
Double jac_h(const DnMat& G, const DnVec& u, Double x, size_t d)
{
    assert(d > 0);

    Double res = 0;

    for (size_t i = 0; i < d; i++) {
        DnVec w = G.col(i) + x * u;
        res += (G.col(i).dot(u) + x) / w.norm();
    }
    res = res - (d - 2);

    return res;
}

// compute hessian of the dual fonction
// In :
//      - const DnMat &G
//      - const DnVec &u
//      - Double x
//      - size_t d
// Out :
//      - Double res : value of the hessian at x
Double hess_h(const DnMat& G, const DnVec& u, Double x, size_t d)
{
    assert(d > 0);

    Double res = 0;

    for (size_t i = 0; i < d; i++) {
        DnVec g = G.col(i);
        DnVec w = g + x * u;
        res += (g.squaredNorm() - pow(g.dot(u), 2)) / pow(w.norm(), 3);
    }

    return res;
}

// Newton method
// In :
//      - const DnMat &G
//      - const DnVec &u
//      - Double x
//      - size_t d
//      - size_t maxiter
//      - size_t &N_it
// Out :
//      - Double x_curr : optimizer of the dual problem
Double newton(const DnMat& G, const DnVec& u, Double x, size_t d, size_t maxiter, size_t& N_it)
{
    assert(d > 0);
    assert(maxiter > 0);

    Double x_curr = x;
    Double jac = jac_h(G, u, x, d);
    Double jac_curr = jac;
    size_t it = 0;
    Double eps = pow(10, -7);

    Double step_size = 0.80;
    Double alpha = 0.60; // correction parameter for the step size

    while (abs(jac) > eps && it < maxiter) {
        x_curr = x_curr - step_size * jac_h(G, u, x_curr, d) / hess_h(G, u, x_curr, d);
        jac = jac_h(G, u, x_curr, d);
        it++;

        if (abs(jac) >= abs(jac_curr)) // step_size has to be updated
        {
            step_size *= alpha; // reduce the step size
            x_curr = x;
            jac = jac_h(G, u, x, d);
            it = 0;
        }

        jac_curr = jac;
    }

    N_it = N_it + it;

    return x_curr;
}

// solve the block optimization problem
// In :
//      - const DnMat &C
//      - DnMat &V
//      - const DnMat &u
//      - const vector<size_t> &domains
//      - Double &delta
//      - size_t k
//      - size_t i
//      - size_t &N_it
void solveBCD(const DnMat& C, DnMat& V, const DnVec& u, const vector<size_t>& domains, Double& delta, size_t k, size_t i, size_t& N_it)
{
    assert(k > 0);
    assert(i >= 0);

    size_t first = domains[i];
    size_t last = domains[i + 1];
    size_t d = last - first;
    size_t maxiter = 100;

    DnMat G = objectiveVec(first, last, C, V, k);
    vector<Double> a = angles(G, u, d);

    DnVec solRelax = relaxedSol(G, u, d, k);
    Double dual = dualInit(G, u, solRelax, k);

    // Dual opti
    Double dualOpti = newton(G, u, dual, d, maxiter, N_it);

    // Build solution vectors
    vector<Double> angleOpti = solAngle(G, a, dualOpti, d);
    DnMat Vb = vecSol(G, angleOpti, a, u, d, k);

    delta += currentDelta(G, V, Vb, first, last);
    solUpdate(V, Vb, first, last);
}

// return the low-rank solution to the SDP using LR-BCD
// In :
//      - WeightedCSP* wcsp
//      - const DnMat &C
//      - Double tol
//      - size_t maxiter
//      - size_t k
// Out :
//      - vector<Double> V : low-rank solution to the SDP
DnMat LR_BCD(WeightedCSP* wcsp, const DnMat& C, Double tol, size_t maxiter, size_t k)
{
    assert(maxiter > 0);
    assert(k > 0);

    vector<size_t> dom = domains(wcsp);
    size_t d = wcsp->getDomainSizeSum();
    size_t BCD_it = 0;
    size_t N_it = 0;

    // Stopping criteria
    bool stop = false;
    Double delta = 0;

    DnVec q = DnVec::Zero(d + 1);
    DnVec g = DnVec::Zero(k);

    // V init
    DnMat V = mixingInit(wcsp, k);
    DnMat V_old = V;
    Double f_val = evalFun(C, V);

    for (size_t i = 0; i < maxiter; i++) {
        DnVec u = V.col(d);
        delta = 0;

        for (size_t j = 0; j < dom.size() - 1; j++) {
            solveBCD(C, V, u, dom, delta, k, j, N_it);
        }

        BCD_it++;

        assert(delta > 0); // delta should always be positive
        DnMat V_diff = V - V_old;
        Double fun_criterion = delta / (1 + fabs(f_val));
        Double step_criterion = V_diff.norm() / (1 + V_old.norm());

        stop = (fun_criterion < tol) || (step_criterion < tol);

        V_old = V;
        f_val -= delta;

        if (stop) {
            break;
        }
    }

    return V;
}

// Apply one opt search to the integer solution
// In :
//      - const DnMat &C
//      - const vector<size_t> &domains
//      - DnVec &x
//      - Double f
//      - size_t n
// Out :
//      - Double min : value of the new solution after 1-opt
Double oneOptSearch(const DnMat& C, const vector<size_t>& domains, DnVec& x, Double f, size_t n)
{
    assert(n > 0);

    bool changed = true;
    DnVec work = DnVec::Zero(n);
    Double min = f;
    Double cont = 0;

    while (changed) {

        changed = false;

        for (size_t i = 0; i < domains.size() - 1; i++) {

            size_t first = domains[i];
            size_t last = domains[i + 1];

            // flip the entry 1 to -1 in the domain of each variable
            size_t ind = first;
            Double current_value = min;

            for (size_t j = first; j < last; j++) {
                if (x(j) == 1) {
                    x(j) = -1;
                    ind = j;
                }
            }
            size_t indMin = ind;

            // remove contribution of x(ind)
            work = C.col(ind);
            work(ind) = 0;

            cont = x.dot(work);
            current_value += 4.0 * cont * x(ind);

            for (size_t j = first; j < last; j++) {
                // add contribution of x(j)
                // remove contribution of x(ind)
                work = C.col(j);
                work(j) = 0;

                cont = x.dot(work);
                current_value -= 4.0 * cont * x(j);

                if (current_value < min) {
                    min = current_value;
                    indMin = j;
                    changed = true;
                }

                current_value += 4.0 * cont * x(j);
            }
            // f = current_min;
            // cout << "\nValues of f: " << f;
            x(indMin) = 1;
        }
    }

    return min;
}

// Do multiple roundings and keep the best solution
// In :
//      - WeightedCSP* wcsp
//      - const DnMat &C
//      - DnMat &V
//      - const vector<size_t> &domains
//      - size_t nbRound
//      - size_t n
//      - size_t k
// Out :
//      - Double min : value of the new solution after multiple 1-opt
//      - DnVec bestSol : best integer solution found
tuple<Double, DnVec> multipleRounding(WeightedCSP* wcsp, const DnMat& C, DnMat& V,
    const vector<size_t> domains, size_t nbRound, size_t n, size_t k)
{
    assert(wcsp);
    assert(nbRound > 0);
    assert(n > 0);
    assert(k > 0);

    DnVec work(n);
    DnVec bestSol = work;
    Double min = DBL_MAX;
    Double trial_sol = min;

    for (size_t i = 0; i < nbRound; i++) {
        if (ToulBar2::interrupted) {
            throw TimeOut();
        }
        work = rounding(V, domains, n, k, i);
        trial_sol = evalFun(C, work) + bias(wcsp);
        trial_sol = oneOptSearch(C, domains, work, trial_sol, n);

        if (trial_sol < min) {
            min = trial_sol;
            bestSol = work;
        }
    }

    return make_tuple(min, bestSol);
}

// LR-BCD
Cost Solver::lrBCD(string cmd, vector<Value>& solution)
{
    size_t maxiter = 5;
    size_t k = getRank(wcsp) / 2;
    size_t nbR = 50;

    istringstream ss_cmd(cmd);
    if (ss_cmd)
        ss_cmd >> maxiter;
    if (ss_cmd)
        ss_cmd >> k;
    if (ss_cmd)
        ss_cmd >> nbR;

    return lrBCD(maxiter, k, nbR, solution);
}

Cost Solver::lrBCD(size_t maxiter, int k, size_t nbR, vector<Value>& solution)
{
    if (ToulBar2::verbose >= 3) {
        cout << " LR-BCD " << maxiter;
        cout << " " << k;
        cout << " " << nbR << endl;
    }
    if (k == 0) {
        cout << "Error: Cannot use LR BCD local search with a zero rank (see -lrbcd option)." << endl;
        throw BadConfiguration();
    }

    // use current time as seed for random generator
    // srand(time(0)); //SdG: already done in tb2main

    // set tolerance for the stopping criterion
    Double tol = 1e-3;

    // set low rank factorization rank
    //    switch (k) {
    //    case -1:
    //        k = getRank(wcsp);
    //        break;
    //    case -2:
    //        k = getRank(wcsp) / 2;
    //        break;
    //    case -4:
    //        k = getRank(wcsp) / 4;
    //        break;
    //    }
    if (k < 0) {
        size_t rank = getRank(wcsp);
        k = rank / -k;
        if (ToulBar2::verbose >= 0) {
            cout << "LR-BCD rank: " << k << endl;
        }
    }
    assert(k > 0);

    size_t nb_cols = wcsp->getDomainSizeSum() + 1;
    DnMat Q = objectiveMatrix(wcsp);
    DnMat V(k, nb_cols);
    vector<size_t> index2var = index2Var(wcsp);
    vector<vector<Value>> values = getValues(wcsp);

    // LR-BCD
    V = LR_BCD(wcsp, Q, tol, maxiter, k);

    // compute lower bound
    Double lb = evalFun(Q, V) + bias(wcsp);
    assert(lb <= wcsp->getUb());
    if (ToulBar2::verbose >= 1) {
        cout << "LR-BCD putative dual bound: " << std::fixed << std::setprecision(ToulBar2::decimalPoint) << wcsp->Cost2ADCost(lb) << std::setprecision(DECIMAL_POINT) << endl;
    }

    // compute upper bound with GW heuristic
    vector<size_t> dom = domains(wcsp);
    DnVec intSol;
    Double ub_oneopt = DBL_MAX;

    // compute upper bound with GW + 1-opt search
    tie(ub_oneopt, intSol) = multipleRounding(wcsp, Q, V, dom, nbR, nb_cols, k);
    // convert Max-Cut solution to wcsp format
    vector<int> sol = mc2tb2(intSol, dom);

    // initial upper bound
    Cost initialUpperBound = wcsp->getUb();

    if (ub_oneopt < wcsp->getUb()) {
        vector<Value> bestsolution(wcsp->numberOfVariables(), 0);
        int depth = Store::getDepth();
        try {
            Store::store();
            vector<int> tabvars(sol.size());
            vector<Value> solution(sol.size());
            for (unsigned i = 0; i < sol.size(); i++) {
                tabvars[i] = index2var[i];
                solution[i] = values[i][sol[i]];
            }
            wcsp->assignLS(tabvars, solution);
            newSolution();
            for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
                bestsolution[i] = wcsp->getValue(i);
                ((WCSP*)wcsp)->setBestValue(i, bestsolution[i]);
            }
        } catch (const Contradiction&) {
            wcsp->whenContradiction();
        }
        Store::restore(depth);
    }

    if (wcsp->getUb() < initialUpperBound) {
        wcsp->enforceUb();
        wcsp->propagate();
    }

    return (wcsp->getUb() < initialUpperBound) ? wcsp->getSolutionCost() : MAX_COST;
}
