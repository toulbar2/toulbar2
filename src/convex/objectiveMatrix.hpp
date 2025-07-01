#ifndef OBJECTIVEMATRIX_HPP_
#define OBJECTIVEMATRIX_HPP_

#include "tb2config.hpp"

#ifdef LR_BCD_BUILD

#include <string>
#include <iostream>
#include <vector>

#include <eigen3/Eigen/Dense>
#include <core/tb2binconstr.hpp>

// typedef for dense vectors and matrices using Double entries
typedef Eigen::Matrix<Double, Eigen::Dynamic, 1> DnVec;
typedef Eigen::Matrix<Double, Eigen::Dynamic, Eigen::Dynamic> DnMat;

// return a vector of int with the domain indices
// In :
//      - WeightedCSP* wcsp
// Out :
//      - vector<size_t> dom
// TODO: do it while iterating on the variables with getEnumDomain
vector<size_t> domains(WeightedCSP* wcsp)
{
    assert(wcsp);

    size_t acc = 0;
    vector<size_t> dom;
    dom.push_back(acc);

    for (size_t i = 0; i < wcsp->numberOfVariables(); i++) {
        if (wcsp->unassigned(i)) {
            acc += wcsp->getEnumDomainAndCost(i).size();
            dom.push_back(acc);
        }
    }

    assert(dom.size() > 0);
    return dom;
}

// return a vector with the updated indices of the
// variables: WCSP -> SDP
// In :
//      - WeightedCSP* wcsp
// Out :
//      - vector<size_t> updatedIndex
vector<size_t> var2Index(WeightedCSP* wcsp)
{
    assert(wcsp);

    size_t index = 0;
    vector<size_t> updatedIndex;

    for (size_t i = 0; i < wcsp->numberOfVariables(); i++) {
        updatedIndex.push_back(index);
        if (wcsp->unassigned(i)) {
            index++;
        }
    }

    assert(updatedIndex.size() > 0);
    return updatedIndex;
}

// return a vector with the indices of the unassigned wcsp variables
// SDP -> WCSP
// In :
//      - WeightedCSP* wcsp
// Out :
//      - vector<size_t> index2var
vector<size_t> index2Var(WeightedCSP* wcsp)
{
    assert(wcsp);

    vector<size_t> index2var;

    for (size_t i = 0; i < wcsp->numberOfVariables(); i++) {
        if (wcsp->unassigned(i)) {
            index2var.push_back(i);
        }
    }

    assert(index2var.size() > 0);
    return index2var;
}

// return a vector containing the values for each unassigned wcsp variable
// In :
//      - WeightedCSP* wcsp
// Out :
//      - vector<vector<Value>> values
vector<vector<Value>> getValues(WeightedCSP* wcsp)
{
    vector<vector<Value>> values;

    for (size_t i = 0; i < wcsp->numberOfVariables(); i++) {
        if (wcsp->unassigned(i)) {
            vector<pair<Value, Cost>> domcosts = wcsp->getEnumDomainAndCost(i); // get variable values and corresponding unary cost
            vector<Value> domvals;

            for (size_t k = 0; k < domcosts.size(); k++) {
                domvals.push_back(domcosts[k].first);
            }

            values.push_back(domvals);
        }
    }

    return values;
}

// compute stacked vector of unary costs from wcsp instance
// In :
//      - WeightedCSP* wcsp
// Out :
//      - DnVec unaryCost
DnVec unaryCost(WeightedCSP* wcsp)
{
    assert(wcsp);

    size_t d = wcsp->getDomainSizeSum();
    vector<size_t> var2index = var2Index(wcsp);
    vector<size_t> dom = domains(wcsp);
    DnVec unaryCost = DnVec::Zero(d);

    // store unary costs
    for (size_t i = 0; i < wcsp->numberOfVariables(); i++) {
        if (wcsp->unassigned(i)) {
            vector<pair<Value, Cost>> domcosts = wcsp->getEnumDomainAndCost(i); // get variable values and corresponding unary cost
            size_t indI = dom[var2index[i]];

            for (size_t k = 0; k < domcosts.size(); k++) {
                unaryCost(indI + k) = (Double)domcosts[k].second;
            }
        }
    }

    return unaryCost;
}

// compute the symmetric matrix for binary costs
// In :
//      - WeightedCSP* wcsp
// Out :
//      - vector<Double> binaryCost
DnMat binaryCost(WeightedCSP* wcsp)
{
    assert(wcsp);

    size_t d = wcsp->getDomainSizeSum();
    vector<size_t> var2index = var2Index(wcsp);
    vector<size_t> dom = domains(wcsp);
    DnMat binaryCost = DnMat::Zero(d, d);

    // store quadratic costs
    for (size_t k = 0; k < wcsp->numberOfConstraints(); k++) {
        auto* ctr = dynamic_cast<BinaryConstraint*>(((WCSP*)wcsp)->getCtr(k));
        if (!ctr)
            continue;
        else if (ctr->connected() && !ctr->isSep() && ctr->isBinary()) {
            vector<Value> domi = wcsp->getEnumDomain(ctr->getVar(0)->wcspIndex);
            vector<Value> domj = wcsp->getEnumDomain(ctr->getVar(1)->wcspIndex);

            size_t indI = dom[var2index[ctr->getVar(0)->wcspIndex]];
            size_t indJ = dom[var2index[ctr->getVar(1)->wcspIndex]];

            for (size_t i = 0; i < domi.size(); i++) {
                for (size_t j = 0; j < domj.size(); j++) {
                    Double cost = (Double)ctr->getCost(domi[i], domj[j]);
                    if (cost != 0) {
                        binaryCost(indI + i, indJ + j) += 0.5 * cost;
                        binaryCost(indJ + j, indI + i) += 0.5 * cost; // binaryCost must be symmetric
                    }
                }
            }
        }
    }

    for (int i = 0; i < ((WCSP*)wcsp)->getElimBinOrder(); i++) {
        BinaryConstraint* ctr = (BinaryConstraint*)((WCSP*)wcsp)->getElimBinCtr(i);
        if (ctr->connected() && !ctr->isSep()) {
            vector<Value> domi = wcsp->getEnumDomain(ctr->getVar(0)->wcspIndex);
            vector<Value> domj = wcsp->getEnumDomain(ctr->getVar(1)->wcspIndex);

            size_t indI = dom[var2index[ctr->getVar(0)->wcspIndex]];
            size_t indJ = dom[var2index[ctr->getVar(1)->wcspIndex]];

            for (size_t i = 0; i < domi.size(); i++) {
                for (size_t j = 0; j < domj.size(); j++) {
                    Double cost = (Double)ctr->getCost(domi[i], domj[j]);
                    if (cost != 0) {
                        binaryCost(indI + i, indJ + j) += 0.5 * cost;
                        binaryCost(indJ + j, indI + i) += 0.5 * cost; // binaryCost must be symmetric
                    }
                }
            }
        }
    }

    return binaryCost;
}

// compute objective matrix for LR-BCD
// from wcsp instance
// In :
//      - WeightedCSP* wcsp
// Out :
//      - DnMat costMatrix
DnMat objectiveMatrix(WeightedCSP* wcsp)
{
    assert(wcsp);

    size_t d = wcsp->getDomainSizeSum();

    DnVec unaryCostVec = unaryCost(wcsp);
    DnMat binaryCostMat = binaryCost(wcsp);
    DnMat costMatrix = DnMat::Zero(d + 1, d + 1);
    DnVec ones_d = DnVec::Ones(d); // vector of ones of size d

    // -1/1 change of variables
    // linear part (homogenization)
    // c = 0.5 * Q * 1_d + 0.5 * c
    unaryCostVec = 0.5 * (binaryCostMat * ones_d + unaryCostVec);
    costMatrix.block(0, d, d, 1) = 0.5 * unaryCostVec;
    costMatrix.block(d, 0, 1, d) = 0.5 * unaryCostVec.transpose();

    // quadratic part
    // Q = 0.25 * Q
    costMatrix.block(0, 0, d, d) = 0.25 * binaryCostMat;

    return costMatrix;
}

#endif // LR_BCD_BUILD

#endif
