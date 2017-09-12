// extreme_value_distribution

#include "tb2wcsp.hpp"
#include <iostream>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>

using namespace std;

vector<TProb> WCSP::Gumbel_noise_vec(vector<TProb> costsProb)
{
    vector<TProb> perturb_costsProb;
    double perturb = 0;
    TProb perturb_prob = 0;
    const double a = -0.577215; // -c --> Euler constant so as to standard the law.
    const double b = 1.0;
    ToulBar2::seed = std::chrono::system_clock::now().time_since_epoch().count(); //initialize the seed to generate different random distrib.
    std::default_random_engine generator(ToulBar2::seed);
    std::extreme_value_distribution<double> distribution(a, b); // Set the distribution between a and b.
    if (ToulBar2::uai > 1) {
        for (unsigned int i = 0; i < costsProb.size(); ++i) {
            perturb = distribution(generator);
            if (costsProb[i] < pow(10, 37)) {
                perturb_prob = costsProb[i] + perturb;
                perturb_costsProb.push_back(perturb_prob);
            } else {
                perturb_costsProb.push_back(costsProb[i]);
            }
        }
    } else {
        for (unsigned int i = 0; i < costsProb.size(); ++i) {
            perturb = distribution(generator);
            if (costsProb[i] != 0) {
                perturb_prob = costsProb[i] + perturb;
                perturb_costsProb.push_back(perturb_prob);
            } else {
                perturb_costsProb.push_back(costsProb[i]);
            }
        }
    }
    return perturb_costsProb;
}

TProb WCSP::Gumbel_noise(TProb prob)
{
    double perturb = 0;
    TProb perturb_prob = 0;
    const double a = -0.577215; // -c --> Euler constant so as to standard the law.
    const double b = 1.0;
    std::default_random_engine generator(ToulBar2::seed);
    std::extreme_value_distribution<double> distribution(a, b);
    perturb = distribution(generator);
    cout << perturb << endl;
    if (prob < pow(10, 37)) {
        perturb_prob = prob + perturb;
    } else {
        perturb_prob = prob;
    }
    return perturb_prob;
}
