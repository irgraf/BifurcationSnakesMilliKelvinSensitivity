//  RandomNumberGenerator.h

#ifndef RandomNumberGenerator_h
#define RandomNumberGenerator_h

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <random>
#include <iostream>
#include <type_traits>

#define short int
typedef double Float;

using namespace std;


extern mt19937 mt;                                          // random number generator defined in main.cpp

inline int BinomialDist(int N, double prob)           // creates a binomially distributed number with number of trials N and success probability prob
{
    binomial_distribution<int> distributionbinom(N,prob);
    return (int) distributionbinom(mt);
}

inline double StandardNormalDist()           // creates a Gaussian distributed number with unit variance and zero mean
{
    normal_distribution<double> distributionnormal(0.,1.);
    return (double) distributionnormal(mt);
}


#endif /* RandomNumberGenerator_h */
