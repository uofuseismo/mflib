#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <random>
#include <numeric>
#include "utilities.hpp"

std::vector<double> generateUniformRandomNumbers(
   const int n, const double lower, const double upper)
{
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<> dist(lower, upper);
    std::vector<double> result(n, 0);
    for (auto &v : result)
    {
        v = dist(rng);
    }
    return result;
}

