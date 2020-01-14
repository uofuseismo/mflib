#ifndef MFLIB_TESTING_UTILITIES_HPP
#define MFLIB_TESTING_UTILITIES_HPP
#include <vector>

/// Generates uniform random numbers
std::vector<double> generateUniformRandomNumbers(
   int n, double lower =-1, double upper=+1);
/// Computes the infinity norm: max(|r|)
template<class T> T infinityNorm(const int n, const T *__restrict__ r)
{
    auto minmax = std::minmax_element(r, r+n);
    auto linf = std::max(std::abs(*minmax.first), std::abs(*minmax.second));
    return linf;
}
/// Computes the infinity norm: max(|x - y|)
template<class T> T infinityNorm(const int n,
                                 const T *__restrict__ x,
                                 const T *__restrict__ y)
{
    T linf = 0;
    #pragma omp simd reduction(max:linf)
    for (int i=0; i<n; ++i)
    {   
        linf = std::max(linf, std::abs(x[i] - y[i]));
    }   
    return linf;
}

#endif
