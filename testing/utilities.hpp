#ifndef MFLIB_TESTING_UTILITIES_HPP
#define MFLIB_TESTING_UTILITIES_HPP
#include <vector>
#include <numeric>
#include <algorithm>
#include <functional>
#include <cmath>

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

template<class T>
std::vector<T> naivePearsonCorrelation(const int nb, const T *__restrict__ b,
                                       const int nx, const T *__restrict__ x)
{
    // Allocate workspace and initialize result
    const T zero = 0;
    std::vector<T> btemp(nb);
    std::vector<T> xtemp(nb);
    std::vector<T> xc(std::max(0, nx-nb+1), zero);
    // Detrend b
    T bmean = std::accumulate(b, b+nb, zero)/static_cast<T> (nb);
    std::transform(b, b+nb, btemp.begin(),
                   [&](auto value){return value - bmean;});
                   //std::bind2nd(std::minus<T>(), bmean));
    // Compute the L2 norm of the template
    T Eb = std::inner_product(btemp.begin(), btemp.end(),
                              btemp.begin(), zero);
    Eb = std::sqrt(Eb);
    for (int i=0; i<nx-nb+1; ++i)
    {   
        // Detrend the signal in this window
        T xmean = std::accumulate(x+i, x+i+nb, zero)/static_cast<T> (nb);
        std::transform(x+i, x+i+nb, xtemp.begin(), 
                       [&](auto value){return value - xmean;});
                       //std::bind2nd(std::minus<T>(), xmean));
        // Compute the L2 norm of the signal
        T Ex = std::inner_product(xtemp.begin(), xtemp.end(),
                                  xtemp.begin(), zero);
        Ex = std::sqrt(Ex);
        Ex = std::max(std::numeric_limits<T>::epsilon(), Ex);
        // Compute the inner product in the numerator then normalize
        xc[i] = std::inner_product(btemp.begin(), btemp.end(),
                                   xtemp.begin(), zero);
        xc[i] = xc[i]/(Eb*Ex);
    }
    return xc;
}

#endif
