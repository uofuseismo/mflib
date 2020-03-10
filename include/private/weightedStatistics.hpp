#ifndef MFLIB_PRIVATE_WEIGHTED_STATISTICS_HPP
#define MFLIB_PRIVATE_WEIGHTED_STATISTICS_HPP
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <numeric>
namespace
{
/*!
 * @brief Computes the median.
 * @param[in] n    The number of elements in x.
 * @param[in] x    The observations (or residuals).  This is array whose
 *                 dimension is [n].
 * @result The median of x.
 */
template<typename T>
T median(const int n, const T x[])
{
    if (n == 0){return 0;} // Driver routine should check
    if (n == 1){return x[0];}
    std::vector<T> xwork(x, x+n);
    std::sort(xwork.begin(), xwork.end());
    if (n%2 == 1)
    {
        return xwork[n/2];
    }
    else
    {
        return (xwork[n/2-1] + xwork[n/2])/2;
    }
}
/*!
 * @brief Computes the weighted median.
 * @param[in] n        The number of elements in x and weights.
 * @param[in] x        The observations (or residuals).  This is array whose
 *                     dimension is [n].
 * @param[in] weights  The weight for each observation (or residual).  This is
 *                     an array whose dimension is [n].
 * @result The weighted median of x.
 */
template<typename T>
T weightedMedian(const int n, const T x[], const T weights[])
{
    if (n == 0){return 0;} // Driver routine should check
    if (n == 1){return x[0];}
    if (weights == nullptr){return median(n, x);}
    const T zero = 0;
    T wnorm = std::accumulate(weights, weights + n, zero);
    // Get sorted weights based on x's
    std::vector<std::pair<T, T>> pairs(n);
    #pragma omp simd
    for (int i=0; i<n; ++i)
    {   
        pairs[i] = std::make_pair(x[i], weights[i]/wnorm);
    }   
    std::sort(pairs.begin(), pairs.end(),
              [](const std::pair<T, T> &a, const std::pair<T, T> &b) 
              {
                 return a.first < b.first;
              }
             );
    // Compute the CDF (of weights)
    std::vector<T> Fx(n);
    Fx[0] = pairs[0].second;
    for (int i=1; i<n; ++i)
    {   
        Fx[i] = Fx[i-1] + pairs[i].second;
    }
    // Early exit? 
    const T half = static_cast<T> (1./2.);
    if (Fx[0] > half){return pairs[0].first;} // Probably can't happen (n == 1)
    // Find the first index that exceeds threshold
    for (int i=1; i<n; ++i)
    {   
        if (Fx[i] >= half)
        {
           auto xl = pairs[i-1].first;
           auto xr = pairs[i].first; 
           if (Fx[i] == half){return xr;} // Exact
           return xl + (xr - xl)*(half - Fx[i-1])/(Fx[i] - Fx[i-1]);
        }
    }
    return pairs[0].first; // Shouldn't reach this point
}
/*!
 * @brief Computes the mean.
 * @param[in] n  The number of elements in x.
 */ 
template<typename T>
T mean(const int n, const T x[])
{
    if (n < 1){return 0;}
    const T zero = 0;
    return std::accumulate(x, x+n, zero)/static_cast<T> (n);
}
/*!
 * @brief Computes the weighted mean.
 * @param[in] n        The number of elements in x and weights.
 * @param[in] x        The observations (or residuals).  This is array whose
 *                     dimension is [n].
 * @param[in] weights  The weight for each observation (or residual).  This is
 *                     an array whose dimension is [n].
 * @result The weighted mean of x.
 */
template<typename T>
T weightedMean(const int n, const T x[], const T weights[])
{
    const T zero = 0;
    if (n == 0){return 0;}
    if (n == 1){return x[0];}
    if (weights == nullptr){return mean(n, x);}
    T wnorm = std::accumulate(weights, weights + n, zero);
    return std::inner_product(x, x+n, weights, zero)/wnorm;
}
}
#endif
