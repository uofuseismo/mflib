#include <cstdio>
#include <algorithm>
#include <cmath>
#include <vector>

/// Estimates the magnitudes using Equation 1 of Schaff and Richards (2014).
/// @param[in] n   The number of samples.
/// @param[in] x   The reference (master) event.  This is an array whose
///                dimension is [n] and is assumed to be demeaned.
/// @param[in] y   The template.  This is an array who dimension is [n]
///                and is assumed to be demeaned.
/// @result The magnitude of y relative to x. 
/// @author Keith Koper and Ben Baker, January 2020.
template<class T>
T computeSandREqn1(const int n,
                   const T *__restrict__ x,
                   const T *__restrict__ y)
{
    const T zero = 0;
    auto num = zero;
    auto den = zero;
    #pragma omp simd reduction(+:num, den)
    for (int i=0; i<n; ++i)
    {
        num = num + x[i]*y[i];
        den = den + x[i]*x[i];
    }
    return std::log10(num/den);
}

/// Estimates the magnitudes using Equation 11 of Schaff and Richards (2014).
/// @param[in] n   The number of samples.
/// @param[in] x   The reference (master) event.  This is an array whose
///                dimension is [n] and is assumed to be demeaned.
/// @param[in] y   The template.  This is an array whose dimension is
///                [n] and is assumed to be demeaned.
/// @result The magnitude of y relative to x. 
/// @author Keith Koper and Ben Baker, January 2020.
template<class T>
T computeSandREqn11(const int n,
                    const T *__restrict__ x,
                    const T *__restrict__ y)
{
    const T zero = 0;
    auto num = zero;
    auto den = zero;
    #pragma omp simd reduction(+:num, den)
    for (int i=0; i<n; ++i)
    {
        num = num + y[i]*y[i];
        den = den + x[i]*x[i];
    }
    // N.B.: Keith does:
    //   log10(sqrt(num/n)/(sqrt(den/n)))
    // = log10(sqrt(num/den))
    // = 1/2*log10(num/den)
    // The normalization by n is strange since that is not how a magnitude is
    // defined. 
    const T half = 0.5;
    return half*std::log10(num/den);
}


int main()
{
    std::vector<double> x(100, 1), y(100, 0.5);
    printf("%lf\n", computeSandREqn1(x.size(), x.data(), y.data()));
    printf("%lf\n", computeSandREqn11(x.size(), x.data(), y.data())); 
}
