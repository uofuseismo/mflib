#include <cstdio>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <functional>
#include <vector>
#if __has_include(<pstl/algorithm>)
#include <pstl/algorithm>
#include <pstl/execution>
#define USE_PSTL 1
#endif
#include <mkl.h>
#include "mflib/channelMagnitude.hpp"

using namespace MFLib;

namespace
{
/// Estimates the magnitudes using Equation 1 of Schaff and Richards (2014).
/// @param[in] n   The number of samples.
/// @param[in] x   The reference (master) event.  This is an array whose
///                dimension is [n] and is assumed to be demeaned.
/// @param[in] y   The observed signal.  This is an array who dimension is [n]
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
/// @param[in] y   The observed signal.  This is an array whose dimension is
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
}

/// Implementation.  Since detections likely will be the exception rather
/// than the rule this class does a slower but more accurate double 
/// precision computation.
template<class T>
class ChannelMagnitude<T>::ChannelMagnitudeImpl
{
    void clear()
    {
        if (mX){MKL_free(mX);}
        if (mY){MKL_free(mY);}
        mX = nullptr;
        mY = nullptr;
    }
public:
    /// The template signal.  This has dimension [mLengthX].
    double *__attribute__((aligned(64))) mX = nullptr;
    /// The observed signal.  This has dimension [mLengthY].
    double *__attribute__((aligned(64))) mY = nullptr;
    int mLengthX = 0;
    int mLengthY = 0;
    bool mHaveTemplateWaveform = false;
    bool mHaveDetectedWaveform = false;    
};

/// Constructor
template<class T>
ChannelMagnitude<T>::ChannelMagnitude() :
    pImpl(ChannelMagnitudeImpl())
{
}

/// Destructor
template<class T>
ChannelMagnitude<T>::~ChannelMagnitude() = default;

/// Sets the template waveform
template<>
void ChannelMagnitude<double>::setTemplateWaveform(
    const int n, const double *__restrict__ x)
{
    pImpl->mHaveTemplateWaveform = false;
    if (n < 1 || x == nullptr)
    {
        if (n < 1){throw std::invalid_argument("No samples in waveform\n");}
        throw std::invalid_argument("Waveform is NULL\n");
    }
    // Reallocate
    if (pImpl->mX && n > pImpl->mLengthX)
    {
        MKL_free(pImpl->mX);
        pImpl->mX = nullptr;
    }
    if (!pImpl->mX)
    {
        auto len = static_cast<size_t> (n)*sizeof(double);
        pImpl->mX = static_cast<double *> (MKL_calloc(len, 1, 64));
    }
    // Remove mean and copy
    constexpr double zero = 0;
    auto mean = std::accumulate(x, x+n, zero)/static_cast<double> (n);
#ifdef USE_PSTL
    std::transform(std::execution::unseq, x, x+n, pImpl->mX,
                   [&](auto value){return value - mean;});
#else
    std::transform(x, x+n, pImpl->mX,
                   [&](auto value){return value - mean;});
#endif
    pImpl->mLengthX = n;
    pImpl->mHaveTemplateWaveform = true;
}

/// Is the template waveform set?
template<class T>
bool ChannelMagnitude<T>::haveTemplateWaveform() const noexcept
{
    return pImpl->mHaveTemplateWaveform;
}

/// Sets the detected waveform
template<>
void ChannelMagnitude<double>::setDetectedWaveform(
    const int n, const double *__restrict__ y)
{
    pImpl->mHaveDetectedWaveform = false;
    if (n < 1 || y == nullptr)
    {
        if (n < 1){throw std::invalid_argument("No samples in waveform\n");}
        throw std::invalid_argument("Waveform is NULL\n");
    }
    // Reallocate
    if (pImpl->mY && n > pImpl->mLengthY)
    {
        MKL_free(pImpl->mY);
        pImpl->mY = nullptr;
    }
    if (!pImpl->mY)
    {
        auto len = static_cast<size_t> (n)*sizeof(double);
        pImpl->mY = static_cast<double *> (MKL_calloc(len, 1, 64));
    }
    // Remove mean and copy
    constexpr double zero = 0;
    auto mean = std::accumulate(y, y+n, zero)/static_cast<double> (n);
#ifdef USE_PSTL
    std::transform(std::execution::unseq, y, y+n, pImpl->mY,
                   [&](auto value){return value - mean;});
#else
    std::transform(y, y+n, pImpl->mY,
                   [&](auto value){return value - mean;});
#endif
    pImpl->mLengthY = n;
    pImpl->mHaveDetectedWaveform = true;        
}

template<>
void ChannelMagnitude<float>::setDetectedWaveform(
    const int n, const float *__restrict__ y)
{
    pImpl->mHaveDetectedWaveform = false;
    if (n < 1 || y == nullptr)
    {
        if (n < 1){throw std::invalid_argument("No samples in waveform\n");}
        throw std::invalid_argument("Waveform is NULL\n");
    }
    // Reallocate
    if (pImpl->mY && n > pImpl->mLengthY)
    {
        MKL_free(pImpl->mY);
        pImpl->mY = nullptr;
    }
    if (!pImpl->mY)
    {
        auto len = static_cast<size_t> (n)*sizeof(double);
        pImpl->mY = static_cast<double *> (MKL_calloc(len, 1, 64));
    }
    // Remove mean and copy
    constexpr float zero = 0;
    auto mean = static_cast<float>
                (std::accumulate(y, y+n, zero)/static_cast<double> (n));
    auto mYPtr = pImpl->mY;
    #pragma omp simd aligned(mYPtr: 64)
    for (int i=0; i<n; ++i)
    {
        mYPtr[i] = static_cast<double> (y[i] - mean);
    }    
    pImpl->mLengthY = n;
    pImpl->mHaveDetectedWaveform = true;    
}

/// Is the detected waveform set?
template<class T>
bool ChannelMagnitude<T>::haveDetectedWaveform() const noexcept
{
    return pImpl->mHaveDetectedWaveform;
}


/*
int main()
{
    std::vector<double> x(100, 1), y(100, 0.5);
    printf("%lf\n", computeSandREqn1(x.size(), x.data(), y.data()));
    printf("%lf\n", computeSandREqn11(x.size(), x.data(), y.data())); 
}
*/
