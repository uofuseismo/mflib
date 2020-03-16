#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <functional>
#include <vector>
#if __has_include(<pstl/algorithm>)
#include <pstl/algorithm>
#include <pstl/numeric>
#include <pstl/execution>
#define USE_PSTL 1
#else
#include <numeric>
#include <algorithm>
#endif
#include <mkl.h>
#include "mflib/waveformTemplate.hpp"
#include "mflib/singleChannel/relativeMagnitude.hpp"

using namespace MFLib::SingleChannel;

namespace
{
/// @brief Estimates the scaling factor from Gibbons and Rigdal, 2006.
/// @param[in] magx2   The magnitude of the template signal squared.
/// @param[in] x       The detrended template event. This is an array whose
///                    dimension is [n].
/// @param[in] y       The observed signal.  This is an array who dimension is [n]
///                and is assumed to be demeaned.
/// @result The magnitude of y relative to x. 
/// @author Keith Koper and Ben Baker, January 2020.
double gibbonsRingdal2006(const int n,
                          const double magx2,
                          const double *__restrict__ x,
                          const double *__restrict__ y)
{
    constexpr double zero = 0;
#ifdef USE_PSTL
    auto xy = std::transform_reduce(std::execution::unseq, x, x+n, y, zero);
#else
    auto xy = std::inner_product(x, x+n, y, zero);
#endif
    return std::abs(xy)/magx2; // abs b/c I don't want a sign flip
}

/// @brief Estimates the scaling factor from Equation 10 of
///        Schaff and Richards (2014).
/// @param[in] magx  The length of the detrended template signal.
/// @param[in] magy  The length of the detrended detected signal.
/// @result The magnitude of y relative to x. 
/// @author Keith Koper and Ben Baker, January 2020.
double schaffRichards2014(const double magx2, const double magy2)
{
    return std::sqrt(magy2/magx2);
}

double demeanCopyComputeNorm2(const int n,
                              const double *__restrict__ y,
                              double *__restrict__ ydemeaned)
{
    // Remove mean and copy
    constexpr double zero = 0;
#ifdef USE_PSTL
    auto mean = std::reduce(std::execution::unseq, y, y+n, zero)
               /static_cast<double> (n);
#else
    auto mean = std::accumulate(y, y+n, zero)
               /static_cast<double> (n);
#endif
    double xsum = 0;
    #pragma omp simd reduction(+:xsum)
    for (int i=0; i<n; ++i)
    {
        ydemeaned[i] = y[i] - mean;
        xsum = xsum + ydemeaned[i]*ydemeaned[i];
    }
    return xsum;
}

double demeanNorm2(const int n, double *__restrict__ y)
{
    // Remove mean and copy
    constexpr double zero = 0;
#ifdef USE_PSTL
    auto mean = std::reduce(std::execution::unseq, y, y+n, zero)
               /static_cast<double> (n);
#else
    auto mean = std::accumulate(y, y+n, zero)
               /static_cast<double> (n);
#endif
    double xsum = 0;
    #pragma omp simd reduction(+:xsum)
    for (int i=0; i<n; ++i)
    {   
        y[i] = y[i] - mean;
        xsum = xsum + y[i]*y[i];
    }   
    return xsum;
}


double demeanCopyComputeNorm2(const int n,
                              const float *__restrict__ y,
                              double *__restrict__ ydemeaned)
{
    // Remove mean and copy
    constexpr float zero = 0;
#ifdef USE_PSTL
    auto sum = std::reduce(std::execution::unseq, y, y+n, zero);
#else
    auto sum = std::accumulate(y, y+n, zero);
#endif
    double mean = static_cast<double> (sum);
    double xsum = 0;
    #pragma omp simd reduction(+:xsum)
    for (int i=0; i<n; ++i)
    {
        ydemeaned[i] = static_cast<double> (y[i]) - mean;
        xsum = xsum + ydemeaned[i]*ydemeaned[i];
    }
    return xsum;
}

}

/// Implementation.  Since detections likely will be the exception rather
/// than the rule this class does a slower but more accurate double 
/// precision computation.
template<class T>
class RelativeMagnitude<T>::RelativeMagnitudeImpl
{
public:
    RelativeMagnitudeImpl() = default;
    ~RelativeMagnitudeImpl()
    {
        clear();
    }
    RelativeMagnitudeImpl(const RelativeMagnitudeImpl &impl)
    {
        *this = impl;
    }
    RelativeMagnitudeImpl& operator=(const RelativeMagnitudeImpl &impl)
    {
        if (&impl == this){return *this;}
        mXL22 = impl.mXL22;
        mYL22 = impl.mYL22;
        mGibbonsRingdalAlpha = impl.mGibbonsRingdalAlpha;
        mSchaffRichardsAlpha = impl.mSchaffRichardsAlpha;
        mSignalLength = impl.mSignalLength;
        mHaveDetectedSignal = impl.mHaveDetectedSignal;
        mHaveGibbonsRingdalAlpha = impl.mHaveGibbonsRingdalAlpha;
        mHaveSchaffRichardsAlpha = impl.mHaveSchaffRichardsAlpha;
        mInitialized = impl.mInitialized;
        if (mSignalLength > 0 && impl.mX)
        {
            auto len = static_cast<size_t> (mSignalLength)*sizeof(double);
            mX = static_cast<double *> (MKL_calloc(len, 1, 64));
            std::copy(impl.mX, impl.mX+mSignalLength, mX);
        }
        if (mSignalLength > 0 && impl.mY)
        {
            auto len = static_cast<size_t> (mSignalLength)*sizeof(double);
            mY = static_cast<double *> (MKL_calloc(len, 1, 64));
            std::copy(impl.mY, impl.mY+mSignalLength, mY);
        }
        return *this;
    }
    void clear() noexcept
    {
        nullifyAlpha();
        if (mX){MKL_free(mX);}
        if (mY){MKL_free(mY);}
        mX = nullptr;
        mY = nullptr;
        mXL22 = 0;
        mYL22 = 0;
        mSignalLength = 0;
        mHaveDetectedSignal = false;
        mInitialized = false;
    }
    void nullifyAlpha() noexcept
    {
        mGibbonsRingdalAlpha = 1;
        mSchaffRichardsAlpha = 1;
        mHaveGibbonsRingdalAlpha = false;
        mHaveSchaffRichardsAlpha = false;
    }
   
    /// The demeaned template signal.  This has dimension [mSignalLength].
    double *__attribute__((aligned(64))) mX = nullptr;
    /// The demeaned observed signal.  This has dimension [mSignalLength].
    double *__attribute__((aligned(64))) mY = nullptr;
    /// The squared L2 norm of mX 
    double mXL22 = 0;
    /// The squared L2 norm of mY
    double mYL22 = 0;
    /// The scaling factor in Gibbons and Ringdal, 2006
    double mGibbonsRingdalAlpha = 1;
    /// The scaling factor in Schaff and Richards, 2006
    double mSchaffRichardsAlpha = 1;
    /// The length of the template and detected signal.
    int mSignalLength = 0;
    bool mHaveDetectedSignal = false;
    bool mHaveGibbonsRingdalAlpha = false;
    bool mHaveSchaffRichardsAlpha = false;
    bool mInitialized = false;
};

/// Constructor
template<class T>
RelativeMagnitude<T>::RelativeMagnitude() :
    pImpl(std::make_unique<RelativeMagnitudeImpl> ())
{
}

/// Copy c'tor
template<class T>
RelativeMagnitude<T>::RelativeMagnitude(const RelativeMagnitude &mag)
{
    *this = mag;
}

/// Move c'tor
template<class T>
RelativeMagnitude<T>::RelativeMagnitude(RelativeMagnitude &&mag) noexcept
{
    *this = std::move(mag);
}

/// Copy assignment operator
template<class T>
RelativeMagnitude<T>&
RelativeMagnitude<T>::operator=(const RelativeMagnitude &mag)
{
    if (&mag == this){return *this;}
    if (pImpl){pImpl->clear();}
    pImpl = std::make_unique<RelativeMagnitude<T>::RelativeMagnitudeImpl>
            (*mag.pImpl);
    return *this;
}

/// Move assignment operator
template<class T>
RelativeMagnitude<T>&
RelativeMagnitude<T>::operator=(RelativeMagnitude &&mag) noexcept
{
    if (&mag == this){return *this;}
    pImpl = std::move(mag.pImpl);
    return *this;
}

/// Destructor
template<class T>
RelativeMagnitude<T>::~RelativeMagnitude() = default;

/// Clears the class
template<class T>
void RelativeMagnitude<T>::clear() noexcept
{
    pImpl->clear();
}

/// Initializes the class
template<class T>
void RelativeMagnitude<T>::initialize(const WaveformTemplate &wt)
{
    clear();
    pImpl->mSignalLength = wt.getSignalLength();
    if (pImpl->mSignalLength < 1)
    {
        throw std::invalid_argument("Template signal length = " 
                                  + std::to_string(pImpl->mSignalLength)
                                  + " must be positive\n");
    }
    // Set space
    auto len = static_cast<size_t> (pImpl->mSignalLength)*sizeof(double);
    pImpl->mX = static_cast<double *> (MKL_calloc(len, 1, 64));
    pImpl->mY = static_cast<double *> (MKL_calloc(len, 1, 64));
    // Copy over the signal
    wt.getSignal(pImpl->mSignalLength, &pImpl->mX);
    // Demean the signal then compute norm
    pImpl->mXL22 = demeanNorm2(pImpl->mSignalLength, pImpl->mX);
    if (pImpl->mXL22 == std::numeric_limits<double>::epsilon()*100)
    {
        if (pImpl->mXL22 == 0)
        {
            throw std::invalid_argument("Template signal cannot be constant\n");
        }
        fprintf(stderr, "RelativeMagnitude: Template signal may be dead\n");
    }
    pImpl->mInitialized = true;
}

/// Gets the signal length
template<class T>
int RelativeMagnitude<T>::getDetectedSignalLength() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized\n");}
    return pImpl->mSignalLength;
}

/// Is the class initialized?
template<class T>
bool RelativeMagnitude<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Sets the detected waveform
template<class T>
void RelativeMagnitude<T>::setDetectedSignal(
    const int n, const T *__restrict__ y)
{
    pImpl->mHaveDetectedSignal = false;
    pImpl->nullifyAlpha();
    auto nref = getDetectedSignalLength();
    if (n != nref || y == nullptr)
    {
        if (n != nref)
        {
            throw std::invalid_argument("Signal length = " + std::to_string(n)
                                      + " must = " + std::to_string(nref) 
                                      + "\n");
        }
        throw std::invalid_argument("Signal is NULL\n");
    }
    // Remove mean, copy, compute squared L2 norm of demeaned signal
    pImpl->mYL22 = demeanCopyComputeNorm2(n, y, pImpl->mY);
    if (pImpl->mYL22 < std::numeric_limits<double>::epsilon()*100)
    {
        if (pImpl->mYL22 == 0)
        {
            throw std::invalid_argument("Detected signal cannot be constant\n");
        }
        fprintf(stderr, "RelativeMagnitude: Detected signal may be dead\n");
    }
    pImpl->mHaveDetectedSignal = true;
}

/// Is the detected waveform set?
template<class T>
bool RelativeMagnitude<T>::haveDetectedSignal() const noexcept
{
    return pImpl->mHaveDetectedSignal;
}

/// Compute the scaling factor
template<class T>
T RelativeMagnitude<T>::computeAmplitudeScalingFactor(
    const MFLib::RelativeMagnitudeType type) const
{
    if (!isInitialized())
    {
        throw std::runtime_error("Class is not initialized\n");
    }
    // Check the waveforms exist
    if (!pImpl->mHaveDetectedSignal)
    {
        throw std::runtime_error("Detected waveform not set\n");
    }
    // Compute the magnitudes
    if (!pImpl->mHaveGibbonsRingdalAlpha)
    {
        pImpl->mGibbonsRingdalAlpha
            = gibbonsRingdal2006(pImpl->mSignalLength,
                                 pImpl->mXL22, pImpl->mX, pImpl->mY);
        pImpl->mHaveGibbonsRingdalAlpha = true;
    }
    if (!pImpl->mHaveSchaffRichardsAlpha)
    {
        pImpl->mSchaffRichardsAlpha
           = schaffRichards2014(pImpl->mXL22, pImpl->mYL22);
        pImpl->mHaveSchaffRichardsAlpha = true;
    }
    // Return what the user wants
    if (type == MFLib::RelativeMagnitudeType::GIBBONS_RINGDAL_2006)
    {
        return static_cast<T> (pImpl->mGibbonsRingdalAlpha);
    }
    else
    {
        return static_cast<T> (pImpl->mSchaffRichardsAlpha);
    }
}

/// Compute the relative magnitude
template<class T>
T RelativeMagnitude<T>::computeMagnitudePerturbation(
    const MFLib::RelativeMagnitudeType type) const
{
    auto alpha = computeAmplitudeScalingFactor(type); 
    return convertAmplitudeScalingFactorToMagnitudePerturbation(alpha);
}

template<>
double MFLib::SingleChannel::convertAmplitudeScalingFactorToMagnitudePerturbation(const double alpha)
{
    if (alpha <= 0)
    { 
        throw std::invalid_argument("Scaling factor = "
                                  + std::to_string(alpha)
                                  + " must be positive\n");
    }
    return std::log10(alpha);
}

template<>
float MFLib::SingleChannel::convertAmplitudeScalingFactorToMagnitudePerturbation(const float alpha)
{
    if (alpha <= 0)
    {   
        throw std::invalid_argument("Scaling factor = "
                                  + std::to_string(alpha)
                                  + " must be positive\n");
    }   
    return std::log10(alpha);
}

/// Template class instantiation
template class MFLib::SingleChannel::RelativeMagnitude<double>;
template class MFLib::SingleChannel::RelativeMagnitude<float>;
