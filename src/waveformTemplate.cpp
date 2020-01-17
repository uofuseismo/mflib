#include <cstdio>
#include <cstring>
#include <string>
#include <algorithm>
#if __has_include(<pstl/algorithm>)
#include <pstl/algorithm>
#include <pstl/execution>
#define USE_PSTL 1
#endif
#include <ipps.h>
#include "mflib/waveformTemplate.hpp"

using namespace MFLib;

class WaveformTemplate::WaveformTemplateImpl
{
public:
    WaveformTemplateImpl() = default;
    /// Copy c'tor
    WaveformTemplateImpl(const WaveformTemplateImpl &tplate)
    {
        *this = tplate;
    }
    /// Deep copy of the template
    WaveformTemplateImpl& operator=(const WaveformTemplateImpl &tplate)
    {
        if (&tplate == this){return *this;}
        mShiftAndStackWeight = tplate.mShiftAndStackWeight;
        mSamplingRate = tplate.mSamplingRate;
        mOnsetTime = tplate.mOnsetTime;
        mTravelTime = tplate.mTravelTime;
        // Copy the waveform
        mSignalLength = tplate.mSignalLength;
        if (mSignalLength > 0)
        {
            mSignal = ippsMalloc_64f(mSignalLength);
#ifdef USE_PSTL
            std::copy(std::execution::unseq,
                      tplate.mSignal, tplate.mSignal+mSignalLength, mSignal);
#else
            std::copy(tplate.mSignal, tplate.mSignal+mSignalLength, mSignal);
#endif
        }
        return *this;
    }
    /// Destructor
    ~WaveformTemplateImpl() 
    {
        clear();
    }
    /// Releases memory
    void clear() noexcept
    {
        if (mSignal){ippsFree(mSignal);}
        mSignal = nullptr;
        mSamplingRate = 0;
        mShiftAndStackWeight = 1;
        mOnsetTime =-1;
        mTravelTime =-1;
        mSignalLength = 0;
    }
    
//private:
    /// The template signal.  This has dimension [mSignalLength].
    double *mSignal = nullptr;
    /// The sampling rate in Hz.
    double mSamplingRate = 0;
    /// The shift and stack weight.
    double mShiftAndStackWeight = 1;
    /// The waveform onset time
    double mOnsetTime =-1;
    /// The travel time
    double mTravelTime =-1;
    /// The number of samples in the template waveform signal.
    int mSignalLength = 0;
};

/// Constructor
WaveformTemplate::WaveformTemplate() :
    pImpl(std::make_unique<WaveformTemplateImpl> ())
{
}

/// Copy constructor
WaveformTemplate::WaveformTemplate(const WaveformTemplate &t)
{
    *this = t;
}

/// Move constructor
WaveformTemplate::WaveformTemplate(WaveformTemplate &&t) noexcept
{
    *this = std::move(t);
}

/// Copy assignment operator
WaveformTemplate& WaveformTemplate::operator=(const WaveformTemplate &t)
{
    if (&t == this){return *this;}
    if (pImpl){pImpl->clear();}
    pImpl = std::make_unique<WaveformTemplateImpl> (*t.pImpl);
    return *this;
}

/// Move assignment operator
WaveformTemplate& WaveformTemplate::operator=(WaveformTemplate &&t) noexcept
{
    if (&t == this){return *this;}
    pImpl = std::move(t.pImpl);
    return *this;
}

/// Destructor
WaveformTemplate::~WaveformTemplate() = default;

/// Clears the template
void WaveformTemplate::clear() noexcept
{
    pImpl->clear();
}

/// Sets the sampling rate
void WaveformTemplate::setSamplingRate(const double samplingRate)
{
    if (samplingRate <= 0)
    {
        throw std::invalid_argument("Sampling rate = " 
                                  + std::to_string(samplingRate)
                                  + " must be positive\n");
    }
    pImpl->mOnsetTime =-1;
    pImpl->mSamplingRate = samplingRate;
}

/// Gets the sampling rate
double WaveformTemplate::getSamplingRate() const
{
    if (!haveSamplingRate())
    {
        throw std::runtime_error("Sampling rate not set\n");
    }
    return pImpl->mSamplingRate;
}

/// Determines if the sampling rate is set
bool WaveformTemplate::haveSamplingRate() const noexcept
{
   if (pImpl->mSamplingRate <= 0){return false;}
   return true;
}

/// Sets the shift and stack weight
void WaveformTemplate::setShiftAndStackWeight(const double weight)
{
    if (weight < 0 || weight > 1)
    {
        throw std::invalid_argument("weight = " + std::to_string(weight)
                                  + " must be in range [0,1]\n");
    }
    pImpl->mShiftAndStackWeight = weight;
}

/// Gets the shift and stack weight
double WaveformTemplate::getShiftAndStackWeight() const noexcept
{
    return pImpl->mShiftAndStackWeight;
}

/// Sets the template signal
void WaveformTemplate::setSignal(const int npts, 
                                 const double x[])
{
    if (npts < 1 || x == nullptr)
    {
        if (npts < 1){throw std::invalid_argument("npts must be positive\n");}
        throw std::invalid_argument("x is NULL");
    }
    pImpl->mOnsetTime =-1;
    pImpl->mSignalLength = npts;
    if (pImpl->mSignal){ippsFree(pImpl->mSignal);}
    pImpl->mSignal = ippsMalloc_64f(pImpl->mSignalLength);
#ifdef USE_PSTL
    std::copy(std::execution::unseq, x, x+npts, pImpl->mSignal);
#else
    std::copy(x, x+npts, pImpl->mSignal);
#endif
}

void WaveformTemplate::setSignal(const int npts,
                                 const float x[])
{
    if (npts < 1 || x == nullptr)
    {   
        if (npts < 1){throw std::invalid_argument("npts must be positive\n");}
        throw std::invalid_argument("x is NULL");
    }
    pImpl->mOnsetTime =-1;
    pImpl->mSignalLength = npts;
    if (pImpl->mSignal){ippsFree(pImpl->mSignal);}
    pImpl->mSignal = ippsMalloc_64f(pImpl->mSignalLength);
    double *mSignal = pImpl->mSignal; 
    #pragma omp simd aligned(mSignal: 64)
    for (int i=0; i<npts; ++i){mSignal[i] = static_cast<double> (x[i]);}
}

/// Gets the signal
void WaveformTemplate::getSignal(const int maxx, double *xIn[]) const
{
    auto npts = getSignalLength(); // Throws
    if (maxx < npts)
    {
        throw std::invalid_argument("maxx = " + std::to_string(maxx)
                                  + " must be at least = "
                                  + std::to_string(npts) + "\n");
    }
    double *x = *xIn;
    if (x == nullptr){throw std::invalid_argument("x is NULL\n");}
    const double *__attribute__((aligned(64))) mSignal = pImpl->mSignal;
#ifdef USE_PSTL
    std::copy(std::execution::unseq, mSignal, mSignal+npts, x);
#else
    std::copy(mSignal, mSignal+npts, x);
#endif
}

/// Gets the signal
void WaveformTemplate::getSignal(const int maxx, float *xIn[]) const
{
    auto npts = getSignalLength(); // Throws
    if (maxx < npts)
    {   
        throw std::invalid_argument("maxx = " + std::to_string(maxx)
                                  + " must be at least = "
                                  + std::to_string(npts) + "\n");
    }   
    float *x = *xIn;
    if (x == nullptr){throw std::invalid_argument("x is NULL\n");}
    const double *mSignal = pImpl->mSignal;
    #pragma omp simd aligned(mSignal: 64)
    for (int i=0; i<npts; ++i){x[i] = static_cast<float> (mSignal[i]);}
}

/// Gets the signal length
int WaveformTemplate::getSignalLength() const
{
    if (!haveSignal()){throw std::runtime_error("Signal not set\n");}
    return pImpl->mSignalLength;
}

/// Check if I have the signal
bool WaveformTemplate::haveSignal() const noexcept
{
    if (pImpl->mSignal && pImpl->mSignalLength > 0){return true;}
    return false; 
}

/// Sets the onset time
void WaveformTemplate::setOnsetTime(const double onsetTime)
{
    if (onsetTime < 0)
    {
        throw std::invalid_argument("onsetTime = " 
                                  + std::to_string(onsetTime)
                                  + " must be positive\n");
    }
    double df = getSamplingRate(); // Throws
    int npts = getSignalLength(); // Throws
    double maxOnsetTime = static_cast<double> (npts - 1)*df;
    if (onsetTime > maxOnsetTime)
    {
        throw std::invalid_argument("onsetTime = " 
                                  + std::to_string(onsetTime)
                                  + " cannot exceed maxOnsetTime = "
                                  + std::to_string(maxOnsetTime) + "\n");
    }
    pImpl->mOnsetTime = onsetTime;
}

/// Gets the onset time
double WaveformTemplate::getOnsetTime() const
{
    if (!haveOnsetTime())
    {
        throw std::runtime_error("Onset time not set\n");
    }
    return pImpl->mOnsetTime;
}

/// Do I have the onset time?
bool WaveformTemplate::haveOnsetTime() const noexcept
{
    if (pImpl->mOnsetTime < 0){return false;}
    return true;
}

/// Sets the travel time
void WaveformTemplate::setTravelTime(const double travelTime)
{
    if (travelTime < 0)
    {
        throw std::invalid_argument("Travel time cannot be negative\n");
    }
    pImpl->mTravelTime = travelTime;
}

/// Gets the travel time
double WaveformTemplate::getTravelTime() const
{
    if (!haveTravelTime())
    {
        throw std::invalid_argument("Travel time not yet set\n");
    }
    return pImpl->mTravelTime;
}

/// Do I have the travel time?
bool WaveformTemplate::haveTravelTime() const noexcept
{
    if (pImpl->mTravelTime < 0){return false;}
    return true;
}
