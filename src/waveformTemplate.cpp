#include <cstdio>
#include <cstring>
#include <string>
#include <algorithm>
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
        mShiftAndStackWeight = tplate.mShiftAndStackWeight;
        mSamplingRate = tplate.mSamplingRate;
        mOnsetTime = tplate.mOnsetTime;
        mTravelTime = tplate.mTravelTime;
        // Copy the waveform
        mSignalLength = tplate.mSignalLength;
        if (mSignalLength > 0)
        {
            mSignal = ippsMalloc_64f(mSignalLength);
            ippsCopy_64f(tplate.mSignal, mSignal, mSignalLength);
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
    ippsCopy_64f(x, pImpl->mSignal, npts);
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
    ippsConvert_32f64f(x, pImpl->mSignal, npts);
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
    ippsCopy_64f(pImpl->mSignal, x, npts);
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
    ippsConvert_64f32f(pImpl->mSignal, x, npts);
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
    int index = static_cast<int> (onsetTime*df + 0.5);
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
