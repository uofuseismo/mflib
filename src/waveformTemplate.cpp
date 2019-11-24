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
    /// Deep copy of the template
    WaveformTemplateImpl& operator=(const WaveformTemplateImpl &tplate)
    {
        mShiftAndStackWeight = tplate.mShiftAndStackWeight;
        mSamplingRate = tplate.mSamplingRate;
        // Copy the waveform
        mSignalLength = tplate.mSignalLength;
        if (tplate.mSignalLength > 0)
        {
            auto len = static_cast<size_t> (mSignalLength);
            mSignal = static_cast<double *> (calloc(len, sizeof(double)));
            std::copy(tplate.mSignal, tplate.mSignal+len, mSignal);
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
        if (mSignal){free(mSignal);}
        mSignal = nullptr;
        mSamplingRate = 0;
        mShiftAndStackWeight = 1;
        mSignalLength = 0;
    }
    
//private:
    /// The template signal.  This has dimension [mSignalLength].
    double *mSignal = nullptr;
    /// The sampling rate in Hz.
    double mSamplingRate = 0;
    /// The shift and stack weight.
    double mShiftAndStackWeight = 1;
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

/// 
