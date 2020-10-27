#include <cstdio>
#include <cstring>
#include <string>
#include <algorithm>
#ifdef HAVE_SOURCE_FILE_INCLUSION
 #if __has_include(<pstl/algorithm>)
  #include <pstl/algorithm>
  #include <pstl/execution>
  #define USE_PSTL 1
 #endif
#endif
#include <mkl.h>
#include "mflib/waveformTemplate.hpp"
#include "mflib/networkStationPhase.hpp"
#include "private/private.hpp"

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
        mPhaseOnsetTime = tplate.mPhaseOnsetTime;
        mTravelTime = tplate.mTravelTime;
        mMagnitude = tplate.mMagnitude;
        mIdentifier = tplate.mIdentifier;
        mPolarity = tplate.mPolarity;
        mHaveIdentifier = tplate.mHaveIdentifier; 
        mHaveMagnitude = tplate.mHaveMagnitude;
        // Copy the waveform
        mSignalLength = tplate.mSignalLength;
        if (mSignalLength > 0)
        {
            mSignal = reinterpret_cast<double *>
                      (mkl_calloc(mSignalLength*sizeof(double), 1, 64));
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
        if (mSignal){mkl_free(mSignal);}
        mSignal = nullptr;
        mSamplingRate = 0;
        mShiftAndStackWeight = 1;
        mPhaseOnsetTime =-1;
        mTravelTime =-1;
        mMagnitude = 0;
        mPolarity = MFLib::Polarity::UNKNOWN;
        mSignalLength = 0;
        mIdentifier.first.clear();
        mIdentifier.second = 0;
        mHaveMagnitude = false;
        mHaveIdentifier = false;
    }
    
//private:
    /// The template signal.  This has dimension [mSignalLength].
    double *mSignal = nullptr;
    /// The sampling rate in Hz.
    double mSamplingRate = 0;
    /// The shift and stack weight.
    double mShiftAndStackWeight = 1;
    /// The waveform onset time.
    double mPhaseOnsetTime =-1;
    /// The travel time.
    double mTravelTime =-1;
    /// The magnitude.
    double mMagnitude = 0;
    /// The waveform identifier.
    std::pair<MFLib::NetworkStationPhase, uint64_t> mIdentifier;
    /// The number of samples in the template waveform signal.
    int mSignalLength = 0;
    /// The onset's polarity.
    MFLib::Polarity mPolarity = MFLib::Polarity::UNKNOWN;
    /// Determines if the waveform identifier was set.
    bool mHaveIdentifier = false;
    /// Determines if the magnitude was set.
    bool mHaveMagnitude = false;
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
    pImpl->mPhaseOnsetTime =-1;
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
    pImpl->mPhaseOnsetTime =-1;
    pImpl->mSignalLength = npts;
    if (pImpl->mSignal){mkl_free(pImpl->mSignal);}
    pImpl->mSignal = reinterpret_cast<double *>
                     (mkl_calloc(pImpl->mSignalLength*sizeof(double), 1, 64));
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
    pImpl->mPhaseOnsetTime =-1;
    pImpl->mSignalLength = npts;
    if (pImpl->mSignal){mkl_free(pImpl->mSignal);}
    pImpl->mSignal = reinterpret_cast<double *>
                     (mkl_calloc(pImpl->mSignalLength*sizeof(double), 1, 64));
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
void WaveformTemplate::setPhaseOnsetTime(const double onsetTime)
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
    pImpl->mPhaseOnsetTime = onsetTime;
}

/// Gets the onset time
double WaveformTemplate::getPhaseOnsetTime() const
{
    if (!havePhaseOnsetTime())
    {
        throw std::runtime_error("Onset time not set\n");
    }
    return pImpl->mPhaseOnsetTime;
}

/// Do I have the onset time?
bool WaveformTemplate::havePhaseOnsetTime() const noexcept
{
    if (pImpl->mPhaseOnsetTime < 0){return false;}
    return true;
}

/// Sets the travel time
void WaveformTemplate::setTravelTime(const double travelTime)
{
    if (travelTime < 0)
    {
        throw std::runtime_error("Travel time cannot be negative\n");
    }
    pImpl->mTravelTime = travelTime;
}

/// Gets the travel time
double WaveformTemplate::getTravelTime() const
{
    if (!haveTravelTime())
    {
        throw std::runtime_error("Travel time not yet set\n");
    }
    return pImpl->mTravelTime;
}

/// Do I have the travel time?
bool WaveformTemplate::haveTravelTime() const noexcept
{
    if (pImpl->mTravelTime < 0){return false;}
    return true;
}

/// Sets the identifier
std::pair<MFLib::NetworkStationPhase, uint64_t> 
WaveformTemplate::getIdentifier() const
{
    if (!haveIdentifier())
    {
        throw std::runtime_error("Waveform identifier not set");
    }
    return pImpl->mIdentifier;
}

/// Gets the identifier
void WaveformTemplate::setIdentifier(
    const std::pair<MFLib::NetworkStationPhase, uint64_t> &id) noexcept
{
    pImpl->mIdentifier = id;
    pImpl->mHaveIdentifier = true;
}

/// Do I have the identifier?
bool WaveformTemplate::haveIdentifier() const noexcept
{
    return pImpl->mHaveIdentifier;
}

/// Sets the magnitude 
double WaveformTemplate::getMagnitude() const
{
    if (!haveMagnitude())
    {
        throw std::runtime_error("Magnitude not set");
    }
    return pImpl->mMagnitude;
}

/// Gets the magnitude 
void WaveformTemplate::setMagnitude(const double mag) noexcept
{
    pImpl->mMagnitude = mag;
    pImpl->mHaveMagnitude = true;
}

/// Do I have the magnitude?
bool WaveformTemplate::haveMagnitude() const noexcept
{
    return pImpl->mHaveMagnitude;
}

/// Sets the polarity
void WaveformTemplate::setPolarity(const MFLib::Polarity polarity) noexcept
{
    pImpl->mPolarity = polarity;
}

/// Get the polarity
MFLib::Polarity WaveformTemplate::getPolarity() const noexcept
{
    return pImpl->mPolarity;
}
