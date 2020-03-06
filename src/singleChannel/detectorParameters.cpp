#include <cstdio>
#include <cstdlib>
#include <string>
#include "mflib/singleChannel/detectorParameters.hpp"

#define DEFAULT_THRESHOLD 0.7
#define DEFAULT_MAXIMUM_POLICY MaximumMatchedFilterPolicy::MAXIMUM
#define DEFAULT_CORRELOGRAM_POLICY DetectionWaveformPolicy::REDUCED_CORRELOGRAMS
#define DEFAULT_MIN_SAMPLE_DISTANCE 0

using namespace MFLib::SingleChannel;

class DetectorParameters::DetectorParametersImpl
{
public:
    double mThreshold = DEFAULT_THRESHOLD;
    MaximumMatchedFilterPolicy mMaxPolicy = DEFAULT_MAXIMUM_POLICY;
    DetectionWaveformPolicy mDetectionWaveformPolicy = DEFAULT_CORRELOGRAM_POLICY;
    int mMinSampleDistance = DEFAULT_MIN_SAMPLE_DISTANCE;
    bool mWantDetectedWaveform = false;
    bool mWantAmplitude = true;
};

/// Constructor
DetectorParameters::DetectorParameters() :
    pImpl(std::make_unique<DetectorParametersImpl> ())
{
}

/// Copy constructor
DetectorParameters::DetectorParameters(const DetectorParameters &parameters)
{
    *this = parameters;
}

/// Move constructor
DetectorParameters::DetectorParameters(DetectorParameters &&parameters) noexcept
{
    *this = std::move(parameters);
}

/// Copy assignment
DetectorParameters&
DetectorParameters::operator=(const DetectorParameters &parameters)
{
    if (&parameters == this){return *this;}
    pImpl = std::make_unique<DetectorParametersImpl> (*parameters.pImpl);
    return *this;
}

/// Move assignment
DetectorParameters&
DetectorParameters::operator=(DetectorParameters &&parameters) noexcept
{
    if (&parameters == this){return *this;}
    pImpl = std::move(parameters.pImpl);
    return *this;
}

/// Destructor
DetectorParameters::~DetectorParameters() = default;

/// Resets the class
void DetectorParameters::clear() noexcept
{
    pImpl->mThreshold = DEFAULT_THRESHOLD;
    pImpl->mMaxPolicy = DEFAULT_MAXIMUM_POLICY;
    pImpl->mMinSampleDistance = DEFAULT_MIN_SAMPLE_DISTANCE;
    pImpl->mDetectionWaveformPolicy = DEFAULT_CORRELOGRAM_POLICY;
    pImpl->mWantDetectedWaveform = false;
    pImpl->mWantAmplitude = true;
}

/// Get/set the threshold
void DetectorParameters::setDetectionThreshold(const double tol)
{
    if (tol <= 0 || tol > 1)
    {
        throw std::invalid_argument("Threshold = "
                                  + std::to_string(tol)
                                  + " must be in range (0,1]\n");
    }
    pImpl->mThreshold = tol;
}

double DetectorParameters::getDetectionThreshold() const noexcept
{
    return pImpl->mThreshold;
}

/// Get/set the distance between detections
void DetectorParameters::setMinimumDetectionSpacing(
    const int nSamples) noexcept
{
    pImpl->mMinSampleDistance = nSamples;
}

int DetectorParameters::getMinimumDetectionSpacing() const noexcept
{
    return pImpl->mMinSampleDistance;
}

/// Get/set policy on detections on correlations
void DetectorParameters::setWaveformPolicy(
    const DetectionWaveformPolicy policy) noexcept
{
    pImpl->mDetectionWaveformPolicy = policy;
}


MFLib::DetectionWaveformPolicy
DetectorParameters::getWaveformPolicy() const noexcept
{
    return pImpl->mDetectionWaveformPolicy;
}


/// Get/set the maximum policy
void DetectorParameters::setMaximaPolicy(
    const MFLib::MaximumMatchedFilterPolicy policy) noexcept
{
    pImpl->mMaxPolicy = policy;
}

MFLib::MaximumMatchedFilterPolicy
DetectorParameters::getMaximaPolicy() const noexcept 
{
    return pImpl->mMaxPolicy;
}

/// Get/set detected waveform
void DetectorParameters::enableSaveDetectedWaveform() noexcept
{
    pImpl->mWantDetectedWaveform = true;
}

void DetectorParameters::disableSaveDetectedWaveform() noexcept
{
    pImpl->mWantDetectedWaveform = false;
}

bool DetectorParameters::wantDetectedWaveform() const noexcept
{
    return pImpl->mWantDetectedWaveform;
}

/// Get/set the amplitudes
void DetectorParameters::enableSaveAmplitudeScalingFactor() noexcept
{
    pImpl->mWantAmplitude = true;
}

void DetectorParameters::disableSaveAmplitudeScalingFactor() noexcept
{
    pImpl->mWantAmplitude = false;
}

bool DetectorParameters::wantAmplitudeScalingFactor() const noexcept
{
    return pImpl->mWantAmplitude;
}
