#include <cstdio>
#include <cstdlib>
#include "mflib/singleChannel/detectorParameters.hpp"
#include "singleChannel/pyDetectorParameters.hpp"

using namespace PBMFLib::SingleChannel;

/// Constructor
DetectorParameters::DetectorParameters() :
    mParameters(std::make_unique<MFLib::SingleChannel::DetectorParameters> ())
{
}

DetectorParameters::DetectorParameters(const DetectorParameters &parms)
{
    *this = parms;
}

DetectorParameters::DetectorParameters(
    const MFLib::SingleChannel::DetectorParameters &parms)
{
    *this = parms;
}

/// Copy assignment operators
DetectorParameters&
DetectorParameters::operator=(const DetectorParameters &parms)
{
    if (&parms == this){return *this;}
    mParameters = std::make_unique<MFLib::SingleChannel::DetectorParameters>
                  (*parms.mParameters);
    return *this;
}

DetectorParameters&
DetectorParameters::operator=(
    const MFLib::SingleChannel::DetectorParameters &parms)
{
    mParameters = std::make_unique<MFLib::SingleChannel::DetectorParameters>
                  (parms);
    return *this;
}

/// Get a native class
MFLib::SingleChannel::DetectorParameters
DetectorParameters::getNativeClass() const
{
    MFLib::SingleChannel::DetectorParameters parms(*mParameters);
    return parms;
}

/// Destructor
DetectorParameters::~DetectorParameters() = default;

/// Clears the class
void DetectorParameters::clear() noexcept
{
    mParameters->clear();
}

/// Max policy
void DetectorParameters::setMaximaPolicy(
    const MFLib::MaximumMatchedFilterPolicy policy) noexcept
{
    mParameters->setMaximaPolicy(policy);
}

MFLib::MaximumMatchedFilterPolicy
DetectorParameters::getMaximaPolicy() const noexcept
{
    return mParameters->getMaximaPolicy();
}

/// Detector reduction policy
void DetectorParameters::setMatchedFilteredSignalDetectorPolicy(
    const MFLib::MatchedFilteredSignalDetectorPolicy policy) noexcept
{
    mParameters->setMatchedFilteredSignalDetectorPolicy(policy);
}

MFLib::MatchedFilteredSignalDetectorPolicy
DetectorParameters::getMatchedFilteredSignalDetectorPolicy() const noexcept
{
    return mParameters->getMatchedFilteredSignalDetectorPolicy();
}


/// Minimum sample spacing
void DetectorParameters::setMinimumDetectionSpacing(const int nSamples) noexcept
{
    mParameters->setMinimumDetectionSpacing(nSamples);
}

int DetectorParameters::getMinimumDetectionSpacing() const noexcept
{
    return mParameters->getMinimumDetectionSpacing();
}

/// Get/set detection threshold
void DetectorParameters::setDetectionThreshold(const double tol) 
{
    mParameters->setDetectionThreshold(tol);
}

double DetectorParameters::getDetectionThreshold() const noexcept
{
    return mParameters->getDetectionThreshold();
}

/// Get/set detected waveform
void DetectorParameters::enableSaveDetectedWaveform() noexcept
{
    mParameters->enableSaveDetectedWaveform();
}

void DetectorParameters::disableSaveDetectedWaveform() noexcept
{
    mParameters->disableSaveDetectedWaveform();
}

bool DetectorParameters::wantDetectedWaveform() const noexcept
{
    return mParameters->wantDetectedWaveform();
}

/// Get/set amplitude scaling factor
void DetectorParameters::enableSaveAmplitudeScalingFactor() noexcept
{
    mParameters->enableSaveAmplitudeScalingFactor();
}

void DetectorParameters::disableSaveAmplitudeScalingFactor() noexcept
{
    mParameters->disableSaveAmplitudeScalingFactor();
}

bool DetectorParameters::wantAmplitudeScalingFactor() const noexcept
{
    return mParameters->wantAmplitudeScalingFactor();
}

void PBMFLib::SingleChannel::initializeDetectorParameters(pybind11::module &m)
{
    pybind11::class_<PBMFLib::SingleChannel::DetectorParameters> parms(m, "DetectorParameters");
    parms.def(pybind11::init<> ());
    parms.doc() = "Defines the parameters which influence the behavior of the detector and the returned values in a detection.";

    parms.def("clear",
              &PBMFLib::SingleChannel::DetectorParameters::clear,
              "Resets the class and releases any allocated memory.");

    parms.def("get_maxima_policy",
              &PBMFLib::SingleChannel::DetectorParameters::getMaximaPolicy,
              "Gets the policy defining how the maxima are computed in the matched filtered signals.");
    parms.def("set_maxima_policy",
              &PBMFLib::SingleChannel::DetectorParameters::setMaximaPolicy,
              "Sets the policy defining how maxima are computed in the matched filtered signals.  For example, this can select the maxima or absolute maxima of the matched filtered signals.");

    parms.def("get_matched_filtered_signal_detector_policy",
              &PBMFLib::SingleChannel::DetectorParameters::getMatchedFilteredSignalDetectorPolicy,
              "Gets the policy for whether detections are computed on invidual matched filtered signals.");
 
    parms.def("set_matched_filtered_signal_detector_policy",
              &PBMFLib::SingleChannel::DetectorParameters::setMatchedFilteredSignalDetectorPolicy,
              "Sets the policy on whether detections are computed on individual or reduced matched filtered signals.");

    parms.def("get_minimum_detection_spacing",
              &PBMFLib::SingleChannel::DetectorParameters::getMinimumDetectionSpacing,
              "Gets the minimum spacing in samples between detections.");
    parms.def("set_minimum_detection_spacing",
              &PBMFLib::SingleChannel::DetectorParameters::setMinimumDetectionSpacing,
              "Sets the minimum spacing in samples between detections.  By setting this to 0 or a negative value this effectively disables this feature."); 

    parms.def("get_detection_threshold",
              &PBMFLib::SingleChannel::DetectorParameters::getDetectionThreshold,
              "Gets the mimimum value of a candidate detection.");
    parms.def("set_detection_threshold",
              &PBMFLib::SingleChannel::DetectorParameters::setDetectionThreshold,
              "Sets the minimum value of a candidate detection.  This must be in the range (0,1].");
              
    parms.def("enable_save_detected_waveform",
              &PBMFLib::SingleChannel::DetectorParameters::enableSaveDetectedWaveform,
              "This will cause the detected waveforms to save a copy of the detected waveform.  While this can be useful for visualation this also comes with memory overhead.");
    parms.def("disable_save_detected_waveform",
              &PBMFLib::SingleChannel::DetectorParameters::disableSaveDetectedWaveform,
              "This will prevent the detections from containing the detected waveform.");
    parms.def("want_detected_waveform",
              &PBMFLib::SingleChannel::DetectorParameters::wantDetectedWaveform,
              "Determines if the detections will contain the corresponding snippet of detected waveform.");

    parms.def("enable_save_amplitude_scaling_factor",
              &PBMFLib::SingleChannel::DetectorParameters::enableSaveAmplitudeScalingFactor,
              "This will cause the detector to compute the amplitude scaling factors and the corresponding relative magnitudes.");
    parms.def("disable_save_amplitude_scaling_factor",
              &PBMFLib::SingleChannel::DetectorParameters::disableSaveAmplitudeScalingFactor,
              "This will prevent the detections from containing the amplitude scaling factors and relative magnitudes.");
    parms.def("want_amplitude_scaling_factor",
              &PBMFLib::SingleChannel::DetectorParameters::wantAmplitudeScalingFactor,
              "Determines if the detections will contain the the amplitude scaling factors and relative magnitudes.");

}
