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

/// Get/set get detected waveform
void DetectorParameters::toggleGetDetectedWaveform(const bool lget) noexcept
{
    mParameters->toggleGetDetectedWaveform(lget);
}

bool DetectorParameters::getDetectedWaveform() const noexcept
{
    return mParameters->getDetectedWaveform();
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
              
    parms.def("toggle_get_detected_waveform",
              &PBMFLib::SingleChannel::DetectorParameters::toggleGetDetectedWaveform,
              "By setting this to true, the detections will contain the snippet of detected waveform.  While this can be useful for visualation this also comes with memory overhead.");
    parms.def("get_detected_waveform",
              &PBMFLib::SingleChannel::DetectorParameters::getDetectedWaveform,
              "Determines if the detections will contain the corresponding snippet of detected waveform.");

    pybind11::enum_<MFLib::MaximumMatchedFilterPolicy> (parms, "MaxmimumMatchedFilterPolicy")
        .value("maximum", MFLib::MaximumMatchedFilterPolicy::MAXIMUM,
               "Selects the maxima of the matched filtered signals.")
        .value("absolute_maximum", MFLib::MaximumMatchedFilterPolicy::ABSOLUTE_MAXIMUM,
               "Selects the maxima of the absolute value of the matched filtered signals."); 
}
