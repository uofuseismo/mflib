#include "mflib/version.hpp"
#include "pyWaveformTemplate.hpp"
#include "pyNetworkStationPhase.hpp"
#include "pyMultiChannelMatchedFilterParameters.hpp"
#include "pyMultiChannelMatchedFilter.hpp"
#include "singleChannel/pyMatchedFilterParameters.hpp"
#include "singleChannel/pyMatchedFilter.hpp"
#include "singleChannel/pyAssociatorParameters.hpp"
#include "singleChannel/pyRelativeMagnitude.hpp"
#include "singleChannel/pyDetection.hpp"
#include "singleChannel/pyDetectorParameters.hpp"
#include "singleChannel/pyDetector.hpp"
#include "singleChannel/pyAssociator.hpp"
#include <pybind11/pybind11.h>

PYBIND11_MODULE(pymflib, m)
{
    // Version number
    m.attr("__version__") = MFLIB_VERSION; 
   
    PBMFLib::initializeWaveformTemplate(m);
    PBMFLib::initializeNetworkStationPhase(m);

    pybind11::enum_<MFLib::RelativeMagnitudeType> (m, "RelativeMagnitudeType")
        .value("gibbons_ringdal_2006", MFLib::RelativeMagnitudeType::GIBBONS_RINGDAL_2006,
               "Computes the relative magnitude using the methodology of Gibbons and Ringdal, 2006.  This method assumes the the template waveform is noise-free.")
        .value("schaff_richards_2014", MFLib::RelativeMagnitudeType::SCHAFF_RICHARDS_2014,
               "Computes the relative magnitude using the methodology of Schaff and Richards, 2014.  This method accomodates noise in both the detected and template waveform.");
    pybind11::enum_<MFLib::MaximumMatchedFilterPolicy> (m, "MaxmimumMatchedFilterPolicy")
        .value("maximum", MFLib::MaximumMatchedFilterPolicy::MAXIMUM,
               "Selects the maxima of the matched filtered signals.")
        .value("absolute_maximum", MFLib::MaximumMatchedFilterPolicy::ABSOLUTE_MAXIMUM,
               "Selects the maxima of the absolute value of the matched filtered signals.");

    pybind11::enum_<MFLib::MatchedFilteredSignalDetectorPolicy> (m, "MatchedFilteredSignalDetectorPolicy")
        .value("reduce", MFLib::MatchedFilteredSignalDetectorPolicy::REDUCED,
               "Combines all matched filtered signals by taking the maximum or absolute maximum then computing detections.")
        .value("single", MFLib::MatchedFilteredSignalDetectorPolicy::SINGLE,
               "Makes detections on individual matched filtered signals.  This may result in multiple detections from a different template occurring simultaneously in a small time window.  In this case, the associator must remove these duplicates."); 

    pybind11::enum_<MFLib::Polarity> (m, "Polarity")
        .value("unknown", MFLib::Polarity::UNKNOWN,
               "The polarity is unknown.")
        .value("compression", MFLib::Polarity::COMPRESSION,
               "The polarity corresponds to the compressional (shaded) quadrant of the beachball and is probably upwards on the seismogram.") 
        .value("dilatation", MFLib::Polarity::DILATATION,
               "The polarity corresopnds to the dilitational (unshaded) quadrant of the beachball and is probably downwards on the seismogram.");

    pybind11::enum_<MFLib::DetectionTimeInterpolationType> (m, "DetectionTimeInterpolationType")
        .value("quadratic", MFLib::DetectionTimeInterpolationType::QUADRATIC,
               "A quadratic polynomial will be fit to the correlation values.  The polynomial will be optimized for an interpolated detection time.  This is fast and relatively accurate.")
        .value("lanczos", MFLib::DetectionTimeInterpolationType::LANCZOS,
               "A truncated sinc (Lanczos) interpolation will be applied to the correlation values.  This is slow but can be more accurate than the quadrtic method.")
        .value("none", MFLib::DetectionTimeInterpolationType::NONE,
               "No detection time interpolation will be performed and the interpolated detection time will equal the observed detection time.");

    pybind11::module scModule = m.def_submodule("SingleChannel");

    PBMFLib::SingleChannel::initializeMatchedFilterParameters(scModule);
    PBMFLib::SingleChannel::initializeMatchedFilter(scModule);
    PBMFLib::SingleChannel::initializeRelativeMagnitude(scModule);
    PBMFLib::SingleChannel::initializeDetection(scModule);
    PBMFLib::SingleChannel::initializeDetectorParameters(scModule);
    PBMFLib::SingleChannel::initializeDetector(scModule);
    PBMFLib::SingleChannel::initializeAssociatorParameters(scModule);
    PBMFLib::SingleChannel::initializeAssociator(scModule);

    pybind11::module mcModule = m.def_submodule("MultiChannel");
    PBMFLib::initializeMultiChannelMatchedFilterParameters(mcModule);
    PBMFLib::initializeMultiChannelMatchedFilter(mcModule);
}
