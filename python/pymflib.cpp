#include "mflib/version.hpp"
#include "pyWaveformTemplate.hpp"
#include "pyMultiChannelMatchedFilterParameters.hpp"
#include "pyMultiChannelMatchedFilter.hpp"
#include "singleChannel/pyMatchedFilterParameters.hpp"
#include "singleChannel/pyMatchedFilter.hpp"
#include "singleChannel/pyRelativeMagnitude.hpp"
#include "singleChannel/pyDetection.hpp"
#include "singleChannel/pyDetectorParameters.hpp"
#include "singleChannel/pyDetector.hpp"
#include <pybind11/pybind11.h>

PYBIND11_MODULE(pymflib, m)
{
    // Version number
    m.attr("__version__") = MFLIB_VERSION; 
   
    PBMFLib::initializeWaveformTemplate(m);

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

    pybind11::module scModule = m.def_submodule("SingleChannel");

    PBMFLib::SingleChannel::initializeMatchedFilterParameters(scModule);
    PBMFLib::SingleChannel::initializeMatchedFilter(scModule);
    PBMFLib::SingleChannel::initializeRelativeMagnitude(scModule);
    PBMFLib::SingleChannel::initializeDetection(scModule);
    PBMFLib::SingleChannel::initializeDetectorParameters(scModule);
    PBMFLib::SingleChannel::initializeDetector(scModule);

    pybind11::module mcModule = m.def_submodule("MultiChannel");
    PBMFLib::initializeMultiChannelMatchedFilterParameters(mcModule);
    PBMFLib::initializeMultiChannelMatchedFilter(mcModule);
}
