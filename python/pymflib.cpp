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
