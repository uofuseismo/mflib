#include "mflib/version.hpp"
#include "pyWaveformTemplate.hpp"
#include "pyMultiChannelMatchedFilterParameters.hpp"
#include "pyMultiChannelMatchedFilter.hpp"
#include "singleChannel/pyMatchedFilterParameters.hpp"
#include "singleChannel/pyMatchedFilter.hpp"
#include "singleChannel/pyRelativeMagnitude.hpp"
#include "singleChannel/pyDetection.hpp"
#include <pybind11/pybind11.h>

PYBIND11_MODULE(pymflib, m)
{
    // Version number
    m.attr("__version__") = MFLIB_VERSION; 
   
    PBMFLib::initializeWaveformTemplate(m);
    PBMFLib::SingleChannel::initializeMatchedFilterParameters(m);
    PBMFLib::SingleChannel::initializeMatchedFilter(m);
    PBMFLib::SingleChannel::initializeRelativeMagnitude(m);
    PBMFLib::SingleChannel::initializeDetection(m);
    PBMFLib::initializeMultiChannelMatchedFilterParameters(m);
    PBMFLib::initializeMultiChannelMatchedFilter(m);
}
