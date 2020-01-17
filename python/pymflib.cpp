#include "mflib/version.hpp"
#include "pyWaveformTemplate.hpp"
#include "pyMultiChannelMatchedFilterParameters.hpp"
#include "pyMultiChannelMatchedFilter.hpp"
#include "pySingleChannelMatchedFilterParameters.hpp"
#include "pySingleChannelMatchedFilter.hpp"
#include <pybind11/pybind11.h>

PYBIND11_MODULE(pymflib, m)
{
    // Version number
    m.attr("__version__") = MFLIB_VERSION; 
   
    PBMFLib::initializeWaveformTemplate(m);
    PBMFLib::SingleChannel::initializeMatchedFilterParameters(m);
    PBMFLib::SingleChannel::initializeMatchedFilter(m);
    PBMFLib::initializeMultiChannelMatchedFilterParameters(m);
    PBMFLib::initializeMultiChannelMatchedFilter(m);
}
