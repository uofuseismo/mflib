#ifndef PYMFLIB_MODULES_HPP
#define PYMFLIB_MODULES_HPP
#include <pybind11/pybind11.h>

/// Initializes the waveform template python interface
void initializeWaveformTemplate(pybind11::module &m);
/// Initializes the matched filter python interface
//void initializeMatchedFilter(pybind11::module &m);

#endif
