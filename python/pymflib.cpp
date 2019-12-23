#include "modules.hpp"
#include <pybind11/pybind11.h>

PYBIND11_MODULE(libpymflib, modules)
{
    pybind11::module m = modules.def_submodule("MF");//WaveformTemplate");
    initializeWaveformTemplate(m);
}
