#include <cstdio>
#include <cstdlib>
#include "mflib/waveformTemplate.hpp"
#include "modules.hpp"
#include "waveformTemplate.hpp"

using namespace PBWaveformTemplate;

/// Constructor
WaveformTemplate::WaveformTemplate() :
    mWaveformTemplate(std::make_unique<MFLib::WaveformTemplate> ())
{
}

/// Destructor
WaveformTemplate::~WaveformTemplate() = default;

/// Setters/getters for signal
void WaveformTemplate::setSignal(
    pybind11::array_t<double, pybind11::array::c_style |
                              pybind11::array::forcecast> &x)
{
    // Use pointers
    pybind11::buffer_info xbuf = x.request();
    auto len = static_cast<int> (xbuf.size);
    const double *xptr = (double *) (xbuf.ptr);
    if (xptr == nullptr)
    {
        throw std::runtime_error("x is null");
    }
    mWaveformTemplate->setSignal(len, xptr);
}

pybind11::array_t<double> WaveformTemplate::getSignal() const
{
    auto npts = mWaveformTemplate->getSignalLength();
    auto y = pybind11::array_t<double, pybind11::array::c_style> (npts);
    pybind11::buffer_info ybuf = y.request();
    auto yptr = static_cast<double *> (ybuf.ptr);
    mWaveformTemplate->getSignal(npts, &yptr);
    return y;
}

/// Setters/getters for sampling rate
void WaveformTemplate::setSamplingRate(const double samplingRate)
{
    if (samplingRate <= 0)
    {
        throw std::invalid_argument("Sampling rate = "
                                  + std::to_string(samplingRate)
                                  + " must be positive\n");
    }
    mWaveformTemplate->setSamplingRate(samplingRate);
}

double WaveformTemplate::getSamplingRate() const
{
    return mWaveformTemplate->getSamplingRate();
}

/// Setters/getters for shift and stack weight
void WaveformTemplate::setShiftAndStackWeight(const double weight)
{
    if (weight < 0 || weight > 1)
    {
        throw std::invalid_argument("Weight = " + std::to_string(weight)
                                 + " must be in range [0,1]\n");
    }
    mWaveformTemplate->setShiftAndStackWeight(weight);
}

double WaveformTemplate::getShiftAndStackWeight() const
{
    return mWaveformTemplate->getShiftAndStackWeight();
}

void initializeWaveformTemplate(pybind11::module &m)
{
    m.doc() = "Defines the waveform template";

    pybind11::class_<PBWaveformTemplate::WaveformTemplate> wt(m, "WaveformTemplate");

    /// The waveform
    wt.def(pybind11::init<> ());
    wt.doc() = "Defines a waveform template to correlate with a signal.";
    wt.def("set_signal", &PBWaveformTemplate::WaveformTemplate::getSignal,
           "Sets the signal defining a waveform template.");
    wt.def("get_signal", &PBWaveformTemplate::WaveformTemplate::setSignal,
           "Gets the signal defining a waveform template.");
    /// Sampling rate (Hz)
    wt.def_property("sampling_rate",
                    &PBWaveformTemplate::WaveformTemplate::getSamplingRate,
                    &PBWaveformTemplate::WaveformTemplate::setSamplingRate,
                    "The sampling rate, in Hz, at which the template was sampled.  This must be positive.");
    /// Shift and stack weight
    wt.def_property("shift_and_stack_weight",
                    &PBWaveformTemplate::WaveformTemplate::getShiftAndStackWeight,
                    &PBWaveformTemplate::WaveformTemplate::setShiftAndStackWeight,
                    "Weight used when shifting and stacking correlations.  This must be in the range of [0,1].");



}
