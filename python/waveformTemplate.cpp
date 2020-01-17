#include <cstdio>
#include <cstdlib>
#include <vector>
#include "mflib/waveformTemplate.hpp"
#include "pyWaveformTemplate.hpp"

using namespace PBMFLib;

/// Constructor
WaveformTemplate::WaveformTemplate() :
    mWaveformTemplate(std::make_unique<MFLib::WaveformTemplate> ())
{
}

WaveformTemplate::WaveformTemplate(const WaveformTemplate &wt)
{
    *this = wt;
}

WaveformTemplate::WaveformTemplate(const MFLib::WaveformTemplate &wt)
{
    *this = wt;
}

WaveformTemplate& WaveformTemplate::operator=(const WaveformTemplate &wt)
{
    if (&wt == this){return *this;}
    mWaveformTemplate
        = std::make_unique<MFLib::WaveformTemplate> (*wt.mWaveformTemplate);
    return *this; 
}

WaveformTemplate& WaveformTemplate::operator=(const MFLib::WaveformTemplate &wt)
{
    mWaveformTemplate = std::make_unique<MFLib::WaveformTemplate> (wt);
    return *this;
}

MFLib::WaveformTemplate WaveformTemplate::getNativeClass() const
{
    MFLib::WaveformTemplate t(*mWaveformTemplate);
    return t;
}

/// Destructor
WaveformTemplate::~WaveformTemplate() = default;

/// Clears the module
void WaveformTemplate::clear() noexcept
{
    mWaveformTemplate->clear();
}

/// Setters/getters for signal
void WaveformTemplate::setSignal(
    const pybind11::array_t<double, pybind11::array::c_style |
                                    pybind11::array::forcecast> &x)
{
    // Use pointers
    pybind11::buffer_info xbuf = x.request();
    auto len = static_cast<int> (xbuf.size);
    const double *xptr = (double *) (xbuf.ptr);
    if (xptr == nullptr)
    {
        throw std::invalid_argument("x is null");
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

/// Setters/getters for onset time
void WaveformTemplate::setOnsetTime(const double onsetTime)
{
    if (onsetTime < 0)
    {
        throw std::invalid_argument("onset_time = " 
                                  + std::to_string(onsetTime)
                                  + " must be positive\n");
    }
    mWaveformTemplate->setOnsetTime(onsetTime);
}

double WaveformTemplate::getOnsetTime() const
{
    return mWaveformTemplate->getOnsetTime();
}

/// Setters/getters for travel time
void WaveformTemplate::setTravelTime(const double travelTime)
{
    if (travelTime < 0)
    {
        throw std::invalid_argument("onset_time = " 
                                  + std::to_string(travelTime)
                                  + " must be positive\n");
    }
    mWaveformTemplate->setTravelTime(travelTime);
}


double WaveformTemplate::getTravelTime() const
{
    return mWaveformTemplate->getTravelTime();
}

void PBMFLib::initializeWaveformTemplate(pybind11::module &m)
{
    pybind11::class_<PBMFLib::WaveformTemplate> wt(m, "WaveformTemplate");
    /// The default constructor:
    wt.def(pybind11::init<> ());
    wt.doc() = "Defines a waveform template to correlate with a signal.";
    /// Waveform template signal:
    wt.def("get_signal",
           &PBMFLib::WaveformTemplate::getSignal,
           "Gets the waveform template signal.");
    wt.def("set_signal",
           &PBMFLib::WaveformTemplate::setSignal,
          "The waveform template signal.  Nominally, this will be a NumPy array representing a time series.  This is required for basic cross-correlation.");
    /// Sampling rate (Hz):
    wt.def("get_sampling_rate",
           &PBMFLib::WaveformTemplate::getSamplingRate,
           "Gets the sampling rate in Hz.");
    wt.def("set_sampling_rate",
           &PBMFLib::WaveformTemplate::setSamplingRate,
           "The sampling rate, in Hz, at which the template was sampled.  This must be positive.  This is required for basic cross-correlation.  Note - this will invalidate the phase_travel_time and onset_time_in_signal.");
    /// Shift and stack weight:
    wt.def("get_shift_and_stack_weight",
           &PBMFLib::WaveformTemplate::getShiftAndStackWeight,
           "The shfit and stack weight.");
    wt.def("set_shift_and_stack_weight",
           &PBMFLib::WaveformTemplate::setShiftAndStackWeight,
           "Weight used when shifting and stacking correlations.  This must be in the range of [0,1].  By default this is unity.  This parameter need only be set if you are shifting and stacking waveforms.");
    /// The arrival onset (relative to the trace start)
    wt.def("get_onset_time_in_signal",
           &PBMFLib::WaveformTemplate::getOnsetTime,
           "The phase onset time of the relative to the signal start in seconds.");
    wt.def("set_onset_time_in_signal",
           &PBMFLib::WaveformTemplate::setOnsetTime,
           "This defines the arrival onset time in seconds relative to the template signal start.  For example, if this is 2, then 2 seconds into the template signal is the expected arrival onset.  This is only necessary if you are shifting and stacking waveforms.  Note - the signal and sampling rate should be set prior to setting this property.");
    /// The travel time 
    wt.def("get_phase_travel_time",
           &PBMFLib::WaveformTemplate::getTravelTime,
           "The phase travel time from the source to the receiver in seconds.");
    wt.def("set_phase_travel_time",
           &PBMFLib::WaveformTemplate::setTravelTime,
           "This defines the phase's travel time in seconds from the origin time.  For example, if this is 9 then it took 9 seconds for the phase to travel from the source to the receiver.  This is only necessary if you are shifting and stacking waveforms.  Note - the signal and sampling rate should be set prior to setting this property.");
    /// Clears the class memory and resets
    wt.def("clear",
           &PBMFLib::WaveformTemplate::clear,
           "Clears the class's memory and resets the class.");
}
