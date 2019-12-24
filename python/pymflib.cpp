#include "modules.hpp"
#include "waveformTemplate.hpp"
#include "matchedFilterOptions.hpp"
#include <pybind11/pybind11.h>

PYBIND11_MODULE(pymflib, m)
{
    //pybind11::module m; // = modules.def_submodule("Classes");
    //m.doc() = "Bindings for performing matched filtering.";

    pybind11::class_<PBMFLib::WaveformTemplate> wt(m, "WaveformTemplate");
    /// The default constructor:
    wt.def(pybind11::init<> ());
    wt.doc() = "Defines a waveform template to correlate with a signal.";
    /// Waveform template signal:
    wt.def_property("signal",
                    &PBMFLib::WaveformTemplate::getSignal,
                    &PBMFLib::WaveformTemplate::setSignal,
                    "The waveform template signal.  Nominally, this will be a NumPy array representing a time series.  This is required for basic cross-correlation.");
    /// Sampling rate (Hz):
    wt.def_property("sampling_rate",
                    &PBMFLib::WaveformTemplate::getSamplingRate,
                    &PBMFLib::WaveformTemplate::setSamplingRate,
                    "The sampling rate, in Hz, at which the template was sampled.  This must be positive.  This is required for basic cross-correlation.  Note - this will invalidate the phase_travel_time and onset_time_in_signal.");
    /// Shift and stack weight:
    wt.def_property("shift_and_stack_weight",
                    &PBMFLib::WaveformTemplate::getShiftAndStackWeight,
                    &PBMFLib::WaveformTemplate::setShiftAndStackWeight,
                    "Weight used when shifting and stacking correlations.  This must be in the range of [0,1].  By default this is unity.  This parameter need only be set if you are shifting and stacking waveforms.");
    /// The arrival onset (relative to the trace start)
    wt.def_property("onset_time_in_signal",
                    &PBMFLib::WaveformTemplate::getOnsetTime,
                    &PBMFLib::WaveformTemplate::setOnsetTime,
                    "This defines the arrival onset time in seconds relative to the template signal start.  For example, if this is 2, then 2 seconds into the template signal is the expected arrival onset.  This is only necessary if you are shifting and stacking waveforms.  Note - the signal and sampling rate should be set prior to setting this property.");
    /// The travel time 
    wt.def_property("phase_travel_time",
                    &PBMFLib::WaveformTemplate::getTravelTime,
                    &PBMFLib::WaveformTemplate::setTravelTime,
                    "This defines the phase's travel time in seconds from the origin time.  For example, if this is 9 then it took 9 seconds for the phase to travel from the source to the receiver.  This is only necessary if you are shifting and stacking waveforms.  Note - the signal and sampling rate should be set prior to setting this property.");
    /// Clears the class memory and resets
    wt.def("clear",
            &PBMFLib::WaveformTemplate::clear,
           "Clears the class's memory and resets the class.");
    //--------------------------------------------------------------------------------------------//
    //                                  Matched Filter Options                                    //
    //--------------------------------------------------------------------------------------------// 
    pybind11::class_<PBMFLib::MatchedFilterOptions> mfOptions(m, "MatchedFilterOptions");
    /// The default constructor:
    mfOptions.def(pybind11::init<> ());
    mfOptions.doc() = "Defines the options class for matched filtering.";
    /// Adds a template to the class
    mfOptions.def("add_template",
                  &PBMFLib::MatchedFilterOptions::addTemplate,
                  "Adds a waveform template.  The waveform template must at minimum have a signal.");
    /// Clear the templates
    mfOptions.def("clear_templates",
                  &PBMFLib::MatchedFilterOptions::clearTemplates,
                  "Clears the existing templates from class's memory.");
    /// Clears the class
    mfOptions.def("clear",
                  &PBMFLib::MatchedFilterOptions::clear,
                  "Clears the class's memory and resets the class.");
}
