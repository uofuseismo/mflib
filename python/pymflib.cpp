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

    PBMFLib::SingleChannel::initializeMatchedFilter(m);
    //--------------------------------------------------------------------------------------------//
    //                                     Waveform Template                                      //
    //--------------------------------------------------------------------------------------------//
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
    //                                  Matched Filter Parameters                                 //
    //--------------------------------------------------------------------------------------------// 
    pybind11::class_<PBMFLib::MultiChannelMatchedFilterParameters>
        mfParameters(m, "MultiChannelMatchedFilterParameters");
    /// The default constructor:
    mfParameters.def(pybind11::init<> ());
    mfParameters.doc() = "Defines the options class for matched filtering.";
    /// Adds a template to the class
    mfParameters.def("add_template",
                     &PBMFLib::MultiChannelMatchedFilterParameters::addTemplate,
                     "Adds a waveform template.  The waveform template must at minimum have a signal.");
    mfParameters.def("get_template",
                     &PBMFLib::MultiChannelMatchedFilterParameters::getTemplate,
                     "Gets the it'th waveform template.");
    mfParameters.def_property("signal_size",
                              &PBMFLib::MultiChannelMatchedFilterParameters::getSignalSize,
                              &PBMFLib::MultiChannelMatchedFilterParameters::setSignalSize,
                              "Defines the size of the signals to which the templates will be applied.  This must be positive and should be called after setting all the templates.");
    mfParameters.def_property_readonly("number_of_templates",
                                       &PBMFLib::MultiChannelMatchedFilterParameters::getNumberOfTemplates,
                                       "This is the number of templates that have been set in the class.");
    /// Tuning parameters
    mfParameters.def_property("fft_length",
                              &PBMFLib::MultiChannelMatchedFilterParameters::getFFTLength,
                              &PBMFLib::MultiChannelMatchedFilterParameters::setFFTLength,
                              "This is a tuning parameter.  From the signal size and template lengths a nominal FFT length will be computed.  While this length will be a power of 2 the underlying Cooley-Tukey Discrete Fourier Transform may be more efficient for smaller FFT lenghts.  This should be set after both the templates and signal size are set."); 
    mfParameters.def("block_length",
                     &PBMFLib::MultiChannelMatchedFilterParameters::getFFTLength,
                     "This is the block length in the overlap-and-add method.  It is computed from the FFT length and the signal size.");
    /// Clear the templates
    mfParameters.def("clear_templates",
                     &PBMFLib::MultiChannelMatchedFilterParameters::clearTemplates,
                     "Clears the existing templates from class's memory.");
    /// Clears the class
    mfParameters.def("clear",
                     &PBMFLib::MultiChannelMatchedFilterParameters::clear,
                     "Clears the class's memory and resets the class.");
    //--------------------------------------------------------------------------------------------//
    //                                   Matched Filtering                                        //
    //--------------------------------------------------------------------------------------------//
    pybind11::class_<PBMFLib::MultiChannelMatchedFilter<double>> dmcmf(m, "MultiChannelMatchedFilter");
    //pybind11::class_<PBMFLib::MultiChannelMatchedFilter<float>> fmcmf(m, "MultiChannelMatchedFilter");
    //pybind11::class_<PBMFLib::MultiChannelMatchedFilter<double>> mcmf(m, "MultiChannelMatchedFilter");
    //pybind11::class_<PBMFLib::MultiChannelMatchedFilter<float>> mcmf(m, "MultiChannelMatchedFilter");
    /// The default constructor:
    dmcmf.def(pybind11::init<> ());
    dmcmf.doc() = "Applies matched filtering to multiple channels where each channel has a specific template.";
    /// Initialize
    dmcmf.def("initialize",
              &PBMFLib::MultiChannelMatchedFilter<double>::initialize,
              "Initializes the multi-channel matched filtering engine.  This is a high-overhead function.");
    /// Sets signals
    dmcmf.def("set_signal",
              &PBMFLib::MultiChannelMatchedFilter<double>::setSignal,
              "Sets the signal corresponding to the it'th template.  This should be called after the class is initialized.  Moreover, the template index must in the range of [0, number_of_templates].");
    dmcmf.def("zero_signal",
              &PBMFLib::MultiChannelMatchedFilter<double>::zeroSignal,
              "Over a time period a station can be down.  This sets the signal corresponding to the it'th template to 0.  The template index must be in the range of [0, number_of_templates].");
    dmcmf.def_property_readonly("number_of_templates", 
                                &PBMFLib::MultiChannelMatchedFilter<double>::getNumberOfTemplates, 
                                "The number of templates in the multi-channel cross-correlation.");
    /// Apply
    dmcmf.def("apply",
              &PBMFLib::MultiChannelMatchedFilter<double>::apply,
              "Applies the templates set during the initialization stage to the signals set by set_signal");
    /// 
    dmcmf.def_property_readonly("have_matched_filtered_signals",
                                &PBMFLib::MultiChannelMatchedFilter<double>::haveMatchedFilteredSignals,
                                "Determines if the matched filtering as been applied and the signals are available for extraction from the class.");
    /// Gets the it'th matched filtered signal
    dmcmf.def("get_matched_filtered_signal",
              &PBMFLib::MultiChannelMatchedFilter<double>::getMatchedFilteredSignal,
              "Gets the it'th matched filtered signal.  The `apply' method must have been called and the template index must be in the range [0,number_of_templates]");
           

}
