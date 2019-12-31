#!/usr/bin/env python3
import pytest
import pymflib
import numpy as np

def test_waveform_template():
    """
    Tests the waveform template.  This is the basic unit of currency in this library.  
    """
    # Some values to set then recover
    template_signal = np.arange(1, 400, 1) # 4 seconds of signal
    sampling_rate = 100
    shift_and_stack_weight = 0.8
    travel_time = 9 # Travel time for this phase is 9 seconds
    onset_time = .5 # Onset is .5 seconds into the signal
    # Figure up the class
    wt = pymflib.WaveformTemplate()

    # Test the signal
    wt.signal = template_signal
    assert np.max(np.abs(wt.signal - template_signal)) == 0, "Signal failed"

    # Test the sampling rate
    wt.sampling_rate = sampling_rate 
    assert wt.sampling_rate == sampling_rate, "Sampling rate failed"

    # Test the shift and stack weight
    wt.shift_and_stack_weight = shift_and_stack_weight
    assert wt.shift_and_stack_weight == shift_and_stack_weight, "Shift and stack weight failed"

    # Test the travel time
    wt.phase_travel_time = travel_time
    assert wt.phase_travel_time == travel_time, "Travel time failed"

    # Test the onset time
    wt.onset_time_in_signal = onset_time
    assert wt.onset_time_in_signal == onset_time, "Onset time failed"

    # Test the copy c'tor
    wt_copy = wt
    assert np.max(np.abs(wt_copy.signal - template_signal)) == 0, "Copy signal failed"
    assert wt_copy.sampling_rate == sampling_rate, "Copy sampling rate failed"
    assert wt_copy.shift_and_stack_weight == shift_and_stack_weight, "Copy shift and stack weight failed"
    assert wt_copy.phase_travel_time == travel_time, "Travel time failed"
    assert wt_copy.onset_time_in_signal == onset_time, "Onset time failed"

    # Clear the module
    wt.clear()
    assert wt.shift_and_stack_weight == 1, "Clear failed" # Should reset to default weight

##########################################################################################################

def test_matched_filter_parameters():
    """
    Performs some basic matched filtering testings.
    """
    wf1 = pymflib.WaveformTemplate()
    wf2 = pymflib.WaveformTemplate()
    mf_parameters = pymflib.MultiChannelMatchedFilterParameters()

    signal_size = 1024  # Signal size
    sampling_rate = 100 # Sampling rate of templates and signal
    n_samples1 = 100    # Length of first template
    n_samples2 = 95     # Length of second template
    t1 = np.arange(1, n_samples1, 1)
    t2 = np.arange(1, n_samples2, 1)

    wf1.signal = t1
    wf2.signal = t2
    wf1.sampling_rate = sampling_rate
    wf2.sampling_rate = sampling_rate
    # Set some other things for a subsequent test
    travel_time = 4
    onset_time = 0.05
    sas_weight = 0.8
    wf1.shift_and_stack_weight = sas_weight
    wf1.onset_time_in_signal = onset_time
    wf1.shift_and_stack_weight = sas_weight
    wf1.phase_travel_time = travel_time

    wf2.shift_and_stack_weight = sas_weight
    wf2.onset_time_in_signal = onset_time
    wf2.shift_and_stack_weight = sas_weight
    wf2.phase_travel_time = travel_time

    mf_parameters.add_template(wf1)
    mf_parameters.add_template(wf2)
    mf_parameters.signal_size = signal_size
    assert mf_parameters.number_of_templates == 2, "n_templates failed"
    # Should default to this
    assert mf_parameters.fft_length == 512, "fft length wrong"
    # Can try changing it
    mf_parameters.fft_length = 550
    assert mf_parameters.fft_length == 550, "fft length change failed"
    # Try to recover one of my templates
    for i in range(mf_parameters.number_of_templates):
        wt_back = mf_parameters.get_template(i)
        assert wt_back.sampling_rate == sampling_rate, "sampling rate copy failed"
        assert wt_back.shift_and_stack_weight == sas_weight, "sas weight copy failed"
        assert wt_back.onset_time_in_signal == onset_time, "onset time copy failed"
        assert wt_back.phase_travel_time == travel_time, "ttime copy failed"
        t_back = wt_back.signal
        if (i == 0):
            t = t1
        else:
            t = t2
        assert np.max(np.abs(t_back - t)) == 0, 'failed to recover template'

    # Dump the templates
    mf_parameters.clear_templates()

    assert mf_parameters.number_of_templates == 0, "clear templates failed"

#############################################################################################################

def dumb_xc(template, signal):
    template = template - np.average(template)
    template_norm = np.linalg.norm(template, ord=2)
    nb = len(template)
    nx = len(signal)
    xc = np.zeros(nx-nb+1)
    for i in range(len(xc)):
        i1 = i
        i2 = i + nb
        y = signal[i1:i2] - np.average(signal[i1:i2])
        y_norm = np.linalg.norm(y, ord=2)
        xc[i] = np.dot(template, y)/(template_norm*y_norm)
    return xc

def test_matched_filtering():
    wf1 = pymflib.WaveformTemplate()
    wf2 = pymflib.WaveformTemplate()
    # Define signals and templates
    signal_size = 3000
    Tmax = 8
    f0 = 4/Tmax
    dt = Tmax/(signal_size - 1)
    sampling_rate = 1/dt
    times = np.linspace(0, Tmax, signal_size, endpoint=False)
    signal1 = np.exp(-0.1*times)*np.sin(2*np.pi*f0*times)
    signal2 = np.exp(-0.15*times)*np.cos(2*np.pi*f0*times) 

    # Define the templates
    Template_window = 2
    nb = int(Template_window/dt + 0.5)
    template_times = np.linspace(0, Template_window, nb, endpoint=False)
    t1 = np.sin(2*np.pi*f0*template_times)
    t2 = np.cos(2*np.pi*f0*template_times)

    wf1 = pymflib.WaveformTemplate()
    wf2 = pymflib.WaveformTemplate()
    wf1.signal = t1
    wf2.signal = t2
    wf1.sampling_rate = sampling_rate
    wf2.sampling_rate = sampling_rate

    # Attached the templates to the parameters
    mf_parameters = pymflib.MultiChannelMatchedFilterParameters()
    mf_parameters.add_template(wf1)
    mf_parameters.add_template(wf2)
    mf_parameters.signal_size = signal_size

    # Initialize the multi-channel matched filter
    mf = pymflib.MultiChannelMatchedFilter()
    mf.initialize(mf_parameters)
    # Attach the signals for this run
    mf.set_signal(0, signal1)
    mf.set_signal(1, signal2)
    # Run it
    mf.apply()
    assert mf.have_matched_filtered_signals, "mf switch not updated" 
    mf1 = mf.get_matched_filtered_signal(0)
    mf2 = mf.get_matched_filtered_signal(1)
    xc1 = dumb_xc(t1, signal1) # Do it the dumb way
    xc2 = dumb_xc(t2, signal2) # Do it the dumb way
    assert np.max(np.abs(mf1 - xc1)) < 1.e-14, 'mfilter 1 failed'
    assert np.max(np.abs(mf2 - xc2)) < 1.e-14, 'mfilter 2 failed'
   

if __name__ == "__main__":
    test_matched_filter_parameters()
    test_matched_filtering()
