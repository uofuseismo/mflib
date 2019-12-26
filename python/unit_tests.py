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
    mfOptions = pymflib.MatchedFilterParameters()

    # Set some templates
    signalSize = 1024
    nSamples1 = 100
    nSamples2 = 95
    t1 = np.arange(1, nSamples1, 1)
    t2 = np.arange(1, nSamples2, 1)

    wf1.signal = t1
    wf2.signal = t2

    mfOptions.addTemplate(wf1)
    mfOptions.addTemplate(wf2)
    