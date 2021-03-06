#!/usr/bin/env python3
import pytest
import pymflib
import time
import pickle
import numpy as np

def test_network_station_phase():
    nsp = pymflib.NetworkStationPhase()
    nsp.set_network("UU")
    nsp.set_station("ALT")
    nsp.set_channel("EHN")
    nsp.set_location_code("01")
    nsp.set_phase("P")

    assert nsp.get_network() == "UU", 'network failed'
    assert nsp.get_station() == "ALT",'station failed'
    assert nsp.get_channel() == "EHN", 'channel failed'
    assert nsp.get_location_code() == "01", 'location failed'
    assert nsp.get_phase() == "P",'phase failed'

    nsp_copy = nsp
    assert nsp_copy.get_network() == "UU", 'copy network failed'
    assert nsp_copy.get_station() == "ALT",'copy station failed'
    assert nsp_copy.get_channel() == "EHN", 'copy channel failed'
    assert nsp_copy.get_location_code() == "01", 'copy location failed'
    assert nsp_copy.get_phase() == "P",'copy phase failed'

    data = pickle.dumps(nsp, 2)
    p_nsp = pickle.loads(data)
    assert p_nsp.get_network() == "UU", 'pickle network failed'
    assert p_nsp.get_station() == "ALT",'pickle station failed'
    assert p_nsp.get_channel() == "EHN", 'pickle channel failed'
    assert p_nsp.get_location_code() == "01", 'pickle location failed'
    assert p_nsp.get_phase() == "P",'pickle phase failed'

def test_waveform_template():
    """
    Tests the waveform template.  This is the basic unit of currency in this library.  
    """
    print("Waveform template test...")
    # Some values to set then recover
    template_signal = np.arange(1, 400, 1) # 4 seconds of signal
    sampling_rate = 100
    shift_and_stack_weight = 0.8
    travel_time = 9 # Travel time for this phase is 9 seconds
    onset_time = .5 # Onset is .5 seconds into the signal
    magnitude = 6
    nsp = pymflib.NetworkStationPhase()
    polarity = pymflib.Polarity.compression
    nsp.set_network("UU")
    nsp.set_station("FORU")
    nsp.set_channel("HHZ")
    nsp.set_location_code("01")
    nsp.set_phase("P")
    identifier = [nsp, 4]
    # Figure up the class
    wt = pymflib.WaveformTemplate()

    # Test the signal
    wt.set_signal(template_signal)
    assert np.max(np.abs(wt.get_signal() - template_signal)) == 0, "Signal failed"

    # Test the sampling rate
    wt.set_sampling_rate(sampling_rate)
    assert wt.get_sampling_rate() == sampling_rate, "Sampling rate failed"

    # Test the shift and stack weight
    wt.set_shift_and_stack_weight(shift_and_stack_weight)
    assert wt.get_shift_and_stack_weight() == shift_and_stack_weight, "Shift and stack weight failed"

    # Test the travel time
    wt.set_phase_travel_time(travel_time)
    assert wt.get_phase_travel_time() == travel_time, "Travel time failed"

    # Test the onset time
    wt.set_phase_onset_time_in_signal(onset_time)
    assert wt.get_phase_onset_time_in_signal() == onset_time, "Onset time failed"
 
    # Test the magnitude
    wt.set_magnitude(magnitude)
    assert wt.get_magnitude() == magnitude, "Magnitude failed"

    # Test the polarity
    wt.set_polarity(polarity)
    assert wt.get_polarity() == polarity, "Polarity failed"

    # Test the identifier
    wt.set_identifier(identifier)
    id_back = wt.get_identifier()
    assert id_back[0].get_network() == identifier[0].get_network(), 'network failed'
    assert id_back[0].get_station() == identifier[0].get_station(), 'station failed'
    assert id_back[0].get_channel() == identifier[0].get_channel(), 'channel failed'
    assert id_back[0].get_location_code() == identifier[0].get_location_code(), 'location failed'
    assert id_back[0].get_phase() == identifier[0].get_phase(), 'phase failed'
    assert id_back[1] == identifier[1], "Identifier failed"

    # Test the copy c'tor
    wt_copy = wt
    assert np.max(np.abs(wt_copy.get_signal() - template_signal)) == 0, "Copy signal failed"
    assert wt_copy.get_sampling_rate() == sampling_rate, "Copy sampling rate failed"
    assert wt_copy.get_shift_and_stack_weight() == shift_and_stack_weight, "Copy shift and stack weight failed"
    assert wt_copy.get_phase_travel_time() == travel_time, "Travel time failed"
    assert wt_copy.get_phase_onset_time_in_signal() == onset_time, "Onset time failed"
    assert wt_copy.get_magnitude() == magnitude, "Magnitude failed"
    id_back = wt_copy.get_identifier()
    assert id_back[0].get_network() == identifier[0].get_network(), 'network failed'
    assert id_back[0].get_station() == identifier[0].get_station(), 'station failed'
    assert id_back[0].get_channel() == identifier[0].get_channel(), 'channel failed'
    assert id_back[0].get_location_code() == identifier[0].get_location_code(), 'location failed'
    assert id_back[0].get_phase() == identifier[0].get_phase(), 'phase failed'
    assert id_back[1] == identifier[1], "Identifier failed"
    assert wt_copy.get_polarity() == polarity, "Polarity copy failed"

    # Clear the module
    wt.clear()
    assert wt.get_shift_and_stack_weight() == 1, "Clear failed" # Should reset to default weight

##########################################################################################################

def test_matched_filter_parameters():
    """
    Performs some basic matched filtering testings.
    """
    wf1 = pymflib.WaveformTemplate()
    wf2 = pymflib.WaveformTemplate()
    mf_parameters = pymflib.MultiChannel.MatchedFilterParameters()

    signal_size = 1024  # Signal size
    sampling_rate = 100 # Sampling rate of templates and signal
    n_samples1 = 100    # Length of first template
    n_samples2 = 95     # Length of second template
    t1 = np.arange(1, n_samples1, 1)
    t2 = np.arange(1, n_samples2, 1)

    wf1.set_signal(t1)
    wf2.set_signal(t2)
    wf1.set_sampling_rate(sampling_rate)
    wf2.set_sampling_rate(sampling_rate)
    # Set some other things for a subsequent test
    travel_time = 4
    onset_time = 0.05
    sas_weight = 0.8
    wf1.set_shift_and_stack_weight(sas_weight)
    wf1.set_phase_onset_time_in_signal(onset_time)
    wf1.set_shift_and_stack_weight(sas_weight)
    wf1.set_phase_travel_time(travel_time)

    wf2.set_shift_and_stack_weight(sas_weight)
    wf2.set_phase_onset_time_in_signal(onset_time)
    wf2.set_shift_and_stack_weight(sas_weight)
    wf2.set_phase_travel_time(travel_time)

    mf_parameters.add_template(wf1)
    mf_parameters.add_template(wf2)
    mf_parameters.set_signal_size( signal_size)
    assert mf_parameters.get_number_of_templates() == 2, "n_templates failed"
    # Should default to this
    assert mf_parameters.get_fft_length() == 512, "fft length wrong"
    # Can try changing it
    mf_parameters.set_fft_length(550)
    assert mf_parameters.get_fft_length() == 550, "fft length change failed"
    # Try to recover one of my templates
    for i in range(mf_parameters.get_number_of_templates()):
        wt_back = mf_parameters.get_template(i)
        assert wt_back.get_sampling_rate() == sampling_rate, "sampling rate copy failed"
        assert wt_back.get_shift_and_stack_weight() == sas_weight, "sas weight copy failed"
        assert wt_back.get_phase_onset_time_in_signal() == onset_time, "onset time copy failed"
        assert wt_back.get_phase_travel_time() == travel_time, "ttime copy failed"
        t_back = wt_back.get_signal()
        if (i == 0):
            t = t1
        else:
            t = t2
        assert np.max(np.abs(t_back - t)) == 0, 'failed to recover template'

    # Dump the templates
    mf_parameters.clear_templates()

    assert mf_parameters.get_number_of_templates() == 0, "clear templates failed"

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
    print("Multi channel matched filtering test...")
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
    wf1.set_signal(t1)
    wf2.set_signal(t2)
    wf1.set_sampling_rate(sampling_rate)
    wf2.set_sampling_rate(sampling_rate)

    # Attached the templates to the parameters
    mf_parameters = pymflib.MultiChannel.MatchedFilterParameters()
    mf_parameters.add_template(wf1)
    mf_parameters.add_template(wf2)
    mf_parameters.set_signal_size(signal_size)

    # Initialize the multi-channel matched filter
    mf = pymflib.MultiChannel.MatchedFilter()
    mf.initialize(mf_parameters)
    # Attach the signals for this run
    mf.set_signal(0, signal1)
    mf.set_signal(1, signal2)
    # Run it
    mf.apply()
    assert mf.have_matched_filtered_signals(), "mf switch not updated" 
    mf1 = mf.get_matched_filtered_signal(0)
    mf2 = mf.get_matched_filtered_signal(1)
    xc1 = dumb_xc(t1, signal1) # Do it the dumb way
    xc2 = dumb_xc(t2, signal2) # Do it the dumb way
    assert np.max(np.abs(mf1 - xc1)) < 1.e-14, 'mfilter 1 failed'
    assert np.max(np.abs(mf2 - xc2)) < 1.e-14, 'mfilter 2 failed'
  
##########################################################################################################

def test_single_channel_matched_filtering():
    print("Single channel matched filtering test...")
    sampling_rate = 100.0
    signal_size = int(sampling_rate*60*3) + 1  # Signal size is 3 minutes
    template_size = int(sampling_rate*10) + 1  # Template size is 10 seconds 
    # Generate a random signal
    time_start = time.time()
    signal = np.random.uniform(-3, 4, signal_size) # Mean is (4 + -3)/2 = 0.5 
    wf1 = pymflib.WaveformTemplate()
    wf2 = pymflib.WaveformTemplate()
    wf3 = pymflib.WaveformTemplate()
    # Extract templates from given pct in signal 
    i1 = int(signal_size*0.25)
    i2 = int(signal_size*0.05)
    i3 = int(signal_size*0.8)
    # Extract parts of the input signal and set those as templates
    wf1.set_sampling_rate(sampling_rate)
    wf2.set_sampling_rate(sampling_rate)
    wf3.set_sampling_rate(sampling_rate)
    wf1.set_signal(signal[i1:i1+template_size])
    wf2.set_signal(signal[i2:i2+template_size])
    wf3.set_signal(signal[i3:i3+template_size])
    # Initialize the parameters
    mf_parms = pymflib.SingleChannel.MatchedFilterParameters() 
    mf_parms.set_signal_size(signal_size)
    mf_parms.add_template(wf1)
    mf_parms.add_template(wf2)
    mf_parms.add_template(wf3)
    assert mf_parms.get_number_of_templates() == 3, 'failed to set templates'
    mf = pymflib.SingleChannel.MatchedFilter()
    mf.initialize(mf_parms)
    mf.set_signal(signal)
    mf.apply()
    xc1 = mf.get_matched_filtered_signal(0)
    xc2 = mf.get_matched_filtered_signal(1)
    xc3 = mf.get_matched_filtered_signal(2)
    time_end = time.time()

    time_start_naive = time.time() 
    xc1_ref = dumb_xc(wf1.get_signal(), signal)
    xc2_ref = dumb_xc(wf2.get_signal(), signal)
    xc3_ref = dumb_xc(wf3.get_signal(), signal)
    time_end_naive = time.time()
    assert np.max(np.abs(xc1 - xc1_ref)) < 1.e-14, 'failed to compute xc1'
    assert np.max(np.abs(xc2 - xc2_ref)) < 1.e-14, 'failed to compute xc2'
    assert np.max(np.abs(xc3 - xc3_ref)) < 1.e-14, 'failed to compute xc3'
    mf_parms.clear_templates()
    assert mf_parms.get_number_of_templates() == 0, 'failed to clear templates'
    print("Library execution time: %f (s)"%(time_end - time_start))
    print("Naive execution time on Python: %f (s)"%(time_end_naive - time_start_naive))

def test_single_channel_relative_magnitude():
    print("Single channel relative magnitude test...")
    signal_size = 355
    sampling_rate = 100.0
    rmag = pymflib.SingleChannel.RelativeMagnitude()
    wt = pymflib.WaveformTemplate()
    xt = np.random.uniform( 5, 7, signal_size)
    yt = np.random.uniform(-2, 1, signal_size)
    wt.set_sampling_rate(sampling_rate)
    wt.set_signal(xt)
 
    rmag.initialize(wt) 
    assert rmag.is_initialized(), 'class should be inititalized'
    assert rmag.get_detected_signal_length() == signal_size, 'signal length is wrong'
    rmag.set_detected_signal(yt)
    assert rmag.have_detected_signal(), 'signal not set'

    mag1 = rmag.compute_magnitude_perturbation(pymflib.RelativeMagnitudeType.gibbons_ringdal_2006)
    mag2 = rmag.compute_magnitude_perturbation(pymflib.RelativeMagnitudeType.schaff_richards_2014)
    # Compute reference solutions
    xtd = xt - np.mean(xt)
    ytd = yt - np.mean(yt)
    mag1Ref = np.log10(np.abs(np.dot(xtd, ytd))/(np.dot(xtd, xtd)))
    mag2Ref = 1./2.*np.log10(np.dot(ytd, ytd)/(np.dot(xtd, xtd)))
    assert np.abs(mag1Ref - mag1) < 1.e-14, 'gibbons-ringdal failed'
    assert np.abs(mag2Ref - mag2) < 1.e-14, 'schaff-richards failed'

def test_single_channel_detection():
    print("Single channel detection test...")
    nsp = pymflib.NetworkStationPhase()
    nsp.set_network("UU")
    nsp.set_station("CMU")
    nsp.set_channel("EHZ")
    nsp.set_location_code("01")
    nsp.set_phase("P")
    onset_time = 4
    template_id = [nsp, 204]
    onset_int_time = 9
    detection_time = 2
    xc_val = 0.3
    det_int_time = 3
    scale_factor1 = 9
    scale_factor2 = 10
    travel_time = 13
    phase_int_time = 6
    signal_size = 400
    xt = np.random.uniform( 5, 7, signal_size) 
    detection = pymflib.SingleChannel.Detection()

    assert not detection.have_template_identifier(), 'template id bool not failed'
    detection.set_template_identifier(template_id)
    assert detection.have_template_identifier(), 'template id bool failed'
    template_back = detection.get_template_identifier()
    assert template_back[0].get_network() == template_id[0].get_network(), 'template id net failed'
    assert template_back[0].get_station() == template_id[0].get_station(), 'template id stat failed'
    assert template_back[0].get_channel() == template_id[0].get_channel(), 'template id chan failed'
    assert template_back[0].get_location_code() == template_id[0].get_location_code(), 'template id loc failed'
    assert template_back[0].get_phase() == template_id[0].get_phase(), 'template id phase failed'
    assert template_back[1] == template_id[1], 'template id failed'

    assert not detection.have_correlation_coefficient(), 'xc coeff bool not failed'
    detection.set_correlation_coefficient(xc_val)
    assert detection.have_correlation_coefficient(), 'onset xc coeff bool failed'
    assert detection.get_correlation_coefficient() == xc_val, 'xc val failed'
 
    assert not detection.have_phase_onset_time(), 'onset time bool not failed'
    detection.set_phase_onset_time(onset_time)
    assert detection.have_phase_onset_time(), 'onset time bool failed'
    assert detection.get_phase_onset_time() == onset_time, 'phase onset time failed'

    assert not detection.have_interpolated_phase_onset_time(), 'int onset time bool not failed'
    detection.set_interpolated_phase_onset_time(onset_int_time)
    assert detection.have_interpolated_phase_onset_time(), 'int onset time bool failed'
    assert detection.get_interpolated_phase_onset_time() == onset_int_time, 'int phase onset time failed'

    assert not detection.have_detection_time(), 'detection bool not failed'
    detection.set_detection_time(detection_time)
    assert detection.have_detection_time(), 'detection time bool failed'
    assert detection.get_detection_time() == detection_time, 'detection time failed'

    assert not detection.have_interpolated_detection_time(), 'interp detection not bool failed'
    detection.set_interpolated_detection_time(det_int_time)
    assert detection.have_interpolated_detection_time(), 'interp detection time bool failed'
    assert detection.get_interpolated_detection_time() == det_int_time, 'interp detection time failed'

    assert not detection.have_travel_time(), 'travel time not bool failed'
    detection.set_travel_time(travel_time)
    assert detection.have_travel_time(), 'travel time bool failed'
    assert detection.get_travel_time() == travel_time, 'travel time failed'

    assert not detection.have_amplitude_scaling_factor(), 'scale factor bool failed'
    rmtype = pymflib.RelativeMagnitudeType
    detection.set_amplitude_scaling_factor(scale_factor1, rmtype.gibbons_ringdal_2006)
    detection.set_amplitude_scaling_factor(scale_factor2, rmtype.schaff_richards_2014) 
    assert detection.get_amplitude_scaling_factor(rmtype.gibbons_ringdal_2006) == scale_factor1, 'gr 1 failed'
    assert detection.get_amplitude_scaling_factor(rmtype.schaff_richards_2014) == scale_factor2, 'sr 1 failed' 

    assert not detection.have_detected_signal(), 'detection signal bool not failed'
    detection.set_detected_signal(xt)
    assert detection.have_detected_signal(), 'detected signal bool failed'
    xback = detection.get_detected_signal()
    assert np.max(np.abs(xback - xt)) < 1.e-14, 'detected signal not recovered'

    # Test copy c'tor
    det_copy = detection
    assert det_copy.have_template_identifier(), 'copy template id bool failed'
    template_back = det_copy.get_template_identifier()
    assert template_back[0].get_network() == template_id[0].get_network(), 'template id net failed'
    assert template_back[0].get_station() == template_id[0].get_station(), 'template id stat failed'
    assert template_back[0].get_channel() == template_id[0].get_channel(), 'template id chan failed'
    assert template_back[0].get_location_code() == template_id[0].get_location_code(), 'template id loc failed'
    assert template_back[0].get_phase() == template_id[0].get_phase(), 'template id phase failed'
    assert template_back[1] == template_id[1], 'template id failed'
    assert det_copy.have_correlation_coefficient(), 'copy xc val bool failed'
    assert det_copy.get_correlation_coefficient(), 'copy xc val failed'
    assert det_copy.have_phase_onset_time(), 'copy onset time bool failed'
    assert det_copy.have_interpolated_phase_onset_time(), 'copy int time bool failed'
    assert det_copy.get_phase_onset_time() == onset_time, 'copy phase onset time failed'
    assert det_copy.get_interpolated_phase_onset_time(), 'copy int phase onset time failed'
    assert det_copy.have_detection_time(), 'copy detection time bool failed'
    assert det_copy.get_detection_time() == detection_time, 'copy detection time failed'
    assert det_copy.get_interpolated_detection_time() == det_int_time, 'copy int detection time failed'
    assert det_copy.get_travel_time() == travel_time, 'copy travel time failed'
    assert det_copy.get_amplitude_scaling_factor(rmtype.gibbons_ringdal_2006) == scale_factor1, 'gr 1 failed'
    assert det_copy.get_amplitude_scaling_factor(rmtype.schaff_richards_2014) == scale_factor2, 'sr 1 failed'
    xback = det_copy.get_detected_signal()
    assert np.max(np.abs(xback - xt)) < 1.e-14, 'copy detected signal not recovered'

    #assert np.max(np.abs(xback - xb)) < 1.e-14, 'detected signal not recovered'

    # Test detection
    data = pickle.dumps(det_copy, 2)
    p_det = pickle.loads(data)
    assert p_det.have_template_identifier(), 'pickle template id bool failed'
    template_back = p_det.get_template_identifier()
    assert template_back[0].get_network() == template_id[0].get_network(), 'pickle template id net failed'
    assert template_back[0].get_station() == template_id[0].get_station(), 'pickle template id stat failed'
    assert template_back[0].get_channel() == template_id[0].get_channel(), 'pickle template id chan failed'
    assert template_back[0].get_location_code() == template_id[0].get_location_code(), 'pickle template id loc failed'
    assert template_back[0].get_phase() == template_id[0].get_phase(), 'pickle template id phase failed'
    assert template_back[1] == template_id[1], 'pickle template id failed'
    assert p_det.have_correlation_coefficient(), 'pickle xc val bool failed'
    assert p_det.get_correlation_coefficient(), 'pickle xc val failed'
    assert p_det.have_phase_onset_time(), 'pickle onset time bool failed'
    assert p_det.have_interpolated_phase_onset_time(), 'pickle int time bool failed'
    assert p_det.get_phase_onset_time() == onset_time, 'pickle phase onset time failed'
    assert p_det.get_interpolated_phase_onset_time(), 'pickle int phase onset time failed'
    assert p_det.have_detection_time(), 'pickle detection time bool failed'
    assert p_det.get_detection_time() == detection_time, 'pickle detection time failed'
    assert p_det.get_interpolated_detection_time() == det_int_time, 'pickle int detection time failed'
    assert p_det.get_travel_time() == travel_time, 'pick travel time failed'
    assert p_det.get_amplitude_scaling_factor(rmtype.gibbons_ringdal_2006) == scale_factor1, 'pickle gr 1 failed'
    assert p_det.get_amplitude_scaling_factor(rmtype.schaff_richards_2014) == scale_factor2, 'pickle sr 1 failed'
    # TODO - I don't know how ot serialize the detected waveform
    #xback = p_det.get_detected_signal()
    #print(xback)


def test_single_channel_detector_parameters():
    print("Single channel detector parameters test...")
    parms = pymflib.SingleChannel.DetectorParameters()
    min_spacing = 5 # samples
    det_tol = 0.5 # detection tolerance
    get_wave = True

    parms.set_maxima_policy(pymflib.MaxmimumMatchedFilterPolicy.absolute_maximum)
    parms.set_matched_filtered_signal_detector_policy(pymflib.MatchedFilteredSignalDetectorPolicy.single)
    parms.set_minimum_detection_spacing(min_spacing)
    parms.set_detection_threshold(det_tol)
    parms.enable_save_detected_waveform()
    parms.disable_save_amplitude_scaling_factor()

    assert parms.get_maxima_policy() == pymflib.MaxmimumMatchedFilterPolicy.absolute_maximum, 'max policy failed'
    assert parms.get_minimum_detection_spacing() == min_spacing, 'min spacing failed'
    assert parms.get_detection_threshold() == det_tol, 'det tol failed'
    assert parms.want_detected_waveform(), 'get det wave failed'
    assert parms.get_matched_filtered_signal_detector_policy() == \
           pymflib.MatchedFilteredSignalDetectorPolicy.single, 'mfs detector policy failed'
    assert not parms.want_amplitude_scaling_factor(), 'get det amp failed'

    # Test copy c'tor
    parms_copy = parms 
    assert parms_copy.get_maxima_policy() == pymflib.MaxmimumMatchedFilterPolicy.absolute_maximum, 'copy max policy failed'
    assert parms_copy.get_minimum_detection_spacing() == min_spacing, 'copy min spacing failed'
    assert parms_copy.get_detection_threshold() == det_tol, 'copy det tol failed'
    assert parms_copy.want_detected_waveform(), 'copy get det wave failed'
    assert parms_copy.get_matched_filtered_signal_detector_policy() == \
           pymflib.MatchedFilteredSignalDetectorPolicy.single, 'copy mfs detector policy failed'
    assert not parms_copy.want_amplitude_scaling_factor(), 'copy get the det amp failed'

def test_single_channel_associator_parameters():
    print("Single channel associator parameters test...")
    parms = pymflib.SingleChannel.AssociatorParameters()
    min_to_cluster = 8
    time_eps = 0.5
    parms.set_origin_time_tolerance(time_eps)
    parms.set_minimum_number_of_picks_in_event(min_to_cluster)
    parms.enable_correlation_coefficient_weighting()

    assert parms.get_origin_time_tolerance() == time_eps, 'ot failed'
    assert parms.get_minimum_number_of_picks_in_event() == min_to_cluster, 'min cluster failed'
    assert parms.use_correlation_coefficient_weighting(), 'toggle weight failed'

    parms_copy = parms
    assert parms_copy.get_origin_time_tolerance() == time_eps, 'ot failed copy'
    assert parms_copy.get_minimum_number_of_picks_in_event() == min_to_cluster, 'min cluster failed copy'
    assert parms_copy.use_correlation_coefficient_weighting(), 'toggle weight failed copy'

if __name__ == "__main__":
    test_network_station_phase()
    test_waveform_template()
    test_matched_filter_parameters()
    test_matched_filtering()
    test_single_channel_matched_filtering()
    test_single_channel_relative_magnitude()
    test_single_channel_detection()
    test_single_channel_detector_parameters()
    test_single_channel_associator_parameters()
