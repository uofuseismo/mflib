#include <iostream>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "mflib/enums.hpp"
#include "mflib/waveformTemplate.hpp"
#include "mflib/networkStationPhase.hpp"
#include "pyWaveformTemplate.hpp"
#include "pyNetworkStationPhase.hpp"

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
    pybind11::array_t<double, pybind11::array::c_style |
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

void WaveformTemplate::setSignalFromVector(const std::vector<double> &signal)
{
    if (signal.empty())
    {
        throw std::invalid_argument("x has not point");
    }
    mWaveformTemplate->setSignal(signal.size(), signal.data());
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

std::vector<double> WaveformTemplate::getSignalAsVector() const
{
     auto npts = mWaveformTemplate->getSignalLength(); 
     std::vector<double> y(npts, 0);
     auto yPtr = y.data();
     mWaveformTemplate->getSignal(npts, &yPtr);
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
void WaveformTemplate::setPhaseOnsetTime(const double onsetTime)
{
    if (onsetTime < 0)
    {
        throw std::invalid_argument("onset_time = " 
                                  + std::to_string(onsetTime)
                                  + " must be positive\n");
    }
    mWaveformTemplate->setPhaseOnsetTime(onsetTime);
}

double WaveformTemplate::getPhaseOnsetTime() const
{
    return mWaveformTemplate->getPhaseOnsetTime();
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

/// Setters/getters for magnitude
void WaveformTemplate::setMagnitude(const double mag) noexcept
{
    mWaveformTemplate->setMagnitude(mag);
}

double WaveformTemplate::getMagnitude() const
{
    return mWaveformTemplate->getMagnitude();
}

/// Setters/getters for polarity
void WaveformTemplate::setPolarity(const MFLib::Polarity polarity) noexcept
{
    mWaveformTemplate->setPolarity(polarity);
}

MFLib::Polarity WaveformTemplate::getPolarity() const noexcept
{
    return mWaveformTemplate->getPolarity();
}

/// Setters/getters for identifier
void WaveformTemplate::setIdentifier(
    const std::pair<PBMFLib::NetworkStationPhase, uint64_t> &idIn) noexcept
{
    std::pair<MFLib::NetworkStationPhase, uint64_t> id;
    id.first = idIn.first.getNativeClass();
    id.second = idIn.second;
    mWaveformTemplate->setIdentifier(id);
}

std::pair<PBMFLib::NetworkStationPhase, uint64_t>
WaveformTemplate::getIdentifier() const
{
    auto id = mWaveformTemplate->getIdentifier();
    PBMFLib::NetworkStationPhase nsp(id.first);
    std::pair<PBMFLib::NetworkStationPhase, uint64_t> idOut(nsp, id.second);
    return id;
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
    wt.def("get_phase_onset_time_in_signal",
           &PBMFLib::WaveformTemplate::getPhaseOnsetTime,
           "The phase onset time of the relative to the signal start in seconds.");
    wt.def("set_phase_onset_time_in_signal",
           &PBMFLib::WaveformTemplate::setPhaseOnsetTime,
           "This defines the arrival onset time in seconds relative to the template signal start.  For example, if this is 2, then 2 seconds into the template signal is the expected arrival onset.  This is only necessary if you are shifting and stacking waveforms.  Note - the signal and sampling rate should be set prior to setting this property.");
    /// The travel time 
    wt.def("get_phase_travel_time",
           &PBMFLib::WaveformTemplate::getTravelTime,
           "The phase travel time from the source to the receiver in seconds.");
    wt.def("set_phase_travel_time",
           &PBMFLib::WaveformTemplate::setTravelTime,
           "This defines the phase's travel time in seconds from the origin time.  For example, if this is 9 then it took 9 seconds for the phase to travel from the source to the receiver.  This is only necessary if you are shifting and stacking waveforms.  Note - the signal and sampling rate should be set prior to setting this property.");
    /// The magnitude
    wt.def("get_magnitude",
           &PBMFLib::WaveformTemplate::getMagnitude,
           "Gets the magnitude associated with this template.");
    wt.def("set_magnitude",
           &PBMFLib::WaveformTemplate::setMagnitude,
           "Sets the magnitude associated with this template.");
    /// The polarity
    wt.def("get_polarity",
           &PBMFLib::WaveformTemplate::getPolarity,
           "Gets the polarity of the onset.");
    wt.def("set_polarity",
           &PBMFLib::WaveformTemplate::setPolarity,
           "Sets the polarity of the onset.");
    /// The waveform identifier
    wt.def("get_identifier",
           &PBMFLib::WaveformTemplate::getIdentifier,
           "Gets the template identifier.");
    wt.def("set_identifier",
           &PBMFLib::WaveformTemplate::setIdentifier,
           "Sets the template identifier.  An identifier consists of the network name, station name, and seismic phase identifier as well as the event identifier to which the arrival was associated in the catalog.  This is essential for association when using the single channel pipeline.");
    /// Clears the class memory and resets
    wt.def("clear",
           &PBMFLib::WaveformTemplate::clear,
           "Clears the class's memory and resets the class.");
    /// Makes this class pickleable
    wt.def("__getstate__",
           [](const PBMFLib::WaveformTemplate &wt)
    {
        // Template identifier
        std::string network = "";
        std::string station = "";
        std::string channel = "";
        std::string location = "";
        std::string phase = "";
        uint64_t waveid = 0;
        try
        {
            auto id = wt.getIdentifier();
            network = id.first.getNetwork();
            station = id.first.getStation();
            channel = id.first.getChannel();
            location = id.first.getLocationCode();
            phase = id.first.getPhase();
            waveid = id.second;
        }
        catch (const std::exception &e)
        {
        }
        double df = 0;
        double travelTime =-1;
        double onsetTime =-1;
        double weight =-1;
        double magnitude = std::numeric_limits<double>::min();
        try
        {
            df = wt.getSamplingRate();
            travelTime = wt.getTravelTime();
            onsetTime = wt.getPhaseOnsetTime();
            weight = wt.getShiftAndStackWeight(); 
            magnitude = wt.getMagnitude();
        }
        catch (const std::exception &e)
        {
        }
        auto polarity = static_cast<int> (wt.getPolarity());
        std::string x;
        try
        {
            auto xWork = wt.getSignalAsVector();
            x.reserve(16*xWork.size());
            for (size_t i=0; i<xWork.size(); ++i)
            {
                 x = x + std::to_string(xWork[i]);
                 if (i < xWork.size() - 1){x = x + ",";}
            }
        }
        catch (const std::exception &e)
        {
        }
        return pybind11::make_tuple(network, station, channel, location, phase, waveid,
                                    df, travelTime, onsetTime, weight, magnitude,
                                    polarity,
                                    x);
    });
    wt.def("__setstate__",
           [](PBMFLib::WaveformTemplate &p,
              pybind11::tuple t)
    {
        if (t.size() != 13)
        {
            std::cerr << "Tuple in invalid state" << std::endl;
        }
        // Call c'tor -> memory leak because of new?
        new (&p) PBMFLib::WaveformTemplate(); 

        PBMFLib::NetworkStationPhase nsp;
        nsp.setNetwork(t[0].cast<std::string> ());
        nsp.setStation(t[1].cast<std::string> ());
        nsp.setChannel(t[2].cast<std::string> ());
        nsp.setLocationCode(t[3].cast<std::string> ());
        nsp.setPhase(t[4].cast<std::string> ());
        uint64_t waveid = t[5].cast<uint64_t> ();
        p.setIdentifier(std::make_pair(nsp, waveid));

        double df = t[6].cast<double> ();
        double travelTime = t[7].cast<double> ();
        double onsetTime = t[8].cast<double> ();
        double weight = t[9].cast<double> ();
        double magnitude = t[10].cast<double> ();
        auto polarity = static_cast<MFLib::Polarity> (t[11].cast<int> ());

        auto csignal = t[12].cast<std::string> ();
        if (csignal.size() > 0)
        {
            std::vector<double> signal;
            signal.reserve(csignal.size()/8); 
            std::stringstream ss(csignal);
            for (double d; ss >> d;)
            {
                signal.push_back(d);
                if (ss.peek() == ','){ss.ignore();}
            }
            p.setSignalFromVector(signal);
 
            if (df > 0){p.setSamplingRate(df);}
            if (travelTime >= 0){p.setTravelTime(travelTime);}
            if (onsetTime >= 0){p.setPhaseOnsetTime(onsetTime);}
            if (weight > 0 && weight < 1){p.setShiftAndStackWeight(weight);}
            if (magnitude > std::numeric_limits<double>::min())
            {
                p.setMagnitude(magnitude);
            }
            p.setPolarity(polarity);
         }
    });
}
