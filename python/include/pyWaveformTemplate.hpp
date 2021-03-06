#ifndef PYMFLIB_PYWAVEFORMTEMPLATE_HPP
#define PYMFLIB_PYWAVEFORMTEMPLATE_HPP
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "mflib/enums.hpp"

namespace MFLib
{
class WaveformTemplate;
}

namespace PBMFLib
{
class NetworkStationPhase;
class WaveformTemplate
{
public:
    /// Constructor
    WaveformTemplate();
    /// Copy constructors
    WaveformTemplate(const WaveformTemplate &wt);
    WaveformTemplate(const MFLib::WaveformTemplate &wt);
    /// Copy assignment operator
    WaveformTemplate& operator=(const WaveformTemplate &wt);
    WaveformTemplate& operator=(const MFLib::WaveformTemplate &wt);
    /// Destructor
    ~WaveformTemplate();
    /// Clears the class
    void clear() noexcept;
    /// Sampling rate
    void setSamplingRate(double samplingRate);
    double getSamplingRate() const;
    /// Onset time
    void setPhaseOnsetTime(double onsetTime);
    double getPhaseOnsetTime() const;
    /// Travel time
    void setTravelTime(double travelTime);
    double getTravelTime() const; 
    /// Sets the magnitude
    void setMagnitude(double mag) noexcept;
    double getMagnitude() const;
    /// Sets the identifier
    void setIdentifier(const std::pair<PBMFLib::NetworkStationPhase, uint64_t> &id) noexcept;
    std::pair<PBMFLib::NetworkStationPhase, uint64_t> getIdentifier() const;
    /// Signal
    void setSignal(pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast> &x);
    void setSignalFromVector(const std::vector<double> &x);
    //void setSignal(const pybind11::array_t<float, pybind11::array::c_style | pybind11::array::forcecast> &x);
    pybind11::array_t<double> getSignal() const;
    std::vector<double> getSignalAsVector() const;
    /// Shift and stack weight
    void setShiftAndStackWeight(double weight);
    double getShiftAndStackWeight() const; 
    /// Polartiy
    void setPolarity(MFLib::Polarity polarity) noexcept;
    MFLib::Polarity getPolarity() const noexcept;
    /// Creates 
    MFLib::WaveformTemplate getNativeClass() const;
private:
    std::unique_ptr<MFLib::WaveformTemplate> mWaveformTemplate;
};
void initializeWaveformTemplate(pybind11::module &m);
}
#endif
