#ifndef PYMFLIB_WAVEFORMTEMPLATE_HPP
#define PYMFLIB_WAVEFORMTEMPLATE_HPP
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace MFLib
{
class WaveformTemplate;
}

namespace PBMFLib
{
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
    void setOnsetTime(double onsetTime);
    double getOnsetTime() const;
    /// Travel time
    void setTravelTime(double travelTime);
    double getTravelTime() const; 
    /// Signal
    void setSignal(pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast> &x);
    //void setSignal(pybind11::array_t<float, pybind11::array::c_style | pybind11::array::forcecast> &x);
    pybind11::array_t<double> getSignal() const;
    /// Shift and stack weight
    void setShiftAndStackWeight(double weight);
    double getShiftAndStackWeight() const; 
    /// Creates 
    MFLib::WaveformTemplate getNativeClass() const;
private:
    std::unique_ptr<MFLib::WaveformTemplate> mWaveformTemplate;
};

}
#endif
