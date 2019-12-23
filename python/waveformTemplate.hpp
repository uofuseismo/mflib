#ifndef PYMFLIB_WAVEFORMTEMPLATE_HPP
#define PYMFLIB_WAVEFORMTEMPLATE_HPP
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace PBWaveformTemplate
{
class WaveformTemplate
{
public:
    WaveformTemplate();
    ~WaveformTemplate();
    /// Sampling rate
    void setSamplingRate(double samplingRate);
    double getSamplingRate() const;
    /// Signal
    void setSignal(pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast> &x);
    //void setSignal(pybind11::array_t<float, pybind11::array::c_style | pybind11::array::forcecast> &x);
    pybind11::array_t<double> getSignal() const;
    /// Shift and stack weight
    void setShiftAndStackWeight(double weight);
    double getShiftAndStackWeight() const; 
private:
    std::unique_ptr<MFLib::WaveformTemplate> mWaveformTemplate;
};

}
#endif
