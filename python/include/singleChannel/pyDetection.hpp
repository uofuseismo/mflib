#ifndef PYMFLIB_SINGLECHANNEL_PYDETECTION_HPP
#define PYMFLIB_SINGLECHANNEL_PYDETECTION_HPP
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "mflib/singleChannel/detection.hpp"
#include "pyWaveformTemplate.hpp"

namespace MFLib
{
namespace SingleChannel
{
template<class T> class Detection;
}
}

namespace PBMFLib
{
namespace SingleChannel
{
template<class T>
class Detection
{
public:
    /// Constructor
    Detection();
    Detection(const Detection &parms);
    Detection(Detection  &&parms) noexcept;
    ~Detection();
    void clear() noexcept;
    /// Operators
    Detection& operator=(const Detection &parms);
    Detection& operator=(Detection &&parms) noexcept;
    //p Waveform identifier
    void setTemplateIdentifier(uint64_t id) noexcept;
    uint64_t getTemplateIdentifier() const;
    bool haveTemplateIdentifier() const noexcept;
    /// Detected signal
    void setDetectedSignal(const pybind11::array_t<double, pybind11::array::c_style |
                                                           pybind11::array::forcecast> &x);
    pybind11::array_t<double> getDetectedSignal() const;
    bool haveDetectedSignal() const noexcept;
    /// Onset time
    void setDetectionTime(double time);
    double getDetectionTime() const;
    bool haveDetectionTime() const noexcept;
    /// Phase onset time
    void setPhaseOnsetTime(double onsetTime);
    double getPhaseOnsetTime() const;
    bool havePhaseOnsetTime() const noexcept; 
    ///
    MFLib::SingleChannel::Detection<T> getNativeClass() const;
private:
    std::unique_ptr<MFLib::SingleChannel::Detection<T>> mDetection;
};
void initializeDetection(pybind11::module &m);
}
}
#endif
