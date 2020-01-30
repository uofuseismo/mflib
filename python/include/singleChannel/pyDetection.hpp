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
    Detection(const Detection &det);
    Detection(const MFLib::SingleChannel::Detection<T> &det);
    Detection(Detection &&det) noexcept;
    ~Detection();
    void clear() noexcept;
    /// Operators
    Detection& operator=(const Detection &det);
    Detection& operator=(const MFLib::SingleChannel::Detection<T> &det);
    Detection& operator=(Detection &&parms) noexcept;
    /// Correlation coefficient
    void setCorrelationCoefficient(double xc);
    double getCorrelationCoefficient() const;
    bool haveCorrelationCoefficient() const noexcept;
    /// Waveform identifier
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
    /// Onset time
    void setInterpolatedDetectionTime(double time);
    double getInterpolatedDetectionTime() const;
    bool haveInterpolatedDetectionTime() const noexcept;
    /// Phase onset time
    void setPhaseOnsetTime(double onsetTime);
    double getPhaseOnsetTime() const;
    bool havePhaseOnsetTime() const noexcept; 
    /// Phase onset time
    void setInterpolatedPhaseOnsetTime(double onsetTime);
    double getInterpolatedPhaseOnsetTime() const;
    bool haveInterpolatedPhaseOnsetTime() const noexcept;
    /// Amplitude sclaing factor
    double getMagnitudePerturbation(MFLib::RelativeMagnitudeType type) const;
    double getAmplitudeScalingFactor(MFLib::RelativeMagnitudeType type) const;
    bool haveAmplitudeScalingFactor() const noexcept;
    ///
    MFLib::SingleChannel::Detection<T> getNativeClass() const;
private:
    std::unique_ptr<MFLib::SingleChannel::Detection<T>> mDetection;
};
void initializeDetection(pybind11::module &m);
}
}
#endif
