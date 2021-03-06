#ifndef PYMFLIB_SINGLECHANNEL_PYDETECTION_HPP
#define PYMFLIB_SINGLECHANNEL_PYDETECTION_HPP
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "mflib/singleChannel/detection.hpp"
#include "mflib/networkStationPhase.hpp"
#include "pyWaveformTemplate.hpp"

/*
namespace MFLib
{
namespace SingleChannel
{
template<class T> class Detection;
}
}
*/

namespace PBMFLib
{
class NetworkStationPhase;
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
    Detection& operator=(MFLib::SingleChannel::Detection<T> &&det) noexcept;
    /// Correlation coefficient
    void setCorrelationCoefficient(double xc);
    double getCorrelationCoefficient() const;
    bool haveCorrelationCoefficient() const noexcept;
    /// Waveform identifier
    void setTemplateIdentifier(const std::pair<PBMFLib::NetworkStationPhase, uint64_t> &id) noexcept;
    std::pair<PBMFLib::NetworkStationPhase, uint64_t> getTemplateIdentifier() const;
    bool haveTemplateIdentifier() const noexcept;
    /// Detected signal
    void setDetectedSignal(pybind11::array_t<double, pybind11::array::c_style |
                                             pybind11::array::forcecast> &x);
    //void setDetectedSignal(pybind11::array_t<float, pybind11::array::c_style |
    //                                         pybind11::array::forcecast> &x);
    pybind11::array_t<double> getDetectedSignal() const;
    bool haveDetectedSignal() const noexcept;
    void setDetectedSignalAsVector(const std::vector<T> &x);
    void getDetectedSignalAsVector(std::vector<T> &x) const;
    /// Onset time
    void setDetectionTime(double time);
    double getDetectionTime() const;
    bool haveDetectionTime() const noexcept;
    /// Inteprolated onset time
    void setInterpolatedDetectionTime(double time);
    double getInterpolatedDetectionTime() const;
    bool haveInterpolatedDetectionTime() const noexcept;
    /// Phase onset time
    void setPhaseOnsetTime(double onsetTime);
    double getPhaseOnsetTime() const;
    bool havePhaseOnsetTime() const noexcept; 
    /// Interpolated phase onset time
    void setInterpolatedPhaseOnsetTime(double onsetTime);
    double getInterpolatedPhaseOnsetTime() const;
    bool haveInterpolatedPhaseOnsetTime() const noexcept;
    /// Travel time
    void setTravelTime(double travelTime);
    double getTravelTime() const;
    bool haveTravelTime() const noexcept;
    /// Amplitude sclaing factor
    void setAmplitudeScalingFactor(double value, MFLib::RelativeMagnitudeType type);
    double getMagnitudePerturbation(MFLib::RelativeMagnitudeType type) const;
    double getAmplitudeScalingFactor(MFLib::RelativeMagnitudeType type) const;
    bool haveAmplitudeScalingFactor() const noexcept;
    /// Detection
    void setPolarity(MFLib::Polarity polarity) noexcept;
    MFLib::Polarity getPolarity() const noexcept;
    ///
    MFLib::SingleChannel::Detection<T> getNativeClass() const;
private:
    std::unique_ptr<MFLib::SingleChannel::Detection<T>> mDetection;
};
void initializeDetection(pybind11::module &m);
}
}
#endif
