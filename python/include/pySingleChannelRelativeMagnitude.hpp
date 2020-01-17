#ifndef PYMFLIB_PYSINGLECHANNEL_RELATIVEMAGNITUDE_HPP
#define PYMFLIB_PYSINGLECHANNEL_RELATIVEMAGNITUDE_HPP
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "mflib/enums.hpp"
namespace MFLib
{
namespace SingleChannel
{
template<class T> class RelativeMagnitude;
}
}

namespace PBMFLib
{
class WaveformTemplate;
namespace SingleChannel
{
template<class T = double>
class RelativeMagnitude
{
public:
    /// Constructor
    RelativeMagnitude();
    /// Copy c'tor
    RelativeMagnitude(const RelativeMagnitude &rmag);
    /// Copy assignment operator
    RelativeMagnitude operator=(const RelativeMagnitude &rmag);
    /// Destructor
    ~RelativeMagnitude();
    /// Initialize
    void initialize(const WaveformTemplate &wt);
    /// Determines if the class is initialized
    bool isInitialized() const noexcept;
    /// Gets the expected length of the detected signal
    int getDetectedSignalLength() const; 
    /// Sets the detected signal
    void setDetectedSignal(const pybind11::array_t<double, pybind11::array::c_style |
                                                           pybind11::array::forcecast> &x);
    /// Determines if the detected signal was set
    bool haveDetectedSignal() const noexcept;
    /// Computes the amplitude scaling factor
    T computeAmplitudeScalingFactor(const MFLib::RelativeMagnitudeType type) const;
    /// Computes the magnitude scaling factor
    T computeMagnitudePerturbation(const MFLib::RelativeMagnitudeType type) const;
private:
    std::unique_ptr<MFLib::SingleChannel::RelativeMagnitude<T>> mMag;
};
void initializeRelativeMagnitude(pybind11::module &m);
}
}
#endif
