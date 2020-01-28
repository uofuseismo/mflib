#ifndef PYMFLIB_SINGLECHANNEL_PYDETECTORPARAMETERS_HPP
#define PYMFLIB_SINGLECHANNEL_PYDETECTORPARAMETERS_HPP
#include <memory>
#include <pybind11/pybind11.h>
#include "mflib/singleChannel/detectorParameters.hpp"

namespace MFLib
{
namespace SingleChannel
{
class DetectorParameters;
}
}

namespace PBMFLib
{
namespace SingleChannel
{
class DetectorParameters
{
public:
    /// Constructors
    DetectorParameters();
    DetectorParameters(const DetectorParameters &parms);
    DetectorParameters(const MFLib::SingleChannel::DetectorParameters &parms);
    /// Operators
    DetectorParameters& operator=(const DetectorParameters &parms);
    DetectorParameters& operator=(const MFLib::SingleChannel::DetectorParameters &parms);
    MFLib::SingleChannel::DetectorParameters getNativeClass() const;
    /// Destructors
    ~DetectorParameters();
    void clear() noexcept;
    /// Behavioral
    void setMaximaPolicy(MFLib::MaximumMatchedFilterPolicy policy) noexcept;
    MFLib::MaximumMatchedFilterPolicy getMaximaPolicy() const noexcept;
    void setMinimumDetectionSpacing(int nSamples) noexcept;
    int getMinimumDetectionSpacing() const noexcept;
    void setDetectionThreshold(double tol);
    double getDetectionThreshold() const noexcept;
    /// Return values
    void toggleGetDetectedWaveform(bool lget) noexcept;
    bool getDetectedWaveform() const noexcept; 
private:
    std::unique_ptr<MFLib::SingleChannel::DetectorParameters> mParameters;
};
void initializeDetectorParameters(pybind11::module &m);
}
}
#endif
