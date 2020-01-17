#ifndef PYMFLIB_SINGLECHANNELFILTERPARAMETERS_HPP
#define PYMFLIB_SINGLECHANNELFILTERPARAMETERS_HPP
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace MFLib
{
class WaveformTemplate;
namespace SingleChannel
{
class MatchedFilterParameters;
}
}

namespace PBMFLib
{
class WaveformTemplate;
namespace SingleChannel
{
class MatchedFilterParameters
{
public:
    /// Constructor
    MatchedFilterParameters();
    MatchedFilterParameters(const MatchedFilterParameters &parms);
    MatchedFilterParameters& operator=(const MatchedFilterParameters &parms);
    /// Destructor
    ~MatchedFilterParameters();
    /// Add template
    void addTemplate(const WaveformTemplate &wt);
    int getNumberOfTemplates() const;
    WaveformTemplate getTemplate(int indx) const;
    /// Clears templates
    void clearTemplates() noexcept;
    /// Sets/gets the signal size
    void setSignalSize(int n);
    int getSignalSize() const;
    /// Sets/gets the FFT length
    void setFFTLength(int fftLength);
    int getFFTLength() const;
    int getBlockLength() const;
    /// Clears the memory
    void clear() noexcept;
    /// Gets the native parameters class
    MFLib::SingleChannel::MatchedFilterParameters getNativeClass() const;
private:
    std::unique_ptr<MFLib::SingleChannel::MatchedFilterParameters> mParameters;
};
}
}
#endif

