#ifndef PYMFLIB_MATCHEDFILTERPARAMETERS_HPP
#define PYMFLIB_MATCHEDFILTERPARAMETERS_HPP
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "waveformTemplate.hpp"

namespace MFLib
{
class WaveformTemplate;
class MatchedFilterParameters;
}

namespace PBMFLib
{
class WaveformTemplate;
class MultiChannelMatchedFilterParameters
{
public:
    /// Constructor
    MultiChannelMatchedFilterParameters();
    MultiChannelMatchedFilterParameters(const MultiChannelMatchedFilterParameters &parms);
    MultiChannelMatchedFilterParameters& operator=(const MultiChannelMatchedFilterParameters &parms);
    /// Destructor
    ~MultiChannelMatchedFilterParameters();
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
    MFLib::MatchedFilterParameters getNativeClass() const; 
private:
    std::unique_ptr<MFLib::MatchedFilterParameters> mParameters;
};
}
#endif
