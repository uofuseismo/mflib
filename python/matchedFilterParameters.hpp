#ifndef PYMFLIB_MATCHEDFILTERPARAMETERS_HPP
#define PYMFLIB_MATCHEDFILTERPARAMETERS_HPP
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace MFLib
{
class WaveformTemplate;
class MatchedFilterParameters;
}

namespace PBMFLib
{
class WaveformTemplate;
class MatchedFilterParameters
{
public:
    /// Constructor
    MatchedFilterParameters();
    /// Destructor
    ~MatchedFilterParameters();
    /// Add template
    void addTemplate(const WaveformTemplate &wt);
    /// Clears templates
    void clearTemplates() noexcept;
    /// Clears the memory
    void clear() noexcept;
private:
    std::unique_ptr<MFLib::MatchedFilterParameters> mOptions;
};
}
#endif
