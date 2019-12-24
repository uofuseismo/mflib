#ifndef PYMFLIB_MATCHEDFILTEROPTIONS_HPP
#define PYMFLIB_MATCHEDFILTEROPTIONS_HPP
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace MFLib
{
class WaveformTemplate;
class MatchedFilterOptions;
}

namespace PBMFLib
{
class WaveformTemplate;
class MatchedFilterOptions
{
public:
    /// Constructor
    MatchedFilterOptions();
    /// Destructor
    ~MatchedFilterOptions();
    /// Add template
    void addTemplate(const WaveformTemplate &wt);
    /// Clears templates
    void clearTemplates() noexcept;
    /// Clears the memory
    void clear() noexcept;
private:
    std::unique_ptr<MFLib::MatchedFilterOptions> mOptions;
};
}
#endif
