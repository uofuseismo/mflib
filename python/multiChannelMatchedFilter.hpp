#ifndef PYMFLIB_MULTICHANNELMATCHEDFILTER_HPP
#define PYMFLIB_MULTICHANNELMATCHEDFILTER_HPP
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace MFLib
{
template<class T> class MatchedFilter;
}

namespace PBMFLib
{
class MultiChannelMatchedFilterParameters;
template<class T = double>
class MultiChannelMatchedFilter
{
public:
    /// Constructors
    MultiChannelMatchedFilter();
    /// Destructors
    ~MultiChannelMatchedFilter();
    /// Set parameters
    void initialize(const MultiChannelMatchedFilterParameters &parms);
    /// Set signal
    void setSignal(int it,
                   const pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast> &x);
    void zeroSignal(int it);
    /// Is this class initiialized
    bool isInitialized() const noexcept; 
    int getNumberOfTemplates() const noexcept;
    /// Apply
    void apply();
    /// Have the matched filtered signals?
    bool haveMatchedFilteredSignals() const noexcept;
    pybind11::array_t<T, pybind11::array::c_style | pybind11::array::forcecast>
    getMatchedFilteredSignal(int it) const;
private:
    std::unique_ptr<MFLib::MatchedFilter<T>> mcmf;
//    MultiChannelMatchedFilter(const MultiChannelMatchedFilter &m) = delete;
//    MultiChannelMatchedFilter& operator=(const MultiChannelMatchedFilter &m) = delete;
};
}
#endif
