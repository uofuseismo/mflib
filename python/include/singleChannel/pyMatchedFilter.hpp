#ifndef PYMFLIB_PYSINGLECHANNEL_MATCHEDFILTER_HPP
#define PYMFLIB_PYSINGLECHANNEL_MATCHEDFILTER_HPP
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace MFLib
{
namespace SingleChannel
{
template<class T> class MatchedFilter;
}
}

namespace PBMFLib
{
namespace SingleChannel
{
class MatchedFilterParameters;
template<class T = double>
class MatchedFilter
{
public:
    /// Constructor
    MatchedFilter();
    /// Destructor
    ~MatchedFilter();
    /// Initialize the class
    void initialize(const PBMFLib::SingleChannel::MatchedFilterParameters &parms);
    /// Set signal
    void setSignal(const pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast> &x);
    /// Is this class initiialized
    bool isInitialized() const noexcept; 
    int getNumberOfTemplates() const noexcept;
    /// Apply
    void apply();
    /// Have the matched filtered signals?
    bool haveMatchedFilteredSignals() const noexcept;
    pybind11::array_t<T, pybind11::array::c_style | pybind11::array::forcecast>
    getMatchedFilteredSignal(int it) const;

    const MFLib::SingleChannel::MatchedFilter<T>* getNativeClassPointer() const;
private:
    std::unique_ptr<MFLib::SingleChannel::MatchedFilter<T>> scmf; 
};
void initializeMatchedFilter(pybind11::module &m);
}
}
#endif
