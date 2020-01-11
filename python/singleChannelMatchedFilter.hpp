#ifndef PYMFLIB_SINGLECHANNEL_MATCHEDFILTER_HPP
#define PYMFLIB_SINGLECHANNEL_MATCHEDFILTER_HPP
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
    MatchedFilter();
    ~MatchedFilter();
private:
    std::unique_ptr<MFLib::SingleChannel::MatchedFilter<T>> scmf; 
};
void initializeMatchedFilter(pybind11::module &m);
}
}
#endif
