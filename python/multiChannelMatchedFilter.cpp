#include <cstdio>
#include <cstdlib>
#include "mflib/matchedFilter.hpp"
#include "mflib/matchedFilterParameters.hpp"
#include "multiChannelMatchedFilter.hpp"
#include "matchedFilterParameters.hpp"

using namespace PBMFLib;

/// Constructor
template<class T>
MultiChannelMatchedFilter<T>::MultiChannelMatchedFilter() :
    mcmf(std::make_unique<MFLib::MatchedFilter<T>> ())
{
}

/// Destructor
template<class T>
MultiChannelMatchedFilter<T>::~MultiChannelMatchedFilter() = default;

/// Initialize
template<class T>
void MultiChannelMatchedFilter<T>::initialize(
    const MatchedFilterParameters &parmsIn)
{
    MFLib::MatchedFilterParameters parms = parmsIn.getNativeClass();
    if (!parms.isValid())
    {
         throw std::invalid_argument("Parameters are invalid\n");
    }
    mcmf->initialize(parms);
}

/// Sets the signal
template<class T>
void MultiChannelMatchedFilter<T>::setSignal(
    const int it,
    const pybind11::array_t<double, pybind11::array::c_style |
                                    pybind11::array::forcecast> &x)
{
    if (!isInitialized())
    {
        throw std::runtime_error("Class is not initialized\n");
    }
    // Use pointers
    pybind11::buffer_info xbuf = x.request();
    auto len = static_cast<int> (xbuf.size);
    const double *xptr = (double *) (xbuf.ptr);
    if (xptr == nullptr)
    {   
        throw std::runtime_error("x is null");
    }
    mcmf->setSignal(it, len, xptr);
}

/// Zeros a signal
template<class T>
void MultiChannelMatchedFilter<T>::zeroSignal(const int it)
{
    if (!isInitialized())
    {
        throw std::runtime_error("Class is not initialized\n");
    }
    mcmf->zeroSignal(it);
}

/// Gets the number of templates
template<class T>
int MultiChannelMatchedFilter<T>::getNumberOfTemplates() const noexcept
{
    if (!isInitialized()){return 0;}
    return mcmf->getNumberOfTemplates();
}

/// Determines if the class is initialized
template<class T>
bool MultiChannelMatchedFilter<T>::isInitialized() const noexcept
{
    return mcmf->isInitialized();
}

/// Instantiation
template class PBMFLib::MultiChannelMatchedFilter<double>;
