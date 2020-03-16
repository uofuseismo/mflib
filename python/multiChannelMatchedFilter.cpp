#include <cstdio>
#include <cstdlib>
#include "mflib/matchedFilter.hpp"
#include "mflib/matchedFilterParameters.hpp"
#include "pyMultiChannelMatchedFilter.hpp"
#include "pyMultiChannelMatchedFilterParameters.hpp"

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
    const MultiChannelMatchedFilterParameters &parmsIn)
{
    MFLib::MatchedFilterParameters parms = parmsIn.getNativeClass();
    if (!parms.isValid())
    {
         if (parms.getSignalSize() < 1)
         {
             fprintf(stderr, "Signal size not set\n");
         }
         if (parms.getNumberOfTemplates() < 1)
         {
             fprintf(stderr, "No templates\n");
         }
         throw std::invalid_argument("Parameters are invalid\n");
    }
    mcmf->initialize(parms);
}

/// Sets the signal
template<>
void MultiChannelMatchedFilter<double>::setSignal(
    const int it,
    pybind11::array_t<double, pybind11::array::c_style |
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
        throw std::invalid_argument("x is null");
    }
    mcmf->setSignal(it, len, xptr);
}

template<>
void MultiChannelMatchedFilter<float>::setSignal(
    const int it,
    pybind11::array_t<double, pybind11::array::c_style |
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
        throw std::invalid_argument("x is null");
    }
    std::vector<float> xIn(len);
    #pragma omp simd
    for (int i=0; i<len; ++i){xIn[i] = static_cast<float> (xptr[i]);}
    mcmf->setSignal(it, len, xIn.data());
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

/// Applies the matched filtering
template<class T>
void MultiChannelMatchedFilter<T>::apply()
{
    if (!isInitialized())
    {
        throw std::runtime_error("Matched filter not initialized\n");
    }
    mcmf->apply();
}

/// Determines if the matched filtering was applied
template<class T>
bool MultiChannelMatchedFilter<T>::haveMatchedFilteredSignals() const noexcept
{
    return mcmf->haveMatchedFilteredSignals();
}

/// Gets the matched filtered signals
template<class T>
pybind11::array_t<T, pybind11::array::c_style | pybind11::array::forcecast>
MultiChannelMatchedFilter<T>::getMatchedFilteredSignal(const int it) const
{
    if (!haveMatchedFilteredSignals())
    {
        throw std::runtime_error("Matched filtering not yet applied\n");
    }
    auto nt = getNumberOfTemplates();
    if (it < 0 || it >= nt)
    {
        throw std::invalid_argument("Template index = " + std::to_string(it)
                                  + " must be in range [0," + std::to_string(nt-1)
                                  + "]\n");
    }
    auto res = mcmf->getMatchedFilteredSignal(it);
    auto npts = res.size();
    auto y = pybind11::array_t<T, pybind11::array::c_style> (npts);
    pybind11::buffer_info ybuf = y.request();
    auto yptr = static_cast<T *> (ybuf.ptr);
    std::copy(res.begin(), res.end(), yptr);
    return y;
}

/// Initialize
void PBMFLib::initializeMultiChannelMatchedFilter(pybind11::module &m)
{
    pybind11::class_<PBMFLib::MultiChannelMatchedFilter<double>> dmcmf(m, "MatchedFilter");
    dmcmf.def(pybind11::init<> ());
    dmcmf.doc() = "Applies matched filtering to multiple channels where each channel has a specific template.";
    /// Initialize
    dmcmf.def("initialize",
              &PBMFLib::MultiChannelMatchedFilter<double>::initialize,
              "Initializes the multi-channel matched filtering engine.  This is a high-overhead function.");
    /// Sets signals
    dmcmf.def("set_signal",
              &PBMFLib::MultiChannelMatchedFilter<double>::setSignal,
              "Sets the signal corresponding to the it'th template.  This should be called after the class is initialized.  Moreover, the template index must in the range of [0, number_of_templates].");
    dmcmf.def("zero_signal",
              &PBMFLib::MultiChannelMatchedFilter<double>::zeroSignal,
              "Over a time period a station can be down.  This sets the signal corresponding to the it'th template to 0.  The template index must be in the range of [0, number_of_templates].");
    dmcmf.def("get_number_of_templates", 
              &PBMFLib::MultiChannelMatchedFilter<double>::getNumberOfTemplates, 
              "The number of templates in the multi-channel cross-correlation.");
    /// Apply
    dmcmf.def("apply",
              &PBMFLib::MultiChannelMatchedFilter<double>::apply,
              "Applies the templates set during the initialization stage to the signals set by set_signal");
    /// 
    dmcmf.def("have_matched_filtered_signals",
              &PBMFLib::MultiChannelMatchedFilter<double>::haveMatchedFilteredSignals,
              "Determines if the matched filtering as been applied and the signals are available for extraction from the class.");
    /// Gets the it'th matched filtered signal
    dmcmf.def("get_matched_filtered_signal",
              &PBMFLib::MultiChannelMatchedFilter<double>::getMatchedFilteredSignal,
              "Gets the it'th matched filtered signal.  The `apply' method must have been called and the template index must be in the range [0,number_of_templates]");
}

/// Instantiation
template class PBMFLib::MultiChannelMatchedFilter<double>;
template class PBMFLib::MultiChannelMatchedFilter<float>;
