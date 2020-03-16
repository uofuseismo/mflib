#include <cstdio>
#include <cstdlib>
#include "mflib/singleChannel/matchedFilter.hpp"
#include "mflib/singleChannel/matchedFilterParameters.hpp"
#include "mflib/enums.hpp"
#include "singleChannel/pyMatchedFilter.hpp"
#include "singleChannel/pyMatchedFilterParameters.hpp"

using namespace PBMFLib::SingleChannel;

/// Constructor
template<class T>
MatchedFilter<T>::MatchedFilter() :
    scmf(std::make_unique<MFLib::SingleChannel::MatchedFilter<T>> ())
{
}

/// Destructor
template<class T>
MatchedFilter<T>::~MatchedFilter() = default;

/// Get handle to underlying class
template<class T>
const MFLib::SingleChannel::MatchedFilter<T>*
MatchedFilter<T>::getNativeClassPointer() const
{
    return scmf.get();
}

/// Initialize
template<class T>
void MatchedFilter<T>::initialize(
    const PBMFLib::SingleChannel::MatchedFilterParameters &parmsIn)
{
    MFLib::SingleChannel::MatchedFilterParameters parms = parmsIn.getNativeClass();
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
    scmf->initialize(parms);
}

/// Sets the signal
template<>
void MatchedFilter<double>::setSignal(
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
        throw std::runtime_error("x is null");
    }
    scmf->setSignal(len, xptr);
}

template<>
void MatchedFilter<float>::setSignal(
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
        throw std::runtime_error("x is null");
    }
    std::vector<float> xIn(len);
    #pragma omp simd
    for (int i=0; i<len; ++i){xIn[i] = static_cast<float> (xptr[i]);}
    scmf->setSignal(len, xIn.data());
}

/// Gets the number of templates
template<class T>
int MatchedFilter<T>::getNumberOfTemplates() const noexcept
{
    if (!isInitialized()){return 0;}
    return scmf->getNumberOfTemplates();
}

/// Determines if the class is initialized
template<class T>
bool MatchedFilter<T>::isInitialized() const noexcept
{
    return scmf->isInitialized();
}

/// Applies the matched filtering
template<class T>
void MatchedFilter<T>::apply()
{
    if (!isInitialized())
    {
        throw std::runtime_error("Matched filter not initialized\n");
    }
    scmf->apply();
}

/// Determines if the matched filtering was applied
template<class T>
bool MatchedFilter<T>::haveMatchedFilteredSignals() const noexcept
{
    return scmf->haveMatchedFilteredSignals();
}

/// Gets the matched filtered signals
template<class T>
pybind11::array_t<T, pybind11::array::c_style | pybind11::array::forcecast>
MatchedFilter<T>::getMatchedFilteredSignal(const int it) const
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
    auto res = scmf->getMatchedFilteredSignal(it);
    auto npts = res.size();
    auto y = pybind11::array_t<T, pybind11::array::c_style> (npts);
    pybind11::buffer_info ybuf = y.request();
    auto yptr = static_cast<T *> (ybuf.ptr);
    std::copy(res.begin(), res.end(), yptr);
    return y;
}

/// Creates the class
void PBMFLib::SingleChannel::initializeMatchedFilter(pybind11::module &m)
{
    pybind11::class_<PBMFLib::SingleChannel::MatchedFilter<double>> scmf(m, "MatchedFilter");
    scmf.def(pybind11::init<> ());
    scmf.doc() = "Utilities for applying multiple templates to a single channel waveform.";
    /// Initialize
    scmf.def("initialize",
             &PBMFLib::SingleChannel::MatchedFilter<double>::initialize,
             "Initializes the single-channel matched filtering engine.  This is a high-overhead function.");
    /// Sets signals
    scmf.def("set_signal",
              &PBMFLib::SingleChannel::MatchedFilter<double>::setSignal,
              "Sets the signal.  This should be called after the class is initialized.");
    //scmf.def("zero_signal",
    //          &PBMFLib::SingleChannel::MatchedFilter<double>::zeroSignal,
    //          "Over a time period a station can be down.  This sets the signal to 0.");
    scmf.def("get_number_of_templates",
             &PBMFLib::SingleChannel::MatchedFilter<double>::getNumberOfTemplates,
             "The number of templates in the single-channel multi-template cross-correlation.");
    /// Apply
    scmf.def("apply",
             &PBMFLib::SingleChannel::MatchedFilter<double>::apply,
             "Applies the templates set during the initialization stage to the signal set by set_signal");
    /// 
    scmf.def("have_matched_filtered_signals",
             &PBMFLib::SingleChannel::MatchedFilter<double>::haveMatchedFilteredSignals,
             "Determines if the matched filtering as been applied and the signals are available for extraction from the class.");
    /// Gets the it'th matched filtered signal
    scmf.def("get_matched_filtered_signal",
             &PBMFLib::SingleChannel::MatchedFilter<double>::getMatchedFilteredSignal,
             "Gets the it'th matched filtered signal.  The `apply' method must have been called and the template index must be in the range [0,number_of_templates]");




}

/// Instantiation
template class PBMFLib::SingleChannel::MatchedFilter<double>;
template class PBMFLib::SingleChannel::MatchedFilter<float>;
