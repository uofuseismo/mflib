#include <cstdio>
#include <cstdlib>
#include "mflib/singleChannel/matchedFilter.hpp"
#include "mflib/singleChannel/matchedFilterParameters.hpp"
#include "mflib/enums.hpp"
#include "pySingleChannelMatchedFilter.hpp"
//#include "pySingleChannelMatchedFilterParameters.hpp"

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

/// Creates the class
void PBMFLib::SingleChannel::initializeMatchedFilter(pybind11::module &m)
{
    //m.doc() = "Utilities for applying multiple templates to a single channel waveform.";
    pybind11::class_<PBMFLib::SingleChannel::MatchedFilter<double>> scmf(m, "SingleChannelMatchedFilter");
    scmf.def(pybind11::init<> ());
    scmf.doc() = "Utilities for applying multiple templates to a single channel waveform.";
}
