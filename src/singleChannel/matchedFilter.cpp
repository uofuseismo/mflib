#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include "mflib/matchedFilterParameters.hpp"
#include "mflib/singleChannel/matchedFilterParameters.hpp"
#include "mflib/singleChannel/matchedFilter.hpp"
#include "mflib/matchedFilter.hpp"
#include "private/private.hpp"

using namespace MFLib::SingleChannel;

template<class T>
class MatchedFilter<T>::MatchedFilterImpl
{
public:
    MFLib::MatchedFilter<T> mMF;
};

/// Constructor
template<class T>
MatchedFilter<T>::MatchedFilter() :
    pImpl(std::make_unique<MatchedFilterImpl> ())
{
}

/// Destructor
template<class T>
MatchedFilter<T>::~MatchedFilter() = default;

/// Resets the class
template<class T>
void MatchedFilter<T>::clear() noexcept
{
    pImpl->mMF.clear();
}

/// Initializes the class
template<class T>
void MatchedFilter<T>::initialize(
    const MFLib::SingleChannel::MatchedFilterParameters &parmsIn)
{
    if (!parmsIn.isValid())
    {
        throw std::invalid_argument("Parameters are invalid\n");
    }
    auto parms = parmsIn.getParameters();
    pImpl->mMF.initialize(parms);    
}

/// Gets a template
template<class T>
MFLib::WaveformTemplate
MatchedFilter<T>::getWaveformTemplate(const int it) const
{
    return pImpl->mMF.getWaveformTemplate(it);
}

/// Sets a signal
template<class T>
void MatchedFilter<T>::setSignal(
    const int nSamples, const T signal[])
{
    constexpr int is = 0;  // Always the first signal
    pImpl->mMF.setSignal(is, nSamples, signal); // Throws
}

/// Gets a pointer to the signal to matched filter
template<class T>
const T* MatchedFilter<T>::getSignalPointer() const
{
    constexpr int is = 0;  // Always the first signal
    return pImpl->mMF.getSignalPointer(is);
}

/// Zeros a signal
template<class T>
void MatchedFilter<T>::zeroSignal()
{
    constexpr int is = 0;  // Always the first signal
    pImpl->mMF.zeroSignal(is); // Throws
}

/// Gets the number of templates
template<class T>
int MatchedFilter<T>::getNumberOfTemplates() const
{
    return pImpl->mMF.getNumberOfTemplates();
}

/// Gets the signal length
template<class T>
int MatchedFilter<T>::getSignalLength() const
{
    return pImpl->mMF.getSignalLength();
}

/// Gets the filtered signal length
template<class T>
int MatchedFilter<T>::getFilteredSignalLength() const
{
    return pImpl->mMF.getFilteredSignalLength();
}

/// Applies the matched filtering
template<class T>
void MatchedFilter<T>::apply()
{
    return pImpl->mMF.apply();
}

/// Gets the matched filtered signal
template<class T>
std::vector<T> MatchedFilter<T>::getMatchedFilteredSignal(
    const int it) const
{
    return pImpl->mMF.getMatchedFilteredSignal(it);
}

/// Gets a pointer to the matched filtered signal
template<class T>
const T* MatchedFilter<T>::getMatchedFilterSignalPointer(
    const int it) const
{
    return pImpl->mMF.getMatchedFilterSignalPointer(it);
}

/// Determines if the class is initialized
template<class T>
bool MatchedFilter<T>::isInitialized() const noexcept
{
    return pImpl->mMF.isInitialized();
}

/// Determines if the matched filter signals have been computed
template<class T>
bool MatchedFilter<T>::haveMatchedFilteredSignals() const noexcept
{
    return pImpl->mMF.haveMatchedFilteredSignals();
}

/// Template instantiation
template class MFLib::SingleChannel::MatchedFilter<double>;
template class MFLib::SingleChannel::MatchedFilter<float>;
