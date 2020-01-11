#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include "mflib/singleChannelMatchedFilter.hpp"
#include "mflib/matchedFilter.hpp"

using namespace MFLib;

template<class T>
class SingleChannelMatchedFilter<T>::SingleChannelMatchedFilterImpl
{
public:
    MatchedFilter<T> mMF;
};

/// Constructor
template<class T>
SingleChannelMatchedFilter<T>::SingleChannelMatchedFilter() :
    pImpl(std::make_unique<SingleChannelMatchedFilterImpl> ())
{
}

/// Destructor
template<class T>
SingleChannelMatchedFilter<T>::~SingleChannelMatchedFilter() = default;

/// Resets the class
template<class T>
void SingleChannelMatchedFilter<T>::clear() noexcept
{
    pImpl->mMF->clear();
}

/// Sets a signal
template<class T>
void SingleChannelMatchedFilter<T>::setSignal(
    const int nSamples, const T signal[])
{
    constexpr int is = 0;  // Always the first signal
    pImpl->mMF->setSignal(is, nSamples, signal); // Throws
}

/// Zeros a signal
template<class T>
void SingleChannelMatchedFilter<T>::zeroSignal()
{
    constexpr int is = 0;  // Always the first signal
    pImpl->mMF->zeroSignal(is); // Throws
}

/// Gets the number of templates
template<class T>
int SingleChannelMatchedFilter<T>::getNumberOfTemplates() const
{
    return pImpl->mMF->getNumberOfTemplates();
}

/// Gets the signal length
template<class T>
int SingleChannelMatchedFilter<T>::getSignalLength() const
{
    return pImpl->mMF->getSignalLength();
}

/// Gets the filtered signal length
template<class T>
int SingleChannelMatchedFilter<T>::getFilteredSignalLength() const
{
    return pImpl->mMF->getFilteredSignalLength();
}

/// Applies the matched filtering
template<class T>
void SingleChannelMatchedFilter<T>::apply()
{
    return pImpl->mMF->apply();
}

/// Gets the matched filtered signal
template<class T>
std::vector<T> SingleChannelMatchedFilter<T>::getMatchedFilteredSignal(
    const int it) const
{
    return pImpl->mMF->getMatchedFilteredSignal(it);
}

/// Gets a pointer to the matched filtered signal
template<class T>
const T* SingleChannelMatchedFilter<T>::getMatchedFilterSignalPointer(
    const int it) const
{
    return pImpl->mMF->getMatchedFilterSignalPointer(it);
}

/// Determines if the class is initialized
template<class T>
bool SingleChannelMatchedFilter<T>::isInitialized() const noexcept
{
    return pImpl->mMF->isInitialized();
}

/// Determines if the matched filter signals have been computed
template<class T>
bool SingleChannelMatchedFilter<T>::haveMatchedFilteredSignals() const noexcept
{
    return pImpl->mMF->haveMatchedFilteredSignals();
}
