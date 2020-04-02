#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <algorithm>
#if __has_include(<pstl/algorithm>)
#include <pstl/algorithm>
#include <pstl/execution>
#define USE_PSTL 1
#endif
#include "mflib/singleChannel/associator.hpp"
#include "mflib/singleChannel/associatorParameters.hpp"
#include "mflib/singleChannel/detection.hpp"
#include "mflib/singleChannel/event.hpp"
#include "mflib/networkStationPhase.hpp"
#include "private/dbscan.hpp"
#include "private/weightedStatistics.hpp"

using namespace MFLib::SingleChannel;

namespace
{

}

template<class T>
class Associator<T>::AssociatorImpl
{
public:
    //MFLib::DBSCAN mDBSCAN;
    MFLib::SingleChannel::AssociatorParameters mParameters;
    std::vector<MFLib::SingleChannel::Detection<T>> mDetections;
    bool mInitialized = false;
};

/// Constructor
template<class T>
Associator<T>::Associator() :
    pImpl(std::make_unique<AssociatorImpl> ())
{
}

/// Copy assignment operator
template<class T>
Associator<T>& Associator<T>::operator=(const Associator &associator)
{
    if (&associator == this){return *this;}
    pImpl = std::make_unique<AssociatorImpl> (*associator.pImpl);
    return *this;
}

/// Move assignment operator
template<class T>
Associator<T>& Associator<T>::operator=(Associator &&associator) noexcept
{
    if (&associator == this){return *this;}
    pImpl= std::move(associator.pImpl);
    return *this;
}

/// Destructor
template<class T>
Associator<T>::~Associator() = default;

/// Releases memory on the class and resets it
template<class T>
void Associator<T>::clear() noexcept
{
    pImpl->mParameters.clear();
    pImpl->mDetections.clear();
    pImpl->mInitialized = false;
}

/// Determines if the class is intitialized
template<class T>
bool Associator<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// 
template<class T>
void Associator<T>::associate()
{
    if (!isInitialized()){throw std::runtime_error("Class is not intialized\n");}
    // Sort the detections based on event IDs
    std::sort(
#ifdef USE_PSTL
              std::execution::par,
#endif
              pImpl->mDetections.begin(), pImpl->mDetections.end(),
              [] (const MFLib::SingleChannel::Detection<T> &a,
                  const MFLib::SingleChannel::Detection<T> &b)
              {
                  return a.getTemplateIdentifier().second <
                         b.getTemplateIdentifier().second;
              });
    // Get the unique event IDs
    std::vector<uint64_t> evids;
    evids.reserve(pImpl->mDetections.size());
    for (const auto &det : pImpl->mDetections)
    {
        evids.push_back(det.getTemplateIdentifier().second);
    }
    auto uniqueEvids = evids;
#ifdef USE_PSTL
    auto it = std::unique(std::execution::par,
                          uniqueEvids.begin(), uniqueEvids.end());
#else
    auto it = std::unique(uniqueEvids.begin(), uniqueEvids.end());
#endif
    uniqueEvids.resize(std::distance(uniqueEvids.begin(), it));
    // Loop over the associated picks
    MFLib::DBSCAN dbscan;
}

/// Gets the number of associated events
/*
template<class T>
int Associator<T>::getNumberOfEvents() const noexcept
{
    return static_cast<int> (pImpl->mEvents.size());
}
*/

/// Implementation
template class MFLib::SingleChannel::Associator<double>;
template class MFLib::SingleChannel::Associator<float>;
