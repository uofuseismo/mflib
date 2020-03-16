#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include "mflib/singleChannel/associator.hpp"
#include "mflib/singleChannel/detection.hpp"
#include "mflib/singleChannel/event.hpp"
#include "mflib/networkStationPhase.hpp"
#include "private/weightedStatistics.hpp"

using namespace MFLib::SingleChannel;

namespace
{

/// Returns the unique event IDs in our list of picks
template<class T>
std::vector<uint64_t> getUniqueEVIDs(
    std::vector<MFLib::SingleChannel::Detection<T>> &detections)
{
    std::vector<uint64_t> evids;
    evids.reserve(detections.size());
    for (int i=0; i<static_cast<int> (detections.size()); ++i)
    {
        auto id = detections[i].getIdentifier();
        auto nsp = id.first;
        auto evid = id.second;
        bool lfound = false;
        for (int j=0; j<evids.static_cast<int> (evids.size()); ++j)
        {
            if (evid == evids[j])
            {
                lfound = true;
                break;
            }
        }
        if (!lfound)
        {
            evids.push_back(evid);
        }
    }
    std::sort(evids.begin(), evids.end());
    return evids;
}

}

class Associator::AssociatorImpl
{
public:
    std::vector<MFLib::SingleChannel::Event> mEvents;
};

/// Constructor
Associator::Associator() :
    pImpl(std::make_unique<AssociatorImpl> ())
{
}

/// Copy assignment operator
Associator& Associator::operator=(const Associator &associator)
{
    if (&associator == this){return *this;}
    pImpl = std::make_unique<AssociatorImpl> (*associator.pImpl);
    return *this;
}

/// Move assignment operator
Associator& Associator::operator=(Associator &&associator) noexcept
{
    if (&associator == this){return *this;}
    pImpl= std::move(associator.pImpl);
    return *this;
}

/// Destructor
Associator::~Associator() = default;

/// 
void Associator::associate()
{
    pImpl->mEvents.clear();
    // Loop over the associated picks
    
}

/// Gets the number of associated events
int Associator::getNumberOfEvents() const noexcept
{
    return static_cast<int> (pImpl->mEvents.size());
}
