#include <cstdio>
#include <cstdlib>
#include <string>
#include <functional>
#include <algorithm>
#include "mflib/networkStationPhase.hpp"

using namespace MFLib;

namespace
{
std::string makeHashString(const std::string &network,
                           const std::string &station,
                           const std::string &phase)
{
    return network + "." + station + "." + phase;
}
}

class NetworkStationPhase::NetworkStationPhaseImpl
{
public:
    std::string mNetwork;
    std::string mStation;
    std::string mPhase;
    size_t mHash = 0;
};

/// Constructor
NetworkStationPhase::NetworkStationPhase() :
    pImpl(std::make_unique<NetworkStationPhaseImpl> ())
{
}

/// Copy c'tor
NetworkStationPhase::NetworkStationPhase(const NetworkStationPhase &nsp)
{
    *this = nsp;
}

/// Move c'tor
NetworkStationPhase::NetworkStationPhase(NetworkStationPhase &&nsp) noexcept
{
    *this = std::move(nsp);
}

/// Destructor
NetworkStationPhase::~NetworkStationPhase() = default;

/// Clears memory and resets the class
void NetworkStationPhase::clear() noexcept
{
    pImpl->mNetwork.clear();
    pImpl->mStation.clear();
    pImpl->mPhase.clear();
    pImpl->mHash = 0; 
}

/// Copy assignment
NetworkStationPhase& 
NetworkStationPhase::operator=(const NetworkStationPhase &nsp)
{
    if (&nsp == this){return *this;}
    pImpl = std::make_unique<NetworkStationPhaseImpl> (*nsp.pImpl);
    return *this;
}

/// Move assignment
NetworkStationPhase&
NetworkStationPhase::operator=(NetworkStationPhase &&nsp) noexcept
{
    if (&nsp == this){return *this;}
    pImpl = std::move(nsp.pImpl);
    return *this; 
}

/// Sets the network
void NetworkStationPhase::setNetwork(const std::string &s) noexcept
{
    // Copy, remove blanks, and convert to uppercase
    pImpl->mNetwork = s;
    pImpl->mNetwork.erase(std::remove(pImpl->mNetwork.begin(),
                                      pImpl->mNetwork.end(), ' '),
                          pImpl->mNetwork.end());
    std::transform(pImpl->mNetwork.begin(), pImpl->mNetwork.end(), 
                   pImpl->mNetwork.begin(), ::toupper);
    auto work = makeHashString(pImpl->mNetwork, pImpl->mStation, pImpl->mPhase);
    pImpl->mHash = std::hash<std::string>()(work);

}

/// Gets the network
std::string NetworkStationPhase::getNetwork() const noexcept
{
    return pImpl->mNetwork;
}

/// Sets the station
void NetworkStationPhase::setStation(const std::string &s) noexcept
{
    // Copy, remove blanks, and convert to upper case
    pImpl->mStation = s;
    pImpl->mStation.erase(std::remove(pImpl->mStation.begin(),
                                      pImpl->mStation.end(), ' '),
                          pImpl->mStation.end()); 
    std::transform(pImpl->mStation.begin(), pImpl->mStation.end(),
                   pImpl->mStation.begin(), ::toupper);
    auto work = makeHashString(pImpl->mNetwork, pImpl->mStation, pImpl->mPhase);
    pImpl->mHash = std::hash<std::string>()(work);
}

/// Gets the station
std::string NetworkStationPhase::getStation() const noexcept
{
    return pImpl->mStation;
}

/// Sets the phase
void NetworkStationPhase::setPhase(const std::string &s) noexcept
{
    // Copy, remove blanks, and convert to upper case
    pImpl->mPhase = s;
    pImpl->mPhase.erase(std::remove(pImpl->mPhase.begin(),
                                    pImpl->mPhase.end(), ' '),
                        pImpl->mPhase.end());
    std::transform(pImpl->mPhase.begin(), pImpl->mPhase.end(), 
                   pImpl->mPhase.begin(), ::toupper);
    auto work = makeHashString(pImpl->mNetwork, pImpl->mStation, pImpl->mPhase);
    pImpl->mHash = std::hash<std::string>()(work);
}

/// Gets the phase name
std::string NetworkStationPhase::getPhase() const noexcept
{
    return pImpl->mPhase;
}

/// Gets the hash
size_t NetworkStationPhase::getHash() const noexcept
{
    return pImpl->mHash;
}

/// Checks for equality
bool MFLib::operator==(const NetworkStationPhase &lhs,
                       const NetworkStationPhase &rhs)
{
    return lhs.getHash() == rhs.getHash();
}

/// Checks for inequality
bool MFLib::operator!=(const NetworkStationPhase &lhs,
                       const NetworkStationPhase &rhs)
{
    return !(lhs == rhs);
}
