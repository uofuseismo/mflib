#include <cstdio>
#include <cstdlib>
#include <string>
#include "pyNetworkStationPhase.hpp"
#include "mflib/networkStationPhase.hpp"

using namespace PBMFLib;

/// C'tor
NetworkStationPhase::NetworkStationPhase() :
    mNetworkStationPhase(std::make_unique<MFLib::NetworkStationPhase> ())
{
}

NetworkStationPhase::NetworkStationPhase(
    const NetworkStationPhase &nsp)
{
    *this = nsp;
}

NetworkStationPhase::NetworkStationPhase(
    const MFLib::NetworkStationPhase &nsp)
{
    *this = nsp;
}

/// Copy assignment
NetworkStationPhase&
NetworkStationPhase::operator=(const NetworkStationPhase &nsp)
{
    if (&nsp == this){return *this;}
    mNetworkStationPhase
        = std::make_unique<MFLib::NetworkStationPhase>
          (*nsp.mNetworkStationPhase);
    return *this;
}

NetworkStationPhase&
NetworkStationPhase::operator=(const MFLib::NetworkStationPhase &nsp)
{
    mNetworkStationPhase = std::make_unique<MFLib::NetworkStationPhase> (nsp); 
    return *this;
}

MFLib::NetworkStationPhase 
NetworkStationPhase::getNativeClass() const
{
    MFLib::NetworkStationPhase nsp(*mNetworkStationPhase);
    return nsp;
}

/// Destructor
NetworkStationPhase::~NetworkStationPhase() = default;

/// Reset the class
void NetworkStationPhase::clear() noexcept
{
    mNetworkStationPhase->clear();
}

/// Set/get network 
void NetworkStationPhase::setNetwork(const std::string &network) noexcept
{
    mNetworkStationPhase->setNetwork(network);
}

std::string NetworkStationPhase::getNetwork() const noexcept
{
    return mNetworkStationPhase->getNetwork();
}

/// Set/get station 
void NetworkStationPhase::setStation(const std::string &station) noexcept
{
    mNetworkStationPhase->setStation(station);
}

std::string NetworkStationPhase::getStation() const noexcept
{
    return mNetworkStationPhase->getStation();
}

/// Set/get phase 
void NetworkStationPhase::setPhase(const std::string &phase) noexcept
{
    mNetworkStationPhase->setPhase(phase);
}

std::string NetworkStationPhase::getPhase() const noexcept
{
    return mNetworkStationPhase->getPhase();
}

/// Get the hash
size_t NetworkStationPhase::getHash() const noexcept
{
    return mNetworkStationPhase->getHash();
}

void PBMFLib::initializeNetworkStationPhase(pybind11::module &m)
{
    pybind11::class_<PBMFLib::NetworkStationPhase>
        nsp(m, "NetworkStationPhase");
    /// The default constructor:
    nsp.def(pybind11::init<> ());
    nsp.doc() = "Defines a network, station, and phase to be used in identifying a template.";
    /// Network name
    nsp.def("set_network",
            &PBMFLib::NetworkStationPhase::setNetwork,
            "Sets the network name.");
    nsp.def("get_network",
            &PBMFLib::NetworkStationPhase::getNetwork,
            "Gets the network name.");
    /// Station name
    nsp.def("set_station",
            &PBMFLib::NetworkStationPhase::setStation,
            "Sets the station name.");
    nsp.def("get_station",
            &PBMFLib::NetworkStationPhase::getStation,
            "Gets the station name."); 
    /// Phase name
    nsp.def("set_phase",
            &PBMFLib::NetworkStationPhase::setPhase,
            "Sets the phase name.");
    nsp.def("get_phase",  
            &PBMFLib::NetworkStationPhase::getPhase,
            "Gets the phase name.");
    /// Clears the class memory and resets
    nsp.def("clear",
            &PBMFLib::NetworkStationPhase::clear,
            "Clears the class's memory and resets the class.");
}
