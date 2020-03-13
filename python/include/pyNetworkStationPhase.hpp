#ifndef PYMFLIB_PYNETWORKSTATIONPHASE_HPP
#define PYMFLIB_PYNETWORKSTATIONPHASE_HPP
#include <memory>
#include <pybind11/pybind11.h>

namespace MFLib
{
class NetworkStationPhase;
}

namespace PBMFLib
{
class NetworkStationPhase
{
public:
    // C'tors
    NetworkStationPhase();
    NetworkStationPhase(const NetworkStationPhase &nsp);
    NetworkStationPhase(const MFLib::NetworkStationPhase &nsp);
    // Operators
    NetworkStationPhase& operator=(const NetworkStationPhase &nsp);
    NetworkStationPhase& operator=(const MFLib::NetworkStationPhase &nsp);
    // Destructors
    ~NetworkStationPhase();
    void clear() noexcept;

    void setNetwork(const std::string &network) noexcept;
    std::string getNetwork() const noexcept;

    void setStation(const std::string &station) noexcept;
    std::string getStation() const noexcept;

    void setPhase(const std::string &phase) noexcept;
    std::string getPhase() const noexcept;
    size_t getHash() const noexcept;
    MFLib::NetworkStationPhase getNativeClass() const;
private:
    std::unique_ptr<MFLib::NetworkStationPhase> mNetworkStationPhase;
};
void initializeNetworkStationPhase(pybind11::module &m);
}
#endif
