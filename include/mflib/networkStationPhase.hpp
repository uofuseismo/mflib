#ifndef MFLIB_NETWORKSTATIONPHASE_HPP
#define MFLIB_NETWORKSTATIONPHASE_HPP
#include <memory>
namespace MFLib
{
class NetworkStationPhase
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor. 
     */
    NetworkStationPhase();
    /*!
     * @brief Copy constructor.
     * @param[in] nsp   The NetworkStationPhase class from which to initialize
     *                  this class.
     */
    NetworkStationPhase(const NetworkStationPhase &nsp);
    /*!
     * @brief Move constructor.
     */
    NetworkStationPhase(NetworkStationPhase &&nsp) noexcept;
    /*! @} */


    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] nsp  The NetworkStationPhase class to copy to this.
     * @result A deep copy of nsp.
     */
    NetworkStationPhase& operator=(const NetworkStationPhase &nsp);    
    /*!
     * @brief Move assignment operator.
     * @param[in,out] nsp  The NetworkStationPhase class whose memory will
     *                     be moved to this.  On exit, nsp's behavior is
     *                     undefined.
     * @result The memory from nsp moved to this.
     */
    NetworkStationPhase& operator=(NetworkStationPhase &&nsp) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor.
     */
    ~NetworkStationPhase();
    /*!
     * @brief Clears the memory and resets the class.
     */ 
    void clear() noexcept;
    /*! @} */

    /*! @name Network
     * @{
     */
    /*!
     * @brief Sets the network name.
     * @param[in] network  The network name.
     */    
    void setNetwork(const std::string &network) noexcept;
    /*!
     * @brief Gets the network name.
     * @result The network name.
     */
    std::string getNetwork() const noexcept;
    /*! @} */

    /*! @name Station
     * @{
     */
    /*!
     * @brief Sets the station name.
     * @param[in] station  The station name.
     */
    void setStation(const std::string &station) noexcept;
    /*!
     * @brief Gets the station name.
     * @result The station name.
     */
    std::string getStation() const noexcept;
    /* @} */

    /*! @name Phase
     * @{
     */
    /*!
     * @brief Sets the seismic phase.
     * @param[in] phase  The name of the seismic phase.
     */
    void setPhase(const std::string &phase) noexcept; 
    /*!
     * @brief Gets the seismic phase.
     * @result The name of the seismic phase.
     */
    std::string getPhase() const noexcept;
    /* @} */
    
    /*!
     * @brief Returns the hash value of the network, station, and phase.
     * @result The hash value of the network, station, and phase.
     */
    size_t getHash() const noexcept;
private:
    class NetworkStationPhaseImpl;
    std::unique_ptr<NetworkStationPhaseImpl> pImpl;
};
bool operator==(const NetworkStationPhase &lhs, const NetworkStationPhase &rhs);
bool operator!=(const NetworkStationPhase &lhs, const NetworkStationPhase &rhs);
}
#endif
