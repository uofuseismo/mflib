#ifndef MFLIB_SINGLECHANNEL_EVENT_HPP
#define MFLIB_SINGLECHANNEL_EVENT_HPP
#include <memory>
namespace MFLib
{
namespace SingleChannel
{
/*!
 * @brief This is a container for an event.  An event is comprised of arrivals.
 */
class Event
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor. 
     */
    Event();
    /*!
     * @brief Copy constructor.
     * @param[in] event  The event class from which to initialize this class.
     */
    Event(const Event &event);
    /*!
     * @brief Move constructor.
     * @param[in,out] event  The event class from which to initialize this
     *                       class.  On exit, event's behavior is undefined.
     */
    Event(Event &&event) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] event  The event class to copy to this.
     * @result A deep copy of event.
     */
    Event& operator=(const Event &event);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] event  The event class whose memory will be moved to this.
     *                       On exit, event's behavior is undefined.
     * @result The memory from event moved to this.
     */
    Event& operator=(Event &&event) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor.
     */
    ~Event();
    /*!
     * @brief Releases memory and resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Event Identifier
     * @{
     */
    /*!
     * @brief Sets the event identifier.
     * @param[in] evid   The event identifier.
     */
    void setIdentifier(uint64_t evid) noexcept;
    /*!
     * @brief Gets the event identifier.
     * @result The event identifier.
     * @throws std::runtime_error if the event identifier was not set.
     * @sa \c haveIdentifier()
     */
    uint64_t getIdentifier() const;
    /*!
     * @brief Determines if the event identifier was set.
     * @result True indicates that the event identifier was set.
     */
    bool haveIdentifier() const noexcept;
    /*! @} */

    /*! @name Magnitude
     * @{
     */
    /*!
     * @brief Sets the magnitude of the event.
     * @param[in] magnitude  The event magnitude.
     */
    void setMagnitude(double magnitude) noexcept;
    /*!
     * @brief Gets the magnitude of the event.
     * @result The magnitude of the event.
     * @throws std::runtime_error if the magnitude was not set.
     * @sa \c haveMagnitude()
     */
    double getMagnitude() const;
    /*!
     * @brief Determines if the magnitude was set.
     * @result True indicates that the magnitude was set.
     */
    bool haveMagnitude() const noexcept;
    /*! @} */

    /*! @name Location
     * @{
     */
    /*!
     * @brief Sets the event latitude.
     * @param[in] latitude  The event's latitude in degrees.  This must be
     *                      in the range [-90,90].
     * @throws std::invalid_argument if the latitude is out of bounds.
     */
    void setLatitude(double latitude);
    /*!
     * @brief Gets the event latitude.
     * @result The event latitude in degrees.
     * @throws std::runtime_error if this was not set.
     * @sa \c haveLatitude()
     */
    double getLatitude() const;
    /*!
     * @brief Determines if the latitude was set.
     * @result True indicates that the latitude was set.
     */
    bool haveLatitude() const noexcept;

    /*!
     * @brief Sets the event longitude.
     * @param[in] longitude  The event's longitude in degrees.  This must be
     *                       in the range (-540,540).
     */
    void setLongitude(double longitude);
    /*!
     * @brief Gets the event longitude. 
     * @result The event longitude in degrees.  This will be in the
     *         range [0,360].
     * @throws std::runtime_error if this was not set.
     * @sa \c haveLongitude()
     */
    double getLongitude() const;
    /*!
     * @brief Determines if the longitude was set.
     * @result True indicates that the longitude was set.
     */
    bool haveLongitude() const noexcept;

    /*!
     * @brief Sets the event depth. 
     * @param[in] depth   The event depth in kilometers.
     */
    void setDepth(double depth);
    /*!
     * @brief Gets the event depth.
     * @result The event depth in kilometers.
     * @throws std::runtime_error if this was not set.
     * @sa \c haveDepth()
     */
    double getDepth() const;
    /*!
     * @brief Determines if the depth was set.
     * @result True indicates that the depth was set.
     */
    bool haveDepth() const noexcept;
    /*! @} */

    /*! @name Location Quality
     * @{
     */
    /*!
     * @brief Sets the solution's root mean squared error.
     * @param[in] rms  The event location's RMS in seconds.
     * @throws std::invalid_argument if this is not positive.
     */
    void setRootMeanSquaredError(double rms);
    /*!
     * @brief Gets the solution's root mean squared error.
     * @result The event location's RMS in seconds.
     */
    double getRootMeanSquaredError() const noexcept;

    /*!
     * @brief Sets the horizontal error.
     * @param[in] erh    The horizontal error in kilometers.
     * @throws std::invalid_argument if the error is negative.
     */
    void setHorizontalError(double erh);
    /*!
     * @brief Gets the horizontal error.
     * @result The horizontal error in kilometers.
     */
    double getHorizontalError() const noexcept;

    /*!
     * @brief Sets the vertical error.
     * @param[in] erz  The vertical error in kilometers.
     * @throws std::invalid_argument if this is not positive.
     */
    void setVerticalError(double erz);
    /*!
     * @brief Gets the vertical error.
     * @result The vertical error in kilometers.
     */
    double getVerticalError() const noexcept; 
    /*! @} */
private:
    class EventImpl;
    std::unique_ptr<EventImpl> pImpl;
};
}
}
#endif
