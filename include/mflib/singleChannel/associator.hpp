#ifndef MFLIB_SINGLECHANNEL_ASSOCIATOR_HPP
#define MFLIB_SINGLECHANNEL_ASSOCIATOR_HPP
#include <memory>
namespace MFLib
{
namespace SingleChannel
{
template<class T> class Detection;
class AssociatorParameters;
/*!
 * @class Associator associator.hpp "mflib/singleChannel/associator.hpp"
 * @brief This class combines single channel matched filter detections into
 *        events.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
template<class T>
class Associator
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor.
     */
    Associator();
    /*!
     * @brief Copy constructor.
     * @param[in] associator  The associator class from which to initialize
     *                        this class. 
     */
    Associator(const Associator &associator);
    /*!
     * @brief Move constructor.
     * @param[in,out] associator  The associator class from which to initialize
     *                            this class.  On exit, associator's behavior
     *                            is undefined.
     */ 
    Associator(Associator &&associator) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] associator  The associator class to copy to this.
     * @result A deep copy of associator.
     */
    Associator& operator=(const Associator &associator);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] associator  The associator class whose memory will be
     *                            moved to this.  On exit, associator's 
     *                            behavior is undefined.
     * @result The memory from associator moved to this.
     */
    Associator& operator=(Associator &&associator) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor. 
     */
    ~Associator();
    /*!
     * @brief Releases all memory and resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Initialization
     * @{
     */
    /*!
     * @brief Initializes the associator class.
     * @param[in] parameters  The parameters from which to initialize
     *                        the associator.
     */
    void initialize(const MFLib::SingleChannel::AssociatorParameters &parameters);
    /*!
     * @brief Determines if the class is initialized.
     * @result True indicates that the class is inititalized.
     */
    bool isInitialized() const noexcept;
    /*! @} */

    /*! @name Detections
     * @{
     */
    /*!
     * @brief Adds a detection. 
     */
    template<typename T>
    void addDetection(const MFLib::SingleChannel::Detection<T> &det);
    /*!
     * @brief Gets the number of detections set on the class.
     * @Result The number of detections.
     */
    int getNumberOfDetections() const noexcept;
    /*!
     * @brief Clears the detections set on the class. 
     */
    void clearDetections() noexcept;
    /*! @} */

    /*!
     * @brief Associates the input detections into events.
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c isInitialized()
     */
    void associate();

    /*!
     * @brief Gets the number of events.
     * @result The number of events created by the associator.
     * @sa \c associate()
     */
    int getNumberOfEvents() const noexcept;

    /*!
     * @brief Returns the detections belonging to the iev'th event index.
     * @param[in] iev  The event index.  This must be in the range
     *                 [0, \c getNumberOfEvents()-1].
     * @result The detections grouped into the iev'th event.
     * @sa \c getNumberOfEvents()
     */
    template<typename T>
    std::vector<MFLib::SingleChannel::Detection<T>> getDetectionsForEvent(int iev) const;
private:
    class AssociatorImpl;
    std::unique_ptr<AssociatorImpl> pImpl;
};
}
}
#endif
