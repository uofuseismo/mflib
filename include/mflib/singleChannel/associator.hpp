#ifndef MFLIB_SINGLECHANNEL_ASSOCIATOR_HPP
#define MFLIB_SINGLECHANNEL_ASSOCIATOR_HPP
#include <memory>
namespace MFLib
{
namespace SingleChannel
{
class AssociatorParameters;
/*!
 * @class Associator associator.hpp "mflib/singleChannel/associator.hpp"
 * @brief This class combines single channel matched filter detections into
 *        events.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
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

    /*!
     * @brief Associates the input detections into events.
     */
    void associate();
    int getNumberOfEvents() const noexcept;
private:
    class AssociatorImpl;
    std::unique_ptr<AssociatorImpl> pImpl;
};
}
}
#endif
