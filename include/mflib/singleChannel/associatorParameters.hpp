#ifndef MFLIB_SINGLECHANNEL_ASSOCIATORPARAMETERS_HPP
#define MFLIB_SINGLECHANNEL_ASSOCIATORPARAMETERS_HPP
#include <memory>
namespace MFLib
{
namespace SingleChannel
{
/*!
 * @name AssociatorParameters "associatorParameters.hpp" "mflib/singleChannel/associatorParameters.hpp"
 * @brief Defines the parameters for the single channel associator.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
class AssociatorParameters
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor.
     */
    AssociatorParameters();
    /*!
     * @brief Copy constructor.
     * @param[in] parameters  The parameters class from which to initialize
     *                        this class.
     */
    AssociatorParameters(const AssociatorParameters &parameters);
    /*!
     * @brief Move constructor.
     * @param[in,out] parameters  The parameters class from which to initialize
     *                            this class.  On exit, parameters's behavior
     *                            is undefined.
     */
    AssociatorParameters(AssociatorParameters &&parameters) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] parameters  The parameters class to copy to this.
     * @result A deep copy of the parameters class.
     */
    AssociatorParameters& operator=(const AssociatorParameters &parameters);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] parameters  The parameters class whose memory will be
     *                            moved to this.  On exit, parameters's beahvior
     *                            is undefined.
     * @result The memory from parameters moved to this.
     */
    AssociatorParameters& operator=(AssociatorParameters &&parameters) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    ~AssociatorParameters();
    /*!
     * @brief Releases any memory and resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*!
     * @brief After the picks are migrated to origin time space this defines
     *        the maximum tolerance that the picks must be within to 
     *        the origin time be associated to an event.
     * @param[in] tolerance  The tolerance in seconds.  This must be
     *                       non-negative. 
     * @throws std::invalid_argument if this is negative.
     */
    void setOriginTimeTolerance(double tolerance);
    /*!
     * @brief Gets the origin time tolerance.
     * @result The origin time tolerance in seconds.
     */
    double getOriginTimeTolerance() const noexcept;
private:
    class AssociatorParametersImpl;
    std::unique_ptr<AssociatorParametersImpl> pImpl;
};
}
}
#endif
