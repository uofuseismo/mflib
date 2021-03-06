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

    /*!
     * @brief Sets the minimum number of picks required to form an event.
     * @param[in] minPicks  The minimum number of picks required to
     *                      form an event.
     * @throws std::invalid_argument if this is not positive.
     */
    void setMinimumNumberOfPicksInEvent(int minPicks);
    /*!
     * @brief Gets the minimum number of picks required to form an event.
     * @result The minimum number of picks requied to form an event.
     */ 
    int getMinimumNumberOfPicksInEvent() const noexcept;

 
    /*!
     * @brief Allows the algorithm to use the correlation coefficients
     *        as weights in the clustering.
     */
    void enableCorrelationCoefficientWeighting() noexcept;
    /*!
     * @brief Prevents the algorithm from using correlation coefficients
     *        as weights in the clustering.  In this case all arrivals have
     *        equal weight.
     */
    void disableCorrelationCoefficientWeighting() noexcept;
    /*!
     * @brief Determines if the clustering algorithm should weight origin
     *        times based on the (absolute value) of the correlation
     *        coefficient.
     * @result True indicates that the clustering should utilize the 
     *         correlation coefficients as weights.
     */
    bool useCorrelationCoefficientWeighting() const noexcept;
private:
    class AssociatorParametersImpl;
    std::unique_ptr<AssociatorParametersImpl> pImpl;
};
}
}
#endif
