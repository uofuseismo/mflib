#ifndef MFLIB_SINGLECHANNEL_DETECTORPARAMETERS_HPP
#define MFLIB_SINGLECHANNEL_DETECTORPARAMETERS_HPP
#include <memory>
#include "mflib/enums.hpp"
namespace MFLib
{
namespace SingleChannel
{
/*!
 * @class "DetectorParameters" "detectorParameters.hpp" "mflib/singleChannel/detectorParameters.hpp"
 * @brief Defines the parameters for the single-channel detector.  This includes
 *        policies on how to interpet the matched filtered signals and
 *        what information the output detections contain.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
class DetectorParameters
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor.
     */
    DetectorParameters();
    /*!
     * @brief Copy constructor.
     * @param[in] parameters  The detector parameters from which to initialize
     *                        this class.
     */
    DetectorParameters(const DetectorParameters &parameters);
    /*!
     * @brief Move constructor.
     * @param[in,out] parameters  The detector parameters from which to
     *                            initialize this class.  On exit, parameters's
     *                            behavior is undefined.
     */
    DetectorParameters(DetectorParameters &&parameters) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] parameters   The parameters class to copy.
     * @result A deep copy of the parameters class.
     */
    DetectorParameters& operator=(const DetectorParameters &parms);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] parameters  The parameters class whose memory will be
     *                            moved to this.  On exit, the behavior of
     *                            parameters is undefined.
     * @result The memory from parameters moved to this.
     */
    DetectorParameters& operator=(DetectorParameters &&parms) noexcept;
    /*! @} */
 
    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor.
     */
    ~DetectorParameters();
    /*!
     * @brief Clears the memory and resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Maximum Policy
     * @{
     */
    /*!
     * @brief Specifies the policy which is used to compute the maxima of
     *        the matched filtered signals.
     * @param[in] policy  The policy for computing the maxima of the
      *                   matched filtered signals.
     */
    void setMaximaPolicy(MaximumMatchedFilterPolicy policy) noexcept;
    /*!
     * @brief Gets the policy for computing the maxima of the matched filtered
     *        signals.
     * @result The policy for computing the maxima of the matched
     *         filtered signals.
     */
    MaximumMatchedFilterPolicy getMaximaPolicy() const noexcept;
    /*! @} */

    /*! @name Minimum Detection Sample Spacing
     * @{
     */
    /*!
     * @brief Sets the minimum number of samples required between detections.
     * @param[in] nSamples  The minimum number of samples between detections. 
     *                      For example, if this is 5, then the nearest  
     *                      detection must be at least 6 samples from this
     *                      detection.  Setting this to zero or a negative value
     *                      disables this behavior.
     */
    void setMinimumDetectionSpacing(int nSamples) noexcept;
    /*!
     * @brief Gets the minimum number of samples required between detections. 
     * @result The minimum number of samples required between detections.
     */
    int getMinimumDetectionSpacing() const noexcept;
    /*! @} */

    /*! @name Detection Threshold
     * @{
     */
    /*!
     * @brief Sets the minimum correlation value for a candidate detection.  
     * @param[in] tol  The minimum correlation-score above which a detection
     *                 is declared.  This must be in the range (0,1].
     * @throws std::invalid_argument if this not postiive or exceeds 1.
     */
    void setDetectionThreshold(double tol);
    /*! 
     * @brief Gets the minimum correlation value of a candidate detection.
     * @result The minimum threshold above which a detection can be declared.
     */
    double getDetectionThreshold() const noexcept; 
    /*! @} */
private:
    class DetectorParametersImpl;
    std::unique_ptr<DetectorParametersImpl> pImpl;
};
}
}
#endif
