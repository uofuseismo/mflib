#ifndef MFLIB_SINGLECHANNEL_DETECTIONTIMEINTERPOLATIONPARAMETERS_HPP
#define MFLIB_SINGLECHANNEL_DETECTIONTIMEINTERPOLATIONPARAMETERS_HPP
#include <memory>
#include "mflib/enums.hpp"
namespace MFLib
{
namespace SingleChannel
{
/*!
 * @class DetectionTimeInterpolationParameters "detectionTimeInterpolationParameters.hpp" "mflib/singleChannel/detectionTimeInterpolationParameters.hpp"
 * @brief Defines the parameters for interpolating the detection time.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
class DetectionTimeInterpolationParameters
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     */
    DetectionTimeInterpolationParameters();
    /*!
     * @brief Copy constructor
     * @param[in] parameters  The detection time interpolation parameters from
     *                        which to initialize this class.
     */
    DetectionTimeInterpolationParameters(
        const DetectionTimeInterpolationParameters &parameters);
    /*!
     * @brief Move constructor
     * @param[in,out] parameters  The detection time interpolation parameters
     *                            from which to initialization this class. 
     *                            On exit, parameters's behavior is undefined.
     */
    DetectionTimeInterpolationParameters(
        DetectionTimeInterpolationParameters &&parameters) noexcept;
    /*! @* */

    /*! @name Destructors
     * @{
     */
    ~DetectionTimeInterpolationParameters();
    /*!
     * @brief Resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] parameters  The detection time interpolation parameters to
     *                        copy to this.
     * @result A deep copy of the detection time parameters.
     */
    DetectionTimeInterpolationParameters&
        operator=(const DetectionTimeInterpolationParameters &parameters);
    /*!
     * @brief Move assignment operators.
     * @param[in,out] parameters  The detection time interpolation parameters 
     *                            whose memory will be moved to this.  On exit,
     *                            parameters's behavior is undefined.
     */
    DetectionTimeInterpolationParameters&
        operator=(DetectionTimeInterpolationParameters &&parameters) noexcept;
    /*! @} */

    /*! @brief Lanczos (Truncated sinc) Interpolation
     * @{
     */
    /*!
     * @brief Enbables Lanczos interpolation.  This is a slow but accurate
     *        interpolation strategy that will grid-search for a refined
     *        detection time.
     */
    void enableLanczos() noexcept;
    /*!
     * @brief Defines the size of the kernel: 
     *        \f$ L(x) = \sinc (x) \sinc \left ( \frac{x}{\alpha} \right ) \f$
     *        for $x \in (-\alpha,\alpha)$.  There will be $2 \alpha - 1$ side
     *        lobes.  As this number is increased the interpolation improves
     *        at the cost of runtime.
     * @param[in] alpha  The width of the window in samples on either side of x.
     * @throws std::invalid_argument if this is not positive.
     */
    void setLanczosAlpha(int alpha);
    /*! @brief
     * @Gets the \f$ \alpha \f$ for the Lanczos interpolation.
     */
    int getLanczosAlpha() const noexcept; 
    /*!
     * @brief Defines the number of interpolation points over the interval
     *        \$
     *           \left (x_{opt}-\frac{dt}/{2}, x_{opt}+\frac{dt}{2} \right )
     *        \$.
     *        The algorithm will interpolate the correlation value at each
     *        point on this interval and select the maximum correlation value.
     * @param[in] npts   The number of interpolation points.  For example, if
     *                   this is $10$ then the result is similar to resampling
     *                   the correlogram \f$ 10 \f$ times.
     * @throws std::invalid_argument if npts is not positive.
     */
    void setLanczosNumberOfInterpolationPoints(int npts);
    /*!
     * @brief Gets the number of samples at which the Lanczos interolate will
     *        estimate a cross-correlation value.
     * @result The number of interpolation samples.
     */
    int getLanczosNumberOfInterpolationPoints() const noexcept;
    /*! @} */

    /*! @name Polynomial Interpolation
     * @{
     */ 
    /*!
     * @brief Enables quadrtic polynomial interpolation.  This is a fast and
     *        fairly accurate interpolation strategy.
     */
    void enableQuadratic() noexcept;
    /*! @} */

    /*! @name No Interpolation
     * @{
     */
    /*!
     * @brief Disables the interpolation of detection times.
     *        In this case, the interolated detection time will equal the
     *        detection time.
     */
    void disable() noexcept;
    /*! @ */

    /*!
     * @brief Gets the detection interpolation type.
     * @brief The type of detection time interpolation to perform.
     */
    DetectionTimeInterpolationType getType() const noexcept; 
private:
    class InterpParametersImpl;
    std::unique_ptr<InterpParametersImpl> pImpl;
};
}
}
#endif
