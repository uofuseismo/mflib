#ifndef MFLIB_SINGLECHANNEL_RELATIVEMAGNITUDE_HPP
#define MFLIB_SINGLECHANNEL_RELATIVEMAGNITUDE_HPP
#include <memory>
#include "mflib/enums.hpp"
namespace MFLib
{
class WaveformTemplate;
}
namespace MFLib
{
namespace SingleChannel
{
/*!
 * @name RelativeMagnitude relativeMagnitude.hpp "mflib/singleChannel/relativeMagnitude.hpp"
 * @brief Given a template and observed waveform snippet, this computes
 *        a relative magnitude.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
template<class T = double>
class RelativeMagnitude
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor.
     */
    RelativeMagnitude();
    /*!
     * @brief Copy constructor.
     * @param[in] magnitude  The magnitude class from which to initialize this
     *                       class.
     */
    RelativeMagnitude(const RelativeMagnitude &magnitude);
    /*!
     * @brief Move constructor.
     * @param[in,out] magnitude  The magnitude class from which to initialize
     *                           this class.  On exit, magnitude's behavior is
     *                           undefined.
     */
    RelativeMagnitude(RelativeMagnitude &&magnitude) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*! 
     * @brief Copy assignment operator.
     * @param[in] magnitude   The relative magnitude class to copy.
     * @result A deep copy of the magnitude class.
     */
    RelativeMagnitude& operator=(const RelativeMagnitude &magnitude);
    /*! 
     * @brief Move assignment operator.
     * @param[in,out] magnitude  The relative magnitude class whose memory will
     *                           be moved to this.  On exit, mag's behavior
     *                           is undefined.
     * @result The memory from magnitude moved to this.
     */
    RelativeMagnitude& operator=(RelativeMagnitude &&magnitude) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*! 
     * @brief Destructor.
     */
    ~RelativeMagnitude();
    /*!
     * @brief Clears the memory and resets variables on the class.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Initialization and Properties
     * @{
     */
    /*!
     * @brief Initializes the class.
     * @param[in] wt   The waveform template.  
     * @throws std::invalid_argument if there is no signal set on the
     *         waveform template or that signal is constant or all zero.
     */
    void initialize(const WaveformTemplate &wt);
    /*!
     * @brief Gets the expected length of the detected signal.
     * @result The expected length of the detected signal.
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c isInitialized()
     */
    int getDetectedSignalLength() const; 
    /*! 
     * @brief Determines if the class is initialized.
     */
    bool isInitialized() const noexcept;
    /*! @} */

    /*! @name Detected Signal
     * @{
     */
    /*! 
     * @brief Sets the detected (slave) waveform.
     * @param[in] n   The length of the detected waveform.  This must
     *                match \c getDetectedSignalLength().
     * @param[in] x   The detected (slave) waveform.  This is an array
     *                whose dimension is [n].
     * @throws std::runtime_error if the class is not initialized.
     * @throws std::invalid_argument if x is NULL or constant or n does not
     *         equal \c getDetectedSignalLength().
     * @sa \c isInitialized(), \c getDetectedSignalLength()
     */
    void setDetectedSignal(int n, const T y[]);
    /*! 
     * @brief Determines whether or not the detected signal has been set.
     * @result True indicates that the detected signal was set.
     */
    bool haveDetectedSignal() const noexcept;
    /*! @} */

    /*! @name Computations
     * @{
     */
    /*!
     * @brief Computes the amplitude scaling factor to, \f$ \alpha \f$,
     *        in \f$ \textbf{y} = \alpha \textbf{x} \f$. 
     *        where $\textbf{y}$ is the detected waveform and $\textbf{x} \f$
     *        is the template waveform. 
     * @param[in] type   Defines the relative magnitude type.
     * @result The scaling factor.
     * @throws std::runtime_error if the class is not initialized or the 
     *         detected waveform has not been set.
     * @sa \c isInitialized(), \c haveDetectedSignal()
     */
    T computeAmplitudeScalingFactor(const MFLib::RelativeMagnitudeType type) const;
    /*! 
     * @brief Computes the magnitude perturbation.  For example, if using
     *        the Gibbons and Ringdal, 2006 method this will compute 
     *        \f$ \delta M = \log_{10} \alpha \f$ where $\f \alpha \f$ is the
     *        scaling factor computed by \c computeAmplitudeScalingFactor().
     * @param[in] type   Defines the relative magnitude type.
     * @result The relative magnitude perturbation.
     * @throws std::runtime_error if the class is not initialized or the 
     *         detected waveform has not been set.
     * @sa \c isInitialized(), \c haveDetectedSignal()
     */
    T computeMagnitudePerturbation(const MFLib::RelativeMagnitudeType type) const;
    /*! @} */
private:
    class RelativeMagnitudeImpl;
    mutable std::unique_ptr<RelativeMagnitudeImpl> pImpl;
};
}
}
#endif
