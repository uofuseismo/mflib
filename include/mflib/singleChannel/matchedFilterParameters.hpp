#ifndef MFLIB_SINGLECHANNEL_FILTERPARAMETERS_HPP
#define MFLIB_SINGLECHANNEL_FILTERPARAMETERS_HPP
#include <memory>
#include <vector>
#include "mflib/matchedFilterParameters.hpp"
#include "mflib/waveformTemplate.hpp"
#include "mflib/enums.hpp"
namespace MFLib
{
namespace SingleChannel
{
/*!
 * @brief Defines the parameters required to perform the matched filtering
 *        for single-station detection.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
class MatchedFilterParameters
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor
     */
    MatchedFilterParameters();
    /*!
     * @brief Copy constructor.
     * @param[in] parameters  The parameters from which to initialize this
     *                        class.
     */
    MatchedFilterParameters(const MatchedFilterParameters &parameters);
    /*!
     * @brief Move constructor.
     * @param[in,out] parameters  The parameters from which to initialize this
     *                            class.  On exit, parameters' behavior is
     *                            undefined.
     */
    MatchedFilterParameters(MatchedFilterParameters &&parameters) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] parameters  The parameters class to copy.
     * @result A deep copy of the parameters.
     */ 
    MatchedFilterParameters&
        operator=(const MatchedFilterParameters &parameters);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] parameters  The parameters class whose memory will be
     *                            moved to this.  On exit, parameters' behavior
     *                            is undefined.
     * @result The memory from parameters moved to this.
     */
    MatchedFilterParameters&
        operator=(MatchedFilterParameters &&parameters) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor
     */
    ~MatchedFilterParameters();
    /*!
     * @brief Clears all memory and resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*!
     * @brief Adds a waveform template signal.
     * @param[in] waveformTemplate  The waveform signal template.
     * @throws std::invalid_argument if nSamples is not positive,
     *         signalTemplate is NULL, or the waveform template's 
     *         sampling rate differs from other first waveform template's
     *         sampling rate.
     * @note This will invalidate the FFT length and the block size.
     */
    void addTemplate(const MFLib::WaveformTemplate &waveformTemplate);
    /*!
     * @brief Clears all templates but leaves other variables untouched.
     */
    void clearTemplates() noexcept;
    /*!
     * @brief Gets the it'th template.
     * @param[in] it   The it'th template to return.
     * @result The it'th template waveform.
     * @throws std::invalid_argument if it is less than 0 or greater than or
     *         equal to the number of templates.
     * @sa \c getNumberOfTemplates()
     */
    WaveformTemplate getTemplate(int it) const;
    /*!
     * @brief Gets the number of signal templates.
     * @result The number of signal templates in the class.
     */
    int getNumberOfTemplates() const noexcept;
    /*!
     * @brief Gets the longest template length.  For example, if given
     *        template with length 300 samples and a second template with
     *        length 400 samples, then this will return 400.
     * @result The maximum length of the output template signals.
     * @throws std::runtime_error if no templates were set.
     */
    int getMaxTemplateLength() const;

    /*!
     * @brief Determines the size of the input signal chunk to process.
     * @param[in] nx   The number of samples in the input signal chunk.
     * @result throws std::invalid_argument if nx is not positive.
     * @note This will invalidate the FFT length and the block size.
     */
    void setSignalSize(int nx);
    /*!
     * @brief Gets the size of the input signal chunk to process.
     * @result The size of the input signal chunk.
     * @throws std::runtime_error if the signal size is not set.
     */
    int getSignalSize() const;

    /*!
     * @brief Sets the FFT length.
     * @param[in] fftLength   The FFT length.
     * @throws std::runtime_error if no templates are set.
     * @throws std::invalid_argument if fftLength < \c getMaxTemplateLength().
     * @sa \c getMaxTemplateLength(), \c getNumberOfTemplates()
     */
    void setFFTLength(int fftLength);
    /*!
     * @brief Gets the FFT length in the overlap and add method.
     * @result The FFT length.
     * @throws std::runtime_error if no templates are set or the signal size is
     *         not set.
     * @sa \c getNumberOfTemplates(), \c getSignalSize()
     */
    int getFFTLength() const;
    /*!
     * @brief Gets the block length in the overlap and add method.
     * @result The block length.
     * @throws std::runtime_error if no templates are set or the signal size is
     *         not set.
     * @sa \c getNumberOfTemplates(), \c getSignalSize()
     */ 
    int getBlockLength() const;

    /*!
     * @brief Gets the matched filter detection mode.
     * @result The detection mode.
     */
    MFLib::MatchedFilterDetectionMode getDetectionMode() const noexcept;
    /*!
     * @brief Defines the matched filter implementation.
     * @param[in] implementation  The matched filter implementation.
     */
    void setMatchedFilterImplementation(
        const MatchedFilterImplementation implementation) noexcept;
    /*!
     * @brief Gets the matched filter implementation.
     * @result The matched filter implementation.
     */
    MatchedFilterImplementation getMatchedFilterImplementation() const noexcept;

    /*!
     * @brief Determines if the matched filtering parameters are well-defined.
     * @result True indicates that this is a valid matched filtering parameters
     *         class.
     */    
    bool isValid() const noexcept;   

    /*!
     * @brief Returns a copy of the matched filter parameters that are suitable
     *        for initializing the matched filtering engines.
     * @result The matched filter parameters suitable for initializing the
     *         matched filtering engine classes.
     * @throws std::runtime_error if the class is not valid.
     */
    MFLib::MatchedFilterParameters getParameters() const;
private:
    class MatchedFilterParametersImpl;
    std::unique_ptr<MatchedFilterParametersImpl> pImpl;
};
}
}
#endif
