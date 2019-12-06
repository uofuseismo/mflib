#ifndef MFLIB_MATCHEDFILTEROPTIONS_HPP
#define MFLIB_MATCHEDFILTEROPTIONS_HPP
#include <memory>
#include <vector>
#include "waveformTemplate.hpp"
#include "enums.hpp"
namespace MFLib
{
class WaveformTemplate;
/*!
 * @brief Defines the options for the matched filtering.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
class MatchedFilterOptions
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor
     */
    MatchedFilterOptions();
    /*!
     * @brief Copy constructor.
     * @param[in] options  The options from which to initialize this class.
     */
    MatchedFilterOptions(const MatchedFilterOptions &options);
    /*!
     * @brief Move constructor.
     * @param[in,out] options  The options from which to initialize this class.
     *                         On exit, options' behavior is undefined.
     */
    MatchedFilterOptions(MatchedFilterOptions &&options) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] options  The options class to copy.
     * @result A deep copy of the options.
     */ 
    MatchedFilterOptions& operator=(const MatchedFilterOptions &options);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] options  The options class whose memory will be moved 
     *                         to this.  On exit, options' behavior is
     *                         undefined.
     * @result The memory from options moved to this.
     */
    MatchedFilterOptions& operator=(MatchedFilterOptions &&options) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor
     */
    ~MatchedFilterOptions();
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
    void addTemplate(const WaveformTemplate &waveformTemplate);
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
     * @brief Gets the longest template length.  For example, if given a
     *        template length 300 samples and 400 samples, this will return
     *        400.
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
     * @throws std::invalid_argument if fftLength < \c getTemplateLength().
     * @sa \c getTemplateLength(), \c getNumberOfTemplates()
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
     * @brief When stacking - stack the absolute value of the matched filters
     *        or the raw values.
     * @param[in] labs   If true then stack the absolute values of the matched
     *                   filtered signals.
     */
    void setStackAbsoluteValues(bool labs) noexcept;
    /*!
     * @brief Determines whether or not to stack the absolute values or raw
     *        values of the matched filtered signals.
     * @result True indicates the absolute values of the matched filtered
     *         signals are to be stacked.
     */
    bool getStackAbsoluteValues() const noexcept;

    /*!
     * @brief Determines if the matched filtering options are well-defined.
     * @result True indicates that this is a valid matched filtering options
     *         class.
     */    
    bool isValid() const noexcept;   
private:
    class MatchedFilterOptionsImpl;
    std::unique_ptr<MatchedFilterOptionsImpl> pImpl;
};

/*!
 * @brief Computes the optimal FFT and block length in the overlap and add 
 *        method.
 * @param[in] nb   The filter length.  This must be positive.
 * @param[in] nx   The signal length.  This must exceed the filter length.
 * @result The optimal FFT length is given by result.first and the optimal
 *         block length is given by result.second. 
 * @note This will make the FFT length a power of 2.  There may exist more 
 *       efficient FFT lengths if using a Cooley-Tukey FFT algorithm.
 */
std::pair<int, int> computeOptimalFFTAndBlockLength(int nb, int nx);
}
#endif
