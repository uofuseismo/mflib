#ifndef MFLIB_MATCHEDFILTER_HPP
#define MFLIB_MATCHEDFILTER_HPP
#include <memory>
#include <complex>
#include <vector>
namespace MFLib
{
class MatchedFilterOptions;
/*!
 * @brief The workhorse algorithm that applies the matched filters.
 * @copyright Ben Baker distributed under the MIT license.
 */
template<class T>
class MatchedFilter
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor.
     */
    MatchedFilter();
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor.
     */
    ~MatchedFilter();
    /*!
     * @brief Resets the class and releases the memory.
     */
    void clear() noexcept;
    /*! @} */

    /*!
     * @brief Initializes the matched filter.
     * @throws std::invalid_argument if options is not valid.
     */
    void initialize(const MatchedFilterOptions &options);

    /*!
     * @brief Sets the signal corresponding to the it'th template.
     * @param[in] it        The signal corresponds to the it'th template.
     * @param[in] nSamples  The number of samples in the signal.
     * @param[in] signal    The signal corresponding to the it'th template
     *                      to set.  This is an array whose dimension is 
     *                      [nSamples].
     * @throws std::runtime_error if the class is not initialized.
     * @throws std::invalid_argument if it is not in the range
     *         [0,\c getNumberOfTemplates()-1], nSamples does not
     *         match \c getSignalLength(), or signal is NULL.
     */
    void setSignal(int it, int nSamples, const T signal[]);
    /*!
     * @brief For the case when the signal corresponding to the it'th
     *        template corresponds to a dead trace this is a convenience
     *        function that allows its values to be set to 0.
     * @param[in] it   The template index.
     * @throws std::runtime_error if the class is not initialized.
     * @throws std::invalid_argument if it is not in the range
     *         [0,\c getNumberOfTemplates()-1].
     */
    void zeroSignal(int it);
    /*!
     * @brief Gets the number of templates.
     * @result The number of templates.
     * @throws std::runtime_error if the class is not initialized.
     */
    int getNumberOfTemplates() const; 
    /*! 
     * @brief Gets the the length of the input signal.
     * @result The number of samples expected in the input signals.
     * @throws std::runtime_error if the class is not initialized.
     */
    int getSignalLength() const;
    /*!
     * @brief Gets the length of the matched filtered signals.
     * @result The length of the matched filtered signals.
     * @throws std::runtime_error if the class is not initialized.
     */
    int getFilteredSignalLength() const;

    /*!
     * @brief Applies the matched filtering to the input signal.
     * @sa \c isInitialized()
     */
    void apply();

    /*!
     * @brief Gets the matched filtered signal corresponding to the it'th
     *        template.
     * @param[in] it   The template index.
     * @result The matched filtered signal corresponding to the it'th template.
     * @throws std::runtime_error if the class is not initialized, or the
     *         matched filtered signals have not been computed.
     * @throws std::invalid_argument if it is not in the range
     *         [0,\c getNumberOfTemplates()-1].
     * @sa \c haveMatchedFilteredSignals()
     */
    std::vector<T> getMatchedFilteredSignal(int it) const;
    /*!
     * @brief Gets a pointer to the matched filtered signal corresponding to
     *        the it'th template.  Notice, this discards the first and last
     *        \c getTemplateLength() samples.
     * @param[in] it   The template index.
     * @result A pointer to the it'th matched filtered signal.  This is an
     *         array whose dimension is [\c getFilteredSignalLength() ]
     * @throws std::runtime_error if the class is not initialized, or the
     *         matched filtered signals have not been computed.
     * @throws std::invalid_argument if it is not in the range
     *         [0,\c getNumberOfTemplates()- \c getTemplateLength() + 1].
     * @sa \c haveMatchedFilteredSignals()
     */  
    const T* getMatchedFilterSignalPointer(int it) const;
    /*!
     * @brief Checks if the matched filtering has been applied to this batch
     *        of signals.
     * @result True indicates that the matched filtering has been applied.
     */
    bool haveMatchedFilteredSignals() const noexcept;

    /*!
     * @brief Determines if the class is initialized.
     * @result True indicates that the class is initialized.
     */
    bool isInitialized() const noexcept;


    /*! @name Debugging Routines
     * @{
     */
    /*!
     * @brief Gets the frequency domain spectrum of the it'th template.
     * @note This is a debugging routine.
     * @param[in] it   The template index.  
     * @result The spectra of the it'th template.
     * @throws std::invalid_error if it is less than 0 or greater than or equal
     *         to the number of templates.
     * @throws std::runtime_error if the class is not initialized.
     */
    std::vector<std::complex<T>> getSpectraOfTemplate(int it) const;
    /*! @} */
private:
    class MatchedFilterImpl;
    std::unique_ptr<MatchedFilterImpl> pImpl;
};
}
#endif
