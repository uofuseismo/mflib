#ifndef MFLIB_SINGLECHANNEL_MATCHEDFILTER_HPP
#define MFLIB_SINGLECHANNEL_MATCHEDFILTER_HPP
#include <memory>
#include <complex>
#include <vector>
namespace MFLib
{
namespace SingleChannel
{
class MatchedFilterParameters;
}
}

namespace MFLib
{
namespace SingleChannel
{
/*!
 * @brief The interface to the single-channel matched filtering method.
 *        This is to be used when the user wishes to apply many templates
 *        to the same signal.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
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
     * @param[in] parameters  The parameters for the matched filtering.
     * @throws std::invalid_argument if the parameters are not valid.
     */
    void initialize(const MFLib::SingleChannel::MatchedFilterParameters &parameters);

    /*!
     * @brief Sets the signal to matched filter.
     * @param[in] nSamples  The number of samples in the signal.
     * @param[in] signal    The signalto matched filter.  This is an array
     *                      whose dimension is [nSamples].
     * @throws std::runtime_error if the class is not initialized.
     * @throws std::invalid_argument if nSamples does not match
     *         \c getSignalLength() or signal is NULL.
     */
    void setSignal(int nSamples, const T signal[]);
    /*! 
     * @brief Gets a pointer to the signal to matched filter.
     * @result A pointer to the signal to matched filter.
     *         This is an array whose dimension is [\c getSignalLength()].
     * @throws std::invalid_argument if the class is not initialized.
     */
    const T* getSignalPointer() const;

    /*!
     * @brief For the case when the signal is all zeros and the user does
     *        not wish to break the processing workflow.
     * @throws std::runtime_error if the class is not initialized.
     */
    void zeroSignal();
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
     * @brief Checks if the matched filtering has been applied to the signal.
     * @result True indicates that the matched filtering has been applied.
     */
    bool haveMatchedFilteredSignals() const noexcept;

    /*!
     * @brief Determines if the class is initialized.
     * @result True indicates that the class is initialized.
     */
    bool isInitialized() const noexcept;
private:
    MatchedFilter(const MatchedFilter &mf) = delete;
    MatchedFilter& operator=(const MatchedFilter &mf) = delete;
    class MatchedFilterImpl;
    std::unique_ptr<MatchedFilterImpl> pImpl;
};
}
}
#endif
