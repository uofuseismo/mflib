#ifndef MFLIB_PEAKFINDER_HPP
#define MFLIB_PEAKFINDER_HPP
#include <memory>
namespace MFLib
{
/*!
 * @name PeakFinder "PeakFinder.hpp" "private/peakFinder.hpp"
 * @brief A utility for finding the maximal peaks in a signal whose
 *        maximal values cannot exceed unity.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
template<class T = double>
class PeakFinder
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor.
     */
    PeakFinder();
    /*!
     * @brief Copy constructor.
     * @param[in] peaks  The peaks class from which to initialize this class.
     */
    PeakFinder(const PeakFinder &peaks);
    /*!
     * @brief Move constructor.
     * @param[in,out] peaks  The peaks class from which to initialize this 
     *                       class.  On exit, peaks's behavior is undefined.
     */
    PeakFinder(PeakFinder &&peaks) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] peaks  The peaks class to copy to this.
     * @result A deep copy of peaks.
     */
    PeakFinder& operator=(const PeakFinder &peaks);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] peaks  The peaks class whose memory will be moved to
     *                       this.  On exit, peaks's behavior is undefined.
     * @result The memory of peaks moved to this.
     */
    PeakFinder& operator=(PeakFinder &&peaks) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor.
     */
    ~PeakFinder();
    /*!
     * @brief Clears the class and resets the memory.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Properties
     * @{
     */
    /*!
     * @brief Sets the minimum distance in samples between candidate
     *        peaks.
     * @param[in] minDistance  The minimum distance in samples that peaks
     *                         must be separated.  If this is 0 or negative
     *                         then this screening is deactiviated.
     * @note This will invalidate any existing computed peaks.
     */ 
    void setMinimumPeakDistance(int minDistance) noexcept;
    /*!
     * @brief Gets the minimum distance in samples between candidate
     *        peaks.
     * @result The minimum distance in samples between peaks.
     */
    int getMinimumPeakDistance() const noexcept; 

    /*!
     * @brief Sets minimal value that a peak can be.
     * @param[in] threshold  The minimum peak threshold.
     *                       This must be in the range [0,1].
     * @throws std::invalid_argument if threshold is not in a valid range.
     * @note This will invalidate any existing computed peaks.
     */
    void setThreshold(T threshhold);
    /*!
     * @brief Gets the threshold above which a peak must exist.
     * @result The threshold.
     */
    T getThreshold() const noexcept;
    /*! @} */

    /*! @brief Usage after initialization
     * @{
     */
    /*!
     * @brief Sets the signal of which to find the peak.
     * @param[in] n    The number of samples in x.  This must be positive.
     * @param[in] x    The signal from which the peaks will be computed.
     * @throws std::invalid_argument if x is NULL or n is not positive.
     */
    void setSignal(const int n, const T x[]);
    /*!
     * @brief Computes the peak finding algorithm.
     */
    void apply();
    /*!
     * @brief Gets the number of peaks.
     * @result The number of peaks in the signal.
     */
    int getNumberOfPeaks() const noexcept;
    /*!
     * @brief Gets the peaks.
     * @param[in] maxSpace       The max number of elements in peaks.  This must
     *                           be at least \c getNumberOfPeaks().
     * @param[out] peakIndices   The indices of the peaks in the signal.
     */
    void getPeaks(int maxSpace, int *peaksIndices[]) const;
    /*! @} */
private:
    class PeakFinderImpl;
    std::unique_ptr<PeakFinderImpl> pImpl;
};
}
#endif
