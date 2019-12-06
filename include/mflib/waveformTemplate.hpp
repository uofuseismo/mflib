#ifndef MFLIB_WAVEFORMTEMPLATE_HPP
#define MFLIB_WAVEFORMTEMPLATE_HPP
#include <memory>
namespace MFLib
{
/*!
 * @brief Defines a waveform template.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
class WaveformTemplate
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor.
     */
    WaveformTemplate(); 
    /*!
     * @brief Copy constructor.
     * @param[in] tplate  The template class from which to initialize this
     *                    template.
     */
    WaveformTemplate(const WaveformTemplate &tplate);
    /*!
     * @brief Move constructor.
     * @param[in,out] tplate  The template class whose memory will be moved
     *                        to this.  On exit, tplate's behavior is undefined.
     */
    WaveformTemplate(WaveformTemplate &&tplate) noexcept;
    /*! } */

    /*!
     * @brief Copy assignment operator.
     * @param[in] tplate   The waveform template class to copy.
     * @result A deep copy of the inupt waveform template.
     */
    WaveformTemplate& operator=(const WaveformTemplate &tplate);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] tplate  The waveform template class whose memory will be
     *                        moved to this class.  On exit, tplate's behavior
     *                        is undefined.
     * @result The memory from tplate moved to this.
     */
    WaveformTemplate& operator=(WaveformTemplate &&tplate) noexcept;

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor.
     */
    ~WaveformTemplate();
    /*!
     * @brief Clears and resets the memory of the class.
     */ 
    void clear() noexcept;
    /*! @} */

    /*! @name Signal (Required)
     * @{
     */
    /*!
     * @brief Sets the template waveform signal.
     * @param[in] npts  The number of samples in the template.
     * @param[in] x     The template waveform signal.  This is an array whose
     *                  dimension is [npts].
     * @throws std::invalid_argument if npts is not positive or x is NULL.
     * @note This will invalidate the onset time.
     */
    void setSignal(int npts, const double x[]);
    /*! @coypdoc setSignal */
    void setSignal(int npts, const float x[]);
    /*!
     * @brief Gets the template waveform signal.
     * @param[in] maxx  The maximum number of samples that x can hold.  This
     *                  must be at least \c getSignalLength().
     * @param[out] x    The waveform template.  This is an array whose dimension
     *                  is [maxx] however only the first \c getSignalLength() 
     *                  samples are accessed.
     * @throws std::invalid_argument if x is NULL or maxx is too small.
     * @throws std::runtime_error if the signal was not set.
     */
    void getSignal(int maxx, double *x[]);
    /*! @copydoc getSignal */
    void getSignal(int maxx, float *x[]);
    /*!
     * @result The length of the template waveform signal.
     * @throws std::runtime_error if the template waveform was not set.
     * @sa \c haveSignal() 
     */
    int getSignalLength() const;
    /*!
     * @brief Determines if the template waveform signal was set.
     * @result True indicates that the template was set.
     */
    bool haveSignal() const noexcept; 
    /*! @} */

    /*! @name Sampling Rate (Required)
     * @{
     */
    /*!
     * @brief Sets the sampling rate.
     * @param[in] samplingRate  The sampling rate in Hz of the template waveform
     *                          signal.
     * @throws std::invalid_argument if this is not positive.
     * @note This will invalidate the onset time.
     */
    void setSamplingRate(double samplingRate); 
    /*!
     * @brief Gets the sampling rate.
     * @result The template's sampling rate.
     * @throws std::runtime_error if the sampling rate was not set.
     */ 
    double getSamplingRate() const;
    /*!
     * @brief Determines if the sampling rate was set.
     * @result True indicates that the sampling rate was set.
     */
    bool haveSamplingRate() const noexcept;
    /*! @} */

    /*! @name Shift and Sum Weight
     * @{
     */
    /*!
     * @brief Defines the template's weight in the shift and stack operation.
     * @param[in] weight  The weight of this template during the shift and sum
     *                    operation.
     * @throws std::invalid_argument if weight is not in the range [0,1].
     */
    void setShiftAndStackWeight(double weight);
    /*!
     * @brief Gets the template's weight during the shift and stack operation.
     * @result The template's weight during the shift and sum operation.
     * @note If \c setShiftStackAndWeight() was not called then this will
     *       be unity.
     */
    double getShiftAndStackWeight() const noexcept;
    /*! @} */

    /*! @name Onset Time (Required for Shifting)
     * @{
     */
    /*!
     * @brief Sets the time in seconds relative to the trace where the onset
     *        occurs.  
     * @param[in] onsetTime  The onset time in seconds where the pick occurs.
     *                       For example, if the pick is 2 seconds into the
     *                       the trace, i.e., there is `noise' 2 seconds prior
     *                       to the pick, then this should be 2.
     * @throws std::runtime_error if the waveform template signal was not set
     *         or the sampling rate was not set.
     * @throws std::invalid_argument if the onset time is not in the trace.
     * @sa \c haveSignal(), \c haveSamplingRate().
     * @note This is invalidated whenever the sampling period or signal is set.
     */ 
    void setOnsetTime(double onsetTime);
    /*!
     * @brief Gets the onset time.
     * @result The time, relative to the trace start, where the pick occurs.
     *         For example, if 2, then 2 seconds into the trace is where the
     *         pick was made.
     * @throws std::runtime_error if the onset time was not set.
     * @sa \c haveOnsetTime()
     */
    double getOnsetTime() const;
    /*!
     * @brief Determines if the onset time was set.
     * @result True indicates that the onset time was set.
     */
    bool haveOnsetTime() const noexcept;
    /*!@ } */

    /*! @brief Travel Time (Required for Shifting)
     * @{
     */
    /*!
     * @brief Sets the observed travel time for the pick.
     * @param[in] travelTime   The observed travel time in seconds for the pick.
     *                         For example, if this is 7, then it took 7 seconds
     *                         for the wave to travel from the source to the
     *                         receiver.
     * @throws std::invalid_argument if travelTime is negative
     * @note The trace commences travelTime - onsetTime seconds after the origin
     *       time.
     */
    void setTravelTime(double travelTime);
    /*!
     * @brief Gets the observed travel time for the pick.
     * @result The observed travel time for this pick in seconds.
     *         For example, if this is 7, then it took 7 seconds for the
     *         wave to travel from the source to the receiver.
     * @throws std::runtime_error if the travel time was not set.
     * @sa \c haveTravelTime() 
     */
    double getTravelTime() const;
    /*!
     * @brief Determines if the travel time was set set.
     */
    bool haveTravelTime() const noexcept;
    /*! @} */
private:
    class WaveformTemplateImpl;
    std::unique_ptr<WaveformTemplateImpl> pImpl;
};
}
#endif
