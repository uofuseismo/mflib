#ifndef MFLIB_ENUMS_HPP
#define MFLIB_ENUMS_HPP 1
namespace MFLib
{
/*!
 * @brief Defines the matched filter implementation.
 */
enum class MatchedFilterImplementation
{
    DIRECT,       /*!< This uses a time-domain filtering strategy. */
    OVERLAP_ADD,  /*!< This uses a frequency domain overlap and add
                       algorithm. */
    AUTO          /*!< This lets the computer decide. */
};

/*!
 * @brief Distinguishes between two types of detection strategies.
 *        The first strategy aims to apply many templates to a single
 *        channel and is necessarily useful for single-station detection.
 *        The other strategy will apply a template to each corresponding
 *        waveform in a network and is necessarily useful for network-based
 *        event detection.
 */ 
enum class MatchedFilterDetectionMode
{
    SINGLE_CHANNEL, /*!< This mode will apply many templates to a single 
                         channel.  This would be useful for single
                         station detection scenarios. */
    MULTI_CHANNEL   /*!< This mode will apply x templates to x corresponding
                         channels.  This would be useful for event
                         detection within a network. */ 
};

/*!
 * @brief Defines the relative magnitude computation strategy.
 */
enum class RelativeMagnitudeType
{
    GIBBONS_RINGDAL_2006 = 0, /*!< This finds the optimal scalar to match the 
                                   amplitude of the template waveform to the
                                   detected waveform. */
    SCHAFF_RICHARDS_2014      /*!< This finds a scalar to match the amplitude
                                  of the template waveform to the detected
                                  waveform whilst including a correction for
                                  the fact that the template and detected event
                                  have a noise contribution. */
};

/*!
 * @brief Defines whether we are seeking detections corresponding to maxima or
 *        absolute maxima of the matched filtered time series.
 */ 
enum class MaximumMatchedFilterPolicy
{
    MAXIMUM = 0,      /*!< Selects the optima of the matched filtered
                           signals by locating the maximum values. */
    ABSOLUTE_MAXIMUM  /*!< Selects the optima of the matched filtered
                           signals by first computing the absolute value
                           and then locating the maximum values. */
};

/*!
 * @brief Defines whether we make detections on all correlograms or a 
 *        reduced version of the correlograms.
 */
enum class DetectionWaveformPolicy
{
    REDUCED_CORRELOGRAMS = 0, /*!< Combines all correlograms by taking
                                   the maximum or absolute maximum
                                   of the correlograms then computing
                                   the detections. */
    SINGLE_CORRELOGRAMS       /*!< Makes detections on individual 
                                   correlograms.  There may however
                                   be multiple detections from 
                                   different templates in a small
                                   time window.  In this case, the
                                   associator must remove these
                                   duplicates. */
};

/*!
 * @brief Defines the interpolation policy for refining detections.
 */
enum class DetectionInterpolationPolicy
{
    LANCZOS = 0,   /*! Performs sinc interpolation (using a Lanczos window) on
                       fine interval then selects a maxima.  This appears to
                       be more stable but more computationally expensive. */
    QUADRATIC,     /*! Fits a quadratic and then optimizes for a delay time
                       which optimizes this quadratic.  This is computationally
                       trivial but may be suboptimal for coarse sampling
                       periods.  */
    NONE           /*!< Do not interpolate the detection time.  This would be
                        appropriate when the template is very long and the 
                        detection is probably more akin to an S + surface wave
                        detection. */
};

}
#endif
