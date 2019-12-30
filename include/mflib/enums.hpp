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

}
#endif
