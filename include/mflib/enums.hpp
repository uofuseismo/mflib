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
}
#endif
