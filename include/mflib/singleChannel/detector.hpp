#ifndef MFLIB_SINGLECHANNEL_DETECTOR_HPP
#define MFLIB_SINGLECHANNEL_DETECTOR_HPP
#include <memory>
namespace MFLib
{
namespace SingleChannel
{
template<class T> class MatchedFilter;
class DetectorParameters;
}
}

namespace MFLib
{
namespace SingleChannel
{
/*!
 * @class Detector "detector.hpp" "mflib/singleChannel/detector.hpp"
 * @brief This class computes detections from time series of correlation
 *        coefficients.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
template<class T = double>
class Detector
{
public:
    /*!
     * @brief Constructor.
     */
    Detector();
    /*!
     * @brief Destructor.
     */
    ~Detector();
    /*!
     * @brief Initializes the detector class.
     */ 
//    void initialize(std::shared_ptr<const MFLib::SingleChannel::MatchedFilter<T>> xc);
    /*!
     * @brief Converts the correlation coefficients from the matched filter
     *        to detections.
     */
    void detect(const MatchedFilter<T> &xc);
    /*!
     * @brief Gets the number of detections.
     */
    int getNumberOfDetections() const noexcept;
    /*!
     * @brief Determines if the class is initialized.
     * @result True indicates that the class is initialized.
     */ 
    bool isInitialized() const noexcept;
private:
    class DetectorImpl;
    std::unique_ptr<DetectorImpl> pImpl;
};
}
}
#endif
