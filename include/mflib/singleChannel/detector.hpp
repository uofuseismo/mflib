#ifndef MFLIB_SINGLECHANNEL_DETECTOR_HPP
#define MFLIB_SINGLECHANNEL_DETECTOR_HPP
#include <memory>
#include <vector>
#include "mflib/singleChannel/detection.hpp"
namespace MFLib
{
namespace SingleChannel
{
template<class T> class MatchedFilter;
template<class T> class Detection;
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
    void clear() noexcept;
    /*!
     * @brief Initializes the detector class.
     * @param[in] parameters  The detector parameters.
     * @param[in] mf          The class containing the matched-filtered signals.
     *                        This must be initialized.
     * @throws std::invalid_argument if the parameters are invalid or the
     *         matched filtering class is not initialized.
     */ 
    void initialize(const MFLib::SingleChannel::DetectorParameters &parameters,
                    const MFLib::SingleChannel::MatchedFilter<T> &mf);
    /*!
     * @brief Converts the correlation coefficients from the matched filter
     *        to detections.
     */
    void detect(const MFLib::SingleChannel::MatchedFilter<T> &xc);
    /*!
     * @brief Gets the number of detections.
     */
    int getNumberOfDetections() const noexcept;
    /*!
     * @brief Determines if the class is initialized.
     * @result True indicates that the class is initialized.
     */ 
    bool isInitialized() const noexcept;

    //std::iterator<MFLib::SingleChannel::Detection> begin();
    //std::iterator end();

    /*!
     * @brief Gets a reference to the first detection.
     */
    //const MFLib::SingleChannel::Detection<T>& begin() const;
    /*!
     * @brief Gets a reference to the last detection.
     */ 
    //const MFLib::SingleChannel::Detection<T>& end() const;
    const MFLib::SingleChannel::Detection<T>& operator[](size_t i) const;
private:
    class DetectorImpl;
    std::unique_ptr<DetectorImpl> pImpl;
};
}
}
#endif
