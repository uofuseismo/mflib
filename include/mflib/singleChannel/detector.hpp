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
     */ 
    void initialize(const MFLib::SingleChannel::DetectorParameters &parameters);

    /*!
     * @brief Gets the number of templates.
     * @returns The number of templates.
     */
    int getNumberOfTemplates() const noexcept;
    /*!
     * @brief Returns the length of the matched filtered signal.
     */
    int getMatchedFilteredSignalLength() const noexcept;

    /*!
     * @brief This will extract the matched filtered signals from the matched
     *        filtering class and reduce these matched filtered signals into one
     *        signal.  The reduction policy is defined in the parameters.
     *        From this reduced signal the detections will be created.
     * @param[in] mf    The matched filtering class.
     * @throws std::invalid_argument if the matched filtered signals are not
     *         yet computed or the number of templates in the matched filtering
     *         class does not match the \c getNumberOfTemplates(). 
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c getNumberOfTemplates(), \c isInitialized()
     */
    void setMatchedFilteredSignals(const MFLib::SingleChannel::MatchedFilter<T> &mf);
    /*!
     * @brief Converts the correlation coefficients from the matched filter
     *        to detections.
     */
    void detect(const MFLib::SingleChannel::MatchedFilter<T> &mf);
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
