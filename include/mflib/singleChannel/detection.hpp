#ifndef MFLIB_SINGLECHANNEL_DETECTION_HPP
#define MFLIB_SINGLECHANNEL_DETECTION_HPP
#include <memory>
namespace MFLib
{
namespace SingleChannel
{
/*!
 * @name Detection "detection.hpp" "mflib/singleChannel/detection.hpp"
 * @brief Defines a waveform detection.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
template<class T = double>
class Detection
{
public: 
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     */
    Detection();
    /*!
     * @brief Copy constructor.
     * @param[in] detection   The detection class from which to initialize
     *                        this class.
     */
    Detection(const Detection &detection);
    /*!
     * @brief Move constructor.
     * @param[in,out] detection  The detection class from which to initialize
     *                           this class.  On exit, detection's behavior
     *                           will be undefined.
     */
    Detection(Detection &&detection) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] detection  The detection class to copy to this.
     * @result A deep copy of of the input detection class.
     */
    Detection& operator=(const Detection &detection);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] detection  The detection class whose memory will be moved
     *                           to this.  On exit, detection's behavior is
     *                           undefined.
     * @result The memory from detection moved to this.
     */
    Detection& operator=(Detection &&detection) noexcept;
    /*! } */

    /*! @name Destrutors
     * @{
     */
    ~Detection();
    /*!
     * @brief Clears the class and resets the memory.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Template Identifier
     * @{
     */
    /*!
     * @brief Sets the waveform template identifier to which this 
     *        detection corresponds.
     * @param[in] id   The template identifier.
     */
    void setTemplateIdentifier(uint64_t id) noexcept;
    /*!
     * @brief The template index.
     * @throws std::runtime_error if the template identifier was not set.
     * @sa \c haveTemplateIdentifier()
     */
    uint64_t getTemplateIdentifier() const;
    /*!
     * @brief Determines if the detection identifier was set.
     */
    bool haveTemplateIdentifier() const noexcept;
    /*! @} */

    /*! @name Detected Signal
     * @{
     */
    /*!
     * @brief Sets the detected waveform.
     * @param[in] n   The number of samples in x.  This must be positive.
     * @param[in] x   The detected signal to set.
     * @throws std::invalid_argument if n is not postiive or x is NULL.
     */
    void setDetectedSignal(int n, const T x[]);
    /*!
     * @brief Gets the deteted signal.
     * @param[in] n   The space allocated to the output.  This must be at
     *                least \c getDetectedSignalLength().
     * @param[out] x  The detected signal.  This is an array whose dimension
     *                is [n] however only the first \c getDetectedSignalLength()
     *                samples are defined.
     * @throws std::runtime_error if the detected signal was not set.
     * @sa \c getDetectedSignalLength()
     */
    void getDetectedSignal(int n, T *x[]) const;
    /*!
     * @brief Gets a pointer to the detected signal.
     * @result This is an array of dimension [\c getDetetedSignalLength()].
     * @throws std::runtime_error if the detected signal is not set.
     */
    const T *getDetectedSignalPointer() const;
    /*!
     * @brief Gets the detected signal length.
     * @Throws std::runtime_error if the signal was not set.
     * @sa \c getDetectedSignalLength() 
     */
    int getDetectedSignalLength() const;
    /*!
     * @brief Determines if the detected signal was set.
     * @result True indicates that the detected signal was set.
     */
    bool haveDetectedSignal() const noexcept;
    /*!} */

    /*! @brief Detection Time
     * @{
     */
    /*!
     * @brief The time corresponding to the detection.
     * @param[in] time  The time of the detection relative to the trace start
     *                  in seconds.
     * @throws std::invalid_argument if the detection time is negative.
     */
    void setDetectionTime(double time);    
    /*!
     * @brief Gets the time corresponding to the detection.
     * @result The time of the detection relative to the start of the trace
     *         in seconds.
     * @throws std::runtime_error if the detection time was not set.
     */
    double getDetectionTime() const;
    /*!
     * @brief Determines if the detection time was set.
     * @result True indicates that the detection time was set.
     * @throws std::runtime_error if the detection time was not set.
     */
    bool haveDetectionTime() const noexcept;
    /*! @} */

    /*! @brief Phase Onset Time
     * @{
     */
    /*!
     * @brief Sets the phase onset time.
     * @param[in] onsetTime   The onset time of the phase relative to the
     *                        start of the trace in seconds.
     * @throws std::invalid_argument if the onset time is negative.
     */
    void setPhaseOnsetTime(double onsetTime);
    /*!
     * @brief Gets the onset time of the phase.
     * @result The phase onset time relative to the start of the trace
     *         in seconds.
     * @throws std::runtime_error if the phase onset time was not set.
     * @sa \c havePhaseOnsetTime()
     */
    double getPhaseOnsetTime() const;
    /*!
     * @brief Determines if the phase onset time was set.
     * @result True indicates that the phase onset time was set.
     */
    bool havePhaseOnsetTime() const noexcept;
    /*! @} */

    /*! 
     * @brief Determines if this class has a detection.
     * @result True indicates that this class is tied to a detection.
     */
    //bool haveDetection() const noexcept;
private:
    class DetectionImpl;
    std::unique_ptr<DetectionImpl> pImpl;
};
}
}
#endif

