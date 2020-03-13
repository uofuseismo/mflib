#include <cstdio>
#include <cstdlib>
#include <string>
#include "mflib/waveformTemplate.hpp"
#include "mflib/singleChannel/detection.hpp"
#include "singleChannel/pyDetection.hpp"
#include "pyNetworkStationPhase.hpp"

using namespace PBMFLib::SingleChannel;

/// Constructor
template<class T>
Detection<T>::Detection() :
    mDetection(std::make_unique<MFLib::SingleChannel::Detection<T>> ())
{
}

/// Copy c'tor
template<class T>
Detection<T>::Detection(const Detection &detection)
{
    *this = detection;
}

template<class T>
Detection<T>::Detection(const MFLib::SingleChannel::Detection<T> &detection)
{
    *this = detection;
}

/// Move c'tor
template<class T>
Detection<T>::Detection(Detection &&detection) noexcept
{
    *this = std::move(detection);
}

/// Copy assignment operator
template<class T>
Detection<T>& Detection<T>::operator=(const Detection &det)
{
    if (&det == this){return *this;}
    mDetection = std::make_unique<MFLib::SingleChannel::Detection<T>>
                 (*det.mDetection);
    return *this;
}

template<class T>
Detection<T>& Detection<T>::operator=(
    const MFLib::SingleChannel::Detection<T> &det)
{
    mDetection = std::make_unique<MFLib::SingleChannel::Detection<T>> (det);
    return *this;
}

/// Move assignment operator
template<class T>
Detection<T>& Detection<T>::operator=(Detection &&det) noexcept
{
    if (&det == this){return *this;}
    mDetection = std::move(det.mDetection);
    return *this;
}

/// Destructor
template<class T>
Detection<T>::~Detection() = default;

/// Clears the class
template<class T>
void Detection<T>::clear() noexcept
{
    return mDetection->clear();
}

/// Waveform identifier
template<class T>
void Detection<T>::setTemplateIdentifier(
    const std::pair<PBMFLib::NetworkStationPhase, uint64_t> &idIn) noexcept
{
    auto nsp = idIn.first.getNativeClass();
    std::pair<MFLib::NetworkStationPhase, uint64_t> id(nsp, idIn.second);
    mDetection->setTemplateIdentifier(id);
}

template<class T> std::pair<PBMFLib::NetworkStationPhase, uint64_t>
Detection<T>::getTemplateIdentifier() const
{
    auto id = mDetection->getTemplateIdentifier();
    PBMFLib::NetworkStationPhase pid(id.first);
    std::pair<PBMFLib::NetworkStationPhase, uint64_t> idOut(pid, id.second);
    return idOut;
}

template<class T>
bool Detection<T>::haveTemplateIdentifier() const noexcept
{
    return mDetection->haveTemplateIdentifier();
}

// Detected signal
template<>
void Detection<double>::setDetectedSignal(
    const pybind11::array_t<double, pybind11::array::c_style |
                                    pybind11::array::forcecast> &x) 
{
    pybind11::buffer_info xbuf = x.request();
    auto len = static_cast<int> (xbuf.size);
    const double *xptr = (double *) (xbuf.ptr);
    if (xptr == nullptr)
    {
        throw std::invalid_argument("x is null");
    }
    mDetection->setDetectedSignal(len, xptr);
}

template<>
void Detection<float>::setDetectedSignal(
    const pybind11::array_t<double, pybind11::array::c_style |
                                    pybind11::array::forcecast> &x) 
{
    pybind11::buffer_info xbuf = x.request();
    auto len = static_cast<int> (xbuf.size);
    const double *xptr = (double *) (xbuf.ptr);
    if (xptr == nullptr)
    {   
        throw std::invalid_argument("x is null");
    }
    std::vector<float> xWrite(len);
    #pragma omp simd
    for (int i=0; i<len; ++i){xWrite[i] = static_cast<float> (xptr[i]);}
    mDetection->setDetectedSignal(len, xWrite.data());
}

template<>
pybind11::array_t<double> Detection<double>::getDetectedSignal() const
{
    auto npts = mDetection->getDetectedSignalLength();
    auto y = pybind11::array_t<double, pybind11::array::c_style> (npts);
    pybind11::buffer_info ybuf = y.request();
    auto yptr = static_cast<double *> (ybuf.ptr);
    mDetection->getDetectedSignal(npts, &yptr);
    return y;
}

template<>
pybind11::array_t<double> Detection<float>::getDetectedSignal() const
{
    auto npts = mDetection->getDetectedSignalLength();
    auto y = pybind11::array_t<double, pybind11::array::c_style> (npts);
    pybind11::buffer_info ybuf = y.request();
    auto yptr = static_cast<double *> (ybuf.ptr);
    auto temp = mDetection->getDetectedSignalPointer();
    #pragma omp simd
    for (int i=0; i<npts; ++i){yptr[i] = static_cast<double> (temp[i]);}
    return y;
}

template<class T>
bool Detection<T>::haveDetectedSignal() const noexcept
{
    return mDetection->haveDetectedSignal();
}

// Detection time
template<class T>
void Detection<T>::setDetectionTime(const double time)
{
    if (time < 0)
    {
        throw std::invalid_argument("Detection time must be positive\n");
    } 
    mDetection->setDetectionTime(time);
}

template<class T>
double Detection<T>::getDetectionTime() const
{
    return mDetection->getDetectionTime();
}

template<class T>
bool Detection<T>::haveDetectionTime() const noexcept
{
    return mDetection->haveDetectionTime();
}

// Interpolated detection time
template<class T>
void Detection<T>::setInterpolatedDetectionTime(const double time)
{
    if (time < 0)
    {
        throw std::invalid_argument("Detection time must be positive\n");
    }   
    mDetection->setInterpolatedDetectionTime(time);
}

template<class T>
double Detection<T>::getInterpolatedDetectionTime() const
{
    return mDetection->getInterpolatedDetectionTime();
}

template<class T>
bool Detection<T>::haveInterpolatedDetectionTime() const noexcept
{
    return mDetection->haveInterpolatedDetectionTime();
}

/// Onset time
template<class T>
void Detection<T>::setPhaseOnsetTime(const double time) 
{
    if (time < 0)
    {
        throw std::invalid_argument("Onset time must be positive\n");
    }
    mDetection->setPhaseOnsetTime(time);
}

template<class T>
double Detection<T>::getPhaseOnsetTime() const
{
    return mDetection->getPhaseOnsetTime();
}

template<class T>
bool Detection<T>::havePhaseOnsetTime() const noexcept
{
    return mDetection->havePhaseOnsetTime();
}

/// Interpolated onset time
template<class T>
void Detection<T>::setInterpolatedPhaseOnsetTime(const double time) 
{
    if (time < 0)
    {
        throw std::invalid_argument("Onset time must be positive\n");
    }
    mDetection->setInterpolatedPhaseOnsetTime(time);
}

template<class T>
double Detection<T>::getInterpolatedPhaseOnsetTime() const
{
    return mDetection->getInterpolatedPhaseOnsetTime();
}

template<class T>
bool Detection<T>::haveInterpolatedPhaseOnsetTime() const noexcept
{
    return mDetection->haveInterpolatedPhaseOnsetTime();
}

/// Detection time
template<class T>
void Detection<T>::setCorrelationCoefficient(const double xc)
{
    mDetection->setCorrelationCoefficient(xc);
}

template<class T>
double Detection<T>::getCorrelationCoefficient() const
{
    return mDetection->getCorrelationCoefficient();
}

template<class T>
bool Detection<T>::haveCorrelationCoefficient() const noexcept
{
    return mDetection->haveCorrelationCoefficient();
}

template<class T>
double Detection<T>::getMagnitudePerturbation(
    const MFLib::RelativeMagnitudeType type) const
{
    return mDetection->getMagnitudePerturbation(type);
}

template<class T>
double Detection<T>::getAmplitudeScalingFactor(
    const MFLib::RelativeMagnitudeType type) const
{
    return mDetection->getAmplitudeScalingFactor(type);
}
 
template<class T>
bool Detection<T>::haveAmplitudeScalingFactor() const noexcept
{

    if (mDetection->haveAmplitudeScalingFactor(MFLib::RelativeMagnitudeType::GIBBONS_RINGDAL_2006) &&
        mDetection->haveAmplitudeScalingFactor(MFLib::RelativeMagnitudeType::SCHAFF_RICHARDS_2014))
    {
        return true;
    }
    return false;
}

/// Gets the native class
template<class T>
MFLib::SingleChannel::Detection<T> 
Detection<T>::getNativeClass() const
{
    auto detOut(*mDetection);
    return detOut;
}

/// Initialize the class
void PBMFLib::SingleChannel::initializeDetection(pybind11::module &m) 
{
    pybind11::class_<PBMFLib::SingleChannel::Detection<double>>
        mDetDouble(m, "Detection");
    /// The default constructor:
    mDetDouble.def(pybind11::init<> ());
    /// Correlation score
    mDetDouble.def("get_correlation_coefficient",
                   &PBMFLib::SingleChannel::Detection<double>::getCorrelationCoefficient,
                   "Gets the correlation coefficient associated with the detection.");
    mDetDouble.def("set_correlation_coefficient",
                   &PBMFLib::SingleChannel::Detection<double>::setCorrelationCoefficient,
                   "Sets the correlation coefficient associated with the detection.  This must be in the range [-1,1]");
    mDetDouble.def("have_correlation_coefficient",
                   &PBMFLib::SingleChannel::Detection<double>::haveCorrelationCoefficient,
                   "Determines if the correlation coefficient was set.");
    /// Template identifier
    mDetDouble.def("get_template_identifier",
                   &PBMFLib::SingleChannel::Detection<double>::getTemplateIdentifier,
                   "Gets the template identifier associated with the detection.");
    mDetDouble.def("set_template_identifier",
                   &PBMFLib::SingleChannel::Detection<double>::setTemplateIdentifier,
                   "Sets the template identifier associated with the detection.");
    mDetDouble.def("have_template_identifier",
                   &PBMFLib::SingleChannel::Detection<double>::haveTemplateIdentifier,
                   "Determines if the template identifier was set.");
    /// Detected signal
    mDetDouble.def("get_detected_signal",
                   &PBMFLib::SingleChannel::Detection<double>::getDetectedSignal,
                   "Gets the detected signal.  The length of the signal should match the length of the waveform template signal.");
    mDetDouble.def("set_detected_signal",
                   &PBMFLib::SingleChannel::Detection<double>::setDetectedSignal,
                   "Sets the detected signal.");
    mDetDouble.def("have_detected_signal",
                   &PBMFLib::SingleChannel::Detection<double>::haveDetectedSignal,
                   "Determines if the detected signal was set.");
    /// Detection time 
    mDetDouble.def("get_detection_time",
                   &PBMFLib::SingleChannel::Detection<double>::getDetectionTime,
                   "Gets the detection time in seconds relative to the trace start.");
    mDetDouble.def("set_detection_time",
                   &PBMFLib::SingleChannel::Detection<double>::setDetectionTime,
                   "Sets the detection time in seconds relative to the trace start.");
    mDetDouble.def("have_detection_time",
                   &PBMFLib::SingleChannel::Detection<double>::haveDetectionTime,
                   "Determines if the detection time was set.");
    /// Interpolated detection time 
    mDetDouble.def("get_interpolated_detection_time",
                   &PBMFLib::SingleChannel::Detection<double>::getInterpolatedDetectionTime,
                   "Gets the interpolated detection time in seconds relative to the trace start.  This can have subsampling-rate resolution.");
    mDetDouble.def("set_interpolated_detection_time",
                   &PBMFLib::SingleChannel::Detection<double>::setInterpolatedDetectionTime,
                   "Sets the detection time in seconds relative to the trace start.  This can ahve subsampling-rate resolution.");
    mDetDouble.def("have_interpolated_detection_time",
                   &PBMFLib::SingleChannel::Detection<double>::haveInterpolatedDetectionTime,
                   "Determines if the interpolated detection time was set.");
    // Phase onset time
    mDetDouble.def("get_phase_onset_time",
                   &PBMFLib::SingleChannel::Detection<double>::getPhaseOnsetTime,
                   "Gets the phase onset time in seconds relative to the trace start.  This likely will only be accessible if the phase onsets were set in the waveform templates.");
    mDetDouble.def("set_phase_onset_time",
                   &PBMFLib::SingleChannel::Detection<double>::setPhaseOnsetTime,
                   "Sets the phase onset time in seconds relative to the trace start.");
    mDetDouble.def("have_phase_onset_time",
                   &PBMFLib::SingleChannel::Detection<double>::havePhaseOnsetTime,
                   "Determines if the phase onset time is set.");
    // Phase onset time
    mDetDouble.def("get_interpolated_phase_onset_time",
                   &PBMFLib::SingleChannel::Detection<double>::getInterpolatedPhaseOnsetTime,
                   "Gets the interpolated phase onset time in seconds relative to the trace start.  This time can have subsampling-rate resolution.  This likely will only be accessible if the phase onsets were set in the waveform templates.");
    mDetDouble.def("set_interpolated_phase_onset_time",
                   &PBMFLib::SingleChannel::Detection<double>::setInterpolatedPhaseOnsetTime,
                   "Sets the phase onset time in seconds relative to the trace start.  This time can have subsampling-rate resolution.");
    mDetDouble.def("have_interpolated_phase_onset_time",
                   &PBMFLib::SingleChannel::Detection<double>::haveInterpolatedPhaseOnsetTime,
                   "Determines if the interpolated phase onset time is set.");
    // 
    mDetDouble.def("get_amplitude_scaling_factor",
                   &PBMFLib::SingleChannel::Detection<double>::getAmplitudeScalingFactor,
                   "Gets the amplitude scaling factor to best match the template to the detected signal.");
    mDetDouble.def("get_magnitude_perturbation",
                   &PBMFLib::SingleChannel::Detection<double>::getMagnitudePerturbation,
                   "GEts the magnitude perturbation that should be added to the template magnitude to obtain the relative magnitude.");
    mDetDouble.def("have_amplitude_scaling_factor",
                   &PBMFLib::SingleChannel::Detection<double>::haveAmplitudeScalingFactor,
                   "Determines if the amplitude scaling factors were computed");
    // Clears
    mDetDouble.def("clear",
                   &PBMFLib::SingleChannel::Detection<double>::clear,
                   "Clears the class and releases memory.");
}

/// Instantiation
template class PBMFLib::SingleChannel::Detection<double>;
template class PBMFLib::SingleChannel::Detection<float>;

