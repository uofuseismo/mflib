#include <cstdio>
#include <cstdlib>
#include <string>
#include "mflib/waveformTemplate.hpp"
#include "mflib/singleChannel/detection.hpp"
#include "singleChannel/pyDetection.hpp"

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
void Detection<T>::setTemplateIdentifier(const uint64_t id) noexcept
{
    mDetection->setTemplateIdentifier(id);
}

template<class T>
uint64_t Detection<T>::getTemplateIdentifier() const
{
    return mDetection->getTemplateIdentifier();
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
        mDetDouble(m, "SingleChannelDetection");
    /// The default constructor:
    mDetDouble.def(pybind11::init<> ());
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
    // Clears
    mDetDouble.def("clear",
                   &PBMFLib::SingleChannel::Detection<double>::clear,
                   "Clears the class and releases memory.");
}

/// Instantiation
template class PBMFLib::SingleChannel::Detection<double>;
template class PBMFLib::SingleChannel::Detection<float>;

