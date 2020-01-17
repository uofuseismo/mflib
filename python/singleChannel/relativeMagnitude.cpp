#include <cstdio>
#include <cstdlib>
#include <string>
#include "mflib/singleChannel/relativeMagnitude.hpp"
#include "mflib/waveformTemplate.hpp"
#include "mflib/enums.hpp"
#include "pyWaveformTemplate.hpp"
#include "pySingleChannelRelativeMagnitude.hpp"

using namespace PBMFLib::SingleChannel;

/// Constructor
template<class T>
RelativeMagnitude<T>::RelativeMagnitude() :
    mMag(std::make_unique<MFLib::SingleChannel::RelativeMagnitude<T>> ())
{
}

/// Copy c'tor
template<class T>
RelativeMagnitude<T>::RelativeMagnitude(const RelativeMagnitude &rmag)
{
    *this = rmag;
}

/// Copy assignment operator
template<class T>
RelativeMagnitude<T>
RelativeMagnitude<T>::operator=(const RelativeMagnitude &rmag)
{
    if (&rmag == this){return *this;}
    mMag = std::make_unique<MFLib::SingleChannel::RelativeMagnitude<T>>
           (*rmag.mMag);
    return *this;
}

/// Destructor
template<class T>
RelativeMagnitude<T>::~RelativeMagnitude() = default;

/// Initialize
template<class T>
void RelativeMagnitude<T>::initialize(const WaveformTemplate &wt)
{
    auto t = wt.getNativeClass();
    mMag->initialize(t);
}

/// Gets the signal length
template<class T>
int RelativeMagnitude<T>::getDetectedSignalLength() const
{
    if (!isInitialized())
    {
        throw std::runtime_error("Class is not initialized\n");
    }
    return mMag->getDetectedSignalLength();
}

/// Determines if the class is initialized
template<class T>
bool RelativeMagnitude<T>::isInitialized() const noexcept
{
    return mMag->isInitialized();
}

/// Set the waveform
template<>
void RelativeMagnitude<double>::setDetectedSignal(
    const pybind11::array_t<double, pybind11::array::c_style |
                                    pybind11::array::forcecast> &x)
{
    auto nref = getDetectedSignalLength();
    pybind11::buffer_info xbuf = x.request();
    auto len = static_cast<int> (xbuf.size);
    if (len != nref)
    {
        throw std::invalid_argument("signal length = "
                                  + std::to_string(len)
                                  + " must = " + std::to_string(nref) + "\n");
    }
    const double *xptr = (double *) (xbuf.ptr);
    if (xptr == nullptr){throw std::invalid_argument("signal is null\n");}
    mMag->setDetectedSignal(len, xptr);
}

template<>
void RelativeMagnitude<float>::setDetectedSignal(
    const pybind11::array_t<double, pybind11::array::c_style |
                                    pybind11::array::forcecast> &x) 
{
    auto nref = getDetectedSignalLength();
    pybind11::buffer_info xbuf = x.request();
    auto len = static_cast<int> (xbuf.size);
    if (len != nref)
    {   
        throw std::invalid_argument("signal length = " 
                                  + std::to_string(len)
                                  + " must = " + std::to_string(nref) + "\n");
    }
    const double *xptr = (double *) (xbuf.ptr);
    if (xptr == nullptr){throw std::invalid_argument("signal is null\n");}
    std::vector<float> xIn(len);
    #pragma omp simd
    for (int i=0; i<len; ++i){xIn[i] = static_cast<float> (xptr[i]);}
    if (xptr == nullptr){throw std::invalid_argument("signal is null\n");}
    mMag->setDetectedSignal(len, xIn.data());
}

template<class T>
T RelativeMagnitude<T>::computeAmplitudeScalingFactor(
    const MFLib::RelativeMagnitudeType type) const
{
    return mMag->computeAmplitudeScalingFactor(type);
}

template<class T>
T RelativeMagnitude<T>::computeMagnitudePerturbation(
    const MFLib::RelativeMagnitudeType type) const
{
    return mMag->computeMagnitudePerturbation(type);
}

template<class T>
bool RelativeMagnitude<T>::haveDetectedSignal() const noexcept
{
    return mMag->haveDetectedSignal();
}

/// Initialize the class
void PBMFLib::SingleChannel::initializeRelativeMagnitude(pybind11::module &m) 
{
    pybind11::class_<PBMFLib::SingleChannel::RelativeMagnitude<double>>
        mRmagDouble(m, "SingleChannelRelativeMagnitude");
    /// The default constructor:
    mRmagDouble.def(pybind11::init<> ());
    mRmagDouble.doc() = "This is class is used for computing relative magnitudes and amplitudes between detected and template events on a channel-by-channel basis.";
    /// Methods
    mRmagDouble.def("initialize",
                    &PBMFLib::SingleChannel::RelativeMagnitude<double>::initialize, 
                    "This is initializes the a relative magnitude class from the given waveform template - i.e., this class will compute magnitudes and amplitudes relative to the given template.  This waveform template must be a non-constant signal."); 
    mRmagDouble.def("is_initialized",
                    &PBMFLib::SingleChannel::RelativeMagnitude<double>::isInitialized,
                    "Determines if the class is initialized or not.");
    mRmagDouble.def("get_detected_signal_length",
                    &PBMFLib::SingleChannel::RelativeMagnitude<double>::getDetectedSignalLength,
                    "Gets the expected length of the detected signal."); 
    mRmagDouble.def("set_detected_signal",
                    &PBMFLib::SingleChannel::RelativeMagnitude<double>::setDetectedSignal,
                    "Sets the signal detected by the template.  The signal is a NumPy array whose length must must equal the value return by get_detected_signal_length.");
    mRmagDouble.def("have_detected_signal",
                    &PBMFLib::SingleChannel::RelativeMagnitude<double>::haveDetectedSignal,
                    "Determines if the detected signal was set.");
    mRmagDouble.def("compute_amplitude_scaling_factor",
                    &PBMFLib::SingleChannel::RelativeMagnitude<double>::computeAmplitudeScalingFactor,
                    "Computes the amplitude scaling factor that scales the template waveform so that it best matches the detected waveform.  The class must be initialized and the detected signal must be set prior to calling this function.");
    mRmagDouble.def("compute_magnitude_perturbation",
                    &PBMFLib::SingleChannel::RelativeMagnitude<double>::computeMagnitudePerturbation,
                    "Computes the magnitude perturbation to apply to the reference magnitude as to obtain the magnitude for the detected event.  The class must be initialized and the detected signal must be set prior to calling this function.");
    pybind11::enum_<MFLib::RelativeMagnitudeType> (mRmagDouble, "RelativeMagnitudeType")
        .value("gibbons_ringdal_2006", MFLib::RelativeMagnitudeType::GIBBONS_RINGDAL_2006,
               "Computes the relative magnitude using the methodology of Gibbons and Ringdal, 2006.  This method assumes the the template waveform is noise-free.")
        .value("schaff_richards_2014", MFLib::RelativeMagnitudeType::SCHAFF_RICHARDS_2014,
               "Computes the relative magnitude using the methodology of Schaff and Richards, 2014.  This method accomodates noise in both the detected and template waveform.");
  
}

/// Instantiation
template class PBMFLib::SingleChannel::RelativeMagnitude<double>;
template class PBMFLib::SingleChannel::RelativeMagnitude<float>;
