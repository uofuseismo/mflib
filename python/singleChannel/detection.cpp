#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
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
PBMFLib::SingleChannel::Detection<T>& 
Detection<T>::operator=(const PBMFLib::SingleChannel::Detection<T> &det)
{
    if (&det == this){return *this;}
    mDetection = std::make_unique<MFLib::SingleChannel::Detection<T>>
                 (*det.mDetection);
    return *this;
}

template<class T>
PBMFLib::SingleChannel::Detection<T>& Detection<T>::operator=(
    const MFLib::SingleChannel::Detection<T> &det)
{
    mDetection = std::make_unique<MFLib::SingleChannel::Detection<T>> (det);
    return *this;
}

/// Move assignment operator
template<class T>
Detection<T>& Detection<T>::operator=(Detection<T> &&det) noexcept
{
    if (&det == this){return *this;}
    mDetection = std::move(det.mDetection);
    return *this;
}

template<class T>
Detection<T>&
Detection<T>::operator=(MFLib::SingleChannel::Detection<T> &&det) noexcept
{
    *mDetection = std::move(det);
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
    pybind11::array_t<double, pybind11::array::c_style |
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
    pybind11::array_t<double, pybind11::array::c_style |
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

template<class T>
void Detection<T>::setDetectedSignalAsVector(const std::vector<T> &x)
{
    mDetection->setDetectedSignal(x.size(), x.data());
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
void Detection<T>::getDetectedSignalAsVector(std::vector<T> &x) const
{
    auto npts = mDetection->getDetectedSignalLength();
    x.resize(npts);
    auto xPtr = x.data();
    mDetection->getDetectedSignal(npts, &xPtr);
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

/// Travel time
template<class T>
void Detection<T>::setTravelTime(const double travelTime)
{
    if (travelTime < 0)
    {
        throw std::invalid_argument("Travel time cannot be negative\n");
    }
    mDetection->setTravelTime(travelTime);
}

template<class T>
double Detection<T>::getTravelTime() const
{
    if (!haveTravelTime()){throw std::runtime_error("Travel time not set\n");}
    return mDetection->getTravelTime();
}

template<class T>
bool Detection<T>::haveTravelTime() const noexcept
{
    return mDetection->haveTravelTime();
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
void Detection<T>::setAmplitudeScalingFactor(
    const double value,
    const MFLib::RelativeMagnitudeType type)
{
    if (value <= 0)
    {
        throw std::invalid_argument("Scaling factor = " + std::to_string(value)
                                  + " must be positive\n");
    }
    mDetection->setAmplitudeScalingFactor(value, type);
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

template<class T>
void Detection<T>::setPolarity(const MFLib::Polarity polarity) noexcept
{
    mDetection->setPolarity(polarity);
}

template<class T>
MFLib::Polarity Detection<T>::getPolarity() const noexcept
{
    return mDetection->getPolarity();
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
    // Travel time
    mDetDouble.def("set_travel_time",
                   &PBMFLib::SingleChannel::Detection<double>::setTravelTime,
                   "Sets the catalog travel time from the source to the receiver in seconds.");
    mDetDouble.def("get_travel_time",
                   &PBMFLib::SingleChannel::Detection<double>::getTravelTime,
                   "Gets the catalog travel time from the source to the receiveri in seconds.");
    mDetDouble.def("have_travel_time",
                   &PBMFLib::SingleChannel::Detection<double>::haveTravelTime,
                   "True indicates the travel time was set."); 
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
    // Amplitude scaling factor
    mDetDouble.def("set_amplitude_scaling_factor",
                   &PBMFLib::SingleChannel::Detection<double>::setAmplitudeScalingFactor,
                   "Sets the amplitude scaling factor to best the template's amplitude to the detected signal's amplitude.");
    mDetDouble.def("get_amplitude_scaling_factor",
                   &PBMFLib::SingleChannel::Detection<double>::getAmplitudeScalingFactor,
                   "Gets the amplitude scaling factor to best match the template's amplitude to the detected signal's amplitude.");
    mDetDouble.def("get_magnitude_perturbation",
                   &PBMFLib::SingleChannel::Detection<double>::getMagnitudePerturbation,
                   "GEts the magnitude perturbation that should be added to the template magnitude to obtain the relative magnitude.");
    mDetDouble.def("have_amplitude_scaling_factor",
                   &PBMFLib::SingleChannel::Detection<double>::haveAmplitudeScalingFactor,
                   "Determines if the amplitude scaling factors were computed");
    // Polarity
    mDetDouble.def("get_polarity",
                   &PBMFLib::SingleChannel::Detection<double>::getPolarity,
                   "Gets the polarity of the detection.");
    mDetDouble.def("set_polarity",
                   &PBMFLib::SingleChannel::Detection<double>::setPolarity,
                   "Sets the polarity of the detection.");
    // Clears
    mDetDouble.def("clear",
                   &PBMFLib::SingleChannel::Detection<double>::clear,
                   "Clears the class and releases memory.");

    /// Makes this class pickleable
    mDetDouble.def("__getstate__",
                   [](const PBMFLib::SingleChannel::Detection<double> &p)
    {
        // XC
        bool haveXC = p.haveCorrelationCoefficient();
        double xc = 0;
        if (haveXC){xc = p.getCorrelationCoefficient();}
        // Template identifier
        bool haveID = p.haveTemplateIdentifier();
        std::string network = "";
        std::string station = "";
        std::string channel = "";
        std::string location = "";
        std::string phase = ""; 
        uint64_t waveid = 0; 
        if (haveID)
        {
            auto id = p.getTemplateIdentifier();
            network = id.first.getNetwork();
            station = id.first.getStation();
            channel = id.first.getChannel();
            location = id.first.getLocationCode();
            phase = id.first.getPhase();
            waveid = id.second; 
        }
        // Detection time
        bool haveDetTime = p.haveDetectionTime();
        double detTime = 0;
        if (haveDetTime){detTime = p.getDetectionTime();}
        // Interpolated detection time
        bool haveIntDetTime = p.haveInterpolatedDetectionTime();
        double intDetTime = 0;
        if (haveIntDetTime){intDetTime = p.getInterpolatedDetectionTime();}
        // Onset time
        bool haveOnsetTime = p.havePhaseOnsetTime();
        double onsetTime = 0;
        if (haveOnsetTime){onsetTime = p.getPhaseOnsetTime();}
        // Interpolated phase onset time
        bool haveIntOnsetTime = p.haveInterpolatedPhaseOnsetTime();
        double intOnsetTime = 0;
        if (haveIntOnsetTime){intOnsetTime = p.getInterpolatedPhaseOnsetTime();} 
        // Travel time
        bool haveTravelTime = p.haveTravelTime();
        double travelTime = 0;
        if (haveTravelTime){travelTime = p.getTravelTime();}
        // Amplitude 1, 2...
        bool haveAmp = p.haveAmplitudeScalingFactor();
        double amp1 = 0;
        double amp2 = 0;
        if (haveAmp)
        {
            amp1 = p.getAmplitudeScalingFactor(
                       MFLib::RelativeMagnitudeType::GIBBONS_RINGDAL_2006);
        }
        if (haveAmp)
        {
            amp2 = p.getAmplitudeScalingFactor(
                       MFLib::RelativeMagnitudeType::SCHAFF_RICHARDS_2014);
        }
        /*
        bool haveSignal = p.haveDetectedSignal();
        std::vector<double> x(1);
        if (haveSignal){p.getDetectedSignalAsVector(x);}
        auto cPtr = reinterpret_cast<char *> (x.data());
        std::string sx(cPtr, cPtr+x.size()*sizeof(double));
        */
        // Polarity
        int polarity = static_cast<int> (p.getPolarity());
        return pybind11::make_tuple(haveXC, xc,
                                    haveID, network, station, channel, location, phase, waveid,
                                    haveDetTime, detTime,
                                    haveIntDetTime, intDetTime,
                                    haveOnsetTime, onsetTime,
                                    haveIntOnsetTime, intOnsetTime,
                                    haveTravelTime, travelTime,
                                    haveAmp, amp1,
                                    haveAmp, amp2,
                                    polarity);
//                                    haveSignal, sx);
/*
        return pybind11::make_tuple(p.getNetwork(),
                                    p.getStation(),
                                    p.getChannel(),
                                    p.getLocationCode(),
                                    p.getPhase());
*/
    });
    mDetDouble.def("__setstate__",
                   [](PBMFLib::SingleChannel::Detection<double> &p,
                      pybind11::tuple t)
    {
        if (t.size() != 24)
        {
            throw std::runtime_error("Tuple in invalid state\n");
        }
        // Call constructor -> memory leak because of new?
        new (&p) PBMFLib::SingleChannel::Detection<double>();

        bool haveXC = t[0].cast<bool> ();
        if (haveXC){p.setCorrelationCoefficient(t[1].cast<double> ());}

        bool haveID = t[2].cast<bool> ();
        if (haveID)
        {
            PBMFLib::NetworkStationPhase nsp;
            nsp.setNetwork(t[3].cast<std::string> ());
            nsp.setStation(t[4].cast<std::string> ());
            nsp.setChannel(t[5].cast<std::string> ());
            nsp.setLocationCode(t[6].cast<std::string> ());
            nsp.setPhase(t[7].cast<std::string> ());
            uint64_t waveid = t[8].cast<uint64_t> (); 
            p.setTemplateIdentifier(std::make_pair(nsp, waveid));
        }

        bool haveDetTime = t[9].cast<bool> ();
        if (haveDetTime){p.setDetectionTime(t[10].cast<double> ());}

        bool haveIntDetTime = t[11].cast<bool> ();
        if (haveIntDetTime)
        {
            p.setInterpolatedDetectionTime(t[12].cast<double> ());
        }

        bool haveOnsetTime = t[13].cast<bool> ();
        if (haveOnsetTime){p.setPhaseOnsetTime(t[14].cast<double> ());}

        bool haveIntOnsetTime = t[15].cast<bool> ();
        if (haveIntOnsetTime)
        {
            p.setInterpolatedPhaseOnsetTime(t[16].cast<double> ());
        }

        bool haveTravelTime = t[17].cast<bool> ();
        if (haveTravelTime)
        {
            p.setTravelTime(t[18].cast<double> ());
        }

        bool haveAmp1 = t[19].cast<bool> ();
        if (haveAmp1)
        {
            p.setAmplitudeScalingFactor(t[20].cast<double> (),
                   MFLib::RelativeMagnitudeType::GIBBONS_RINGDAL_2006);
        }
        bool haveAmp2 = t[21].cast<bool> ();
        if (haveAmp2)
        {
            p.setAmplitudeScalingFactor(t[22].cast<double> (),
                   MFLib::RelativeMagnitudeType::SCHAFF_RICHARDS_2014);
        }

        auto polarity = static_cast<MFLib::Polarity> (t[23].cast<int> ());
        p.setPolarity(polarity);

        /*
        auto haveSignal = t[22].cast<bool> ();
        if (haveSignal)
        {
           std::string cx = t[23].cast<std::string> ();
           auto xPtr = reinterpret_cast<const double *> (cx.data());
           auto npts = cx.size()/sizeof(double);
           std::vector<double> x(xPtr, xPtr+npts);
printf("%lf, %lf\n", x[0], x[1]);
           p.setDetectedSignalAsVector(x);
        }
        */
    });
}

/// Instantiation
template class PBMFLib::SingleChannel::Detection<double>;
template class PBMFLib::SingleChannel::Detection<float>;

