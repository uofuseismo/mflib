#include <cstdio>
#include <cstdlib>
#include <vector>
#include <pybind11/stl.h>
#include "singleChannel/pyAssociator.hpp"
#include "singleChannel/pyAssociatorParameters.hpp"
#include "singleChannel/pyDetector.hpp"
#include "mflib/singleChannel/associator.hpp"
#include "mflib/singleChannel/associatorParameters.hpp"
#include "mflib/singleChannel/detection.hpp"

using namespace PBMFLib::SingleChannel;

/// Constructors
template<class T>
Associator<T>::Associator() :
    mAssociator(std::make_unique<MFLib::SingleChannel::Associator<T>> ())
{
}

template<class T>
Associator<T>::Associator(const Associator &associator)
{
    *this = associator;
}

template<class T>
Associator<T>::Associator(const MFLib::SingleChannel::Associator<T> &associator)
{
    *this = associator;
}

/// Assignment operators
template<class T>
Associator<T>& Associator<T>::operator=(const Associator &associator)
{
    if (&associator == this){return *this;}
    mAssociator = std::make_unique<MFLib::SingleChannel::Associator<T>>
                  (*associator.mAssociator);
    return *this;
}

template<class T>
Associator<T>& Associator<T>::operator=(
    const MFLib::SingleChannel::Associator<T> &associator)
{
    mAssociator = std::make_unique<MFLib::SingleChannel::Associator<T>>
                  (associator);
    return *this;
}

/// Destructors
template<class T>
Associator<T>::~Associator() = default;

template<class T>
void Associator<T>::clear() noexcept
{
    mAssociator->clear();
}

/// Initialize
template<class T>
void Associator<T>::initialize(const AssociatorParameters &parameters)
{
    auto parms = parameters.getNativeClass();
    mAssociator->initialize(parms);
}

/// Initialized?
template<class T>
bool Associator<T>::isInitialized() const noexcept
{
    return mAssociator->isInitialized();
}

/// Add a detection
template<class T>
void Associator<T>::addDetection(
    const PBMFLib::SingleChannel::Detection<T> &detection)
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized\n");}
    auto det = detection.getNativeClass();
    if (!det.haveTravelTime())
    {
        throw std::invalid_argument("The detection must have a traveltime\n");
    }
    if (!det.haveTemplateIdentifier())
    {
        throw std::invalid_argument("Template identifier must be set\n");
    }
    mAssociator->addDetection(det);
}

template<class T>
void Associator<T>::clearDetections() noexcept
{
    mAssociator->clearDetections();
}

/// Gets the number of detections 
template<class T>
int Associator<T>::getNumberOfDetections() const noexcept
{
    return mAssociator->getNumberOfDetections();
}

/// Associate
template<class T>
void Associator<T>::associate()
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized\n");}
    mAssociator->associate();
} 

/// Get number of events
template<class T>
int Associator<T>::getNumberOfEvents() const noexcept
{
    return mAssociator->getNumberOfEvents();
}

/// Get the associated detections in an event
template<class T>
std::vector<PBMFLib::SingleChannel::Detection<T>>
Associator<T>::getDetectionsInEvent(const int iev) const
{
    auto nev = getNumberOfEvents();
    std::vector<PBMFLib::SingleChannel::Detection<T>> result;
    if (nev < 1)
    {
        fprintf(stderr, "No events\n");
        return result;
    }
    if (iev < 0 || iev >= nev)
    {
        throw std::invalid_argument("Event index = " + std::to_string(iev)
                                  + " must be in range [0,"
                                  + std::to_string(nev-1) + "]\n");
    }
    auto work = mAssociator->getDetectionsInEvent(iev);
    result.resize(work.size());
    for (int i=0; i<static_cast<int> (work.size()); ++i)
    {
        result[i] = work[i];
    }
    return result;

}

/// Initialize the class
void PBMFLib::SingleChannel::initializeAssociator(pybind11::module &m)
{
    pybind11::class_<PBMFLib::SingleChannel::Associator<double>>
        mAssociatorDouble(m, "Associator");
    /// The default constructor:
    mAssociatorDouble.def(pybind11::init<> ());
    /// Initialize
    mAssociatorDouble.def("initialize",
                          &PBMFLib::SingleChannel::Associator<double>::initialize,
                          "Initializes the associator.");
    mAssociatorDouble.def("is_initialized",
                          &PBMFLib::SingleChannel::Associator<double>::isInitialized,
                          "True indicates that the class is initialized.");
    /// Clears the detections
    mAssociatorDouble.def("add_detection",
                          &PBMFLib::SingleChannel::Associator<double>::addDetection,
                          "Adds a candidate detection to associate.  This class must be initialized prior to calling this.");
    mAssociatorDouble.def("get_number_of_detections",
                          &PBMFLib::SingleChannel::Associator<double>::getNumberOfDetections,
                          "Gets the number of detections set on the class.");
    mAssociatorDouble.def("clear_detections",
                          &PBMFLib::SingleChannel::Associator<double>::clearDetections, 
                          "Clears the detections but leaves the class initialized.");
    /// Associate
    mAssociatorDouble.def("associate",
                          &PBMFLib::SingleChannel::Associator<double>::associate,
                          "After the class is initialized and the detections are set this function will attempt to associate the detections into events with DBSCAN."); 
    mAssociatorDouble.def("get_number_of_events",
                          &PBMFLib::SingleChannel::Associator<double>::getNumberOfEvents,
                          "Following the association phase this is the number of events built from the associations.");
    mAssociatorDouble.def("get_detections_in_event",
                          &PBMFLib::SingleChannel::Associator<double>::getDetectionsInEvent,
                          "Provided there were events created, this function will return the detections that have been associated into the event with index iev.  The index iev should be in the range [0,get_number_of_events()-1].  If get_number_of_events() is zero then this function should not be called.");
}
/// Instantiation
template class PBMFLib::SingleChannel::Associator<double>;
//template class PBMFLib::SingleChannel::Associator<float>;
