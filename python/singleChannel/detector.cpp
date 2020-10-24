#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include "mflib/singleChannel/detector.hpp"
#include "mflib/singleChannel/detectorParameters.hpp"
#include "mflib/singleChannel/matchedFilter.hpp"
#include "singleChannel/pyMatchedFilter.hpp"
#include "singleChannel/pyDetectorParameters.hpp"
#include "singleChannel/pyDetector.hpp"

using namespace PBMFLib::SingleChannel;

/// Constructor
template<class T>
Detector<T>::Detector() :
    mDetector(std::make_unique<MFLib::SingleChannel::Detector<T>> ())
{
}

/// Destructor
template<class T>
Detector<T>::~Detector() = default;

/// Initialize
template<class T>
void Detector<T>::initialize(const DetectorParameters &parmsIn)
{
    auto parms = parmsIn.getNativeClass();
    mDetector->initialize(parms);
}

/// Perform detections 
template<class T>
void Detector<T>::detect(const MatchedFilter<T> &mf)
{
    if (!mf.haveMatchedFilteredSignals())
    {
        throw std::invalid_argument("Matched filtered signals not yet computed\n");
    }
    auto mfPtr = mf.getNativeClassPointer();
    mDetector->detect(*mfPtr);
    mfPtr = nullptr;
}

/// Get the number of detections
template<class T>
int Detector<T>::getNumberOfDetections() const noexcept
{
    return mDetector->getNumberOfDetections();
}

/// Get a detection
template<class T>
PBMFLib::SingleChannel::Detection<T> Detector<T>::getDetection(const int it) const
{
    auto nt = getNumberOfDetections();
    if (nt == 0){throw std::invalid_argument("No detections\n");}
    if (it < 0 || it >= getNumberOfDetections())
    {
        throw std::invalid_argument("template = " + std::to_string(it)
                                  + " must be in range [0,"
                                  + std::to_string(nt-1) + "]\n"); 
    }
    Detection<T> result(mDetector->getDetection(it));
    return result;
}

void PBMFLib::SingleChannel::initializeDetector(pybind11::module &m)
{
    pybind11::class_<PBMFLib::SingleChannel::Detector<double>> ddetector(m, "Detector");
    ddetector.def(pybind11::init<> ());
    ddetector.doc() = "Utilities for computing detections from the matched filtered signals.";
    /// Initialize
    ddetector.def("initialize",
                  &PBMFLib::SingleChannel::Detector<double>::initialize,
                  "Defines the detector behavior and requested metrics associated with a detector.\n");
    ddetector.def("detect",
                  &PBMFLib::SingleChannel::Detector<double>::detect,
                  "Computes the detections from the matched filtered signals.  The matched filtered class must have already computed the matched filtered signals.");
    ddetector.def("get_number_of_detections",
                  &PBMFLib::SingleChannel::Detector<double>::getNumberOfDetections, 
                  "Gets the number of detections found by the detector.");
    ddetector.def("get_detection",
                  &PBMFLib::SingleChannel::Detector<double>::getDetection,
                  "Gets the id'th detection computed by the detector.");
}

/// Template instatiation
template class PBMFLib::SingleChannel::Detector<double>; 
template class PBMFLib::SingleChannel::Detector<float>;
