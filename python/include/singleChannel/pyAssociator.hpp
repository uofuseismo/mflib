#ifndef PYMFLIB_SINGLECHANNEL_PYASSOCIATOR_HPP
#define PYMFLIB_SINGLECHANNEL_PYASSOCIATOR_HPP
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "mflib/singleChannel/associator.hpp"

namespace MFLib
{
namespace SingleChannel
{
template<class T> class Associator;
}
}

namespace PBMFLib
{
namespace SingleChannel
{
template<class T> class Detection;
class AssociatorParameters;
template<class T>
class Associator
{
public:
    /// Constructors
    Associator();
    Associator(const Associator &associator);
    Associator(const MFLib::SingleChannel::Associator<T> &associator);
    /// Operators
    Associator& operator=(const Associator &associator);
    Associator& operator=(const MFLib::SingleChannel::Associator<T> &associator);
    /// Destructors
    ~Associator();
    void clear() noexcept;
    /// Inititialize
    void initialize(const AssociatorParameters &parameters);
    bool isInitialized() const noexcept;
    /// Adds a detection
    void addDetection(const Detection<T> &detection); 
    int getNumberOfDetections() const noexcept;
    void clearDetections() noexcept;
    /// Does the association
    void associate();
    /// Gets the number of events
    int getNumberOfEvents() const noexcept;
    /// Gets the detections for ane vent
    std::vector<PBMFLib::SingleChannel::Detection<T>>
        getDetectionsInEvent(int iev) const; 
private:
    
    std::unique_ptr<MFLib::SingleChannel::Associator<T>> mAssociator;
};
void initializeAssociator(pybind11::module &m);;
}
}
#endif
