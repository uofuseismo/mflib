#ifndef PYMFLIB_SINGLECHANNEL_PYASSOCIATORPARAMETERS_HPP
#define PYMFLIB_SINGLECHANNEL_PYASSOCIATORPARAMETERS_HPP
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "mflib/singleChannel/associatorParameters.hpp"

namespace MFLib
{
namespace SingleChannel
{
class AssociatorParameters;
}
}

namespace PBMFLib
{
namespace SingleChannel
{
class AssociatorParameters
{
public:
    /// Constructor
    AssociatorParameters();
    AssociatorParameters(const AssociatorParameters &parms);
    AssociatorParameters& operator=(const AssociatorParameters &parms);
    /// Destructor
    ~AssociatorParameters();
    void clear() noexcept;
    // Origin time tolerance.
    void setOriginTimeTolerance(double tolerance);
    double getOriginTimeTolerance() const noexcept;
    // Minimum number of picks to create an event
    void setMinimumNumberOfPicksInEvent(int minPicks);
    int getMinimumNumberOfPicksInEvent() const noexcept;
    // Weighting
    void enableCorrelationCoefficientWeighting() noexcept;
    void disableCorrelationCoefficientWeighting() noexcept;
    bool useCorrelationCoefficientWeighting() const noexcept;
    /// Gets the native parameters class
    MFLib::SingleChannel::AssociatorParameters getNativeClass() const;
private:
    std::unique_ptr<MFLib::SingleChannel::AssociatorParameters> mParameters;
};
void initializeAssociatorParameters(pybind11::module &m);
}
}
#endif
