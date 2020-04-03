#include <cstdio>
#include <cstdlib>
#include "mflib/singleChannel/associatorParameters.hpp"
#include "singleChannel/pyAssociatorParameters.hpp"

using namespace PBMFLib::SingleChannel;

/// Constructor
AssociatorParameters::AssociatorParameters() :
    mParameters(std::make_unique<MFLib::SingleChannel::AssociatorParameters> ())
{
}

AssociatorParameters::AssociatorParameters(
    const AssociatorParameters &parms)
{
    *this = parms;
}

AssociatorParameters& AssociatorParameters::operator=(
    const AssociatorParameters &parms)
{
    if (&parms == this){return *this;}
    mParameters
        = std::make_unique<MFLib::SingleChannel::AssociatorParameters>
          (*parms.mParameters);
    return *this;
}

/// Destructor
AssociatorParameters::~AssociatorParameters() = default;

/// Origin time tolerance
void AssociatorParameters::setOriginTimeTolerance(const double tolerance)
{
    if (tolerance < 0)
    {
        throw std::invalid_argument("tolerance = " + std::to_string(tolerance)
                                  + " must be positive\n");
    }
    mParameters->setOriginTimeTolerance(tolerance);
}

double AssociatorParameters::getOriginTimeTolerance() const noexcept
{
    return mParameters->getOriginTimeTolerance();
}

/// Minimum number of picks to create an event
void AssociatorParameters::setMinimumNumberOfPicksInEvent(const int minPicks)
{
    if (minPicks < 1)
    {
        throw std::invalid_argument("Min number of picks must be positive\n");
    }
    mParameters->setMinimumNumberOfPicksInEvent(minPicks);
}

int AssociatorParameters::getMinimumNumberOfPicksInEvent() const noexcept
{
    return mParameters->getMinimumNumberOfPicksInEvent();
}

/// Weighting
void AssociatorParameters::enableCorrelationCoefficientWeighting() noexcept
{
    mParameters->enableCorrelationCoefficientWeighting();
}

void AssociatorParameters::disableCorrelationCoefficientWeighting() noexcept
{
    mParameters->disableCorrelationCoefficientWeighting();
}

bool AssociatorParameters::useCorrelationCoefficientWeighting() const noexcept
{
   return mParameters->useCorrelationCoefficientWeighting();
}

/// Clears the class
void AssociatorParameters::clear() noexcept
{
    mParameters->clear();
}

/// Get the native class
MFLib::SingleChannel::AssociatorParameters
AssociatorParameters::getNativeClass() const
{
    auto parmsOut(*mParameters);
    return parmsOut;
}

/// Initialize the class
void PBMFLib::SingleChannel::initializeAssociatorParameters(pybind11::module &m)
{
    pybind11::class_<PBMFLib::SingleChannel::AssociatorParameters>
        ascParameters(m, "AssociatorParameters");
    /// The default constructor:
    ascParameters.def(pybind11::init<> ());
    ascParameters.doc() = "Defines the options class for DBSCAN-based association.";
    // Origin time tolerance
    ascParameters.def("set_origin_time_tolerance",
                      &PBMFLib::SingleChannel::AssociatorParameters::setOriginTimeTolerance,
                         "This defines the tolerance in seconds that all picks must be within to be considered a cluster.");
    ascParameters.def("get_origin_time_tolerance",
                      &PBMFLib::SingleChannel::AssociatorParameters::getOriginTimeTolerance,
                      "The origin time tolerance for the clustering in seconds.");
    // Min picks
    ascParameters.def("set_minimum_number_of_picks_in_event",
                      &PBMFLib::SingleChannel::AssociatorParameters::setMinimumNumberOfPicksInEvent,
                      "Sets the minimum number of picks required when creating a cluster origin times.");
    ascParameters.def("get_minimum_number_of_picks_in_event",
                      &PBMFLib::SingleChannel::AssociatorParameters::getMinimumNumberOfPicksInEvent,
                     "Gets the minimum number of picks required when creating a cluster origin times.");
    // Weighting
    ascParameters.def("enable_correlation_coefficient_weighting",
                      &PBMFLib::SingleChannel::AssociatorParameters::enableCorrelationCoefficientWeighting,
                      "This will cause the algorithm to use the correlation coefficients as weights in the clustering.");
    ascParameters.def("disable_correlation_coefficient_weighting",
                      &PBMFLib::SingleChannel::AssociatorParameters::disableCorrelationCoefficientWeighting,
                      "This will cause the algorithm to use treat all origin times derived from detections equally.");
    ascParameters.def("use_correlation_coefficient_weighting",
                      &PBMFLib::SingleChannel::AssociatorParameters::useCorrelationCoefficientWeighting,
                      "True indicates that the algorithm will use correlation-based weights during the clustering.");
    /// Clears the class
    ascParameters.def("clear",
                      &PBMFLib::SingleChannel::AssociatorParameters::clear,
                      "Clears the class's memory and resets the class.");
}
