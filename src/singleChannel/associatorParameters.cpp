#include <cstdio>
#include <cstdlib>
#include <string>
#include "mflib/singleChannel/associatorParameters.hpp"

using namespace MFLib::SingleChannel;

#define DEFAULT_OT_TOLERANCE 2
#define DEFAULT_MIN_PICKS 5

class AssociatorParameters::AssociatorParametersImpl
{
public:
    double mOriginTimeTolerance = DEFAULT_OT_TOLERANCE;
    int mMinNumberOfPicks = DEFAULT_MIN_PICKS;
};

/// Default constructor
AssociatorParameters::AssociatorParameters() :
    pImpl(std::make_unique<AssociatorParametersImpl> ())
{
}

/// Copy c'tor
AssociatorParameters::AssociatorParameters(
    const AssociatorParameters &parameters)
{
    *this = parameters;
}

/// Move c'tor
AssociatorParameters::AssociatorParameters(
    AssociatorParameters &&parameters) noexcept
{
    *this = std::move(parameters);
}

/// Copy assignment operator
AssociatorParameters& 
AssociatorParameters::operator=(const AssociatorParameters &parameters)
{
    if (&parameters == this){return *this;}
    pImpl = std::make_unique<AssociatorParametersImpl> (*parameters.pImpl);
    return *this;
}

/// Move assignment operator
AssociatorParameters&
AssociatorParameters::operator=(AssociatorParameters &&parameters) noexcept
{
    if (&parameters == this){return *this;}
    pImpl = std::move(parameters.pImpl);
    return *this;
}

/// Destructor
AssociatorParameters::~AssociatorParameters() = default;

/// Release memory and clear class
void AssociatorParameters::clear() noexcept
{
    pImpl->mOriginTimeTolerance = DEFAULT_OT_TOLERANCE;
    pImpl->mMinNumberOfPicks = DEFAULT_MIN_PICKS;
}

/// Sets the origin time tolerance
void AssociatorParameters::setOriginTimeTolerance(const double tol)
{
    if (tol < 0)
    {
        throw std::invalid_argument("tolerance = " + std::to_string(tol)
                                  + " must be positive\n");
    }
    pImpl->mOriginTimeTolerance = tol;
}

double AssociatorParameters::getOriginTimeTolerance() const noexcept
{
    return pImpl->mOriginTimeTolerance;
}

/// Sets the minimum number of picks required to associate
void AssociatorParameters::setMinimumNumberOfPicksInEvent(const int minPicks)
{
    if (minPicks < 1)
    {
        throw std::invalid_argument("minPicks = " + std::to_string(minPicks)
                                  + " must be positive\n");
    } 
    pImpl->mMinNumberOfPicks = minPicks;
}

int AssociatorParameters::getMinimumNumberOfPicksInEvent() const noexcept
{
    return pImpl->mMinNumberOfPicks;
}

/*
bool AssociatorParameters::haveOriginTimeTolerance() const noexcept
{
    return pImpl->mOriginTimeTolerance >= 0 ? true : false;
}
*/
