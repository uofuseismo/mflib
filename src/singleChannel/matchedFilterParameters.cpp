#include <cstdio>
#include <cstdlib>
#include <string>
#include "mflib/singleChannelMatchedFilterParameters.hpp"
#include "mflib/matchedFilterParameters.hpp"
#include "mflib/enums.hpp"

using namespace MFLib::SingleChannel;

class MatchedFilterParameters::MatchedFilterParametersImpl
{
public:
    MFLib::MatchedFilterParameters mParameters;
};

/// Constructor
MatchedFilterParameters::MatchedFilterParameters() :
    pImpl(std::make_unique<MatchedFilterParametersImpl> ())
{
    pImpl->mParameters.setDetectionMode(
        MatchedFilterDetectionMode::SINGLE_CHANNEL);
}

/// Copy constructor
MatchedFilterParameters::MatchedFilterParameters(
    const MatchedFilterParameters &parms)
{
    *this = parms;
}

/// Move constructor
MatchedFilterParameters::MatchedFilterParameters(
    MatchedFilterParameters &&parms) noexcept
{
    *this = parms;
}

/// Copy assignment
MatchedFilterParameters& MatchedFilterParameters::operator=(
    const MatchedFilterParameters &parms)
{
    if (&parms == this){return *this;}
    pImpl = std::make_unique<MatchedFilterParametersImpl> (*parms.pImpl);
    return *this;
} 

/// Move assignment
MatchedFilterParameters& MatchedFilterParameters::operator=(
    MatchedFilterParameters &&parms) noexcept
{
    if (&parms == this){return *this;}
    pImpl = std::move(parms.pImpl);
    return *this;
}

/// Destructor
MatchedFilterParameters::~MatchedFilterParameters() = default;

/// Clear
void MatchedFilterParameters::clear() noexcept
{
    pImpl->mParameters.clear();
    pImpl->mParameters.setDetectionMode(
        MatchedFilterDetectionMode::SINGLE_CHANNEL);
}

