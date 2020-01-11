#include <cstdio>
#include <cstdlib>
#include <string>
#include "mflib/singleChannel/matchedFilterParameters.hpp"
#include "mflib/matchedFilterParameters.hpp"
#include "mflib/waveformTemplate.hpp"
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

/// Adds a template
void MatchedFilterParameters::addTemplate(
    const MFLib::WaveformTemplate &waveformTemplate)
{
    pImpl->mParameters.addTemplate(waveformTemplate);
}

/// Clears the templates
void MatchedFilterParameters::clearTemplates() noexcept
{
    pImpl->mParameters.clearTemplates();
}

/// Gets the it'th waveform template
MFLib::WaveformTemplate MatchedFilterParameters::getTemplate(const int it) const
{
    return pImpl->mParameters.getTemplate(it);
}

/// Gets the number of templates
int MatchedFilterParameters::getNumberOfTemplates() const noexcept
{
    return pImpl->mParameters.getNumberOfTemplates();
}

/// Gets the max template length
int MatchedFilterParameters::getMaxTemplateLength() const
{
    return pImpl->mParameters.getMaxTemplateLength();
}

/// Set/get signal size
void MatchedFilterParameters::setSignalSize(const int nx)
{
    pImpl->mParameters.setSignalSize(nx);
}

int MatchedFilterParameters::getSignalSize() const
{
    return pImpl->mParameters.getSignalSize();
}

/// Get/set FFT length
void MatchedFilterParameters::setFFTLength(const int fftLength)
{
    return pImpl->mParameters.setFFTLength(fftLength);
}

int MatchedFilterParameters::getFFTLength() const
{
    return pImpl->mParameters.getFFTLength();
}

int MatchedFilterParameters::getBlockLength() const
{
    return pImpl->mParameters.getBlockLength();
}

/// Get the detection mode
MFLib::MatchedFilterDetectionMode
MatchedFilterParameters::getDetectionMode() const noexcept
{
    return pImpl->mParameters.getDetectionMode();
}

/// Defines the matched filter implementation
void MatchedFilterParameters::setMatchedFilterImplementation(
    const MFLib::MatchedFilterImplementation implementation) noexcept
{
    pImpl->mParameters.setMatchedFilterImplementation(implementation);
}

MFLib::MatchedFilterImplementation
MatchedFilterParameters::getMatchedFilterImplementation() const noexcept
{
    return pImpl->mParameters.getMatchedFilterImplementation();
}

/// Determines if the class is valid
bool MatchedFilterParameters::isValid() const noexcept
{
    return pImpl->mParameters.isValid();
}

/// Gets a copy of the underlying matched filters parameter class
MFLib::MatchedFilterParameters
MatchedFilterParameters::getParameters() const
{
    return pImpl->mParameters;
}
