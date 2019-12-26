#include <cstdio>
#include <cstdlib>
#include "mflib/waveformTemplate.hpp"
#include "mflib/matchedFilterParameters.hpp"
#include "modules.hpp"
#include "waveformTemplate.hpp"
#include "matchedFilterParameters.hpp"

using namespace PBMFLib;

/// Constructor
MatchedFilterParameters::MatchedFilterParameters() :
    mParameters(std::make_unique<MFLib::MatchedFilterParameters> ())
{
}

/// Destructor
MatchedFilterParameters::~MatchedFilterParameters() = default;

/// Adds a template
void MatchedFilterParameters::addTemplate(const WaveformTemplate &wt)
{
    auto t = wt.getNativeClass();
    mParameters->addTemplate(t);
}

int MatchedFilterParameters::getNumberOfTemplates() const
{
    return mParameters->getNumberOfTemplates();
}

/*
WaveformTemplate MatchedFilterParameters::getTemplate(const int it)
{
    
}
*/

/// FFT Length 
void MatchedFilterParameters::setFFTLength(const int fftLength)
{
    auto nt = mParameters->getNumberOfTemplates();
    if (nt < 1)
    {
        throw std::runtime_error("Templates must be set prior to calling this\n");
    }
    mParameters->setFFTLength(fftLength);
}

int MatchedFilterParameters::getFFTLength() const
{
    auto nt = getNumberOfTemplates(); //mParameters->getNumberOfTemplates();
    if (nt < 1){return 0;}
    return mParameters->getFFTLength();
}

int MatchedFilterParameters::getBlockLength() const
{
    auto nt = getNumberOfTemplates(); //mParameters->getNumberOfTemplates();
    if (nt < 1){return 0;}
    return mParameters->getBlockLength();
}

/// Clears the templates from the class
void MatchedFilterParameters::clearTemplates() noexcept
{
    mParameters->clearTemplates();
}

/// Clears the class
void MatchedFilterParameters::clear() noexcept
{
    mParameters->clear();
}

