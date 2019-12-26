#include <cstdio>
#include <cstdlib>
#include "mflib/waveformTemplate.hpp"
#include "mflib/matchedFilterParameters.hpp"
#include "waveformTemplate.hpp"
#include "matchedFilterParameters.hpp"

using namespace PBMFLib;

/// Constructor
MatchedFilterParameters::MatchedFilterParameters() :
    mParameters(std::make_unique<MFLib::MatchedFilterParameters> ())
{
}

MatchedFilterParameters::MatchedFilterParameters(
    const MatchedFilterParameters &parms)
{
    *this = parms;
}

MatchedFilterParameters& 
MatchedFilterParameters::operator=(const MatchedFilterParameters &parms)
{
    if (&parms == this){return *this;}
    mParameters
        = std::make_unique<MFLib::MatchedFilterParameters> (*parms.mParameters);
    return *this;
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

WaveformTemplate MatchedFilterParameters::getTemplate(const int it) const
{
    if (it < 0 || it >= getNumberOfTemplates())
    {
        throw std::invalid_argument("template index = "
                                  + std::to_string(it)
                                  + " must be in range [0,"
                                  + std::to_string(getNumberOfTemplates())
                                  + "]\n");
    }
    auto wt = mParameters->getTemplate(it);
    WaveformTemplate wtOut(wt);
    return wtOut;
}

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

/// Sets/gets the signals size
void MatchedFilterParameters::setSignalSize(const int n)
{
    if (n < 1)
    {
        throw std::invalid_argument("Signal size must be positive\n");
    }
    mParameters->setSignalSize(n);     
}

int MatchedFilterParameters::getSignalSize() const
{
    if (getNumberOfTemplates() < 1){return 0;}
    return mParameters->getSignalSize();
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

