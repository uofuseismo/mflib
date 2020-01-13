#include <cstdio>
#include <cstdlib>
#include "mflib/enums.hpp"
#include "mflib/waveformTemplate.hpp"
#include "mflib/matchedFilterParameters.hpp"
#include "pyWaveformTemplate.hpp"
#include "pyMultiChannelMatchedFilterParameters.hpp"

using namespace PBMFLib;

/// Constructor
MultiChannelMatchedFilterParameters::MultiChannelMatchedFilterParameters() :
    mParameters(std::make_unique<MFLib::MatchedFilterParameters> ())
{
    mParameters->setDetectionMode(
        MFLib::MatchedFilterDetectionMode::MULTI_CHANNEL);
}

MultiChannelMatchedFilterParameters::MultiChannelMatchedFilterParameters(
    const MultiChannelMatchedFilterParameters &parms)
{
    *this = parms;
}

MultiChannelMatchedFilterParameters& 
MultiChannelMatchedFilterParameters::operator=(
    const MultiChannelMatchedFilterParameters &parms)
{
    if (&parms == this){return *this;}
    mParameters
        = std::make_unique<MFLib::MatchedFilterParameters> (*parms.mParameters);
    return *this;
}

/// Destructor
MultiChannelMatchedFilterParameters::~MultiChannelMatchedFilterParameters()
    = default;

/// Adds a template
void MultiChannelMatchedFilterParameters::addTemplate(
    const WaveformTemplate &wt)
{
    auto t = wt.getNativeClass();
    mParameters->addTemplate(t);
}

int MultiChannelMatchedFilterParameters::getNumberOfTemplates() const
{
    return mParameters->getNumberOfTemplates();
}

WaveformTemplate 
MultiChannelMatchedFilterParameters::getTemplate(const int it) const
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
void MultiChannelMatchedFilterParameters::setFFTLength(const int fftLength)
{
    auto nt = mParameters->getNumberOfTemplates();
    if (nt < 1)
    {
        throw std::runtime_error("Templates must be set prior to calling this\n");
    }
    mParameters->setFFTLength(fftLength);
}

int MultiChannelMatchedFilterParameters::getFFTLength() const
{
    auto nt = getNumberOfTemplates(); //mParameters->getNumberOfTemplates();
    if (nt < 1){return 0;}
    return mParameters->getFFTLength();
}

int MultiChannelMatchedFilterParameters::getBlockLength() const
{
    auto nt = getNumberOfTemplates(); //mParameters->getNumberOfTemplates();
    if (nt < 1){return 0;}
    return mParameters->getBlockLength();
}

/// Sets/gets the signals size
void MultiChannelMatchedFilterParameters::setSignalSize(const int n)
{
    if (n < 1)
    {
        throw std::invalid_argument("Signal size must be positive\n");
    }
    mParameters->setSignalSize(n);     
}

int MultiChannelMatchedFilterParameters::getSignalSize() const
{
    if (getNumberOfTemplates() < 1){return 0;}
    return mParameters->getSignalSize();
}

/// Clears the templates from the class
void MultiChannelMatchedFilterParameters::clearTemplates() noexcept
{
    mParameters->clearTemplates();
}

/// Clears the class
void MultiChannelMatchedFilterParameters::clear() noexcept
{
    mParameters->clear();
}

/// Get the native class
MFLib::MatchedFilterParameters
MultiChannelMatchedFilterParameters::getNativeClass() const
{
    auto parmsOut(*mParameters);
    return parmsOut;
}
