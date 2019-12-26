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
    mOptions(std::make_unique<MFLib::MatchedFilterParameters> ())
{
}

/// Destructor
MatchedFilterParameters::~MatchedFilterParameters() = default;

/// Adds a template
void MatchedFilterParameters::addTemplate(const WaveformTemplate &wt)
{
    auto t = wt.getNativeClass();
    mOptions->addTemplate(t);
}

/// Gets a template


/// Clears the templates from the class
void MatchedFilterParameters::clearTemplates() noexcept
{
    mOptions->clearTemplates();
}

/// Clears the class
void MatchedFilterParameters::clear() noexcept
{
    mOptions->clear();
}

