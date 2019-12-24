#include <cstdio>
#include <cstdlib>
#include "mflib/waveformTemplate.hpp"
#include "mflib/matchedFilterOptions.hpp"
#include "modules.hpp"
#include "waveformTemplate.hpp"
#include "matchedFilterOptions.hpp"

using namespace PBMFLib;

/// Constructor
MatchedFilterOptions::MatchedFilterOptions() :
    mOptions(std::make_unique<MFLib::MatchedFilterOptions> ())
{
}

/// Destructor
MatchedFilterOptions::~MatchedFilterOptions() = default;

/// Adds a template
void MatchedFilterOptions::addTemplate(const WaveformTemplate &wt)
{
    auto t = wt.getNativeClass();
    mOptions->addTemplate(t);
}

/// Gets a template


/// Clears the templates from the class
void MatchedFilterOptions::clearTemplates() noexcept
{
    mOptions->clearTemplates();
}

/// Clears the class
void MatchedFilterOptions::clear() noexcept
{
    mOptions->clear();
}

