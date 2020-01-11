#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <random>
#include <omp.h>
#include <chrono>
#include "mflib/waveformTemplate.hpp"
#include "mflib/singleChannel/matchedFilter.hpp"
#include "mflib/singleChannel/matchedFilterParameters.hpp"
#include <gtest/gtest.h>

namespace
{
using namespace MFLib;

TEST(singleChannelMatchedFilter, parameters)
{
    SingleChannel::MatchedFilterParameters parms;
    EXPECT_EQ(parms.getDetectionMode(),
              MatchedFilterDetectionMode::SINGLE_CHANNEL);
}

TEST(singleChannelMatchedFilter, testDouble)
{
    SingleChannel::MatchedFilterParameters parms;
    SingleChannel::MatchedFilter<double> mf;
}


}
