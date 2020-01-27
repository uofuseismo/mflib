#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include "mflib/singleChannel/detectorParameters.hpp"
#include <gtest/gtest.h>

namespace
{
using namespace MFLib::SingleChannel;

TEST(singleChannelDetection, detectorParameters)
{
    DetectorParameters parms;
    int minSamples = 50;
    double thresh = 0.8;
    auto policy = MFLib::MaximumMatchedFilterPolicy::ABSOLUTE_MAXIMUM; 
    EXPECT_NO_THROW(parms.setMinimumDetectionSpacing(minSamples));
    EXPECT_NO_THROW(parms.setDetectionThreshold(thresh));
    EXPECT_NO_THROW(parms.setMaximaPolicy(policy));

    EXPECT_EQ(parms.getMinimumDetectionSpacing(), minSamples);
    EXPECT_NEAR(parms.getDetectionThreshold(), thresh, 1.e-14);
    EXPECT_EQ(parms.getMaximaPolicy(), policy);

    DetectorParameters parmsCopy(parms);
    EXPECT_EQ(parmsCopy.getMinimumDetectionSpacing(), minSamples);
    EXPECT_NEAR(parmsCopy.getDetectionThreshold(), thresh, 1.e-14);
    EXPECT_EQ(parmsCopy.getMaximaPolicy(), policy);

}

}
