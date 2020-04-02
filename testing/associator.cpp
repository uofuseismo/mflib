#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <random>
#include <vector>
#include <chrono>
#include <random>
#include "mflib/waveformTemplate.hpp"
#include "mflib/singleChannel/associatorParameters.hpp"
#include "utilities.hpp"
#include "private/weightedStatistics.hpp"
#include <gtest/gtest.h>

namespace
{
using namespace MFLib;

TEST(singleChannelAssociator, parameters)
{
    SingleChannel::AssociatorParameters parms;
    
    EXPECT_NO_THROW(parms.setOriginTimeTolerance(4));
    EXPECT_NO_THROW(parms.setMinimumNumberOfPicksInEvent(8));
    EXPECT_NEAR(parms.getOriginTimeTolerance(), 4, 1.e-14);
    EXPECT_EQ(parms.getMinimumNumberOfPicksInEvent(), 8);

    SingleChannel::AssociatorParameters parmsCopy(parms);
    EXPECT_NEAR(parmsCopy.getOriginTimeTolerance(), 4, 1.e-14);
    EXPECT_EQ(parmsCopy.getMinimumNumberOfPicksInEvent(), 8);
}

TEST(singleChannelAssociator, weightedStatistics)
{
    // https://rdrr.io/cran/spatstat/man/weighted.median.html
    // library(spatstat)
    // x <- c(1.1, 5.3, 3.7, 2.1, 7.0, 9.9)
    // y <- c(1.1, 0.4, 2.1, 3.5, 1.2, 0.8)
    // weighted.median(x, y)
    std::vector<double> x({1.1, 5.3, 3.7, 2.1, 7.0, 9.9});
    std::vector<double> wts({1.1, 0.4, 2.1, 3.5, 1.2, 0.8});
    auto wm = weightedMedian(x.size(), x.data(), wts.data());
    EXPECT_NEAR(wm, 2.085714, 1.e-6);
    // xpost <- c(0.1,0.35,0.05,0.1,0.15,0.05,0.2)
    // weighted.median(xpost, xpost) 
    std::vector<double> xpost({0.1,0.35,0.05,0.1,0.15,0.05,0.2});
    std::vector<double> wpost(xpost);
    wm = weightedMedian(xpost.size(), xpost.data(), wpost.data());
    EXPECT_NEAR(wm, 0.1625, 1.e-6);
    // Do an even length example
    wm = weightedMedian(xpost.size() - 1, xpost.data(), wpost.data());
    EXPECT_NEAR(wm, 0.1333333, 1.e-6);
    // And test an if statement conditional
    std::vector<double> x2({1, 2});
    std::vector<double> w2({1, 4});
    wm = weightedMedian(x2.size(), x2.data(), w2.data());
    EXPECT_NEAR(wm, 1.375, 1.e-6);
    // And an edge case of n=1
    EXPECT_NEAR(x[0], weightedMedian(1, x.data(), wts.data()), 1.e-6); 
    //------------------------------------------------------------------------//
    // Regular median
    wm = median(x.size(), x.data());
    EXPECT_NEAR(wm, 4.5, 1.e-6);
    wm = median(x.size()-1, x.data());
    EXPECT_NEAR(wm, 3.7, 1.e-6);
    EXPECT_NEAR(x[0], median(1, x.data()), 1.e-12);
    //------------------------------------------------------------------------//
    // Weighted mean and mean
    wm = weightedMean(xpost.size(), xpost.data(), xpost.data());
    EXPECT_NEAR(wm, 0.21, 1.e-6);
    EXPECT_NEAR(weightedMean(1, x.data(), x.data()), x[0], 1.e-12);
    wm = mean(xpost.size(), xpost.data());
    EXPECT_NEAR(wm, 0.1428571, 1.e-6);
    EXPECT_NEAR(mean(1, x.data()), x[0], 1.e-12); 
}

TEST(singleChannelAssociator, associator)
{

}

}
