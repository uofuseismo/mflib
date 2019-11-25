#include <cstdio>
#include <cstdlib>
#include <array>
#include "mflib/waveformTemplate.hpp"
#include <gtest/gtest.h>
namespace
{
using namespace MFLib;
TEST(waveformTemplate, basic)
{
    WaveformTemplate tplate;
    double df = 100;
    EXPECT_NO_THROW(tplate.setSamplingRate(df));
    EXPECT_EQ(tplate.getSamplingRate(), df);

    double weight = 0.5;
    EXPECT_EQ(tplate.getShiftAndStackWeight(), 1); // Test default
    EXPECT_NO_THROW(tplate.setShiftAndStackWeight(weight));
    EXPECT_EQ(tplate.getShiftAndStackWeight(), weight);

    EXPECT_FALSE(tplate.haveSignal());
    std::array<double, 10> x8({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    int npts = static_cast<int> (x8.size());
    EXPECT_NO_THROW(tplate.setSignal(npts, x8.data()));
    EXPECT_TRUE(tplate.haveSignal());
    EXPECT_EQ(tplate.getSignalLength(), npts);
    // Recover the signals
    std::fill(x8.begin(), x8.end(), 0);
    double *x8Ptr = x8.data();
    EXPECT_NO_THROW(tplate.getSignal(npts, &x8Ptr));
    for (int i=0; i<x8.size(); ++i)
    {
        EXPECT_NEAR(x8[i], i+1, 1.e-14);
    }

    // Repeat for float precision
    std::array<float, 10> x4({-1, -2, -3, -4, -5, -6, -7, -8, -9, -10});
    EXPECT_NO_THROW(tplate.setSignal(npts, x4.data()));
    std::fill(x4.begin(), x4.end(), 0);
    float *x4Ptr = x4.data();
    EXPECT_NO_THROW(tplate.getSignal(npts, &x4Ptr));
    for (int i=0; i<x4.size(); ++i)
    {
        EXPECT_NEAR(x4[i], -(i+1), 1.e-7);
    }
   
    // Test copy c'tor 
    WaveformTemplate tplateCopy(tplate); 
    EXPECT_EQ(tplateCopy.getShiftAndStackWeight(), weight);
    EXPECT_EQ(tplateCopy.getSamplingRate(), df);
    EXPECT_EQ(tplateCopy.getSignalLength(), npts);
    EXPECT_NO_THROW(tplateCopy.getSignal(npts, &x8Ptr)); 
    for (int i=0; i<x8.size(); ++i)
    {   
        EXPECT_NEAR(x8[i], -(i+1), 1.e-14);
    }
}

}
