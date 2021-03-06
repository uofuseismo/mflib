#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
#include "mflib/waveformTemplate.hpp"
#include "mflib/matchedFilterParameters.hpp"
#include "utilities.hpp"
#include <gtest/gtest.h>

namespace
{
using namespace MFLib;

TEST(MatchedFilter, blockAndFFTLength)
{
    int nx = 75000;
    int nb = 4500;
    std::pair<int, int> result;
    EXPECT_NO_THROW(result = computeOptimalFFTAndBlockLength(nb, nx));
    int fftRef = 32768;
    int lRef = 28269;
    EXPECT_EQ(result.first,  fftRef);
    EXPECT_EQ(result.second, lRef);

    nx = 86400*100;
    nb = 4500;
    EXPECT_NO_THROW(result = computeOptimalFFTAndBlockLength(nb, nx));
    fftRef = 65536;
    lRef = 61037;
    EXPECT_EQ(result.first,  fftRef);
    EXPECT_EQ(result.second, lRef);
    //printf("%d, %d\n", result.first, result.second);
}
 
TEST(MatchedFilter, options)
{
    MatchedFilterParameters options;
    int signalSize = 1024;
    int nSamples = 100;
    double error = 0;
    std::vector<double> tp1(nSamples);
    std::vector<double> tp2(nSamples+1);
    std::vector<double> work;
    WaveformTemplate tplate1, tplate2;
    for (int i=0; i<static_cast<int> (tp1.size()); ++i)
    {
        tp1[i] = i;
    }
    for (int i=0; i<static_cast<int> (tp2.size()); ++i)
    {
        tp2[i] = -i;
    }
    tplate1.setSignal(tp1.size(), tp1.data());
    tplate2.setSignal(tp2.size(), tp2.data());
    tplate1.setSamplingRate(1);
    tplate2.setSamplingRate(1);
    options.setMatchedFilterImplementation(
       MatchedFilterImplementation::AUTO);
    EXPECT_EQ(options.getMatchedFilterImplementation(),
              MatchedFilterImplementation::AUTO);
    options.setSignalSize(signalSize);
    EXPECT_EQ(options.getSignalSize(), signalSize);

    options.setStackAbsoluteValues(false);
    EXPECT_FALSE(options.getStackAbsoluteValues());

    options.setDetectionMode(MatchedFilterDetectionMode::MULTI_CHANNEL);
    EXPECT_EQ(options.getDetectionMode(),
              MatchedFilterDetectionMode::MULTI_CHANNEL);
    options.setDetectionMode(MatchedFilterDetectionMode::SINGLE_CHANNEL);
 
    EXPECT_FALSE(options.isValid());
    EXPECT_NO_THROW(options.addTemplate(tplate1)); //tp1.size(), tp1.data()));
    EXPECT_EQ(options.getMaxTemplateLength(), nSamples);
    EXPECT_EQ(options.getNumberOfTemplates(), 1);
    EXPECT_TRUE(options.isValid());

    EXPECT_NO_THROW(options.addTemplate(tplate2)); //tp2.size(), tp2.data()));
    EXPECT_EQ(options.getMaxTemplateLength(), nSamples+1);
    EXPECT_EQ(options.getNumberOfTemplates(), 2);

    int fftLenRef = 512;
    int fftLen = 0;
    EXPECT_NO_THROW(fftLen = options.getFFTLength());
    int blockLenRef = 412;
    int blockLen = 0;
    EXPECT_NO_THROW(blockLen = options.getBlockLength());
    EXPECT_EQ(fftLen, fftLenRef);
    EXPECT_EQ(blockLen, blockLenRef);

    // Get the templates
    WaveformTemplate tplateOut;
    EXPECT_NO_THROW(tplateOut = options.getTemplate(0));
    std::vector<double> work8(options.getMaxTemplateLength());
    double *work8Ptr = work8.data();
    tplateOut.getSignal(work8.size(), &work8Ptr);
    EXPECT_EQ(static_cast<int> (tp1.size()), tplateOut.getSignalLength());
    error = infinityNorm(tp1.size(), tp1.data(), work8.data());
    EXPECT_NEAR(error, 0, 1.e-14);

    EXPECT_NO_THROW(tplateOut = options.getTemplate(1));
    tplateOut.getSignal(work8.size(), &work8Ptr);
    EXPECT_EQ(static_cast<int> (tp2.size()), tplateOut.getSignalLength());
    error = infinityNorm(tp2.size(), tp2.data(), work8.data());
    EXPECT_NEAR(error, 0, 1.e-14);

    // Test copy c'tor
    MatchedFilterParameters optionsCopy(options); 
    EXPECT_EQ(optionsCopy.getMatchedFilterImplementation(),
              MatchedFilterImplementation::AUTO);
    EXPECT_FALSE(optionsCopy.getStackAbsoluteValues());
    EXPECT_EQ(optionsCopy.getSignalSize(), signalSize);
    EXPECT_EQ(optionsCopy.getDetectionMode(),
              MatchedFilterDetectionMode::SINGLE_CHANNEL);
    EXPECT_TRUE(optionsCopy.isValid());
    EXPECT_NO_THROW(tplateOut = optionsCopy.getTemplate(0));
    tplateOut.getSignal(work8.size(), &work8Ptr);
    error = 0;
    for (int i=0; i<static_cast<int> (tp1.size()); ++i)
    {   
        EXPECT_NEAR(tp1[i], work8[i], 1.e-14);
    }   
    //EXPECT_NEAR(work[nSamples], 0, 1.e-14);

    EXPECT_NO_THROW(tplateOut = options.getTemplate(1));
    tplateOut.getSignal(work8.size(), &work8Ptr);
    error = 0;
    for (int i=0; i<static_cast<int> (tp2.size()); ++i)
    {
        EXPECT_NEAR(tp2[i], work8[i], 1.e-14);
    }

    // Test clear
    EXPECT_NO_THROW(options.clearTemplates());
    EXPECT_EQ(options.getNumberOfTemplates(), 0);
    options.clear();
    EXPECT_FALSE(options.isValid());
}

}
