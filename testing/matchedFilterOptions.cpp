#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
#include "mflib/matchedFilterOptions.hpp"
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
    MatchedFilterOptions<double> options;
    int signalSize = 1024;
    int nSamples = 100;
    double error = 0;
    std::vector<double> tp1(nSamples);
    std::vector<double> tp2(nSamples+1);
    std::vector<double> work;
    for (int i=0; i<static_cast<int> (tp1.size()); ++i)
    {
        tp1[i] = i;
    }
    for (int i=0; i<static_cast<int> (tp2.size()); ++i)
    {
        tp2[i] = i;
    }
    options.setMatchedFilterImplementation(
       MatchedFilterImplementation::AUTO);
    EXPECT_EQ(options.getMatchedFilterImplementation(),
              MatchedFilterImplementation::AUTO);
    options.setSignalSize(signalSize);
    EXPECT_EQ(options.getSignalSize(), signalSize);
 
    EXPECT_FALSE(options.isValid());
    EXPECT_NO_THROW(options.addTemplate(tp1.size(), tp1.data()));
    EXPECT_EQ(options.getTemplateLength(), nSamples);
    EXPECT_EQ(options.getNumberOfTemplates(), 1);
    EXPECT_TRUE(options.isValid());

    EXPECT_NO_THROW(options.addTemplate(tp2.size(), tp2.data()));
    EXPECT_EQ(options.getTemplateLength(), nSamples+1);
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
    EXPECT_NO_THROW(work = options.getTemplate(0));
    error = 0;
    for (int i=0; i<static_cast<int> (tp1.size()); ++i)
    {
        EXPECT_NEAR(tp1[i], work[i], 1.e-14);
    }
    EXPECT_NEAR(work[nSamples], 0, 1.e-14);

    EXPECT_NO_THROW(work = options.getTemplate(1));
    error = 0;
    for (int i=0; i<static_cast<int> (tp2.size()); ++i)
    {   
        EXPECT_NEAR(tp2[i], work[i], 1.e-14);
    }
 
    // Test copy c'tor
    MatchedFilterOptions<double> optionsCopy(options); 
    EXPECT_EQ(optionsCopy.getMatchedFilterImplementation(),
              MatchedFilterImplementation::AUTO);
    EXPECT_EQ(optionsCopy.getSignalSize(), signalSize);
    EXPECT_TRUE(optionsCopy.isValid());
    EXPECT_NO_THROW(work = optionsCopy.getTemplate(0));
    error = 0;
    for (int i=0; i<static_cast<int> (tp1.size()); ++i)
    {   
        EXPECT_NEAR(tp1[i], work[i], 1.e-14);
    }   
    EXPECT_NEAR(work[nSamples], 0, 1.e-14);

    EXPECT_NO_THROW(work = options.getTemplate(1));
    error = 0;
    for (int i=0; i<static_cast<int> (tp2.size()); ++i)
    {
        EXPECT_NEAR(tp2[i], work[i], 1.e-14);
    }

    // Test clear
    EXPECT_NO_THROW(options.clearTemplates());
    EXPECT_EQ(options.getNumberOfTemplates(), 0);
    options.clear();
    EXPECT_FALSE(options.isValid());
}

}
