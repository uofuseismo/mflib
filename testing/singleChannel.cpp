#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <random>
#include <vector>
#include <omp.h>
#include <chrono>
#include <random>
#include "mflib/waveformTemplate.hpp"
#include "mflib/singleChannel/matchedFilter.hpp"
#include "mflib/singleChannel/matchedFilterParameters.hpp"
#include <gtest/gtest.h>

namespace
{
using namespace MFLib;

/// Generates uniform random numbers
std::vector<double> generateUniformRandomNumbers(
   int n, double lower =-1, double upper=+1);
/// Computes the infinity norm: max(|r|)
template<class T> T infinityNorm(int n, const T *r);
/// computes the inifnity norm: max(|x - y|)
template<class T> T infinityNorm(int n, const T *x, const T *y);


TEST(singleChannelMatchedFilter, parameters)
{
    int signalSize = 3823;
    int nSamples = 410;
    // Initialize the random number generator
    SingleChannel::MatchedFilterParameters parms;

    auto tp1 = generateUniformRandomNumbers(nSamples);
    auto tp2 = generateUniformRandomNumbers(nSamples+1);
    std::vector<double> work;
    WaveformTemplate tplate1, tplate2;

    tplate1.setSignal(tp1.size(), tp1.data());
    tplate2.setSignal(tp2.size(), tp2.data());
    tplate1.setSamplingRate(1);
    tplate2.setSamplingRate(1);
    parms.setMatchedFilterImplementation(
       MatchedFilterImplementation::AUTO);
    EXPECT_EQ(parms.getMatchedFilterImplementation(),
              MatchedFilterImplementation::AUTO);
    parms.setSignalSize(signalSize);
    EXPECT_EQ(parms.getSignalSize(), signalSize);
    EXPECT_EQ(parms.getDetectionMode(),
              MatchedFilterDetectionMode::SINGLE_CHANNEL);

    EXPECT_FALSE(parms.isValid());
    EXPECT_NO_THROW(parms.addTemplate(tplate1)); //tp1.size(), tp1.data()));
    EXPECT_EQ(parms.getMaxTemplateLength(), nSamples);
    EXPECT_EQ(parms.getNumberOfTemplates(), 1);
    EXPECT_TRUE(parms.isValid());

    EXPECT_NO_THROW(parms.addTemplate(tplate2)); //tp1.size(), tp1.data()));
    EXPECT_EQ(parms.getMaxTemplateLength(), nSamples+1);
    EXPECT_EQ(parms.getNumberOfTemplates(), 2); 

    // Get the templates
    WaveformTemplate tplateOut;
    EXPECT_NO_THROW(tplateOut = parms.getTemplate(0));
    std::vector<double> work8(parms.getMaxTemplateLength());
    double *work8Ptr = work8.data();
    tplateOut.getSignal(work8.size(), &work8Ptr);
    EXPECT_EQ(static_cast<int> (tp1.size()), tplateOut.getSignalLength());
    auto error = infinityNorm(tp1.size(), tp1.data(), work8Ptr);
    EXPECT_NEAR(error, 0, 1.e-14);

    EXPECT_NO_THROW(tplateOut = parms.getTemplate(1));
    tplateOut.getSignal(work8.size(), &work8Ptr);
    EXPECT_EQ(static_cast<int> (tp2.size()), tplateOut.getSignalLength());
    error = infinityNorm(tp2.size(), tp2.data(), work8Ptr);
    EXPECT_NEAR(error, 0, 1.e-14);

    // Test copy c'tor
    SingleChannel::MatchedFilterParameters parmsCopy(parms);
    EXPECT_EQ(parmsCopy.getMatchedFilterImplementation(),
              MatchedFilterImplementation::AUTO);
    EXPECT_EQ(parmsCopy.getSignalSize(), signalSize);
    EXPECT_EQ(parmsCopy.getDetectionMode(),
              MatchedFilterDetectionMode::SINGLE_CHANNEL);
    EXPECT_TRUE(parmsCopy.isValid());
    EXPECT_NO_THROW(tplateOut = parmsCopy.getTemplate(0));
    tplateOut.getSignal(work8.size(), &work8Ptr);
    error = infinityNorm(tp1.size(), tp1.data(), work8.data());
    EXPECT_NEAR(error, 0, 1.e-14);

    EXPECT_NO_THROW(tplateOut = parms.getTemplate(1));
    tplateOut.getSignal(work8.size(), &work8Ptr);
    error = infinityNorm(tp2.size(), tp2.data(), work8.data());
    EXPECT_NEAR(error, 0, 1.e-14);

    // Test clear
    EXPECT_NO_THROW(parms.clearTemplates());
    EXPECT_EQ(parms.getNumberOfTemplates(), 0);
    parms.clear();
    EXPECT_FALSE(parms.isValid());
}

TEST(singleChannelMatchedFilter, testDouble)
{
    SingleChannel::MatchedFilterParameters parms;
    SingleChannel::MatchedFilter<double> mf;
}



std::vector<double> generateUniformRandomNumbers(
   const int n, const double lower, const double upper)
{
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<> dist(lower, upper);
    std::vector<double> result(n, 0); 
    for (auto &v : result)
    {
        v = dist(rng);
    }   
    return result;
}

template<class T> T infinityNorm(const int n, const T *__restrict__ r)
{
    auto minmax = std::minmax_element(r, r+n);
    auto linf = std::max(std::abs(*minmax.first), std::abs(*minmax.second));
    return linf;    
}

template<class T> T infinityNorm(const int n,
                                 const T *__restrict__ x,
                                 const T *__restrict__ y)
{
    T linf = 0;
    #pragma omp simd reduction(max:linf)
    for (int i=0; i<n; ++i)
    {
        linf = std::max(linf, std::abs(x[i] - y[i]));
    } 
    return linf;
}


}