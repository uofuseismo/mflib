#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <random>
#include <vector>
#include <chrono>
#include <random>
#include "mflib/waveformTemplate.hpp"
#include "mflib/singleChannel/relativeMagnitude.hpp"
#include "utilities.hpp"
#include <gtest/gtest.h>

namespace
{
using namespace MFLib;

double get_dmag_1(double *x,double *y,int np);
double get_dmag_2(double *x,double *y,int np);
std::vector<double> demean(const std::vector<double> &x);

TEST(singleChannelRelativeMagnitude, relativeMagnitude)
{
    int signalSize = 400;
   
    WaveformTemplate wt;
    MFLib::SingleChannel::RelativeMagnitude mag;
    // Make a random signal to test this out
    auto randomTemplate = generateUniformRandomNumbers(signalSize, -2, 1);
    wt.setSignal(randomTemplate.size(), randomTemplate.data());
    auto testSignal1 = generateUniformRandomNumbers(signalSize, -3, 8);
    auto testSignal2 = generateUniformRandomNumbers(signalSize, -7, 1); 
    EXPECT_NO_THROW(mag.initialize(wt));
    EXPECT_EQ(mag.getDetectedSignalLength(), signalSize);
    EXPECT_TRUE(mag.isInitialized());
    EXPECT_FALSE(mag.haveDetectedSignal());
    // Do two computations
    mag.setDetectedSignal(testSignal1.size(), testSignal1.data());
    EXPECT_TRUE(mag.haveDetectedSignal());
    // Compute
    auto x = demean(randomTemplate);
    auto y = demean(testSignal1);
    double dmagRef1 = get_dmag_1(x.data(), y.data(), x.size());
    double dmagRef2 = get_dmag_2(x.data(), y.data(), x.size());
    auto dmag1 = mag.computeMagnitudePerturbation(
                    RelativeMagnitudeType::GIBBONS_RINGDAL_2006); 
    EXPECT_NEAR(dmag1, dmagRef1, 1.e-12);
    auto dmag2 = mag.computeMagnitudePerturbation(
                    RelativeMagnitudeType::SCHAFF_RICHARDS_2014);
    EXPECT_NEAR(dmag2, dmagRef2, 1.e-12);
    // Test copy c'tor by copying second signal onto mag and computing
    mag.setDetectedSignal(testSignal1.size(), testSignal2.data());
    y = demean(testSignal2);
    dmagRef1 = get_dmag_1(x.data(), y.data(), x.size());
    dmagRef2 = get_dmag_2(x.data(), y.data(), x.size());
    MFLib::SingleChannel::RelativeMagnitude magCopy(mag);
    dmag1 = magCopy.computeMagnitudePerturbation(
              RelativeMagnitudeType::GIBBONS_RINGDAL_2006);
    EXPECT_NEAR(dmag1, dmagRef1, 1.e-12);
    dmag2 = magCopy.computeMagnitudePerturbation(
              RelativeMagnitudeType::SCHAFF_RICHARDS_2014);
    EXPECT_NEAR(dmag2, dmagRef2, 1.e-12);
}

std::vector<double> demean(const std::vector<double> &x)
{
    std::vector<double> y(x.size());
    auto mean = std::accumulate(x.begin(), x.end(), 0.0)
               /static_cast<double> (x.size());
    std::transform(x.begin(), x.end(), y.begin(),
                   [&](auto value){return value - mean;}); 
    return y;
}
/*
*  Assuming the master event is x, this gives the 
*  the delta_mag for y relative to x, using the typical
*  Gibbons and Ringdal (2006) least squares solution,
*  as shown in equation 1 of Schaff & Richards (2014).
*  @author Keith Koper
*/
double get_dmag_1(double *x,double *y,int np)
{
    int i;
    double num,denom;
    num=0.0;
    denom=0.0;
    for (i=0;i<np;i++)
    {
        num+=x[i]*y[i];
        denom+=x[i]*x[i];
    }
    num = std::abs(num);
    return(log10(num/denom));
}

/*
*  Assuming the master event is x, this gives the 
*  the delta_mag for y relative to x, using ratios
*  of stdvs  as advocated by Derrick and
*  equation 11 of Schaff & Richards (2014)
 * @author Keith Koper
*/
double get_dmag_2(double *x,double *y,int np)
{
    int i;
    double num,denom;
    num=0.0;
    denom=0.0;
    for (i=0;i<np;i++)
    {
        num+=y[i]*y[i];
        denom+=x[i]*x[i];
    }
    num=sqrt(num/(double)np);
    denom=sqrt(denom/(double)np);
    return(log10(num/denom));
}
}
