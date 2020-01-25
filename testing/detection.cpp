#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "mflib/singleChannel/detection.hpp"
#include "utilities.hpp"
#include <gtest/gtest.h>

namespace
{
using namespace MFLib::SingleChannel;

TEST(singleChannelDetection, detection)
{
    Detection det;    
    int npts = 10;
    std::vector<double> signal(npts), sigTest(npts);
    for (int i=1; i<npts+1; ++i){signal[i-1] = i;}
    double error;
    double *sigPtr = sigTest.data();
    double detectionTime = 1;
    double intDetectionTime = 2;
    double onsetTime = 3;
    double intOnsetTime = 4; 

    EXPECT_FALSE(det.haveDetectedSignal());
    EXPECT_FALSE(det.haveDetectionTime());
    EXPECT_FALSE(det.haveInterpolatedDetectionTime());
    EXPECT_FALSE(det.havePhaseOnsetTime());
    EXPECT_FALSE(det.haveInterpolatedPhaseOnsetTime());

    EXPECT_NO_THROW(det.setDetectedSignal(signal.size(), signal.data()));
    EXPECT_NO_THROW(det.setDetectionTime(detectionTime));
    EXPECT_NO_THROW(det.setInterpolatedDetectionTime(intDetectionTime));
    EXPECT_NO_THROW(det.setPhaseOnsetTime(onsetTime));
    EXPECT_NO_THROW(det.setInterpolatedPhaseOnsetTime(intOnsetTime));
 
    EXPECT_EQ(det.getDetectedSignalLength(), npts);
    EXPECT_NO_THROW(det.getDetectedSignal(npts, &sigPtr));
    error = infinityNorm(npts, sigPtr, signal.data());
    EXPECT_NEAR(error, 0, 1.e-14);
    EXPECT_NEAR(det.getDetectionTime(), detectionTime, 1.e-14);
    EXPECT_NEAR(det.getInterpolatedDetectionTime(), intDetectionTime, 1.e-14);
    EXPECT_NEAR(det.getPhaseOnsetTime(), onsetTime, 1.e-14);
    EXPECT_NEAR(det.getInterpolatedPhaseOnsetTime(), intOnsetTime, 1.e-14); 

    Detection detCopy(det);
    EXPECT_EQ(detCopy.getDetectedSignalLength(), npts);
    EXPECT_NO_THROW(detCopy.getDetectedSignal(npts, &sigPtr));
    error = infinityNorm(npts, sigPtr, signal.data());
    EXPECT_NEAR(error, 0, 1.e-14);
    EXPECT_NEAR(detCopy.getDetectionTime(), detectionTime, 1.e-14);
    EXPECT_NEAR(detCopy.getInterpolatedDetectionTime(),
                intDetectionTime, 1.e-14);
    EXPECT_NEAR(detCopy.getPhaseOnsetTime(), onsetTime, 1.e-14);
    EXPECT_NEAR(detCopy.getInterpolatedPhaseOnsetTime(), intOnsetTime, 1.e-14);

    det.clear();
    EXPECT_FALSE(det.haveDetectedSignal());
    EXPECT_FALSE(det.haveDetectionTime());
    EXPECT_FALSE(det.haveInterpolatedDetectionTime());
    EXPECT_FALSE(det.havePhaseOnsetTime());
    EXPECT_FALSE(det.haveInterpolatedPhaseOnsetTime());
}

}