#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "mflib/networkStationPhase.hpp"
#include "mflib/singleChannel/detection.hpp"
#include "mflib/singleChannel/relativeMagnitude.hpp"
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
    double similarity = 0.8;
    double detectionTime = 1;
    double intDetectionTime = 2;
    double onsetTime = 3;
    double intOnsetTime = 4; 
    double alphaGR = 10;  // Scaling factor
    double alphaSR = 100; 
    MFLib::Polarity polarity = MFLib::Polarity::DILATATION;
    auto magGR = convertAmplitudeScalingFactorToMagnitudePerturbation(alphaGR);
    auto magSR = convertAmplitudeScalingFactorToMagnitudePerturbation(alphaSR);
    
    MFLib::NetworkStationPhase nsp;
    nsp.setNetwork("PB");
    nsp.setStation("B206");
    nsp.setPhase("S"); 
    uint64_t eventID = 20;
    std::pair<MFLib::NetworkStationPhase, uint64_t> id(nsp, eventID);

    EXPECT_FALSE(det.haveCorrelationCoefficient());
    EXPECT_FALSE(det.haveTemplateIdentifier());
    EXPECT_FALSE(det.haveDetectedSignal());
    EXPECT_FALSE(det.haveDetectionTime());
    EXPECT_FALSE(det.haveInterpolatedDetectionTime());
    EXPECT_FALSE(det.havePhaseOnsetTime());
    EXPECT_FALSE(det.haveInterpolatedPhaseOnsetTime());
    EXPECT_FALSE(det.haveAmplitudeScalingFactor(
        MFLib::RelativeMagnitudeType::GIBBONS_RINGDAL_2006));
    EXPECT_FALSE(det.haveAmplitudeScalingFactor(
        MFLib::RelativeMagnitudeType::SCHAFF_RICHARDS_2014));
    EXPECT_EQ(det.getPolarity(), MFLib::Polarity::UNKNOWN);

    EXPECT_NO_THROW(det.setCorrelationCoefficient(similarity));
    EXPECT_NO_THROW(det.setTemplateIdentifier(id));
    EXPECT_NO_THROW(det.setDetectedSignal(signal.size(), signal.data()));
    EXPECT_NO_THROW(det.setDetectionTime(detectionTime));
    EXPECT_NO_THROW(det.setInterpolatedDetectionTime(intDetectionTime));
    EXPECT_NO_THROW(det.setPhaseOnsetTime(onsetTime));
    EXPECT_NO_THROW(det.setInterpolatedPhaseOnsetTime(intOnsetTime));
    EXPECT_NO_THROW(det.setAmplitudeScalingFactor(alphaGR,
                    MFLib::RelativeMagnitudeType::GIBBONS_RINGDAL_2006));
    EXPECT_NO_THROW(det.setAmplitudeScalingFactor(alphaSR,
                    MFLib::RelativeMagnitudeType::SCHAFF_RICHARDS_2014));
    det.setPolarity(polarity);

    EXPECT_NEAR(det.getCorrelationCoefficient(), similarity, 1.e-14);
    EXPECT_EQ(det.getTemplateIdentifier().first, nsp);
    EXPECT_EQ(det.getTemplateIdentifier().second, eventID);
    EXPECT_EQ(det.getDetectedSignalLength(), npts);
    EXPECT_NO_THROW(det.getDetectedSignal(npts, &sigPtr));
    error = infinityNorm(npts, sigPtr, signal.data());
    EXPECT_NEAR(error, 0, 1.e-14);
    EXPECT_NEAR(det.getDetectionTime(), detectionTime, 1.e-14);
    EXPECT_NEAR(det.getInterpolatedDetectionTime(), intDetectionTime, 1.e-14);
    EXPECT_NEAR(det.getPhaseOnsetTime(), onsetTime, 1.e-14);
    EXPECT_NEAR(det.getInterpolatedPhaseOnsetTime(), intOnsetTime, 1.e-14); 
    EXPECT_NEAR(det.getAmplitudeScalingFactor(
                MFLib::RelativeMagnitudeType::GIBBONS_RINGDAL_2006),
                alphaGR, 1.e-14);
    EXPECT_NEAR(det.getMagnitudePerturbation(
                MFLib::RelativeMagnitudeType::GIBBONS_RINGDAL_2006),
                magGR, 1.e-14);
    EXPECT_NEAR(det.getAmplitudeScalingFactor(
                MFLib::RelativeMagnitudeType::SCHAFF_RICHARDS_2014),
                alphaSR, 1.e-14);
    EXPECT_NEAR(det.getMagnitudePerturbation(
                MFLib::RelativeMagnitudeType::SCHAFF_RICHARDS_2014),
                magSR, 1.e-14);
    EXPECT_EQ(det.getPolarity(), polarity);

    Detection detCopy(det);
    EXPECT_NEAR(detCopy.getCorrelationCoefficient(), similarity, 1.e-14);
    EXPECT_EQ(detCopy.getTemplateIdentifier().first, nsp);
    EXPECT_EQ(detCopy.getTemplateIdentifier().second, eventID);
    EXPECT_EQ(detCopy.getDetectedSignalLength(), npts);
    EXPECT_NO_THROW(detCopy.getDetectedSignal(npts, &sigPtr));
    error = infinityNorm(npts, sigPtr, signal.data());
    EXPECT_NEAR(error, 0, 1.e-14);
    EXPECT_NEAR(detCopy.getDetectionTime(), detectionTime, 1.e-14);
    EXPECT_NEAR(detCopy.getInterpolatedDetectionTime(),
                intDetectionTime, 1.e-14);
    EXPECT_NEAR(detCopy.getPhaseOnsetTime(), onsetTime, 1.e-14);
    EXPECT_NEAR(detCopy.getInterpolatedPhaseOnsetTime(), intOnsetTime, 1.e-14);
    EXPECT_NEAR(detCopy.getAmplitudeScalingFactor(
                MFLib::RelativeMagnitudeType::GIBBONS_RINGDAL_2006),
                alphaGR, 1.e-14);
    EXPECT_NEAR(detCopy.getAmplitudeScalingFactor(
                MFLib::RelativeMagnitudeType::SCHAFF_RICHARDS_2014),
                alphaSR, 1.e-14);
    EXPECT_NEAR(detCopy.getMagnitudePerturbation(
                MFLib::RelativeMagnitudeType::GIBBONS_RINGDAL_2006),
                magGR, 1.e-14);
    EXPECT_NEAR(detCopy.getMagnitudePerturbation(
                MFLib::RelativeMagnitudeType::SCHAFF_RICHARDS_2014),
                magSR, 1.e-14);
    EXPECT_EQ(detCopy.getPolarity(), polarity);

    det.clear();
    EXPECT_FALSE(det.haveCorrelationCoefficient());
    EXPECT_FALSE(det.haveTemplateIdentifier());
    EXPECT_FALSE(det.haveDetectedSignal());
    EXPECT_FALSE(det.haveDetectionTime());
    EXPECT_FALSE(det.haveInterpolatedDetectionTime());
    EXPECT_FALSE(det.havePhaseOnsetTime());
    EXPECT_FALSE(det.haveInterpolatedPhaseOnsetTime());
    EXPECT_FALSE(det.haveAmplitudeScalingFactor(
        MFLib::RelativeMagnitudeType::GIBBONS_RINGDAL_2006));
    EXPECT_FALSE(det.haveAmplitudeScalingFactor(
        MFLib::RelativeMagnitudeType::SCHAFF_RICHARDS_2014));
    EXPECT_EQ(det.getPolarity(), MFLib::Polarity::UNKNOWN);
}

}
