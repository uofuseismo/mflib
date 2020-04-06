#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <array>
#include <algorithm>
#include "utilities.hpp"
#include "mflib/waveformTemplate.hpp"
#include "mflib/networkStationPhase.hpp"
#include "mflib/singleChannel/detector.hpp"
#include "mflib/singleChannel/detectorParameters.hpp"
#include "mflib/singleChannel/matchedFilter.hpp"
#include "mflib/singleChannel/matchedFilterParameters.hpp"
#include <gtest/gtest.h>

namespace
{
using namespace MFLib::SingleChannel;

std::vector<double> createTemplateSignal(
    const int npts,
    const double amplitude,
    const double dt,
    const double freq,
    const double decay,
    const double phaseShift);

TEST(singleChannelDetection, interpParameters)
{
    DetectionTimeInterpolationParameters parms;
    int nalpha = 21;
    int npts = 50;
    EXPECT_EQ(parms.getType(),
              MFLib::DetectionTimeInterpolationType::QUADRATIC);

    parms.enableLanczos();
    parms.setLanczosAlpha(nalpha);
    parms.setLanczosNumberOfInterpolationPoints(npts);
    EXPECT_EQ(parms.getType(),
              MFLib::DetectionTimeInterpolationType::LANCZOS);
    EXPECT_EQ(parms.getLanczosAlpha(), nalpha);
    EXPECT_EQ(parms.getLanczosNumberOfInterpolationPoints(), npts);

    DetectionTimeInterpolationParameters parmsCopy(parms);
    EXPECT_EQ(parmsCopy.getType(),
              MFLib::DetectionTimeInterpolationType::LANCZOS);
    EXPECT_EQ(parmsCopy.getLanczosAlpha(), nalpha);
    EXPECT_EQ(parmsCopy.getLanczosNumberOfInterpolationPoints(), npts);
 
}

TEST(singleChannelDetection, detectorParameters)
{
    DetectorParameters parms;
    DetectionTimeInterpolationParameters interpParms;
    interpParms.disable();
    int minSamples = 50;
    double thresh = 0.8;
    auto policy = MFLib::MaximumMatchedFilterPolicy::ABSOLUTE_MAXIMUM; 
    EXPECT_NO_THROW(parms.setMinimumDetectionSpacing(minSamples));
    EXPECT_NO_THROW(parms.setDetectionThreshold(thresh));
    EXPECT_NO_THROW(parms.setMaximaPolicy(policy));
    parms.setDetectionTimeInterpolationParameters(interpParms);
    auto detPolicy = MFLib::MatchedFilteredSignalDetectorPolicy::SINGLE;
    parms.setMatchedFilteredSignalDetectorPolicy(detPolicy);

    parms.disableSaveDetectedWaveform();
    EXPECT_FALSE(parms.wantDetectedWaveform());
    parms.enableSaveDetectedWaveform();

    parms.disableSaveAmplitudeScalingFactor();
    EXPECT_FALSE(parms.wantAmplitudeScalingFactor());
    parms.enableSaveAmplitudeScalingFactor();

    EXPECT_EQ(parms.getMinimumDetectionSpacing(), minSamples);
    EXPECT_NEAR(parms.getDetectionThreshold(), thresh, 1.e-14);
    EXPECT_EQ(parms.getMaximaPolicy(), policy);
    EXPECT_TRUE(parms.wantDetectedWaveform());
    EXPECT_EQ(parms.getMatchedFilteredSignalDetectorPolicy(), detPolicy);
    auto interpParmsBack = parms.getDetectionTimeInterpolationParameters();
    EXPECT_EQ(interpParmsBack.getType(),
              MFLib::DetectionTimeInterpolationType::NONE);    

    DetectorParameters parmsCopy(parms);
    EXPECT_EQ(parmsCopy.getMinimumDetectionSpacing(), minSamples);
    EXPECT_NEAR(parmsCopy.getDetectionThreshold(), thresh, 1.e-14);
    EXPECT_EQ(parmsCopy.getMaximaPolicy(), policy);
    EXPECT_TRUE(parmsCopy.wantDetectedWaveform());
    EXPECT_EQ(parmsCopy.wantDetectedWaveform(), parms.wantDetectedWaveform());
    EXPECT_EQ(parmsCopy.wantAmplitudeScalingFactor(),
              parms.wantAmplitudeScalingFactor());
    EXPECT_EQ(parmsCopy.getMatchedFilteredSignalDetectorPolicy(),
              parms.getMatchedFilteredSignalDetectorPolicy()); 
    interpParmsBack = parmsCopy.getDetectionTimeInterpolationParameters();
    EXPECT_EQ(interpParmsBack.getType(),
              MFLib::DetectionTimeInterpolationType::NONE);
}

TEST(singleChannelDetection, detector)
{
    // Set some parameters
    double df = 100;
    double dt = 1/df;
    double templateDur = 4;
    double signalDur = 120;
    int minSpacing = static_cast<int> (templateDur/2/dt);
    int signalSize = static_cast<int> (signalDur/dt + 0.5) + 1;
    int templateSize = static_cast<int> (templateDur/dt + 0.5) + 1;
    // Embed a few templates
    const int nt = 5;
    std::array<double, nt>   amp({  1.0,   2.0,   0.5,  0.75,  3.0});
    std::array<double, nt>  freq({ 10.0,   5.0,   7.5,  12.5,  2.5});
    std::array<double, nt> decay({  2.8,   2.9,   2.6,   2.2,  2.3});
    std::array<double, nt>   phi({0.001,  0.002, 0.004,-0.003, -0.002});
    std::vector<MFLib::WaveformTemplate> wts(nt);
    // Make this for a particular station
    MFLib::NetworkStationPhase nsp;
    nsp.setNetwork("UU");
    nsp.setStation("CCUT");
    nsp.setPhase("P");
    // Data needs to be noisy otherwise we just get a constant decay of spikes
    auto signal = generateUniformRandomNumbers(signalSize, -0.05, 0.05);
    std::vector<double> tonsets;
    std::vector<double> tonsetsShifted;
    std::vector<double> ampScaled;
    std::vector<std::pair<double, int>> perm;
    std::vector<std::pair<MFLib::NetworkStationPhase, uint64_t>> ids;
    int nc = 0;
    //signal(signalSize, 0);
    for (int i=0; i<nt; ++i)
    {
        auto t = createTemplateSignal(templateSize,
                                      amp[i], dt, freq[i],
                                      decay[i], 0);
        wts[i].setSamplingRate(df);
        std::pair<MFLib::NetworkStationPhase, uint64_t> id(nsp, i);
        wts[i].setIdentifier(id);
        wts[i].setSignal(t.size(), t.data());
        double xscal = 2;
        if (i%2 == 0){xscal = 4;}
        auto scaleWave = 1./xscal;
        t = createTemplateSignal(templateSize,
                                 scaleWave*amp[i], dt, freq[i],
                                 decay[i], phi[i]);
        std::vector<double> tinsert;
        if (i == 0)
        {
            tinsert.push_back(4);
            tinsert.push_back(28);
            tinsert.push_back(90);
        }
        else if (i == 1)
        {
            tinsert.push_back(8);
            tinsert.push_back(24);
            tinsert.push_back(100);
        }
        else if (i == 2)
        {
            tinsert.push_back(33);
            tinsert.push_back(66);
            tinsert.push_back(86); // Causes cycle skip problem
        }
        else if (i == 4)
        {
            tinsert.push_back(16);
            tinsert.push_back(44);
            tinsert.push_back(106);
        }
        else // Intentionally skip 3 
        {
            continue;
        }
        for (const auto &toff : tinsert)
        {
            auto index = static_cast<int> (toff/dt + 0.5); 
            for (int k=0; k<static_cast<int> (t.size()); ++k)
            {
                signal[k+index] = t[k];
            }
            tonsets.push_back(toff);
            tonsetsShifted.push_back(toff + phi[i]);
            ampScaled.push_back(scaleWave);
            perm.push_back(std::pair(toff, nc));
            ids.push_back(id);
            nc = nc + 1;
        }
    }
/*
FILE *fout = fopen("templates.txt", "w");
for (int i=0; i<wts.size(); ++i)
{
std::string tfile = "templates" + std::to_string(i+1) + ".txt";
FILE *fout = fopen(tfile.c_str(), "w");
std::vector<double> temp(wts[i].getSignalLength());
auto tempPtr = temp.data();
wts[i].getSignal(temp.size(), &tempPtr);
for (int j=0; j<wts[i].getSignalLength(); ++j)
{
fprintf(fout, "%lf, %lf\n", j*dt, tempPtr[j]);
}
fclose(fout);
}
*/
    // Set the matched filtered detector parameters
    MatchedFilterParameters mfParms;
    for (int i=0; i<nt; ++i)
    {
        mfParms.addTemplate(wts[i]);
    }
    mfParms.setSignalSize(signalSize);
   
    // Compute the detections
    MatchedFilter<double> mf;
    mf.initialize(mfParms);
    mf.setSignal(signal.size(), signal.data());
    mf.apply();
    // Now compute the detections
    DetectorParameters detParms;
    DetectionTimeInterpolationParameters interpParms;
    interpParms.enableLanczos();
    detParms.setDetectionTimeInterpolationParameters(interpParms);
    detParms.setMinimumDetectionSpacing(minSpacing);
    detParms.setMatchedFilteredSignalDetectorPolicy(
        MFLib::MatchedFilteredSignalDetectorPolicy::REDUCED);
    detParms.enableSaveDetectedWaveform();
    detParms.enableSaveAmplitudeScalingFactor();
    Detector<double> detector;
    EXPECT_NO_THROW(detector.initialize(detParms));
    std::vector<Detection<double>> detections;
    EXPECT_NO_THROW(detector.detect(mf));
    auto ndet = detector.getNumberOfDetections();
    EXPECT_EQ(ndet, static_cast<int> (tonsets.size())); 
    std::sort(tonsets.begin(), tonsets.end());
    std::sort(tonsetsShifted.begin(), tonsetsShifted.end());
    std::sort(perm.begin(), perm.end(), 
              [](const std::pair<double, int> &a, std::pair<double, int> &b)
              {
                  return a.first < b.first;
              });
    for (int i=0; i<ndet; ++i)
    {
        auto detection = detector[i];
        EXPECT_NEAR(detection.getDetectionTime(),
                    tonsets[i], 1.e-10);
        auto id = detection.getTemplateIdentifier();
        auto j = perm[i].second;
        EXPECT_EQ(ids[j].first, id.first);
        EXPECT_EQ(ids[j].second, id.second);
        auto intDetection = detection.getInterpolatedDetectionTime();
        EXPECT_NEAR(std::abs(tonsetsShifted[i] - intDetection), 0.0, dt/2);
        auto alpha1 = detection.getAmplitudeScalingFactor(
            MFLib::RelativeMagnitudeType::GIBBONS_RINGDAL_2006);
        auto alpha2 = detection.getAmplitudeScalingFactor(
            MFLib::RelativeMagnitudeType::SCHAFF_RICHARDS_2014); 
        EXPECT_NEAR(ampScaled[j], alpha1, 1.e-2);
        EXPECT_NEAR(ampScaled[j], alpha2, 1.e-2);
 printf("%d, %d, %lf, %lf, %lf, %lf, %lf, %lf\n",
        static_cast<int> (ids[j].second), static_cast<int> (id.second),
        tonsetsShifted[i]- intDetection,
        tonsetsShifted[i], intDetection, alpha1, alpha2, ampScaled[j]);
    }
}

std::vector<double> createTemplateSignal(
    const int npts,
    const double amplitude,
    const double dt,
    const double freq,
    const double decay,
    const double phaseShift)
{
    std::vector<double> result(npts, 0);
    auto omega = 2*M_PI*freq;
    for (int i=0; i<npts; ++i)
    {
        auto t = i*dt;
        result[i] = amplitude*std::exp(-t*decay)
                   *std::sin(omega*(t - phaseShift));
    }
    return result;
}

}
