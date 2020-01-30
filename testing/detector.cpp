#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <array>
#include <algorithm>
#include "utilities.hpp"
#include "mflib/waveformTemplate.hpp"
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

TEST(singleChannelDetection, detectorParameters)
{
    DetectorParameters parms;
    int minSamples = 50;
    double thresh = 0.8;
    auto policy = MFLib::MaximumMatchedFilterPolicy::ABSOLUTE_MAXIMUM; 
    EXPECT_NO_THROW(parms.setMinimumDetectionSpacing(minSamples));
    EXPECT_NO_THROW(parms.setDetectionThreshold(thresh));
    EXPECT_NO_THROW(parms.setMaximaPolicy(policy));
 
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

    DetectorParameters parmsCopy(parms);
    EXPECT_EQ(parmsCopy.getMinimumDetectionSpacing(), minSamples);
    EXPECT_NEAR(parmsCopy.getDetectionThreshold(), thresh, 1.e-14);
    EXPECT_EQ(parmsCopy.getMaximaPolicy(), policy);
    EXPECT_TRUE(parmsCopy.wantDetectedWaveform());
    EXPECT_EQ(parmsCopy.wantDetectedWaveform(), parms.wantDetectedWaveform());
    EXPECT_EQ(parmsCopy.wantAmplitudeScalingFactor(),
              parms.wantAmplitudeScalingFactor());
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
    // Data needs to be noisy otherwise we just get a constant decay of spikes
    auto signal = generateUniformRandomNumbers(signalSize, -0.05, 0.05);
    std::vector<double> tonsets;
    std::vector<double> tonsetsShifted;
    //signal(signalSize, 0);
    for (int i=0; i<nt; ++i)
    {
        auto t = createTemplateSignal(templateSize,
                                      amp[i], dt, freq[i],
                                      decay[i], 0);
        wts[i].setSamplingRate(df);
        wts[i].setIdentifier(i);
        wts[i].setSignal(t.size(), t.data());
        double xscal = 2;
        if (i%2 == 0){xscal = 4;}
        t = createTemplateSignal(templateSize,
                                 1./xscal*amp[i], dt, freq[i],
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
    detParms.setMinimumDetectionSpacing(minSpacing);
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
    for (int i=0; i<ndet; ++i)
    {
        auto detection = detector[i];
        EXPECT_NEAR(detection.getDetectionTime(),
                    tonsets[i], 1.e-10);
        int id = detection.getTemplateIdentifier();
        auto intDetection = detection.getInterpolatedDetectionTime();
        EXPECT_NEAR(std::abs(tonsetsShifted[i] - intDetection), 0.0, dt/2);
        auto alpha1 = detection.getAmplitudeScalingFactor(
            MFLib::RelativeMagnitudeType::GIBBONS_RINGDAL_2006);
        auto alpha2 = detection.getAmplitudeScalingFactor(
            MFLib::RelativeMagnitudeType::SCHAFF_RICHARDS_2014); 
 printf("%d, %lf, %lf, %lf, %lf, %lf\n", id,
        tonsetsShifted[i]- intDetection,
        tonsetsShifted[i], intDetection, alpha1, alpha2);
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
