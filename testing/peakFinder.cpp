#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "private/peakFinder.hpp"
#include <gtest/gtest.h>

namespace
{
using namespace MFLib;

template<class T>
std::vector<T> gaussianSine(const T tdur, const T dt, 
                            const T sigma, const T sineFreq);

TEST(peakFinder, peakFinder)
{
    PeakFinder<double> peaks;
    int signalSize = 4000;
    const double tdur = 4;
    const double dt = 0.01;
    const double sigma = tdur*0.3;
    const double sineFreq = 12;
    double tol = 0.8;
    int minPeakDistance = static_cast<int> (0.5/dt + 1) + 1;
    peaks.setThreshold(tol);
    peaks.setMinimumPeakDistance(minPeakDistance);
    EXPECT_NEAR(peaks.getThreshold(), tol, 1.e-14);
    EXPECT_EQ(peaks.getMinimumPeakDistance(), minPeakDistance);

    // Let's do an easy case of a constant below the threshold
    std::vector<double> signal(signalSize, 0.4);
    EXPECT_NO_THROW(peaks.setSignal(signalSize, signal.data()));
    EXPECT_NO_THROW(peaks.apply());
    EXPECT_EQ(peaks.getNumberOfPeaks(), 0);

    // Make some sines modulated by a Gaussians
    auto tplate = gaussianSine(tdur, dt, sigma, sineFreq);
    // Embed this template at a few points in the signal
    std::fill(signal.begin(), signal.end(), 0);
    auto ihalf = static_cast<int> (tplate.size()/2);
    std::copy(tplate.begin()+ihalf, tplate.end(), signal.begin());

    int i1 = signalSize*0.25 - ihalf;
    std::copy(tplate.begin(), tplate.end(), signal.begin()+i1);
 
    i1 = signalSize*0.5 - ihalf;
    std::copy(tplate.begin(), tplate.end(), signal.begin()+i1);
   
    i1 = signalSize*0.75 - ihalf;
    std::copy(tplate.begin(), tplate.end(), signal.begin()+i1);
 
    i1 = signal.size() - ihalf;
    std::copy(tplate.begin()+1, tplate.begin()+ihalf+1, signal.begin()+i1);

    // Set the signal
    EXPECT_NO_THROW(peaks.setSignal(signal.size(), signal.data()));
    peaks.apply();
    auto npeaks = peaks.getNumberOfPeaks();
    std::vector<int> peakIndices(signalSize);
    auto peakPtr = peakIndices.data();
    // Should have 3*3 + 2*2 = 13 peaks
    EXPECT_EQ(npeaks, 13);
    peaks.getPeaks(npeaks, &peakPtr);
    // The first and last peak should be at the end points.
    // The other peaks should be the main and second side lobes
    std::vector<int> ref1(13);
    ref1[0] = 0;
    ref1[1] = 67; 

    ref1[2] = signalSize*0.25 - 67;
    ref1[3] = signalSize*0.25;
    ref1[4] = signalSize*0.25 + 67;

    ref1[5] = signalSize*0.5 - 67; 
    ref1[6] = signalSize*0.5;
    ref1[7] = signalSize*0.5 + 67; 

    ref1[8] = signalSize*0.75 - 67; 
    ref1[9] = signalSize*0.75;
    ref1[10]= signalSize*0.75 + 67; 

    ref1[11] = signalSize - 1 - 67;
    ref1[12] = signalSize - 1;
    // The other peaks should be the main and second side lobes
    for (int i=0; i<npeaks; ++i)
    {
        EXPECT_EQ(peakPtr[i], ref1[i]);
    }

    // Now let's keep those near-in peaks
    EXPECT_NO_THROW(peaks.setThreshold(0.8));
    EXPECT_NO_THROW(peaks.setMinimumPeakDistance(0));
    EXPECT_EQ(peaks.getNumberOfPeaks(), 0);
    peaks.setSignal(signal.size(), signal.data());
    peaks.apply();
    // 5*3 + 2*3 = 21 peaks
    EXPECT_EQ(peaks.getNumberOfPeaks(), 21);
    peaks.getPeaks(peaks.getNumberOfPeaks(), &peakPtr);
    std::vector<int> ref2(21);
    ref2[0] = 0;
    ref2[1] = 33;
    ref2[2] = 67; 

    ref2[3] = signalSize*0.25 - 67; 
    ref2[4] = signalSize*0.25 - 33;
    ref2[5] = signalSize*0.25;
    ref2[6] = signalSize*0.25 + 33;
    ref2[7] = signalSize*0.25 + 67; 

    ref2[8]  = signalSize*0.5 - 67; 
    ref2[9]  = signalSize*0.5 - 33;
    ref2[10] = signalSize*0.5;
    ref2[11] = signalSize*0.5 + 33; 
    ref2[12]  = signalSize*0.5 + 67; 

    ref2[13] = signalSize*0.75 - 67; 
    ref2[14] = signalSize*0.75 - 33;
    ref2[15] = signalSize*0.75;
    ref2[16] = signalSize*0.75 + 33;
    ref2[17] = signalSize*0.75 + 67; 

    ref2[18] = signalSize - 1 - 67;
    ref2[19] = signalSize - 1 - 33;
    ref2[20] = signalSize - 1;     
    for (int i=0; i<ref2.size(); ++i)
    {
        EXPECT_EQ(ref2[i], peakPtr[i]);
        //printf("%d, %d\n", ref2[i], peakPtr[i]);
    }
/*
    FILE *fout = fopen("gaussSine.txt", "w");
    //for (int i=0; i<tplate.size(); ++i)
    //{
    //fprintf(fout, "%lf, %lf\n", dt*i, tplate[i]);
    //}
    for (int i=0; i<signal.size(); ++i)
    {
    fprintf(fout, "%lf, %lf\n", dt*i, signal[i]);
    }
    fclose(fout);
*/
}

template<class T>
std::vector<T> gaussianSine(const T tdur, const T dt,
                            const T sigma, const T sineFreq)
{
    auto t0 =-tdur/2; 
    auto t1 = tdur/2;
    auto npts = static_cast<int> ((t1 - t0)/dt + 0.5) + 1;
    std::vector<T> result(npts);
    T max = 0;
    for (int i=0; i<npts; ++i)
    {
        auto t = t0 + static_cast<T> (i)*dt;
        auto expArg = (t/sigma)*(t/sigma);
        auto sinArg = 2*M_PI*((t - t0))/tdur*sineFreq;
        result[i] = std::exp(-(expArg*expArg))*std::cos(sinArg);
        max = std::max(std::abs(result[i]), max);
    }
    // Normalize
    for (int i=0; i<npts; ++i)
    {
        result[i] = result[i]/max;
    }
    result[0] = 0;
    result[result.size()-1] = 0;
    return result;
}

}
