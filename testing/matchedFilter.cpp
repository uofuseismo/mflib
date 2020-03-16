#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <random>
#include <chrono>
#include <ctime>
#include <ipps.h>
#include "mflib/waveformTemplate.hpp"
#include "mflib/matchedFilterParameters.hpp"
#include "mflib/matchedFilter.hpp"
#include "utilities.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace MFLib;

/*
void matchedFilter(const int nb, const double b[],
                   const int nx, const double x[],
                   double y[]);
*/
void dumbXC(const int nb, const double b[],
            const int nx, const double x[],
            double xc[]);
void dumbXC(const int nb, const float b[],
            const int nx, const float x[],
            float xc[]);

TEST(matchedFilter, basicTestDouble)
{
    MatchedFilterParameters options;
    MatchedFilter<double> mf;

    int signalSize = 1000;
    double Tmax = 8;
    double dt = Tmax/(signalSize - 1);
    double f0 = 4.0/Tmax; // 4 cycles in desired output 
    double amp = 2;
    std::vector<double> x(signalSize);
    std::vector<double> x2(signalSize);
    for (int i=0; i<signalSize; ++i)
    {
        auto t = i*dt;
        auto y = amp*std::exp(-0.1*t)*std::sin(2*M_PI*f0*t);
        x[i] = y;
        x2[i] = amp*std::exp(-0.15*t)*std::cos(2*M_PI*f0*t);
    }
    // Make a template - have a sine/cosine wave look for the sine/cosine
    // wave in an exponentially damped wave.  This works because of the 
    // scaling in the denominator.
    double Twin = 2;
    int nb = static_cast<int> (Twin/dt + 0.5) + 1;
    std::vector<double> b(nb);
    std::vector<double> b2(nb);
    for (int i=0; i<nb; ++i)
    {
        auto t = i*dt;
        auto y = std::sin(2*M_PI*f0*t);
        b[i] = y;
        b2[i] = std::cos(2*M_PI*f0*t);
    }
    // Create a reference
    std::vector<double> yRef(signalSize, 0), yRef2(signalSize, 0);
    //matchedFilter(nb, b.data(), signalSize, x.data(), yRef.data());
    dumbXC(nb, b.data(),  signalSize, x.data(),  yRef.data()); 
    dumbXC(nb, b2.data(), signalSize, x2.data(), yRef2.data());
    // Set the templates
    options.setMatchedFilterImplementation(
        MatchedFilterImplementation::AUTO);
    options.setSignalSize(signalSize);
    WaveformTemplate bt, bt2;
    bt.setSignal(b.size(), b.data());
    bt2.setSignal(b2.size(), b2.data());
    bt.setSamplingRate(1);
    bt2.setSamplingRate(1);
    bt.setTravelTime(0);
    bt2.setTravelTime(0);
    bt.setPhaseOnsetTime(0);
    bt2.setPhaseOnsetTime(0);
    options.addTemplate(bt); //b.size(), b.data());
    EXPECT_EQ(options.getNumberOfTemplates(), 1);
    // Initialize the matched filter
    mf.initialize(options);
    auto templateSpectra = mf.getSpectraOfTemplate(0);
    EXPECT_NO_THROW(mf.setSignal(0, signalSize, x.data()));
    const double *waveDataPtr = nullptr;
    EXPECT_NO_THROW(waveDataPtr = mf.getSignalPointer(0));
    auto error = infinityNorm(x.size(), x.data(), waveDataPtr); 
    EXPECT_LT(error, 1.e-14);

    EXPECT_FALSE(mf.haveMatchedFilteredSignals());

    EXPECT_NO_THROW(mf.apply());
    EXPECT_TRUE(mf.haveMatchedFilteredSignals());
    std::vector<double> result;
    EXPECT_NO_THROW(result = mf.getMatchedFilteredSignal(0));
    ippsNormDiff_Inf_64f(yRef.data(), result.data(), result.size(), &error);
    EXPECT_LT(error, 1.e-13);
    //printf("%e\n", error);

    // Try again with four templates
    options.clearTemplates();
    options.addTemplate(bt); //b.size(),  b.data());
    options.addTemplate(bt2); //b2.size(), b2.data());
    options.addTemplate(bt); //b.size(),  b.data());
    options.addTemplate(bt2); //b2.size(), b2.data());
    EXPECT_EQ(options.getNumberOfTemplates(), 4);
    EXPECT_NO_THROW(mf.initialize(options));
    EXPECT_NO_THROW(mf.setSignal(0, signalSize, x.data()));
    EXPECT_NO_THROW(mf.setSignal(1, signalSize, x2.data()));
    EXPECT_NO_THROW(mf.setSignal(2, signalSize, x.data()));
    EXPECT_NO_THROW(mf.setSignal(3, signalSize, x2.data()));
    EXPECT_NO_THROW(waveDataPtr = mf.getSignalPointer(3));
    error = infinityNorm(signalSize, x2.data(), waveDataPtr); 
    EXPECT_LT(error, 1.e-14);
    EXPECT_EQ(mf.getNumberOfTemplates(), 4);
    EXPECT_NO_THROW(mf.apply());
    EXPECT_NO_THROW(result = mf.getMatchedFilteredSignal(0));
    ippsNormDiff_Inf_64f(yRef.data(), result.data(), result.size(), &error);
    EXPECT_LT(error, 1.e-13);
    //printf("e %e\n", error);

    EXPECT_NO_THROW(result = mf.getMatchedFilteredSignal(1));
    ippsNormDiff_Inf_64f(yRef2.data(), result.data(), result.size(), &error);
    EXPECT_LT(error, 1.e-13);
    //printf("e %e\n", error);

    EXPECT_NO_THROW(result = mf.getMatchedFilteredSignal(2));
    ippsNormDiff_Inf_64f(yRef.data(), result.data(), result.size(), &error);
    EXPECT_LT(error, 1.e-13);
    //printf("e %e\n", error);

    EXPECT_NO_THROW(result = mf.getMatchedFilteredSignal(3));
    ippsNormDiff_Inf_64f(yRef2.data(), result.data(), result.size(), &error);
    EXPECT_LT(error, 1.e-13);
    //printf("e %e\n", error);

    // Put this in a loop
}

//----------------------------------------------------------------------------//

TEST(matchedFilter, basicTestFloat)
{
    MatchedFilterParameters options;
    MatchedFilter<float> mf; 

    int signalSize = 1000;
    double Tmax = 8;
    double dt = Tmax/(signalSize - 1); 
    double f0 = 4.0/Tmax; // 4 cycles in desired output 
    double amp = 2;
    std::vector<float> x(signalSize);
    std::vector<float> x2(signalSize);
    for (int i=0; i<signalSize; ++i)
    {
        auto t = i*dt;
        auto y = amp*std::exp(-0.1*t)*std::sin(2*M_PI*f0*t);
        x[i] = y;
        x2[i] = amp*std::exp(-0.15*t)*std::cos(2*M_PI*f0*t);
    }   
    // Make a template - have a sine/cosine wave look for the sine/cosine
    // wave in an exponentially damped wave.  This works because of the 
    // scaling in the denominator.
    double Twin = 2;
    int nb = static_cast<int> (Twin/dt + 0.5) + 1;
    std::vector<float> b(nb);
    std::vector<float> b2(nb);
    for (int i=0; i<nb; ++i)
    {
        auto t = i*dt;
        auto y = std::sin(2*M_PI*f0*t);
        b[i] = y;
        b2[i] = std::cos(2*M_PI*f0*t);
    }   
    // Create a reference
    std::vector<float> yRef(signalSize, 0), yRef2(signalSize, 0); 
    //matchedFilter(nb, b.data(), signalSize, x.data(), yRef.data());
    dumbXC(nb, b.data(),  signalSize, x.data(),  yRef.data()); 
    dumbXC(nb, b2.data(), signalSize, x2.data(), yRef2.data());
    // Set the templates
    options.setMatchedFilterImplementation(
        MatchedFilterImplementation::AUTO);
    options.setSignalSize(signalSize);
    WaveformTemplate bt, bt2;
    bt.setSignal(b.size(), b.data());
    bt2.setSignal(b2.size(), b2.data());
    bt.setSamplingRate(1);
    bt2.setSamplingRate(1);
    bt.setTravelTime(0);
    bt2.setTravelTime(0);
    bt.setPhaseOnsetTime(0);
    bt2.setPhaseOnsetTime(0);
    options.addTemplate(bt); //b.size(), b.data());
    EXPECT_EQ(options.getNumberOfTemplates(), 1);
    // Initialize the matched filter
    mf.initialize(options);
    auto templateSpectra = mf.getSpectraOfTemplate(0);
    EXPECT_NO_THROW(mf.setSignal(0, signalSize, x.data()));
    EXPECT_FALSE(mf.haveMatchedFilteredSignals());

    EXPECT_NO_THROW(mf.apply());
    EXPECT_TRUE(mf.haveMatchedFilteredSignals());
    std::vector<float> result;
    EXPECT_NO_THROW(result = mf.getMatchedFilteredSignal(0));
    float error;
    ippsNormDiff_Inf_32f(yRef.data(), result.data(), result.size(), &error);
    EXPECT_LT(error, 2.e-6);

    // Try again with four templates
    options.clearTemplates();
    options.addTemplate(bt); //b.size(),  b.data());
    options.addTemplate(bt2); //b2.size(), b2.data());
    options.addTemplate(bt); //b.size(),  b.data());
    options.addTemplate(bt2); //b2.size(), b2.data());
    EXPECT_EQ(options.getNumberOfTemplates(), 4); 
    EXPECT_NO_THROW(mf.initialize(options));
    EXPECT_NO_THROW(mf.setSignal(0, signalSize, x.data()));
    EXPECT_NO_THROW(mf.setSignal(1, signalSize, x2.data()));
    EXPECT_NO_THROW(mf.setSignal(2, signalSize, x.data()));
    EXPECT_NO_THROW(mf.setSignal(3, signalSize, x2.data()));
    EXPECT_EQ(mf.getNumberOfTemplates(), 4); 
    EXPECT_NO_THROW(mf.apply());
    EXPECT_NO_THROW(result = mf.getMatchedFilteredSignal(0));
    ippsNormDiff_Inf_32f(yRef.data(), result.data(), result.size(), &error);
    EXPECT_LT(error, 2.e-6);
    //printf("e %e\n", error);

    EXPECT_NO_THROW(result = mf.getMatchedFilteredSignal(1));
    ippsNormDiff_Inf_32f(yRef2.data(), result.data(), result.size(), &error);
    EXPECT_LT(error, 2.e-6);
    //printf("e %e\n", error);

    EXPECT_NO_THROW(result = mf.getMatchedFilteredSignal(2));
    ippsNormDiff_Inf_32f(yRef.data(), result.data(), result.size(), &error);
    EXPECT_LT(error, 2.e-6);
    //printf("e %e\n", error);

    EXPECT_NO_THROW(result = mf.getMatchedFilteredSignal(3));
    ippsNormDiff_Inf_32f(yRef2.data(), result.data(), result.size(), &error);
    EXPECT_LT(error, 2.e-6);
}

//----------------------------------------------------------------------------//

TEST(matchedFilter, shiftAndStack)
{
    double Tmax = 120;
    double df = 100;
    auto dt = 1/df;
    int signalSize = static_cast<int> (Tmax/dt + 0.5) + 1;
    EXPECT_NEAR(Tmax, (signalSize-1)*dt, 1.e-14);
    double templateDuration = 4; // + 2;
    int templateSize = static_cast<int> (templateDuration/dt + 0.5) + 1;
    EXPECT_NEAR(templateDuration, (templateSize-1)*dt, 1.e-14);
    double toff1 = 1; // Arrival offset into template: 1 s 
    double toff2 = 2; //2; // 2;
    double toff3 = 1.5; //1.5; // 2.5;
    int nextra1 = static_cast<int> (toff1/dt + 0.5);// - 10;
    int nextra2 = static_cast<int> (toff2/dt + 0.5);// + 10;
    int nextra3 = static_cast<int> (toff3/dt + 0.5);// + 20;
    // Templates of different lengths are padded to same size.
    int nalloc = std::max(templateSize+nextra1, templateSize+nextra2);
    nalloc = std::max(nalloc, templateSize+nextra3);
    std::vector<double> t1(nalloc, 0); // + 1 s
    std::vector<double> t2(nalloc, 0); // + 2 s
    std::vector<double> t3(nalloc, 0); // + 2.5 s
/*
    std::vector<double> t1(templateSize+nextra1, 0); // + 1 s
    std::vector<double> t2(templateSize+nextra2, 0); // + 2 s
    std::vector<double> t3(templateSize+nextra3, 0); // + 2.5 s
*/

    std::vector<double> signal1(signalSize, 0); 
    std::vector<double> signal2(signalSize, 0);
    std::vector<double> signal3(signalSize, 0);
    double originTime = 8;
    double tt1 = 20; // Travel time 20 s
    double tt2 = 40;
    double tt3 = 80;
    int id1 = static_cast<int> ((originTime + tt1)/dt + 0.5);
    int id2 = static_cast<int> ((originTime + tt2)/dt + 0.5);
    int id3 = static_cast<int> ((originTime + tt3)/dt + 0.5);
    for (int i=0; i<static_cast<int> (signal1.size()); ++i)
    {
        double s1 = 0.001;
        if (i%2 == 0){s1 =-s1;}
        signal1[i] = s1;
        signal2[i] = s1;
        signal3[i] = s1;
    }
    for (int i=0; i<templateSize; ++i)
    {
        auto t = i*dt;
        auto s1 = 2.0*std::exp(-t)*std::sin(2*M_PI*t/templateDuration);
        auto s2 = 1.5*std::exp(-1.5*t)*std::sin(2*M_PI*t/templateDuration);
        auto s3 = 1.2*std::exp(-2.5*t)*std::sin(2*M_PI*t/templateDuration);
        if (i == 0 || i == templateSize - 1)
        {
            s1 = 0;
            s2 = 0;
            s3 = 0;
        }
        t1[i+nextra1] = s1;
        t2[i+nextra2] = s2;
        t3[i+nextra3] = s3;
        signal1[id1+i] = s1;
        signal2[id2+i] = s2;
        signal3[id3+i] = s3;
    }
    // Set the matched filter
    MatchedFilterParameters options;
    MatchedFilter<double> mf; 
    WaveformTemplate bt1, bt2, bt3;

/*
    bt1.setSignal(t1.size(), t1.data());
    bt2.setSignal(t2.size(), t2.data());
    bt3.setSignal(t3.size(), t3.data());
*/
    bt1.setSignal(templateSize+nextra1, t1.data());
    bt2.setSignal(templateSize+nextra2, t2.data());
    bt3.setSignal(templateSize+nextra3, t3.data());

    bt1.setSamplingRate(df);
    bt2.setSamplingRate(df);
    bt3.setSamplingRate(df);
    bt1.setTravelTime(tt1);
    bt2.setTravelTime(tt2);
    bt3.setTravelTime(tt3);
    bt1.setPhaseOnsetTime(toff1);
    bt2.setPhaseOnsetTime(toff2);
    bt3.setPhaseOnsetTime(toff3);

    options.setMatchedFilterImplementation(MatchedFilterImplementation::AUTO);
    options.setSignalSize(signalSize);
    options.addTemplate(bt1);
    options.addTemplate(bt2);
    options.addTemplate(bt3);
    
    mf.initialize(options);
    mf.setSignal(0, signal1.size(), signal1.data());
    mf.setSignal(1, signal2.size(), signal2.data());
    mf.setSignal(2, signal3.size(), signal3.data());
    mf.apply();

    auto pc1 = mf.getMatchedFilteredSignal(0);
    auto pc2 = mf.getMatchedFilteredSignal(1);
    auto pc3 = mf.getMatchedFilteredSignal(2);
/*
FILE *xcFile = fopen("mfs.est.txt", "w");
for (int i=0; i<pc1.size(); ++i)
{
fprintf(xcFile, "%lf, %lf, %lf, %lf\n", i*dt, pc1[i], pc2[i], pc3[i]); 
}
fclose(xcFile);
*/
    // Compute a reference solution
    std::vector<double> yRef1(signalSize, 0);
    std::vector<double> yRef2(signalSize, 0);
    std::vector<double> yRef3(signalSize, 0);
//printf("%ld, %ld, %ld\n", t1.size(), t2.size(), t3.size());
printf("%d, %d, %d\n", templateSize+nextra1, templateSize+nextra2, templateSize+nextra3);
 
/*
    dumbXC(t1.size(), t1.data(), signalSize, signal1.data(), yRef1.data());
    dumbXC(t2.size(), t2.data(), signalSize, signal2.data(), yRef2.data());
    dumbXC(t3.size(), t3.data(), signalSize, signal3.data(), yRef3.data());
*/
    dumbXC(templateSize+nextra1, t1.data(), signalSize, signal1.data(), yRef1.data());
    dumbXC(templateSize+nextra2, t2.data(), signalSize, signal2.data(), yRef2.data());
    dumbXC(templateSize+nextra3, t3.data(), signalSize, signal3.data(), yRef3.data());
/*
xcFile = fopen("mfs.txt", "w");
for (int i=0; i<pc1.size(); ++i)
{
fprintf(xcFile, "%lf, %lf, %lf, %lf\n", i*dt, yRef1[i], yRef2[i], yRef3[i]); 
}
fclose(xcFile);
*/

    double error = 0;
/*
    ippsNormDiff_Inf_64f(yRef1.data()+100, pc1.data(), pc1.size(), &error);
printf("%e\n", error);
    ippsNormDiff_Inf_64f(yRef2.data(), pc2.data(), pc2.size(), &error);
printf("%e\n", error);
    ippsNormDiff_Inf_64f(yRef3.data()+50, pc3.data(), pc3.size(), &error);
printf("%e\n", error);
for (int i=pc1.size()-50; i<pc1.size(); ++i)
{
printf("%d %lf, %lf, %lf\n", i, yRef1[i], pc1[i], yRef1[i] - pc1[i]);
}
*/

    ippsNormDiff_Inf_64f(yRef1.data(), pc1.data(), pc1.size(), &error);
    EXPECT_LE(error, 1.e-11);
printf("%e\n", error);
    ippsNormDiff_Inf_64f(yRef2.data(), pc2.data(), pc2.size(), &error);
    EXPECT_LE(error, 1.e-11);
printf("%e\n", error);
    ippsNormDiff_Inf_64f(yRef3.data(), pc3.data(), pc3.size(), &error);
    EXPECT_LE(error, 1.e-11);
printf("%e\n", error);
//getchar();

    // Should migrate to the origin time
    auto stack = mf.shiftAndStack();
    double xmax;
    int imax = 0;
    ippsMaxIndx_64f(stack.data(), stack.size(), &xmax, &imax);
    auto tori = imax/df;
    EXPECT_NEAR(xmax, 1.0, 1.e-3); // With weighting 3 perfect stacks sum to 1
    EXPECT_NEAR(tori, originTime, 1/df/2);  // Within a sample
    //printf("%lf, %lf\n", tori, xmax);
/*
FILE *ofl1 = fopen("testShift.ref.txt", "w");
FILE *ofl2 = fopen("testShift.txt", "w");
for (int i=0; i<signal1.size(); ++i)
{
fprintf(ofl1, "%lf, %lf, %lf, %lf\n", i*dt, yRef1[i], yRef2[i], yRef3[i]);
fprintf(ofl2, "%lf, %lf, %lf, %lf\n", i*dt, pc1[i], pc2[i], pc3[i]);
}
fclose(ofl1);
fclose(ofl2);
*/
/*
FILE *ofl3 = fopen("stack.txt", "w");
    for (int i=0; i<stack.size(); ++i)
    {
fprintf(ofl3, "%lf, %lf\n", i*dt, stack[i]); 
    }
fclose(ofl3);
*/

}


/*
TEST(matchedFilter, stressTest)
{
    MatchedFilterParameters<double> options;
    MatchedFilter<double> mf;

    //int signalSize = 86400*100; //8640000; //360000; //1000;
    int signalSize = 372000;
    double Tmax = 8;
    double dt = Tmax/(signalSize - 1);
    double f0 = 4.0/Tmax; // 4 cycles in desired output 
    double amp = 2;
    std::vector<double> x(signalSize);
    std::vector<double> x2(signalSize);
    #pragma omp simd
    for (int i=0; i<signalSize; ++i)
    {
        auto t = i*dt;
        auto y = amp*std::exp(-0.1*t)*std::sin(2*M_PI*f0*t);
        x[i] = y;
        x2[i] = amp*std::exp(-0.15*t)*std::cos(2*M_PI*f0*t);
    }
    // Make a template - have a sine/cosine wave look for the sine/cosine
    // wave in an exponentially damped wave.  This works because of the 
    // scaling in the .
    double Twin = 2;
    int nb = static_cast<int> (Twin/dt + 0.5) + 1;
    std::vector<double> b(nb);
    std::vector<double> b2(nb);
    #pragma omp simd
    for (int i=0; i<nb; ++i)
    {
        auto t = i*dt;
        auto y = std::sin(2*M_PI*f0*t);
        b[i] = y;
        b2[i] = std::cos(2*M_PI*f0*t);
    }
    // Set the templates
    options.setMatchedFilterImplementation(
        MatchedFilterImplementation::AUTO);
    options.setSignalSize(signalSize);
    options.addTemplate(b.size(), b.data());
    options.setFFTLength(8640004);
    EXPECT_EQ(options.getNumberOfTemplates(), 1); 
    // Initialize the matched filter
    mf.initialize(options);
    auto templateSpectra = mf.getSpectraOfTemplate(0);
printf("setting signal\n");
    EXPECT_NO_THROW(mf.setSignal(0, signalSize, x.data()));
    EXPECT_FALSE(mf.haveMatchedFilteredSignals());
printf("working...\n");
    EXPECT_NO_THROW(mf.apply());
}

TEST(matchedFilter, stressTestFloat)
{
    MatchedFilterParameters<float> options;
    MatchedFilter<float> mf; 

    int signalSize = 86400*100; //8640000; //360000; //1000;
    double Tmax = 8;
    double dt = Tmax/(signalSize - 1); 
    double f0 = 4.0/Tmax; // 4 cycles in desired output 
    double amp = 2;
    std::vector<float> x(signalSize);
    std::vector<float> x2(signalSize);
    #pragma omp simd
    for (int i=0; i<signalSize; ++i)
    {   
        auto t = i*dt;
        auto y = amp*std::exp(-0.1*t)*std::sin(2*M_PI*f0*t);
        x[i] = y;
        x2[i] = amp*std::exp(-0.15*t)*std::cos(2*M_PI*f0*t);
    }   
    // Make a template - have a sine/cosine wave look for the sine/cosine
    // wave in an exponentially damped wave.  This works because of the 
    // scaling in the .
    double Twin = 2;
    int nb = static_cast<int> (Twin/dt + 0.5) + 1;
    std::vector<float> b(nb);
    std::vector<float> b2(nb);
    #pragma omp simd
    for (int i=0; i<nb; ++i)
    {   
        auto t = i*dt;
        auto y = std::sin(2*M_PI*f0*t);
        b[i] = y;
        b2[i] = std::cos(2*M_PI*f0*t);
    }   
    // Set the templates
    options.setMatchedFilterImplementation(
        MatchedFilterImplementation::AUTO);
    options.setSignalSize(signalSize);
    options.addTemplate(b.size(), b.data());
    options.setFFTLength(8640004);
    EXPECT_EQ(options.getNumberOfTemplates(), 1); 
    // Initialize the matched filter
    mf.initialize(options);
    auto templateSpectra = mf.getSpectraOfTemplate(0);
printf("setting signal\n");
    EXPECT_NO_THROW(mf.setSignal(0, signalSize, x.data()));
    EXPECT_FALSE(mf.haveMatchedFilteredSignals());
printf("working...\n");
    EXPECT_NO_THROW(mf.apply());
}
*/

/*
void matchedFilter(const int nb, const double b[], // Template
                   const int nx, const double x[], // Signal
                   double y[])
{
    // Add a little extra to template to account for start-up
    Ipp64f *xwork = ippsMalloc_64f(nb+nx);
    ippsZero_64f(xwork, nb+nx);
    ippsCopy_64f(x, xwork, nx);
    Ipp64f *ywork = ippsMalloc_64f(nb+nx);

    // Reverse the template so that the convolution is a correlation
    Ipp64f *bflip = ippsMalloc_64f(nb);
    std::reverse_copy(b, b+nb, bflip);
    double norm;
    ippsNorm_L2_64f(bflip, nb, &norm);
    ippsDivC_64f_I(norm, bflip, nb);

    int specSize, bufferSize;
    ippsFIRSRGetSize(nb, ipp64f, &specSize, &bufferSize);
    auto pSpec = reinterpret_cast<IppsFIRSpec_64f *> (ippsMalloc_8u(specSize));
    Ipp8u *pBuf = ippsMalloc_8u(bufferSize);
    ippsFIRSRInit_64f(bflip, nb, ippAlgDirect, pSpec); 
    ippsFIRSR_64f(xwork, ywork, nx+nb, pSpec, NULL, NULL, pBuf);
    ippsCopy_64f(&ywork[nb-1], y, nx);
    ippsFree(pSpec);
    ippsFree(pBuf);
    ippsFree(bflip);
    ippsFree(xwork);
    ippsFree(ywork);
}
*/

/// Implementation
void dumbXC(const int nb, const double b[],
            const int nx, const double x[],
            double xc[])
{
    double *btemp = ippsMalloc_64f(nb);
    double *xtemp = ippsMalloc_64f(nb);
    double mean = 0;
    ippsMean_64f(b, nb, &mean);
    ippsSubC_64f(b, mean, btemp, nb);
    double Eb;
    ippsNorm_L2_64f(btemp, nb, &Eb);
    for (int i=0; i<nx-nb+1; ++i)
    {
        double Ex;
        ippsMean_64f(&x[i], nb, &mean);
        ippsSubC_64f(&x[i], mean, xtemp, nb);
        ippsNorm_L2_64f(xtemp, nb, &Ex);
        Ex = std::max(std::numeric_limits<double>::epsilon(), Ex);
        //ippsDotProd_64f(btemp, xtemp, nb, &xc[i]);
        xc[i] = std::inner_product(btemp, btemp+nb, xtemp, 0.0); //ippsDotProd_64f(btemp, xtemp, nb, &xc[i]);
        xc[i] = xc[i]/(Eb*Ex);
    }
    ippsFree(btemp);
    ippsFree(xtemp);
}

void dumbXC(const int nb, const float b[],
            const int nx, const float x[],
            float xc[])
{
    float *btemp = ippsMalloc_32f(nb);
    float *xtemp = ippsMalloc_32f(nb);
    float mean = 0;
    ippsMean_32f(b, nb, &mean, ippAlgHintAccurate);
    ippsSubC_32f(b, mean, btemp, nb);
    float Eb; 
    ippsNorm_L2_32f(btemp, nb, &Eb);
    for (int i=0; i<nx-nb+1; ++i)
    {   
        float Ex; 
        ippsMean_32f(&x[i], nb, &mean, ippAlgHintAccurate);
        ippsSubC_32f(&x[i], mean, xtemp, nb);
        ippsNorm_L2_32f(xtemp, nb, &Ex);
        Ex = std::max(std::numeric_limits<float>::epsilon(), Ex);
        ippsDotProd_32f(btemp, xtemp, nb, &xc[i]);
        xc[i] = xc[i]/(Eb*Ex);
    }   
    ippsFree(btemp);
    ippsFree(xtemp);
}

/*
#include "advisor-annotate.h"

template<typename T>
void demeanSignal(const size_t n, const size_t nWin,
                  const T *__restrict__ y,
                  T *__restrict__ yMean,
                  T *__restrict__ yStd)
{
    // Set space on all processes
omp_set_num_threads(2);
ANNOTATE_SITE_BEGIN(demeanSignal)
    #pragma omp parallel \
     shared(y, yMean, yStd) \
     firstprivate(n, nWin) \
     default(none)
    {
    auto meanWork2 = static_cast<T *> (std::aligned_alloc(64, nWin*sizeof(T)));
    auto y2 = static_cast<T *> (std::aligned_alloc(64, 2*nWin*sizeof(T)));
    auto xnorm = 1/static_cast<T> (nWin);
    // Parallel loop on the waveform chunks
    #pragma omp for schedule(static)
    for (size_t i=0; i<n-nWin; i=i+nWin)
    {
        const T *yp = &y[i];
        T *mean = &yMean[i];
        T *ys = &yStd[i];
        // Get the mean
        T s = 0;
        T s2 = 0;
        #pragma omp simd reduction(+:s, s2) aligned(y2: 64)
        for (size_t j=0; j<nWin; ++j)
        {
            s  = s + yp[j];
            y2[j] = yp[j]*yp[j];
            s2 = s2 + y2[j];
        }
        #pragma omp simd
        for (size_t j=nWin; j<2*nWin; ++j){y2[j] = yp[j]*yp[j];}
        // Now we calculate the mean for an updated sample
        // \mu_n     = \frac{1}{N} \sum_{i=1}^N y_i
        // \mu_{n+1} = \frac{1}{N} \sum_{i=2}^{N+1} y_i
        //           = \frac{1}{N} \sum_{i=1}^N y_i
        //            + \frac{y_{n+1} - y_1}{N}
        //           = \mu_n + \frac{1}{N} (y_{n+1} - y_1) 
        mean[0] = xnorm*s;
        for (size_t j=1; j<nWin; ++j)
        {
            mean[j] = mean[j-1] + xnorm*(yp[nWin-1+j] - yp[j-1]);
        }
        // Next, put E[y]^2 into the workspace
        #pragma omp simd aligned(meanWork2: 64)
        for (size_t j=0; j<nWin; ++j)
        {
            meanWork2[j] = mean[j]*mean[j];
        }
        // Now, we need the variance signal
        //   \var{y_n,y_n} = E[y^2] - E[y_n]^2
        // Basically, we play the same trick as with the mean:
        //   S_0 = \frac{1}{N} \sum_{i=1}^N y_i^2 - \mu_0^2
        //   S_1 = \frac{1}{N} \sum_{i=2}^{N+2} y_i^2 - \mu_1^2
        //       = S_0 + \mu_0^2
        //       + \frac{1}{N} (y_{n+1}^2 - y_n)
        //       - \mu_1^2
        ys[0] = xnorm*s2 - meanWork2[0];
        for (size_t j=1; j<nWin; ++j)
        {
            ys[j] = ys[j-1] + meanWork2[j-1]
                  + xnorm*(y2[nWin-1+j] - y2[j-1])
                  - meanWork2[j];
        }
        #pragma omp simd
        for (size_t j=0; j<nWin; ++j){ys[j] = std::sqrt(ys[j]);}
    }
    free(meanWork2);
    }
ANNOTATE_SITE_END();
}

TEST(testDemean, demean)
{
    size_t nWin = 901;
    size_t n = 118192;
    double *y = static_cast<double *> (std::aligned_alloc(64, n*sizeof(double)));
    double *yStd = static_cast<double *> (std::aligned_alloc(64, n*sizeof(double)));
    double *yMean = static_cast<double *> (std::aligned_alloc(64, n*sizeof(double)));
    double *yStdRef = static_cast<double *> (std::aligned_alloc(64, n*sizeof(double)));
    double *yMeanRef = static_cast<double *> (std::aligned_alloc(64, n*sizeof(double)));
    std::random_device device;
    std::mt19937 rng(device());
    std::uniform_real_distribution<> dist(1.0,6.0); 
    for (size_t i=0; i<n; ++i){y[i] = dist(rng);}
    // Demean the signal
    //std::chrono::high_resolution_clock::time_point t1
    int maxit = 10;
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int iter=0; iter<maxit; ++iter)
    {
        demeanSignal(n, nWin, y, yMean, yStd);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    auto deltaT = std::chrono::duration_cast<std::chrono::duration<double>> (t1 - t0);
    double fastTime = deltaT.count();
    printf("Average optimized time: %e\n", fastTime/maxit);
    // Compute the reference
    t0 = std::chrono::high_resolution_clock::now();
    for (int iter=0; iter<maxit; ++iter)
    {
        for (int i=0; i<n-nWin; ++i)
        {
            double xmean = 0;
            #pragma omp simd reduction(+:xmean)
            for (int j=0; j<nWin; ++j){xmean = xmean + y[i+j];}
            xmean = xmean/static_cast<double> (nWin);
            double var = 0;
            #pragma omp simd reduction(+:var)
            for (int j=0; j<nWin; ++j)
            {
                double res = y[i+j] - xmean;
                var = var + res*res;
            }
            yMeanRef[i] = xmean;
            yStdRef[i] = std::sqrt(var/static_cast<double> (nWin));
         }
    }
    t1 = std::chrono::high_resolution_clock::now();
    deltaT = std::chrono::duration_cast<std::chrono::duration<double>> (t1 - t0);
    double slowTime = deltaT.count();
    double emaxStd = 0;
    double emaxMean = 0;
    for (int i=0; i<n-nWin; ++i)
    {
         emaxMean = std::max(emaxMean, std::abs(yMeanRef[i] - yMean[i]));
         emaxStd  = std::max(emaxStd,  std::abs(yStdRef[i]  - yStd[i]));
    }
    printf("Average naive time: %e\n", deltaT.count()/maxit); 
    printf("Speedup: %lf\n", slowTime/fastTime);
    EXPECT_LE(emaxMean, 1.e-12);
    EXPECT_LE(emaxStd,  1.e-12);
printf("%e, %e\n", emaxMean, emaxStd);
    free(y);
    free(yStd);
    free(yMean);
    free(yStdRef);
    free(yMeanRef);
}
*/

}
