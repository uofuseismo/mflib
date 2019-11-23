#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <random>
#include <omp.h>
#include <chrono>
#include <ctime>
#include <ipps.h>
#include "mflib/matchedFilterOptions.hpp"
#include "mflib/matchedFilter.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace MFLib;

void matchedFilter(const int nb, const double b[],
                   const int nx, const double x[],
                   double y[]);
void dumbXC(const int nb, const double b[],
            const int nx, const double x[],
            double xc[]);

TEST(matchedFilter, basicTest)
{
    MatchedFilterOptions<double> options;
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
    // scaling in the .
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
/*
    FILE *filterRef = fopen("filterRef.txt", "w");
    for (int i=0; i<yRef.size(); ++i)
    {
        fprintf(filterRef, "%lf, %lf\n", i*dt, yRef[i]);
    }
    fclose(filterRef);
*/
    // Set the templates
    options.setMatchedFilterImplementation(
        MatchedFilterImplementation::AUTO);
    options.setSignalSize(signalSize);
    options.addTemplate(b.size(), b.data());
    EXPECT_EQ(options.getNumberOfTemplates(), 1);
    // Initialize the matched filter
    mf.initialize(options);
    auto templateSpectra = mf.getSpectraOfTemplate(0);
    EXPECT_NO_THROW(mf.setSignal(0, signalSize, x.data()));
    EXPECT_FALSE(mf.haveMatchedFilteredSignals());

    EXPECT_NO_THROW(mf.apply());
    EXPECT_TRUE(mf.haveMatchedFilteredSignals());
    std::vector<double> result;
    EXPECT_NO_THROW(result = mf.getMatchedFilteredSignal(0));
    double error;
    ippsNormDiff_Inf_64f(yRef.data(), result.data(), result.size(), &error);
    EXPECT_LT(error, 5.e-3);
    //printf("%e\n", error);

    // Try again with four templates
    options.clearTemplates();
    options.addTemplate(b.size(),  b.data());
    options.addTemplate(b2.size(), b2.data());
    options.addTemplate(b.size(),  b.data());
    options.addTemplate(b2.size(), b2.data());
    EXPECT_EQ(options.getNumberOfTemplates(), 4);
    EXPECT_NO_THROW(mf.initialize(options));
    EXPECT_NO_THROW(mf.setSignal(0, signalSize, x.data()));
    EXPECT_NO_THROW(mf.setSignal(1, signalSize, x2.data()));
    EXPECT_NO_THROW(mf.setSignal(2, signalSize, x.data()));
    EXPECT_NO_THROW(mf.setSignal(3, signalSize, x2.data()));
    EXPECT_EQ(mf.getNumberOfTemplates(), 4);
    EXPECT_NO_THROW(mf.apply());
    EXPECT_NO_THROW(result = mf.getMatchedFilteredSignal(0));
    ippsNormDiff_Inf_64f(yRef.data(), result.data(), result.size(), &error);
    EXPECT_LT(error, 5.e-3);
    //printf("%e\n", error);

    EXPECT_NO_THROW(result = mf.getMatchedFilteredSignal(1));
    ippsNormDiff_Inf_64f(yRef2.data(), result.data(), result.size(), &error);
    EXPECT_LT(error, 5.e-3);
    //printf("%e\n", error);

    EXPECT_NO_THROW(result = mf.getMatchedFilteredSignal(2));
    ippsNormDiff_Inf_64f(yRef.data(), result.data(), result.size(), &error);
    EXPECT_LT(error, 5.e-3);
    //printf("%e\n", error);

    EXPECT_NO_THROW(result = mf.getMatchedFilteredSignal(3));
    ippsNormDiff_Inf_64f(yRef2.data(), result.data(), result.size(), &error);
    EXPECT_LT(error, 5.e-3);
    //printf("%e\n", error);
}

/*
TEST(matchedFilter, stressTest)
{
    MatchedFilterOptions<double> options;
    MatchedFilter<double> mf;

    int signalSize = 86400*100; //8640000; //360000; //1000;
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
*/

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

/// Implementation
void dumbXC(const int nb, const double b[],
            const int nx, const double x[],
            double xc[])
{
    double Eb;
    ippsDotProd_64f(b, b, nb, &Eb);
    for (int i=0; i<nx-nb+1; ++i)
    {
        double Ex;
        ippsDotProd_64f(&x[i], &x[i], nb, &Ex);
        ippsDotProd_64f(b, &x[i], nb, &xc[i]); 
        xc[i] = xc[i]/std::sqrt(Eb*Ex);
    }
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
