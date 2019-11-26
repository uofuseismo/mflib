#include <cstdio>
#include <cstdlib>
#include <complex>
#include <array>
#include <vector>
#include <algorithm>
#include <numeric>
#include <mkl.h>
#include <mkl_cblas.h>
#include <ipps.h>
#include <fftw/fftw3.h>
#if __has_include(<pstl/algorithm>)
#include <pstl/algorithm>
#include <pstl/execution>
#define USE_PSTL 1
#endif
#ifdef USE_CUDA
#include <cufft.h>
#endif
#include "mflib/waveformTemplate.hpp"
#include "mflib/matchedFilter.hpp"
#include "mflib/matchedFilterOptions.hpp"

using namespace MFLib;

namespace
{

/// @brief Computes the matrix padding length so that we get bit aligned memory.
/// @param[in] n              The array length.  This must be positive.
/// @param[in] precisionSize  The size of the precision, ex: sizeof(float) = 4.
/// @param[in] alignment      The byte alignment.  This should be a power of 2.
/// @result The padded array length.  This will be greater than or equal to n.
inline int padLength(const int n, 
                     const size_t precisionSize = sizeof(double),
                     const int alignment=64)
{
    auto size = static_cast<int> (precisionSize);
    int padLength = 0;
    auto xmod = (n*size)%alignment;
    if (xmod != 0){padLength = (alignment - xmod)/size;}
    auto nptsPadded = n + padLength;
    return nptsPadded;
}

/// @brief Applies the sliding mean and sliding standard deviation to the
///        filtered signal.
/// @param[in] n         The length of the input signal.
/// @param[in] nWin      The window length.
/// @param[in] y         The signal from which to compute the sliding mean
///                      and standard deviation.  This is an array with
///                      dimension [n].
/// @param[in,out] yNum  On input, this is the numerator signal.  
///                      On exit, this is the normalized signal which will
///                      be the Pearson correlation coefficient at each sample.
/// @param[in] tol       If the denominator falls below tol, the Pearson 
///                      correlation coefficient will be set to 0.  This avoids
///                      division by 0.
/// @note I need to lash this int something like:
///       https://dbs.ifi.uni-heidelberg.de/files/Team/eschubert/publications/SSDBM18-covariance-authorcopy.pdf
///       to get a greater speedup.
template<typename T>
void normalizeSignal(const int n, const int nWin,
                     const T *__restrict__ y,
                     T *__restrict__ yNum,
                     const T tol = 1.e-12)
{
    #pragma omp parallel \
     shared(n, nWin, y, yNum) \
     firstprivate(tol) \
     default(none)
    {
    auto s  = static_cast<T *> (std::aligned_alloc(64, nWin*sizeof(T)));
    auto s2 = static_cast<T *> (std::aligned_alloc(64, nWin*sizeof(T)));
    T nWinInv = 1/static_cast<T> (nWin);
    // Parallel loop on the waveform chunks
    #pragma omp for schedule(static)
    for (int i=0; i<n-nWin; i=i+nWin)
    {
        const T *yp = &y[i];
        T *yn = &yNum[i];
        // Get the initial sums 
        T sum  = 0;
        T sum2 = 0;
        #pragma omp simd reduction(+:sum, sum2)
        for (int j=0; j<nWin; ++j)
        {
            sum  = sum  + yp[j];
            sum2 = sum2 + yp[j]*yp[j];
        }
        // Now we calculate the sum for an updated sample
        // S_0 = \sum_{i=1}^N y_i
        // S_1 = \sum_{i=2}^{N+1} y_i
        //     = \sum_{i=1}^N y_i + y_{N+1} - y_1
        //     = S_0 + y_{N+1} - y_1
        s[0] = sum;
        s2[0] = sum2;
        for (int j=1; j<nWin; ++j)
        {
            auto yn_m_y0 = yp[nWin-1+j] - yp[j-1];
            auto yn_p_y0  = yp[nWin-1+j] + yp[j-1];
            // Update sum
            s[j] = s[j-1] + yn_m_y0;
            // Update sum^2.  Note:
            // (y_n^2 - y_0^2) = (y_n + y_0)*(y_n - y_0)
            // requires 2 adds and one multiplication.  However, we already
            // performed an add.  So this costs one add and one multiply.
            // y_n^2 - y_0^2 requires 2 multiplications and one add.
            s2[j] = s2[j-1] + yn_p_y0*yn_m_y0;
        }
        // Compute the denominator
        #pragma omp simd aligned(s, s2: 64)
        for (int j=0; j<nWin; ++j)
        {
            T den = std::max(tol, std::sqrt(s2[j] - (s[j]*s[j])*nWinInv));
            yn[j] = yn[j]/den;
            if (den < tol){yn[j] = 0;}
        }
    } // End loop
    // Release space on all processes
    free(s);
    free(s2);
    } // End parallel
}

void demeanNormalizeAndReverseTemplate(const int nb, const double b[],
                                       double bOut[])
{
    double mean;
    ippsMean_64f(b, nb, &mean);       // \bar{b}
    ippsSubC_64f(b, mean, bOut, nb);  // bOut = b - Mean
    double norm;
    ippsNorm_L2_64f(bOut, nb, &norm); // Compute L2 norm
    ippsDivC_64f_I(norm, bOut, nb);   // Normalize bOut = bOut/||bOut||_2
    ippsFlip_64f_I(bOut, nb);         // Reverse
}

void demeanNormalizeAndReverseTemplate(const int nb, const float b[],
                                       float bOut[])
{
    float mean;
    ippsMean_32f(b, nb, &mean, ippAlgHintAccurate); // \bar{b}
    ippsSubC_32f(b, mean, bOut, nb);  // bOut = b - Mean
    float norm;
    ippsNorm_L2_32f(bOut, nb, &norm); // Compute L2 norm
    ippsDivC_32f_I(norm, bOut, nb);   // Normalize bOut = bOut/||bOut||_2
    ippsFlip_32f_I(bOut, nb);         // Reverse
}

}

template<>
class MatchedFilter<double>::MatchedFilterImpl
{
public:
    /// Destructor
    ~MatchedFilterImpl()
    {
        clear(); 
    }
    /// Resets the class and clears the memory
    void clear()
    {
        if (mHaveFFTwPlans)
        {
            fftw_destroy_plan(mForwardPlan);
            fftw_destroy_plan(mInversePlan);
        }
        if (mB){fftw_free(mB);}
        //if (mCorrelogramSpectra){MKL_free(mCorrelogramSpectra);}
        if (mSignalSegment){MKL_free(mSignalSegment);}
        if (mSegmentSpectra){MKL_free(mSegmentSpectra);}
        //if (mOVACorrelograms){MKL_free(mOVACorrelograms);}
        if (mInputSignals){MKL_free(mInputSignals);}
        if (mFilteredSignals){MKL_free(mFilteredSignals);}
        if (mDenominator){MKL_free(mDenominator);}
        if (mSkipZeroSignal){MKL_free(mSkipZeroSignal);}
        mB = nullptr;
        mBPtr = nullptr;
        //mCorrelogramSpectra = nullptr;
        //mCorrelogramSpectraPtr = nullptr;
        mSignalSegment = nullptr;
        mSegmentSpectra = nullptr;
        mSegmentSpectraPtr = nullptr;
        mInputSignals = nullptr;
        mFilteredSignals = nullptr;
        mDenominator = nullptr;
        mSkipZeroSignal = nullptr;
        /// Set everything else to 0
        mSpectraLength = 0;
        mFFTLength = 0;
        mConvolutionLeadingDimension = 0;
        mSpectraLeadingDimension = 0;
        mSamplesLeadingDimension = 0;
        mSamplesExtra = 0;
        mSamples = 0;
        mFilterLength = 0;
        mTemplates = 0;
        mL = 0;
        mInitialized = false;
        mHaveFFTwPlans = false;
        mHaveMatchedFilters = false;
    }
//private:
    /// Holds the options
    MatchedFilterOptions mOptions;
    /// FFT forward plan - maps the input signal to the frequency domain.
    fftw_plan mForwardPlan;
    /// FFT inverse plan - brings block of templates convolved with 
    /// with input signal back to the time domain.
    fftw_plan mInversePlan; 
    /// The FFT of the time reversed, demeaned, and normalized templates.
    /// This has dimension: [mTemplates x mSpectraLeadingDimension] and is
    /// stored in row major format.
    fftw_complex *mB = nullptr;
    /// This a pointer to mB which is useful when computing X*B in the 
    /// frequency domain.
    std::complex<double> *mBPtr = nullptr;
    /// The spectra of a segment of the input time domain signals.  This has
    /// dimension [mSpectraLeadingDimension x mTemplates].
    fftw_complex *mSegmentSpectra = nullptr;
    std::complex<double> *mSegmentSpectraPtr = nullptr;
    /// The matched filtered signals.  This has dimension
    /// [mSamplesLeadingDimension x mTemplates].
    double *mFilteredSignals = nullptr;
    /// The input signals.  This has dimension
    /// [mSamplesLeadingDimension x mTemplates].
    double *mInputSignals = nullptr;
    /// The denominator in the Pearson correlation.  This has dimension
    /// [mSamplesLeadingDimension x mTemplates]
    double *mDenominator = nullptr;
    /// The time domain signal to Fourier transform.  This has dimension
    /// [mSamplesLeadingDimension x mTemplates].
    double *mSignalSegment = nullptr;
    /// Determines if I should skip this signal because it is 0.
    /// This has dimension [mTemplates].
    bool *mSkipZeroSignal = nullptr;
    /// The length of the Fourier transforms.  This should equal 
    /// mFFTLength/2 + 1.
    int mSpectraLength = 0;
    /// The FFT length.  This is the number of zero-padded points in the input
    /// time domain signals and output correlograms.
    int mFFTLength = 0;
    /// The leading dimension of the output convolutions in the window.
    size_t mConvolutionLeadingDimension = 0;
    /// The leading dimension of the spectra.
    size_t mSpectraLeadingDimension = 0;
    /// The leading dimension of the numerator/denominator/filtered signals.
    size_t mSamplesLeadingDimension = 0;
    /// The filter needs to `start-up'.  This requires [nb-1] samples which
    /// we disregard when extracting.
    int mSamplesExtra = 0;
    /// The length of the input signals.
    int mSamples = 0;
    /// The length of the templates (filters).
    int mFilterLength = 0;
    /// The number of templates.
    int mTemplates = 0;
    /// The window length.
    int mL = 0;
    /// Determines if the class is initialized.
    bool mInitialized = false;
    /// Determines if the FFTw plan was made
    bool mHaveFFTwPlans = false;
    /// Determines if the matched filtering has been applied.
    bool mHaveMatchedFilters = false;
};

//--------------------------Float Implementation------------------------------//
template<>
class MatchedFilter<float>::MatchedFilterImpl
{
public:
    /// Destructor
    ~MatchedFilterImpl()
    {
        clear();
    }
    /// Resets the class and clears the memory
    void clear()
    {
        if (mHaveFFTwPlans)
        {
            fftwf_destroy_plan(mForwardPlan);
            fftwf_destroy_plan(mInversePlan);
        }
        if (mB){fftw_free(mB);}
        //if (mCorrelogramSpectra){MKL_free(mCorrelogramSpectra);}
        if (mSignalSegment){MKL_free(mSignalSegment);}
        if (mSegmentSpectra){MKL_free(mSegmentSpectra);}
        //if (mOVACorrelograms){MKL_free(mOVACorrelograms);}
        if (mInputSignals){MKL_free(mInputSignals);}
        if (mFilteredSignals){MKL_free(mFilteredSignals);}
        if (mDenominator){MKL_free(mDenominator);}
        if (mSkipZeroSignal){MKL_free(mSkipZeroSignal);}
        mB = nullptr;
        mBPtr = nullptr;
        //mCorrelogramSpectra = nullptr;
        //mCorrelogramSpectraPtr = nullptr;
        mSignalSegment = nullptr;
        mSegmentSpectra = nullptr;
        mSegmentSpectraPtr = nullptr;
        mInputSignals = nullptr;
        mFilteredSignals = nullptr;
        mDenominator = nullptr;
        mSkipZeroSignal = nullptr;
        /// Set everything else to 0
        mSpectraLength = 0;
        mFFTLength = 0;
        mConvolutionLeadingDimension = 0;
        mSpectraLeadingDimension = 0;
        mSamplesLeadingDimension = 0;
        mSamplesExtra = 0;
        mSamples = 0;
        mFilterLength = 0;
        mTemplates = 0;
        mL = 0;
        mInitialized = false;
        mHaveFFTwPlans = false;
        mHaveMatchedFilters = false;
    }
//private:
    /// Holds the options
    MatchedFilterOptions mOptions;
    /// FFT forward plan - maps the input signal to the frequency domain.
    fftwf_plan mForwardPlan;
    /// FFT inverse plan - brings block of templates convolved with 
    /// with input signal back to the time domain.
    fftwf_plan mInversePlan;
    /// The FFT of the time reversed, demeaned, and normalized templates.
    /// This has dimension: [mTemplates x mSpectraLeadingDimension] and is
    /// stored in row major format.
    fftwf_complex *mB = nullptr;
    /// This a pointer to mB which is useful when computing X*B in the 
    /// frequency domain.
    std::complex<float> *mBPtr = nullptr;
    /// The spectra of a segment of the input time domain signals.  This has
    /// dimension [mSpectraLeadingDimension x mTemplates].
    fftwf_complex *mSegmentSpectra = nullptr;
    std::complex<float> *mSegmentSpectraPtr = nullptr;
    /// The matched filtered signals.  This has dimension
    /// [mSamplesLeadingDimension x mTemplates].
    float *mFilteredSignals = nullptr;
    /// The input signals.  This has dimension
    /// [mSamplesLeadingDimension x mTemplates].
    float *mInputSignals = nullptr;
    /// The denominator in the Pearson correlation.  This has dimension
    /// [mSamplesLeadingDimension x mTemplates]
    float *mDenominator = nullptr;
    /// The time domain signal to Fourier transform.  This has dimension
    /// [mSamplesLeadingDimension x mTemplates].
    float *mSignalSegment = nullptr;
    /// Determines if I should skip this signal because it is 0.
    /// This has dimension [mTemplates].
    bool *mSkipZeroSignal = nullptr;
    /// The length of the Fourier transforms.  This should equal 
    /// mFFTLength/2 + 1.
    int mSpectraLength = 0;
    /// The FFT length.  This is the number of zero-padded points in the input
    /// time domain signals and output correlograms.
    int mFFTLength = 0;
    /// The leading dimension of the output convolutions in the window.
    size_t mConvolutionLeadingDimension = 0;
    /// The leading dimension of the spectra.
    size_t mSpectraLeadingDimension = 0;
    /// The leading dimension of the numerator/denominator/filtered signals.
    size_t mSamplesLeadingDimension = 0;
    /// The filter needs to `start-up'.  This requires [nb-1] samples which
    /// we disregard when extracting.
    int mSamplesExtra = 0;
    /// The length of the input signals.
    int mSamples = 0;
    /// The length of the templates (filters).
    int mFilterLength = 0;
    /// The number of templates.
    int mTemplates = 0;
    /// The window length.
    int mL = 0;
    /// Determines if the class is initialized.
    bool mInitialized = false;
    /// Determines if the FFTw plan was made
    bool mHaveFFTwPlans = false;
    /// Determines if the matched filtering has been applied.
    bool mHaveMatchedFilters = false;
};

//----------------------------------------------------------------------------//
//                              End Implementations                           //
//----------------------------------------------------------------------------//
/// Constructor
template<class T>
MatchedFilter<T>::MatchedFilter() :
    pImpl(std::make_unique<MatchedFilterImpl> ())
{
}

/// Destructor
template<class T>
MatchedFilter<T>::~MatchedFilter() = default;

/// Resets class
template<class T> void MatchedFilter<T>::clear() noexcept
{
    pImpl->clear();
}

/// Initializes the FFTs for the double precision matched filtering
template<>
void MatchedFilter<double>::initialize(
    const MatchedFilterOptions &options)
{
    clear();
    // Set the templates
    if (!options.isValid())
    {
        throw std::invalid_argument("Options is invalid\n");
    }
    // Figure out the window length
    pImpl->mSamples = options.getSignalSize();
    pImpl->mFilterLength = options.getMaxTemplateLength();
    pImpl->mSamplesExtra = pImpl->mSamples + pImpl->mFilterLength - 1;
    auto result = computeOptimalFFTAndBlockLength(pImpl->mFilterLength,
                                                  pImpl->mSamplesExtra);
    pImpl->mFFTLength = result.first; //options.getFFTLength();
    pImpl->mL = result.second; //options.getBlockLength();
    // Note that Matlab/Octave assume the signal can be complex.  Hence, they
    // use the FFT.  That's unecessary since our signals and templates are
    // real-valued.  Recall, for real signals the complex spectra 
    // at the negative frequencies equals the conjugate of the complex spectra
    // at the positive frequencies.  Next, using the fact that convolution
    // in the time domain is multiplication in the frequency domain, to handle
    // the complex spectra we have to demonstrate that for the negative 
    // frequencies we implicitly obtain
    // \f$ \bar{X}(\omega) \bar{B}(\omega) = \bar{ X(\omega) B(\omega) } \f$.
    // This says that after multiplying the complex spectra at the positive
    // frequencies we implicitly have the complex spectra at the negative
    // frequencies.  For this to be true we require that
    //   \bar{ (x_r + i x_i) } \bar{ (b_r + i b_i) }
    // = (x_r - i x_i)(b_r - i b_i) 
    // = x_r b_r + x_i b_i - i(x_i b_r -  x_r b_i)
    // = \bar{ (x_r + i x_i) (b_r + i b_i) }
    // which is exactly what we want.
    pImpl->mSpectraLength = pImpl->mFFTLength/2 + 1;
    // Our goal is to correlate.  However, we are performing convolutions.
    // Hence, we have to reverse our templates.
    pImpl->mTemplates = options.getNumberOfTemplates();
    pImpl->mConvolutionLeadingDimension = padLength(pImpl->mFFTLength,
                                                    sizeof(double));
    auto len = pImpl->mConvolutionLeadingDimension
              *static_cast<size_t> (pImpl->mTemplates)
              *sizeof(double);
    // Calloc is important for signal padding
    pImpl->mSignalSegment = static_cast<double *> (MKL_calloc(len, 1, 64));
    auto b = pImpl->mSignalSegment;
    std::vector<double> tData(pImpl->mFilterLength);
    for (int it=0; it<pImpl->mTemplates; ++it)
    {
        auto t = options.getTemplate(it);
        auto offset = static_cast<size_t> (it)
                     *pImpl->mConvolutionLeadingDimension;
        double *tPtr = tData.data();
        int ntlen = t.getSignalLength();
        t.getSignal(pImpl->mFilterLength, &tPtr);
        demeanNormalizeAndReverseTemplate(ntlen, tPtr, &b[offset]);
    }
    // Transform
    constexpr int rank = 1;
    constexpr int istride = 1;
    constexpr int ostride = 1;
    constexpr int inembed[1] = {0};
    constexpr int onembed[1] = {0};
    int nf[1] = {pImpl->mFFTLength};
    pImpl->mSpectraLeadingDimension = padLength(pImpl->mSpectraLength,
                                                sizeof(fftw_complex)); 
    len = pImpl->mSpectraLeadingDimension
         *static_cast<size_t> (pImpl->mTemplates)
         *sizeof(fftw_complex);
    pImpl->mB = static_cast<fftw_complex *> (MKL_calloc(len, 1, 64));
    pImpl->mBPtr = reinterpret_cast<std::complex<double> *> (pImpl->mB);
    pImpl->mSegmentSpectra
        = static_cast<fftw_complex *> (MKL_calloc(len, 1, 64));
    pImpl->mSegmentSpectraPtr
        = reinterpret_cast<std::complex<double> *> (pImpl->mSegmentSpectra);
    // Create the forward plan
    pImpl->mForwardPlan
        = fftw_plan_many_dft_r2c(rank, nf, pImpl->mTemplates,
                                 pImpl->mSignalSegment, inembed,
                                 istride, pImpl->mConvolutionLeadingDimension,
                                 pImpl->mSegmentSpectra, onembed,
                                 ostride, pImpl->mSpectraLeadingDimension,
                                 FFTW_PATIENT);
    fftw_execute_dft_r2c(pImpl->mForwardPlan, b, pImpl->mB);
    // Create the inverse plan
    int ni[1] = {pImpl->mFFTLength};
    pImpl->mInversePlan
        = fftw_plan_many_dft_c2r(rank, ni, pImpl->mTemplates,
                                 pImpl->mSegmentSpectra, inembed,
                                 istride, pImpl->mSpectraLeadingDimension,
                                 pImpl->mSignalSegment, onembed, 
                                 ostride, pImpl->mConvolutionLeadingDimension,
                                 FFTW_PATIENT);
    // We used the space for the input signals as b's workspace - zero it
    len = pImpl->mConvolutionLeadingDimension
         *static_cast<size_t> (pImpl->mTemplates);
#ifdef USE_PSTL
    std::fill(std::execution::unseq,
              pImpl->mSignalSegment, pImpl->mSignalSegment+len, 0);
#else
    std::fill(pImpl->mSignalSegment, pImpl->mSignalSegment+len, 0);
#endif
    // Allocate space for denominator signal and the filtered (numerator)
    // signal, and the filtered signal.
    pImpl->mSamplesLeadingDimension = padLength(pImpl->mSamplesExtra,
                                                sizeof(double));
    len = pImpl->mSamplesLeadingDimension
         *static_cast<size_t> (pImpl->mTemplates)
         *sizeof(double);
    pImpl->mInputSignals = static_cast<double *> (MKL_calloc(len, 1, 64));
    pImpl->mFilteredSignals = static_cast<double *> (MKL_calloc(len, 1, 64));
    pImpl->mDenominator = static_cast<double *> (MKL_calloc(len, 1, 64));
    len = static_cast<size_t> (pImpl->mTemplates)*sizeof(bool);
    pImpl->mSkipZeroSignal = static_cast<bool *> (MKL_calloc(len, 1, 64));
    // Finish up initialization 
    pImpl->mHaveFFTwPlans = true;
    pImpl->mInitialized = true;
}

/// Initializes the FFTs for the float precision matched filtering
template<>
void MatchedFilter<float>::initialize(
    const MatchedFilterOptions &options)
{
    clear();
    // Set the templates
    if (!options.isValid())
    {   
        throw std::invalid_argument("Options is invalid\n");
    }   
    // Figure out the window length
    pImpl->mSamples = options.getSignalSize();
    pImpl->mFilterLength = options.getMaxTemplateLength();
    pImpl->mSamplesExtra = pImpl->mSamples + pImpl->mFilterLength - 1;
    auto result = computeOptimalFFTAndBlockLength(pImpl->mFilterLength,
                                                  pImpl->mSamplesExtra);
    pImpl->mFFTLength = result.first; //options.getFFTLength();
    pImpl->mL = result.second;
    pImpl->mSpectraLength = pImpl->mFFTLength/2 + 1;
    // See double precision for explanation
    pImpl->mTemplates = options.getNumberOfTemplates();
    pImpl->mConvolutionLeadingDimension = padLength(pImpl->mFFTLength,
                                                    sizeof(float));
    auto len = pImpl->mConvolutionLeadingDimension
              *static_cast<size_t> (pImpl->mTemplates)
              *sizeof(float);
    // Calloc is important for padding
    pImpl->mSignalSegment = static_cast<float *> (MKL_calloc(len, 1, 64));
    auto b = pImpl->mSignalSegment;
    std::vector<float> tData(pImpl->mFilterLength);
    for (int it=0; it<pImpl->mTemplates; ++it)
    {
        auto t = options.getTemplate(it);
        auto offset = static_cast<size_t> (it)
                     *pImpl->mConvolutionLeadingDimension;
        float *tPtr = tData.data();
        int ntlen = t.getSignalLength();
        t.getSignal(pImpl->mFilterLength, &tPtr);
        demeanNormalizeAndReverseTemplate(ntlen, tPtr, &b[offset]);
    }
    // Transform
    constexpr int rank = 1;
    constexpr int istride = 1;
    constexpr int ostride = 1;
    constexpr int inembed[1] = {0};
    constexpr int onembed[1] = {0};
    int nf[1] = {pImpl->mFFTLength};
    pImpl->mSpectraLeadingDimension = padLength(pImpl->mSpectraLength,
                                                sizeof(fftwf_complex));
    len = pImpl->mSpectraLeadingDimension
         *static_cast<size_t> (pImpl->mTemplates)
         *sizeof(fftwf_complex);
    pImpl->mB = static_cast<fftwf_complex *> (MKL_calloc(len, 1, 64));
    pImpl->mBPtr = reinterpret_cast<std::complex<float> *> (pImpl->mB);
    pImpl->mSegmentSpectra
        = static_cast<fftwf_complex *> (MKL_calloc(len, 1, 64));
    pImpl->mSegmentSpectraPtr
        = reinterpret_cast<std::complex<float> *> (pImpl->mSegmentSpectra);
    // Create the forward plan
    pImpl->mForwardPlan
        = fftwf_plan_many_dft_r2c(rank, nf, pImpl->mTemplates,
                                  pImpl->mSignalSegment, inembed,
                                  istride, pImpl->mConvolutionLeadingDimension,
                                  pImpl->mSegmentSpectra, onembed,
                                  ostride, pImpl->mSpectraLeadingDimension,
                                  FFTW_PATIENT);
    fftwf_execute_dft_r2c(pImpl->mForwardPlan, b, pImpl->mB);
    // Create the inverse plan
    int ni[1] = {pImpl->mFFTLength};
    pImpl->mInversePlan
        = fftwf_plan_many_dft_c2r(rank, ni, pImpl->mTemplates,
                                  pImpl->mSegmentSpectra, inembed,
                                  istride, pImpl->mSpectraLeadingDimension,
                                  pImpl->mSignalSegment, onembed, 
                                  ostride, pImpl->mConvolutionLeadingDimension,
                                  FFTW_PATIENT);
    // We used the space for the input signals as b's workspace - zero it
    len = pImpl->mConvolutionLeadingDimension
         *static_cast<size_t> (pImpl->mTemplates);
#ifdef USE_PSTL
    std::fill(std::execution::unseq,
              pImpl->mSignalSegment, pImpl->mSignalSegment+len, 0); 
#else
    std::fill(pImpl->mSignalSegment, pImpl->mSignalSegment+len, 0); 
#endif
    // Allocate space for denominator signal and the filtered (numerator)
    // signal, and the filtered signal.
    pImpl->mSamplesLeadingDimension = padLength(pImpl->mSamplesExtra,
                                                sizeof(float));
    len = pImpl->mSamplesLeadingDimension
         *static_cast<size_t> (pImpl->mTemplates)
         *sizeof(float);
    pImpl->mInputSignals = static_cast<float *> (MKL_calloc(len, 1, 64));
    pImpl->mFilteredSignals = static_cast<float *> (MKL_calloc(len, 1, 64));
    pImpl->mDenominator = static_cast<float *> (MKL_calloc(len, 1, 64));
    len = static_cast<size_t> (pImpl->mTemplates)*sizeof(bool);
    pImpl->mSkipZeroSignal = static_cast<bool *> (MKL_calloc(len, 1, 64));
    // Finish up initialization 
    pImpl->mHaveFFTwPlans = true;
    pImpl->mInitialized = true;
}

/// Zeros out the it'th signal
template<class T>
void MatchedFilter<T>::zeroSignal(const int it)
{
    pImpl->mHaveMatchedFilters = false;
    auto nt = getNumberOfTemplates(); // throws
    if (it < 0 || it >= nt) 
    {
        throw std::invalid_argument("it = " + std::to_string(it)
                                  + " must be in the range [0,"
                                  + std::to_string(nt-1) + "]\n");
    }
    auto ns = getSignalLength(); // throws
    auto offset = pImpl->mSamplesLeadingDimension
                 *static_cast<size_t> (it);
    T *__attribute__((aligned(64))) ptr = &pImpl->mInputSignals[offset];
#ifdef USE_PSTL
    std::fill(std::execution::unseq, ptr, ptr+ns, 0);
#else
    std::fill(ptr, ptr+ns, 0);
#endif
    pImpl->mSkipZeroSignal[it] = true;
}

/// Sets the it'th signal
template<class T>
void MatchedFilter<T>::setSignal(const int it, const int nSamples,
                                 const T signal[])
{
    pImpl->mHaveMatchedFilters = false;
    auto nt = getNumberOfTemplates(); // throws
    if (it < 0 || it >= nt)
    {
        throw std::invalid_argument("it = " + std::to_string(it)
                                  + " must be in the range [0,"
                                  + std::to_string(nt-1) + "]\n");
    }
    auto ns = getSignalLength(); // throws
    if (nSamples != ns)
    {
        throw std::invalid_argument("Number of samples = "
                                  + std::to_string(nSamples)
                                  + " must equal" + std::to_string(ns) + "\n");
    }
    auto offset = pImpl->mSamplesLeadingDimension
                 *static_cast<size_t> (it);
    T *__attribute__((aligned(64))) ptr = &pImpl->mInputSignals[offset];
#ifdef USE_PSTL
    std::copy(std::execution::unseq,
              signal, signal+nSamples, ptr);
#else
    std::copy(signal, signal+nSamples, ptr);
#endif
    pImpl->mSkipZeroSignal[it] = false;
}

/// Gets it'th matched filtered signal
template<class T>
std::vector<T> MatchedFilter<T>::getMatchedFilteredSignal(const int it) const
{
    auto ptr = getMatchedFilterSignalPointer(it); // Throws
    auto npts = getFilteredSignalLength();
    std::vector<T> result(ptr, ptr+npts);
    return result;
}

/// Gets a pointer to the it'th matched filtered signal
template<class T>
const T* MatchedFilter<T>::getMatchedFilterSignalPointer(const int it) const
{
    auto nt = getNumberOfTemplates(); // throws
    if (it < 0 || it >= nt) 
    {
        throw std::invalid_argument("it = " + std::to_string(it)
                                  + " must be in the range [0,"
                                  + std::to_string(nt-1) + "]\n");
    }
    auto offset = pImpl->mSamplesLeadingDimension
                 *static_cast<size_t> (it)
                + static_cast<size_t> (pImpl->mFilterLength) - 1;
    const T* ptr = &pImpl->mFilteredSignals[offset];
    return ptr;
}

/// Applies the matched template filtering for double precision
template<>
void MatchedFilter<double>::apply()
{
    pImpl->mHaveMatchedFilters = false;
    if (!isInitialized())
    {
        throw std::runtime_error("Class is not initialized\n");
    }
    size_t nzero = pImpl->mSamplesLeadingDimension
                  *static_cast<size_t> (pImpl->mTemplates);
    std::fill(pImpl->mFilteredSignals, pImpl->mFilteredSignals+nzero, 0);
    // Get the sizes 
    auto L = pImpl->mL; 
    auto nb = pImpl->mFilterLength;
    auto nx = pImpl->mSamplesExtra;
    auto nfft = pImpl->mFFTLength;
    auto nTemplates = pImpl->mTemplates;
    auto spectraLength = pImpl->mSpectraLength;
    auto xnorm = 1/static_cast<double> (pImpl->mFFTLength);

    std::complex<double> *B = pImpl->mBPtr;
    std::complex<double> *X = pImpl->mSegmentSpectraPtr;
    // Loop on the windows - parallelizing requires buffering
    for (int istart=0; istart<nx; istart=istart+L)
    {
        for (int it=0; it<nTemplates; ++it)
        {
            auto ioff = static_cast<size_t> (it)
                       *static_cast<size_t> (pImpl->mSamplesLeadingDimension);
            const double *signalPtr = &pImpl->mInputSignals[ioff+istart];
            ioff = static_cast<size_t> (it)
                  *static_cast<size_t> (pImpl->mConvolutionLeadingDimension); 
            double *__attribute__((aligned(64))) dest
                = &pImpl->mSignalSegment[ioff];
            int iend = std::min(nx, istart + L); // Exclusive
            int ncopy = iend - istart;
#ifdef USE_PSTL
            // Extract the padded signal
            std::copy(std::execution::unseq,
                      signalPtr, signalPtr+ncopy, dest);
            // Zero out signal until end
            std::fill(std::execution::unseq, dest+ncopy, dest+nfft, 0);
#else
            std::copy(signalPtr, signalPtr+ncopy, dest);
            std::fill(dest+ncopy, dest+nfft, 0);
#endif
        }
        // Fourier transform
        fftw_execute_dft_r2c(pImpl->mForwardPlan,
                             pImpl->mSignalSegment,
                             pImpl->mSegmentSpectra);
        // Convolve by multiplying spectra.  Note, that here we have removed
        // the mean from the template and normalized.  Hence, the numerator
        // looks like \tilde{X}*Y.
        for (int it=0; it<nTemplates; ++it)
        {
            auto ioff = static_cast<size_t> (it)
                       *static_cast<size_t> (pImpl->mSpectraLeadingDimension);
            std::complex<double> __attribute__((aligned(64))) *bRow = &B[ioff];
            std::complex<double> __attribute__((aligned(64))) *xcRow = &X[ioff];
            #pragma omp simd
            for (int w=0; w<pImpl->mSpectraLength; ++w)
            {
                xcRow[w] = bRow[w]*xcRow[w];
            }
        }
        // Inverse transform
        fftw_execute_dft_c2r(pImpl->mInversePlan,
                             pImpl->mSegmentSpectra,
                             pImpl->mSignalSegment);
        // Apply overlap and add method
        int yend = std::min(nx, istart+nfft); // Exclusive
        int nupdate = yend - istart;
        for (int it=0; it<nTemplates; ++it)
        {
            // Filtered signal destination
            auto idst = pImpl->mSamplesLeadingDimension
                       *static_cast<size_t> (it)
                      + static_cast<size_t> (istart);
            double *ydst = &pImpl->mFilteredSignals[idst];
            // Contribution from segment
            auto isrc = pImpl->mConvolutionLeadingDimension
                       *static_cast<size_t> (it);
            double *__attribute__((aligned(64))) ysrc
                = &pImpl->mSignalSegment[isrc];
            #pragma omp simd
            for (int i=0; i<nupdate; ++i)
            {
                ydst[i] = ydst[i] + ysrc[i]*xnorm;
            }
        }
    } // Loop on windows
    // Compute the normalization 
    #pragma omp parallel \
     firstprivate(nb) \
     default(none)
    {
    auto nxUnpadded = pImpl->mSamples; // Normalization uses unpadded input pts
    constexpr double tol = 1.e-14;
    #pragma omp for
    for (int it=0; it<pImpl->mTemplates; ++it)
    {
        // Avoid division by zero for obviously dead traces.
        // Since the numerators were zero'd on entry simply skip it.
        if (pImpl->mSkipZeroSignal[it]){continue;}
        // Signal index
        auto isrc = pImpl->mSamplesLeadingDimension
                   *static_cast<size_t> (it);
        // Skip the filter start-up
        auto idst = isrc + static_cast<size_t> (nb) - 1;
        const double *y = &pImpl->mInputSignals[isrc];
        double *yNum = &pImpl->mFilteredSignals[idst];
        normalizeSignal(nxUnpadded, nb, y, yNum, tol);
    }
    } // End parallel
    pImpl->mHaveMatchedFilters = true;
}

/// Applies the matched template filtering for float precision
template<>
void MatchedFilter<float>::apply()
{
    pImpl->mHaveMatchedFilters = false;
    if (!isInitialized())
    {
        throw std::runtime_error("Class is not initialized\n");
    }
    size_t nzero = pImpl->mSamplesLeadingDimension
                  *static_cast<size_t> (pImpl->mTemplates);
    std::fill(pImpl->mFilteredSignals, pImpl->mFilteredSignals+nzero, 0);
    // Get the sizes 
    auto L = pImpl->mL;
    auto nb = pImpl->mFilterLength;
    auto nx = pImpl->mSamplesExtra;
    auto nfft = pImpl->mFFTLength;
    auto nTemplates = pImpl->mTemplates;
    auto spectraLength = pImpl->mSpectraLength;
    auto xnorm = static_cast<float> (1/static_cast<double> (pImpl->mFFTLength));

    std::complex<float> *B = pImpl->mBPtr;
    std::complex<float> *X = pImpl->mSegmentSpectraPtr;
    // Loop on the windows - parallelizing requires buffering
    for (int istart=0; istart<nx; istart=istart+L)
    {
        for (int it=0; it<nTemplates; ++it)
        {
            auto ioff = static_cast<size_t> (it)
                       *static_cast<size_t> (pImpl->mSamplesLeadingDimension);
            const float *signalPtr = &pImpl->mInputSignals[ioff+istart];
            ioff = static_cast<size_t> (it)
                  *static_cast<size_t> (pImpl->mConvolutionLeadingDimension);
            float *__attribute__((aligned(64))) dest
                = &pImpl->mSignalSegment[ioff];
            int iend = std::min(nx, istart + L); // Exclusive
            int ncopy = iend - istart;
#ifdef USE_PSTL
            // Extract the padded signal
            std::copy(std::execution::unseq,
                      signalPtr, signalPtr+ncopy, dest);
            // Zero out signal until end
            std::fill(std::execution::unseq, dest+ncopy, dest+nfft, 0);
#else
            std::copy(signalPtr, signalPtr+ncopy, dest);
            std::fill(dest+ncopy, dest+nfft, 0);
#endif
        }
        // Fourier transform
        fftwf_execute_dft_r2c(pImpl->mForwardPlan,
                              pImpl->mSignalSegment,
                              pImpl->mSegmentSpectra);
        // Convolve by multiplying spectra.  Note, that here we have removed
        // the mean from the template and normalized.  Hence, the numerator
        // looks like \tilde{X}*Y.
        for (int it=0; it<nTemplates; ++it)
        {
            auto ioff = static_cast<size_t> (it)
                       *static_cast<size_t> (pImpl->mSpectraLeadingDimension);
            std::complex<float> __attribute__((aligned(64))) *bRow = &B[ioff];
            std::complex<float> __attribute__((aligned(64))) *xcRow = &X[ioff];
            #pragma omp simd
            for (int w=0; w<pImpl->mSpectraLength; ++w)
            {
                xcRow[w] = bRow[w]*xcRow[w];
            }
        }
        // Inverse transform
        fftwf_execute_dft_c2r(pImpl->mInversePlan,
                              pImpl->mSegmentSpectra,
                              pImpl->mSignalSegment);
        // Apply overlap and add method
        int yend = std::min(nx, istart+nfft); // Exclusive
        int nupdate = yend - istart;
        for (int it=0; it<nTemplates; ++it)
        {
            // Filtered signal destination
            auto idst = pImpl->mSamplesLeadingDimension
                       *static_cast<size_t> (it)
                      + static_cast<size_t> (istart);
            float *ydst = &pImpl->mFilteredSignals[idst];
            // Contribution from segment
            auto isrc = pImpl->mConvolutionLeadingDimension
                       *static_cast<size_t> (it);
            float *__attribute__((aligned(64))) ysrc
                = &pImpl->mSignalSegment[isrc];
            #pragma omp simd
            for (int i=0; i<nupdate; ++i)
            {
                ydst[i] = ydst[i] + ysrc[i]*xnorm;
            }
        }
    } // Loop on windows
    // Compute the normalization 
    #pragma omp parallel \
     firstprivate(nb) \
     default(none)
    {
    auto nxUnpadded = pImpl->mSamples; // Normalization uses unpadded input pts
    constexpr float tol = 1.e-6;
    #pragma omp for
    for (int it=0; it<pImpl->mTemplates; ++it)
    {
        // Avoid division by zero for obviously dead traces.
        // Since the numerators were zero'd on entry simply skip it.
        if (pImpl->mSkipZeroSignal[it]){continue;}
        // Signal index
        auto isrc = pImpl->mSamplesLeadingDimension
                   *static_cast<size_t> (it);
        /// Skip the filter start-up
        auto idst = isrc + static_cast<size_t> (nb) - 1;
        const float *y = &pImpl->mInputSignals[isrc];
        float *yNum = &pImpl->mFilteredSignals[idst];
        normalizeSignal(nxUnpadded, nb, y, yNum, tol);
    }
    } // End parallel
    pImpl->mHaveMatchedFilters = true;
}

/// Determine if the class is initialized
template<class T>
bool MatchedFilter<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Get number of templates
template<class T>
int MatchedFilter<T>::getNumberOfTemplates() const
{
    if (!pImpl->mInitialized)
    {   
        throw std::runtime_error("Class not initialized\n");
    }   
    return pImpl->mTemplates;
}

/// Get signal length
template<class T>
int MatchedFilter<T>::getSignalLength() const
{
    if (!pImpl->mInitialized)
    {   
        throw std::runtime_error("Class not initialized\n");
    }   
    return pImpl->mSamples;
}

/// Get the filtered signal length
template<class T>
int MatchedFilter<T>::getFilteredSignalLength() const
{
    int npts = getSignalLength(); // throws
    int nb = pImpl->mFilterLength;
    return npts - nb + 1; 
}

/// Checks if I have the matched filtered signals 
template<class T>
bool MatchedFilter<T>::haveMatchedFilteredSignals() const noexcept
{
    return pImpl->mHaveMatchedFilters;
}

///-----------------------------Debugging Routines---------------------------///
template<class T>
std::vector<std::complex<T>>
MatchedFilter<T>::getSpectraOfTemplate(const int it) const
{
    if (!isInitialized())
    {
        throw std::runtime_error("Class not initialized\n");
    }
    if (it < 0 || it >= pImpl->mTemplates)
    {
        throw std::invalid_argument("Template = " + std::to_string(it)
                                  + " must be in range [0,"
                                  + std::to_string(pImpl->mTemplates) + "]\n");
    }
    std::vector<std::complex<T>> wTemplate(pImpl->mSpectraLength);
    int i0 = pImpl->mSpectraLeadingDimension*it;
    int i1 = i0 + pImpl->mSpectraLength;
    std::copy(pImpl->mBPtr+i0, pImpl->mBPtr+i1, wTemplate.begin());
    return wTemplate;
}

/// Template class instantiation
template class MatchedFilter<double>;
template class MatchedFilter<float>;
