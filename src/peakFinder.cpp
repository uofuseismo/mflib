#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include <limits>
#include <vector>
#include <functional>
#if __has_include(<pstl/algorithm>)
#include <pstl/algorithm>
#include <pstl/execution>
#include <pstl/numeric>
#define USE_PSTL 1
#else
#include <algorithm>
#endif
#include <numeric>
#include <mkl.h>
#include "private/peakFinder.hpp"

/// TODO: Implement this as at some point?
/// https://www.mathworks.com/matlabcentral/fileexchange/25500-peakfinder-x0-sel-thresh-extrema-includeendpoints-interpolate
#define DEFAULT_THRESH 0.7

using namespace MFLib;

namespace 
{

/// Zeros an array
template<class T>
void zero(const int n, T *x)
{
    const T zero = 0;
#ifdef USE_PSTL
    std::fill(std::execution::unseq, x, x+n, zero);
#else
    std::fill(x, x+n, zero);
#endif
}

/// Copies an array
template<class T>
void copy(const int n, const T x[], T y[])
{
#ifdef USE_PSTL
    std::copy(std::execution::unseq, x, x+n, y);
#else
    std::copy(x, x+n, y);
#endif
}

template<class T>
void signFn(const int n, const T x[], int y[])
{
    const T zero = 0;
    #pragma omp simd
    for (int i=0; i<n; ++i)
    {
        y[i] = 0;
        if (x[i] > zero){y[i] = 1;}
        if (x[i] < zero){y[i] =-1;}
    }
}

/// Differentiates an array x
template<class T>
void diff(const int n, const T *x, T *y)
{
    #pragma omp simd
    for (int i=0; i<n-1; ++i)
    {
        y[i] = x[i+1] - x[i];
    }
}

}

/// Holds the variables
template<class T>
class PeakFinder<T>::PeakFinderImpl
{
public:
    PeakFinderImpl() = default;
    /// Copy c'tor
    PeakFinderImpl(const PeakFinderImpl &impl)
    {
        *this = impl;
    }
    /// Copy
    PeakFinderImpl& operator=(const PeakFinderImpl &impl)
    {
        mThreshold = impl.mThreshold;
        mMinPeakDistance = impl.mMinPeakDistance;
        mLength = impl.mLength;
        mNumberOfPeaks = impl.mNumberOfPeaks;
        if (mLength > 0)
        {
            size_t nbytes = 0;
            if (impl.mX)
            {
                nbytes = sizeof(T)*static_cast<size_t> (mLength);
                mX = static_cast<T *> (MKL_calloc(nbytes, 1, 64));
                copy(mLength, impl.mX, mX);
            }
            if (impl.mDiffX)
            {
                nbytes = sizeof(T)*static_cast<size_t> (mLength);
                mDiffX = static_cast<T *> (MKL_calloc(nbytes, 1, 64));
                copy(mLength, impl.mDiffX, mDiffX);
            }
            if (impl.mPeakIndices)
            {
                nbytes = sizeof(int)*static_cast<size_t> (mLength);
                mPeakIndices = static_cast<int *> (MKL_calloc(nbytes, 1, 64));
                copy(mLength, impl.mPeakIndices, mPeakIndices);
            }
            if (impl.mSign)
            {
                nbytes = sizeof(int)*static_cast<size_t> (mLength);
                mSign = static_cast<int *> (MKL_calloc(nbytes, 1, 64));
                copy(mLength, impl.mSign, mSign);
            }
        }
        return *this;
    }
    /// Destructor
    ~PeakFinderImpl()
    {
        clear();
    }
    /// Releases memory and clears class
    void clear() noexcept
    {
        if (mX){MKL_free(mX);}
        if (mDiffX){MKL_free(mDiffX);}
        if (mSign){MKL_free(mSign);}
        if (mPeakIndices){MKL_free(mPeakIndices);}
        mX = nullptr;
        mPeakIndices = nullptr;
        mThreshold = DEFAULT_THRESH;
        mMinPeakDistance = 0;
        mLength = 0;
        mNumberOfPeaks = 0;
    }
    /// The time series of which we are to find the peaks.  This is an 
    /// array of dimension [mLength].
    T *__attribute__((aligned(64))) mX = nullptr;
    /// The derivative of the time series.  This uses a backwards difference
    /// on all elements with the exception of the first element which uses
    /// a forwards difference. This is an array whose dimension is [mLength].
    T *__attribute__((aligned(64))) mDiffX = nullptr;
    /// The indices of the peaks in the signal.  This is an array of 
    /// dimension [mLength] however only the first mNumberOfPeaks samples
    /// is defined.
    int *__attribute__((aligned(64))) mPeakIndices = nullptr;
    /// Holds the sign of the derivative.  This is an array whose dimension
    /// is [mLength].
    int *__attribute__((aligned(64))) mSign = nullptr;
    /// The peak finding threshold
    T mThreshold = DEFAULT_THRESH;
    /// The peaks must be at least this many samples apart to be considered.
    /// Otherwise, the policy will be to use the largest peak in the window.
    int mMinPeakDistance = 0;
    /// Length of x 
    int mLength = 0;
    /// The number of peaks
    int mNumberOfPeaks = 0;
};

/// Constructor
template<class T>
PeakFinder<T>::PeakFinder() : 
    pImpl(std::make_unique<PeakFinderImpl> ())
{
}

/// Copy c'tor
template<class T>
PeakFinder<T>::PeakFinder(const PeakFinder &peaks)
{
    *this = peaks;
}

/// Move c'tor
template<class T>
PeakFinder<T>::PeakFinder(PeakFinder &&peaks) noexcept
{
    *this = std::move(peaks);
}

/// Copy assignment operator
template<class T>
PeakFinder<T>& PeakFinder<T>::operator=(const PeakFinder &peaks)
{
    if (&peaks == this){return *this;}
    pImpl = std::make_unique<PeakFinderImpl> (*peaks.pImpl);
    return *this;
}

/// Move assignment operator
template<class T>
PeakFinder<T>& PeakFinder<T>::operator=(PeakFinder &&peaks) noexcept
{
    if (&peaks == this){return *this;}
    pImpl = std::move(peaks.pImpl);
    return *this;
}

/// Destructor
template<class T>
PeakFinder<T>::~PeakFinder() = default;

/// Sets the threshold
template<class T>
void PeakFinder<T>::setThreshold(const T thresh) 
{
    pImpl->mNumberOfPeaks = 0;
    if (thresh < 0 || thresh > 1)
    {
        throw std::invalid_argument("Threshold = " + std::to_string(thresh)
                                  + " must be in the range [0,1]\n");
    }
    pImpl->mThreshold = thresh;
}

/// Gets the threshold 
template<class T>
T PeakFinder<T>::getThreshold() const noexcept
{
    return pImpl->mThreshold;
}

/// Gets the minimim window size
template<class T>
int PeakFinder<T>::getMinimumPeakDistance() const noexcept
{
    return pImpl->mMinPeakDistance;
}

template<class T>
void PeakFinder<T>::setMinimumPeakDistance(const int minDistance) noexcept
{
    pImpl->mNumberOfPeaks = 0; 
    pImpl->mMinPeakDistance = minDistance;
}

/// Sets the signal
template<class T>
void PeakFinder<T>::setSignal(const int n, const T x[])
{
    // Check inputs
    pImpl->mNumberOfPeaks = 0;
    if (n < 1 || x == nullptr)
    {
        if (n < 1){throw std::invalid_argument("No points in x\n");}
        throw std::invalid_argument("x is NULL\n");
    }
    // Resize
    if (n > pImpl->mLength)
    {
        if (pImpl->mX){MKL_free(pImpl->mX);}
        if (pImpl->mPeakIndices){MKL_free(pImpl->mPeakIndices);}
        size_t nbytes = sizeof(T)*static_cast<size_t> (n);
        pImpl->mX = static_cast<T *> (MKL_calloc(nbytes, 1, 64));
        pImpl->mDiffX = static_cast<T *> (MKL_calloc(nbytes, 1, 64));
        nbytes = sizeof(int)*static_cast<size_t> (n);
        pImpl->mPeakIndices = static_cast<int *> (MKL_calloc(nbytes, 1, 64));
        pImpl->mSign = static_cast<int *> (MKL_calloc(nbytes, 1, 64));
    }
    pImpl->mLength = n;
    copy(n, x, pImpl->mX);
}

/// Get the number of peaks
template<class T>
int PeakFinder<T>::getNumberOfPeaks() const noexcept 
{
    return pImpl->mNumberOfPeaks;
}

/// Get the peaks
template<class T>
void PeakFinder<T>::getPeaks(const int maxx, int peaks[]) const
{
    auto nPeaks = getNumberOfPeaks();
    if (maxx < nPeaks || peaks == nullptr)
    {
        if (maxx < nPeaks)
        {
            throw std::invalid_argument("maxx = " + std::to_string(maxx)
                                      + " must be at least = "
                                      + std::to_string(nPeaks) + "\n");
        }
        throw std::invalid_argument("peaks is NULL\n");
    }
    if (nPeaks < 1){return;} // Nothing to do
    copy(nPeaks, pImpl->mPeakIndices, peaks);
}

/// Does the hardwork of computing the peaks.
/// This works by:
///   (1) Finding the local maxima that exceed some threshold in the
///       detection signal.
///   (2) Sorting the detections in descending order.
///   (3) Saving a detection when it is the largest trigger that is not
///       within some pre-defined distance of an existing trigger.
template<class T>
void PeakFinder<T>::apply()
{
    // Initialize to zero
    auto peakIndices = pImpl->mPeakIndices;
    pImpl->mNumberOfPeaks = 0;
    if (peakIndices){zero(pImpl->mLength, peakIndices);}
    auto n = pImpl->mLength;
    // No data - i.e., nothing to do
    if (n < 1){return;}
    // Possibly no triggers
    auto x = pImpl->mX;
    auto maxElement = std::max_element(x, x+n);
    if (*maxElement < pImpl->mThreshold){return;}
    // Compute derivatives of the input signal
    auto nbytes = static_cast<size_t> (n)*sizeof(T);
    auto dx = pImpl->mDiffX;
    // Differentiate x.
#ifdef USE_PSTL
    std::adjacent_difference(std::execution::unseq, x, x+n, dx);
#else
    std::adjacent_difference(x, x+n, dx);
#endif
    // Maxima must exceed a tolerance
    auto tol = pImpl->mThreshold;
    // Take advantage of the superfluous first element and compute a forwards 
    // difference
    dx[0] = x[1] - x[0];
    // The sign of the first derivative transitioning from a positive to a 
    // negative gives us a local optima.
    auto sign = pImpl->mSign;
    signFn(n, dx, sign);
    // First element: If it is decreasing then it is a candidate for an optimum.
    std::vector<std::pair<T, int>> maxes;
    maxes.reserve(n);
    if (sign[0] < 0 && x[0] > tol)
    {
        maxes.push_back(std::make_pair(x[0], 0));
    }
    // General case: Sign goes from positive to negative.
    for (int i=1; i<n-1; ++i)
    {
        if (sign[i] >= 0 && sign[i+1] < 0)
        {
            if (x[i] > tol)
            {
                maxes.push_back(std::make_pair(x[i], i));
            }
        }
    }
    // Last element:  If it is increasing then it is a candidate for an optimum.
    if (sign[n-1] > 0 && x[n-1] > tol)
    {
        maxes.push_back(std::make_pair(x[n-1], n-1));
    }
    // If there are no maxes then exit early
    if (maxes.size() == 0)
    {
        pImpl->mNumberOfPeaks = 0;
        return;
    }
    // Handle edge case of one max 
    if (maxes.size() == 1)
    {
        pImpl->mNumberOfPeaks = 1;
        peakIndices[0] = maxes[0].second;
    }
    // Proceed with general argument
    std::sort(maxes.begin(), maxes.end(), 
              [](const std::pair<T, int> &a, const std::pair<T, int> &b)
              {
                 return a.first > b.first; // Sort descending order
              });
    // This is the pruning stage
    auto minPeakDistance = pImpl->mMinPeakDistance;
    if (minPeakDistance > 0)
    { 
        // First max is definitely a keeper
        int nPeaks = 1;
        peakIndices[0] = maxes[0].second;
        // Now we only keep maxes within a winner
        for (auto &p : maxes)
        {
            auto index = p.second;
            // Look for the closest trigger in my expanding trigger list
            auto dmin = std::abs(peakIndices[0] - index);
            #pragma omp simd reduction(min: dmin)
            for (int i=1; i<nPeaks; ++i)
            {
                dmin = std::min(index, std::abs(peakIndices[i] - index));
            }
            // If this peak is sufficiently far from the nearest trigger
            // then keep it
            if (dmin > minPeakDistance)
            {
                peakIndices[nPeaks] = index;
                nPeaks = nPeaks + 1;
            }
        }
        pImpl->mNumberOfPeaks = nPeaks;
    }
    else
    {
        // No pruning - this is a straight copy of the peaks
        pImpl->mNumberOfPeaks = static_cast<int> (maxes.size());
        #pragma omp simd
        for (int i=0; i<pImpl->mNumberOfPeaks; ++i)
        {
            peakIndices[i] = maxes[i].second;
        }
    }
}
/*
template<class T>
void PeakFinder<T>::apply()
{
    // Initialize to zero
    pImpl->mNumberOfPeaks = 0;
    if (pImpl->mPeakIndices){zero(pImpl->mLength, pImpl->mPeakIndices);}
    auto len = pImpl->mLength;
    if (len < 1){return;} // Nothing to do
    // Define threshold
    T sign = 1;
    if (!pImpl->mFindMaxima){sign =-1;}
    auto thresh = pImpl->mThreshold;
    // Handle the length-one edge case
    if (len == 1)
    {
        if (pImpl->mX[0] > thresh && pImpl->mIncludeEndPoints)
        {
            pImpl->mPeakIndices[0] = 0;
            pImpl->mNumberOfPeaks = 1;
        }
        return;
    }
    // Handle the length-two edge case
    if (len == 2)
    {
        if (pImpl->mIncludeEndPoints)
        {
            if (pImpl->mX[0] > thresh || pImpl->mX[1] > thresh)
            {
                pImpl->mPeakIndices[0] = 0;
                if (pImpl->mX[1] > pImpl->mX[0]){pImpl->mPeakIndices[0] = 1;}
                pImpl->mNumberOfPeaks = 1;
            }
        }
        return;
    }
    auto pair = std::minmax_element(pImpl->mX, pImpl->mX+len);
    T sel = (*pair.second - *pair.first)/4;
    // Allocate space
    std::vector<T> dx0;
    std::vector<int> ind;
    dx0.reserve(len-1);
    ind.reserve(len+2);
    // Do the general algorithm: Step 1 - Compute gradient 
    T *x0 = pImpl->mX;
    diff(len, x0, dx0.data());
    // If the gradient is constant
    #pragma omp simd
    for (int i=0; i<len; ++i)
    {
        if (dx0[i] == 0){dx0[i] =-std::numeric_limits<T>::epsilon();}
    }
    // Find candidate peaks, valleys, and possible end points
    if (pImpl->mIncludeEndPoints)
    {
        ind.push_back(0);
    }
    // Find sign changes to the derivatives
    for (int i=0; i<len-1; ++i)
    {
        // Product of differing sign denotes sign change 
        if (dx0[i]*dx0[i+1] < 0)
        {
            ind.push_back(i+1);
        }
    }
    if (pImpl->mIncludeEndPoints)
    {
        ind.push_back(len-1);
    }
    // Copy the candidate peaks, valleys, and possible end points
    len = static_cast<int> (ind.size());
    auto minMag = std::numeric_limits<T>::max();
    std::vector<T> x(len);
    #pragma omp simd reduction(min: minMag)
    for (int i=0; i<len; ++i)
    {
        x[i] = x0[ind[i]];
        minMag = std::min(x[i], minMag); 
    }
    auto leftMin = minMag;

    std::vector<T> peakMags;
    std::vector<int> peakInds;
    peakMags.reserve(len);
    peakInds.reserve(len);
    // Find the peaks from the candidate list
    if (len > 2)
    {
        // Set initial parameters for loop
        auto tempMag = minMag;
        auto foundPeak = false;
        // Since we tacked the first point on, we need take care with the 
        // derivative.
        if (pImpl->mIncludeEndPoints)
        {
            T dx3[3], signDx[3];
            diff(3, x.data(), dx3);
            signFn(3, dx3, signDx);
            if (signDx[0] <= 0)
            {
                // Want alternating signs
                if (signDx[1] == signDx[2]) 
                {
                    x.erase(x.begin()+2);
                    ind.erase(ind.begin()+2);
                    len = len - 1;
                }
            }
            else
            {
                if (signDx[1] == signDx[2])
                {
                    x.erase(x.begin()+1);
                    ind.erase(ind.begin()+1);
                    len = len - 1;
                }
            }
        }
        // Skip the first point if it is smaller than the min
        int i1 = 1;
        if (x[0] >= x[1]){i1 = 0;}
        // Preallocate max number of maxima
        auto maxPeaks = static_cast<int> (std::ceil(len/2.0)); 
        

        // Create output
    }
    else
    {
        int xInd = std::distance(x.begin(),
                                 std::max_element(x.begin(), x.end()));
        auto peakMag = x[xInd]; 
        if (pImpl->mIncludeEndPoints && peakMag > minMag + sel)
        {
            peakMags.push_back(peakMag);
            peakInds.push_back(ind[xInd]);
        }
    }
    // Apply thresholding?
}
*/

/// Template instantiation
template class MFLib::PeakFinder<double>;
template class MFLib::PeakFinder<float>;
