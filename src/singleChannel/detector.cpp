#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <numeric>
#include <mkl.h>
#include "mflib/waveformTemplate.hpp"
#include "mflib/singleChannel/detector.hpp"
#include "mflib/singleChannel/detection.hpp"
#include "mflib/singleChannel/matchedFilter.hpp"
#include "private/peakFinder.hpp"

using namespace MFLib::SingleChannel;

namespace
{
/// This function fits a second order polynomial through 3 regularly spaced
/// points and finds the optimum.
template<class T>
T quadraticRefinement(const T y1, const T y2, const T y3, const T dt)
{
    auto num = dt*(y1 - y3);
    auto den = 2*(y1 - 2*y2 + y3);
    return num/den;
}

template<class T>
T quadraticRefinement(const int n, 
                      const T xc[],
                      const int optIndex,
                      const double dt,
                      const bool lAbs)
{
    // Typical interpolation
    if (optIndex > 0 && optIndex < n - 1)
    {
        double y1 = xc[optIndex-1];
        double y2 = xc[optIndex];
        double y3 = xc[optIndex+1];
        if (lAbs)
        {
            y1 = std::abs(y1);
            y2 = std::abs(y2);
            y3 = std::abs(y3);
        } 
        // If this isn't the set's maximum then return no perturbation
        if (y2 < y1 || y2 < y3)
        {
            // This is the set's minimum.  Recast as optimization by negating.
            if (y2 < y1 && y2 < y3)
            {
                y1 =-y1;
                y2 =-y2;
                y3 =-y3;
            }
            else
            {
                return 0; 
            }
        }
        // Avoid a division by zero
        y1 = std::min(y2-100*std::numeric_limits<double>::epsilon(), y1);
        y3 = std::min(y2-100*std::numeric_limits<double>::epsilon(), y3);
        // Now interpolate
        return quadraticRefinement(y1, y2, y3, dt);
    }
    // If extrema is an endpoint then simply return the endpoint
    else
    {
        return 0;
    }
} 
                      

}

template<class T>
class Detector<T>::DetectorImpl
{
public:
    /// Destructor
    ~DetectorImpl()
    {
        clear();
    }
    /// Release memory
    void clear() noexcept
    {
        if (mMaxXC){MKL_free(mMaxXC);}
        if (mMaxTemplate){MKL_free(mMaxTemplate);}
        mMaxXC = nullptr;
        mMaxTemplate = nullptr; 
        mMaxSizeXC = 0;
        mPeakFinder.clear();
        mUseAbsoluteValue = false;
        mInitialized = false;
    }

    T *__attribute__((aligned(64))) mMaxXC = nullptr;
    T *__attribute__((aligned(64))) mMaxTemplate = nullptr;
    std::vector<MFLib::SingleChannel::Detection<T>> mDetections;
    MFLib::PeakFinder<T> mPeakFinder;
    int mMaxSizeXC = 0;
    bool mUseAbsoluteValue = false;
    bool mInitialized = false;
};

/// Constructor
template<class T>
Detector<T>::Detector() :
    pImpl(std::make_unique<Detector<T>::DetectorImpl> ())
{
}

/// Destructor
template<class T>
Detector<T>::~Detector() = default;

/// Clear the class
template<class T>
void Detector<T>::clear() noexcept
{
    pImpl->clear();
}

template<class T>
const MFLib::SingleChannel::Detection<T>& Detector<T>::operator[]
    (const size_t i) const
{
    if (i >= pImpl->mDetections.size())
    {
        throw std::invalid_argument("i is out of bounds\n");
    }
    return pImpl->mDetections[i];
}

/// Initialize

/// Compute the detections
template<class T>
void Detector<T>::detect(const MFLib::SingleChannel::MatchedFilter<T> &xc)
{
    pImpl->mDetections.clear();
    if (!xc.haveMatchedFilteredSignals())    
    {
        throw std::invalid_argument("Matched filtered signals not computed\n");
    }
    // Compute the detections
    auto detectionLength = xc.getFilteredSignalLength();
    auto nt = xc.getNumberOfTemplates();
    // Set space for results
    auto nbytes = sizeof(T)*static_cast<size_t> (detectionLength);
    T *__attribute__((aligned(64))) det
        = static_cast<T *> (MKL_calloc(nbytes, 1, 64));
    nbytes = sizeof(int)*static_cast<size_t> (detectionLength);
    int *__attribute__((aligned(64))) id
        = static_cast<int *> (MKL_calloc(nbytes, 1, 64));
    // For detections we only care about the top performers of all the templates
    for (int it=0; it<nt; ++it)
    {
        const T *__attribute__((aligned(64))) 
        xcPtr = xc.getMatchedFilterSignalPointer(it);
        if (pImpl->mUseAbsoluteValue)
        {
            #pragma omp simd aligned(xcPtr, det, id: 64)
            for (int i=0; i<detectionLength; ++i) 
            {
                auto absXC = std::abs(xcPtr[i]);
                if (absXC > det[i])
                {
                    det[i] = absXC;
                    id[i] = it;
                }
            }
        }
        else
        {
            #pragma omp simd aligned(xcPtr, det, id: 64)
            for (int i=0; i<detectionLength; ++i) 
            {
                if (xcPtr[i] > det[i])
                {
                    det[i] = xcPtr[i];
                    id[i] = it;
                }
            }
        }
    }
    // Compute the peaks which are the detections
    pImpl->mPeakFinder.setSignal(detectionLength, det);
    pImpl->mPeakFinder.apply();
    auto nDetections = pImpl->mPeakFinder.getNumberOfPeaks();
    if (nDetections < 1){return;} // Swing and a miss
    // Copy the templates
    auto nTemplates = xc.getNumberOfTemplates();
    std::vector<MFLib::WaveformTemplate> templates;
    for (int it=0; it<nTemplates; ++it)
    {
        templates[it] = xc.getWaveformTemplate(it);
    }
    // Get the detection indices
    auto peaksIndices = pImpl->mPeakFinder.getPeakIndicesPointer();
    // Create detections
    pImpl->mDetections.resize(nDetections); 
    auto signalPtr = xc.getSignalPointer();
    auto signalLength = xc.getSignalLength();
    for (int i=0; i<nDetections; ++i)
    {
        // Release stale memory
        pImpl->mDetections[i].clear();
        // Get the detection index in the signal and the corresponding template
        int peakIndex = peaksIndices[i];
        int it = id[peakIndex]; 
        // Get the corresponding chunk of signal 
        int templateLength = templates[it].getSignalLength();
        // Stay in bounds
        if (peakIndex + templateLength > signalLength)
        {
            templateLength = signalLength - peakIndex;
        }
        pImpl->mDetections[i].setDetectedSignal(templateLength,
                                                signalPtr + peakIndex);
        // Compute the onset time
        auto dt = 1/templates[it].getSamplingRate(); // Required
        auto detectionTime = peakIndex*dt;
        pImpl->mDetections[i].setDetectionTime(detectionTime);
        // Compute the interpolated onset time
        auto xcPtr = xc.getMatchedFilterSignalPointer(it);
        auto shift = quadraticRefinement(detectionLength,
                                         xcPtr,
                                         peakIndex,
                                         dt,
                                         pImpl->mUseAbsoluteValue);
        auto intDetTime = detectionTime + shift;
        pImpl->mDetections[i].setDetectionTime(detectionTime);
        pImpl->mDetections[i].setInterpolatedDetectionTime(intDetTime);
        // Try getting the phase onset time
        if (templates[it].havePhaseOnsetTime())
        {
            auto pickTime = detectionTime + templates[it].getPhaseOnsetTime();
            auto intPickTime = intDetTime
                             + templates[it].getPhaseOnsetTime();
            pImpl->mDetections[i].setPhaseOnsetTime(pickTime);
            pImpl->mDetections[i].setInterpolatedPhaseOnsetTime(intPickTime);
        }
    }
    // Cleanup
    MKL_free(det);
    MKL_free(id);
}

/// Template instantation
template class MFLib::SingleChannel::Detector<double>;
template class MFLib::SingleChannel::Detector<float>;
