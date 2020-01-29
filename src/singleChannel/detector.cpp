#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include <set>
#include <mkl.h>
#if __has_include(<pstl/algorithm>)
#include <pstl/algorithm>
#include <pstl/execution>
#define USE_PSTL
#endif
#include "mflib/waveformTemplate.hpp"
#include "mflib/singleChannel/detector.hpp"
#include "mflib/singleChannel/detectorParameters.hpp"
#include "mflib/singleChannel/detection.hpp"
#include "mflib/singleChannel/matchedFilter.hpp"
#include "mflib/singleChannel/relativeMagnitude.hpp"
#include "private/peakFinder.hpp"

using namespace MFLib::SingleChannel;

namespace
{

struct Amplitude
{
    bool mHave = false;
};

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

std::vector<int> getUniqueTemplateIDs(const int nDetections,
                                      const int detectionIndex[],
                                      const int nSamples,
                                      const int ids[])
{
    std::vector<int> uniqueIDs(nDetections+1);
    #pragma omp simd
    for (int i=0; i<nDetections; ++i)
    {
        int index = detectionIndex[i];
        uniqueIDs[i] = ids[index];
    }
    // Sort
#ifdef USE_PSTL
    std::sort(std::execution::unseq, uniqueIDs.begin(), uniqueIDs.end());
#else
    std::sort(uniqueIDs.begin(), uniqueIDs.end());
#endif
    uniqueIDs[nDetections] = uniqueIDs[0] - 1; 
    // Count the unique indices
    int n = 0;
    for (int i=0; i<nDetections; ++i)
    {
        if (uniqueIDs[i+1] != uniqueIDs[i]){n = n + 1;}
    }
    // Make a set of unique IDs
    std::vector<int> result(n);
    int j = 0;
    for (int i=0; i<nDetections; ++i)
    {
        if (uniqueIDs[i+1] != uniqueIDs[i])
        {
            result[j] = uniqueIDs[i];
            j = j + 1;
        }
    }
    return result;
}

int getTemplateIndex(const int val, const std::vector<int> &v)  
{
    auto low = std::lower_bound(v.begin(), v.end(), val);
    int index = std::distance(v.begin(), low); 
    if (v[index] != val)
    {
        fprintf(stderr, "val = %d not found in v\n", val);
        return -1; 
    }
    return index;
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
        mDetections.clear();
        mPeakFinder.clear(); 
        mParameters.clear();
        mDetectionLength = 0;
        mUseAbsoluteValue = false;
        mHaveMatchedFilteredSignals = false;
        mInitialized = false;
    }

    std::vector<MFLib::SingleChannel::Detection<T>> mDetections;
    MFLib::PeakFinder<T> mPeakFinder;
    MFLib::SingleChannel::DetectorParameters mParameters;
 
    int mDetectionLength = 0;
    int mSignalLength = 0;
    bool mUseAbsoluteValue = false;
    bool mHaveMatchedFilteredSignals = false;
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
template<class T>
void Detector<T>::initialize(
    const MFLib::SingleChannel::DetectorParameters &parameters)
{
    clear();
    // Copy the parameters
    pImpl->mParameters = parameters;
    pImpl->mUseAbsoluteValue = true;
    if (pImpl->mParameters.getMaximaPolicy() ==
        MaximumMatchedFilterPolicy::MAXIMUM)
    {
        pImpl->mUseAbsoluteValue = false;
    }
    // Set the tolerance on detections
    auto tol = pImpl->mParameters.getDetectionThreshold();
    if (tol < 0){tol = 0;}
    if (tol > 1){tol = 1;}
    pImpl->mPeakFinder.setThreshold(tol);
    // Set minimum spacing between detections
    int minSpacing = pImpl->mParameters.getMinimumDetectionSpacing();
    pImpl->mPeakFinder.setMinimumPeakDistance(minSpacing);
/*
    // Copy the templates
    auto nTemplates = mf.getNumberOfTemplates();
    if (nTemplates < 1)
    {
        throw std::invalid_argument("No templates on detector\n");
    }
    pImpl->mTemplates.resize(nTemplates);
    for (int it=0; it<nTemplates; ++it)
    {
        pImpl->mTemplates[it] = mf.getWaveformTemplate(it);
    }
    // Initialize the corresponding relative amplitude classes
    pImpl->mMagnitudes.resize(nTemplates);
    for (int it=0; it<nTemplates; ++it)
    {
        pImpl->mMagnitudes[it].initialize(pImpl->mTemplates[it]);
    } 
*/
    // Class is ready to go
    pImpl->mInitialized = true;
}

/// Sets the matched filtered signal
/*
template<class T>
void Detector<T>::setMatchedFilteredSignals(
    const MFLib::SingleChannel::MatchedFilter<T> &mf)
{
    pImpl->mHaveMatchedFilteredSignals = false;
    pImpl->mDetections.clear();
    if (!isInitialized())
    {
        throw std::invalid_argument("Detector class not initialized\n");
    }
    if (!mf.haveMatchedFilteredSignals())    
    {
        throw std::invalid_argument("Matched filtered signals not computed\n");
    }
    auto nt = mf.getNumberOfTemplates();
    auto ntRef = static_cast<int> (pImpl->mTemplates.size());
    if (nt != ntRef)
    {
        throw std::invalid_argument("Number of templates in matched filter = "
                                  + std::to_string(nt)
                                  + " must = " + std::to_string(ntRef) + "\n");
    }
    // Copy the signal to filter
    auto signalLength = mf.getSignalLength();
    if (signalLength > pImpl->mSignalLength)
    {
        if (pImpl->mSignal){MKL_free(pImpl->mSignal);}
        auto nbytes = sizeof(T)*static_cast<size_t> (signalLength);
        pImpl->mSignal = static_cast<T *> (MKL_calloc(nbytes, 1, 64));
    } 
    auto signalToCopy = mf.getSignalPointer();
    std::copy(signalToCopy, signalToCopy+signalLength, pImpl->mSignal);
    pImpl->mSignalLength = signalLength;
    // Resize
    auto detectionLength = mf.getFilteredSignalLength();
    if (detectionLength < 1)
    {
        throw std::invalid_argument("Matched filtered signals are length 0\n");
    }
    if (detectionLength > pImpl->mDetectionLength)
    {
        if (pImpl->mMaxMF){MKL_free(pImpl->mMaxMF);}
        if (pImpl->mMaxTemplate){MKL_free(pImpl->mMaxTemplate);}
        auto nbytes = sizeof(T)*static_cast<size_t> (detectionLength);
        pImpl->mMaxMF = static_cast<T *> (MKL_calloc(nbytes, 1, 64));
        nbytes = sizeof(int)*static_cast<size_t> (detectionLength);
        pImpl->mMaxTemplate = static_cast<int *> (MKL_calloc(nbytes, 1, 64));
    }
    pImpl->mDetectionLength = detectionLength;
    auto det = pImpl->mMaxMF;
    auto id = pImpl->mMaxTemplate;
    // For detections we only care about the top performers of all the templates
    for (int it=0; it<nt; ++it)
    {
        const T *__attribute__((aligned(64))) 
        mfPtr = mf.getMatchedFilterSignalPointer(it);
        if (pImpl->mUseAbsoluteValue)
        {
            #pragma omp simd aligned(mfPtr, det, id: 64)
            for (int i=0; i<detectionLength; ++i) 
            {
                auto absXC = std::abs(mfPtr[i]);
                if (absXC > det[i])
                {
                    det[i] = absXC;
                    id[i] = it;
                }
            }
        }
        else
        {
            #pragma omp simd aligned(mfPtr, det, id: 64)
            for (int i=0; i<detectionLength; ++i) 
            {
                if (mfPtr[i] > det[i])
                {
                    det[i] = mfPtr[i];
                    id[i] = it;
                }
            }
        }
    }
    pImpl->mHaveMatchedFilteredSignals = true;
}
*/

/// Compute the detections
template<class T>
void Detector<T>::detect(const MFLib::SingleChannel::MatchedFilter<T> &mf)
{
    // Clear the output
    pImpl->mDetections.clear();
    // Check state of class and inputs
    if (!isInitialized())
    {
        throw std::invalid_argument("Detector class not initialized\n");
    }
    if (!mf.haveMatchedFilteredSignals())    
    {
        throw std::invalid_argument("Matched filtered signals not computed\n");
    }
    auto nt = mf.getNumberOfTemplates();
    if (nt < 1)
    {
        throw std::invalid_argument("No templates on matched filter class\n");
    }
    // Reduce the detections (by taking the maximum)
    auto detectionLength = pImpl->mDetectionLength;
    auto nbytes = sizeof(T)*static_cast<size_t> (detectionLength);
    auto det = static_cast<T *> (MKL_calloc(nbytes, 1, 64));
    nbytes = sizeof(int)*static_cast<size_t> (detectionLength);
    auto id = static_cast<int *> (MKL_calloc(nbytes, 1, 64));
    // Initialize
    const T *__attribute__((aligned(64))) 
    mfPtr = mf.getMatchedFilterSignalPointer(0);
    if (pImpl->mUseAbsoluteValue)
    {
        #pragma omp simd aligned(mfPtr, det: 64)
        for (int i=0; i<detectionLength; ++i){det[i] = std::abs(mfPtr[i]);}
    }
    else
    {
#ifdef USE_PSTL
        std::copy(std::execution::unseq, mfPtr, mfPtr+detectionLength, det);
#else
        std::copy(mfPtr, mfPtr+detectionLength, det);
#endif
    }
    // Now reduce
    for (int it=1; it<nt; ++it)
    {
        if (pImpl->mUseAbsoluteValue)
        {
            #pragma omp simd aligned(mfPtr, det, id: 64)
            for (int i=0; i<detectionLength; ++i) 
            {
                auto absXC = std::abs(mfPtr[i]);
                if (absXC > det[i])
                {
                    det[i] = absXC;
                    id[i] = it;
                }
            }
        }
        else
        {
            #pragma omp simd aligned(mfPtr, det, id: 64)
            for (int i=0; i<detectionLength; ++i) 
            {
                if (mfPtr[i] > det[i])
                {
                    det[i] = mfPtr[i];
                    id[i] = it;
                }
            }
        }
    } // Loop on templates
    // Compute the peaks which are the detections
    pImpl->mPeakFinder.setSignal(detectionLength, det);
    pImpl->mPeakFinder.apply();
    auto nDetections = pImpl->mPeakFinder.getNumberOfPeaks();
    // Swing and a miss
    if (nDetections < 1)
    {
        MKL_free(det);
        MKL_free(id);
        return;
    }
    // Get the detection indices
    auto peaksIndices = pImpl->mPeakFinder.getPeakIndicesPointer();
    // Figure out which templates I'll need and copy them
    auto uniqueTemplateIDs
        = getUniqueTemplateIDs(nDetections, peaksIndices, detectionLength, id);
    std::vector<MFLib::WaveformTemplate> templates(uniqueTemplateIDs.size());
    std::vector<MFLib::SingleChannel::RelativeMagnitude<T>> magnitudes;
    std::vector<T> templateSignal;
    for (int i=0; i<static_cast<int> (uniqueTemplateIDs.size()); ++i)
    {
        auto it = uniqueTemplateIDs[i];
        templates[i] = mf.getWaveformTemplate(it);
        // Ensure we can pair this template to something
        if (!templates[i].haveIdentifier()){templates[i].setIdentifier(it);}
        magnitudes[i].initialize(templates[i]);
    }
    // Create detections (i.e., make the products)
    auto getWaveforms = pImpl->mParameters.getDetectedWaveform();
    pImpl->mDetections.resize(0);
    pImpl->mDetections.reserve(nDetections); 
    auto signalLength = mf.getSignalLength();
    auto signalPtr = mf.getSignalPointer();
    // Extract the detections
    for (int i=0; i<nDetections; ++i)
    {
        // Get the detection index in the signal and the corresponding template
        int peakIndex = peaksIndices[i];
        int it = id[peakIndex];
        // Get the template ID
        auto jt = getTemplateIndex(it, uniqueTemplateIDs);
        if (jt < 0){continue;} // Should never happen
        // Set the basics
        pImpl->mDetections[i].setTemplateIdentifier(
            templates[jt].getIdentifier());
        pImpl->mDetections[i].setCorrelationCoefficient(det[peakIndex]);
        // Get the detected chunk of signal
        int templateLength = templates[jt].getSignalLength();
        bool computeAmplitude = true;
        if (peakIndex + templateLength > signalLength)
        {
            fprintf(stderr, "Cannot compute amplitude for detection %d\n", i+1);
            templateLength = signalLength - peakIndex;
            computeAmplitude = false;
        }
        const T *detectedSignalPtr = signalPtr + peakIndex;
        if (getWaveforms)
        {
            pImpl->mDetections[i].setDetectedSignal(templateLength,
                                                    detectedSignalPtr);
        }
        // Compute the onset time
        auto dt = 1/templates[jt].getSamplingRate(); // Required
        auto detectionTime = peakIndex*dt;
        pImpl->mDetections[i].setDetectionTime(detectionTime);
        // Compute the interpolated onset time
        auto mfPtr = mf.getMatchedFilterSignalPointer(it);
        auto shift = quadraticRefinement(detectionLength,
                                         mfPtr,
                                         peakIndex,
                                         dt,
                                         pImpl->mUseAbsoluteValue);
        auto intDetTime = detectionTime + shift;
        pImpl->mDetections[i].setDetectionTime(detectionTime);
        pImpl->mDetections[i].setInterpolatedDetectionTime(intDetTime);
        // Try getting the phase onset time
        if (templates[jt].havePhaseOnsetTime())
        {
            auto pickTime = detectionTime
                          + templates[jt].getPhaseOnsetTime();
            auto intPickTime = intDetTime
                             + templates[jt].getPhaseOnsetTime();
            pImpl->mDetections[i].setPhaseOnsetTime(pickTime);
            pImpl->mDetections[i].setInterpolatedPhaseOnsetTime(intPickTime);
        }
        // Compute the relative amplitudes / magnitudes
        if (computeAmplitude)
        {
            magnitudes[jt].setDetectedSignal(templateLength, detectedSignalPtr);
            auto magType = MFLib::RelativeMagnitudeType::GIBBONS_RINGDAL_2006;
            T alpha = magnitudes[jt].computeAmplitudeScalingFactor(magType);
            pImpl->mDetections[it].setAmplitudeScalingFactor(alpha, magType);

            magType = MFLib::RelativeMagnitudeType::SCHAFF_RICHARDS_2014;
            alpha = magnitudes[jt].computeAmplitudeScalingFactor(magType);
            pImpl->mDetections[it].setAmplitudeScalingFactor(alpha, magType);
        }
    }
}

/// Determines if class is initialized
template<class T>
bool Detector<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Template instantation
template class MFLib::SingleChannel::Detector<double>;
template class MFLib::SingleChannel::Detector<float>;
