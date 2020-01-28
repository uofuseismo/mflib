#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <numeric>
#include <mkl.h>
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
        if (mMaxMF){MKL_free(mMaxMF);}
        if (mSignal){MKL_free(mSignal);}
        if (mMaxTemplate){MKL_free(mMaxTemplate);}
        mMaxMF = nullptr;
        mSignal = nullptr;
        mMaxTemplate = nullptr; 
        mDetections.clear();
        mMagnitudes.clear();
        mTemplates.clear();
        mPeakFinder.clear(); 
        mParameters.clear();
        mDetectionLength = 0;
        mUseAbsoluteValue = false;
        mHaveMatchedFilteredSignals = false;
        mInitialized = false;
    }

    /// Holds the max of the matched filtered signals.
    /// This is an array whose dimension is [mDetectionLength].
    T *__attribute__((aligned(64))) mMaxMF = nullptr;
    /// This holds the signal which was matched filtered.  
    /// This is an array whose dimension is [mSignalLength].
    T *__attribute__((aligned(64))) mSignal = nullptr;
    /// Holds the template index at each sample corresponding to mMaxMF.
    /// This is an array whose dimension is [mDetectionLength].
    T *__attribute__((aligned(64))) mMaxTemplate = nullptr;
    std::vector<MFLib::SingleChannel::Detection<T>> mDetections;
    std::vector<MFLib::SingleChannel::RelativeMagnitude<T>> mMagnitudes;
    std::vector<MFLib::WaveformTemplate> mTemplates;
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
        pImpl->mMaxTemplate = static_cast<T *> (MKL_calloc(nbytes, 1, 64));
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

/// Compute the detections
template<class T>
void Detector<T>::detect()
{
    pImpl->mDetections.clear();
    if (!isInitialized())
    {
        throw std::invalid_argument("Detector class not initialized\n");
    }
    auto nt = static_cast<int> (pImpl->mTemplates.size());
    auto detectionLength = pImpl->mDetectionLength;
    const T *det = pImpl->mMaxMF;
    const T *id = pImpl->mMaxTemplate;
    // Compute the peaks which are the detections
    pImpl->mPeakFinder.setSignal(detectionLength, det);
    pImpl->mPeakFinder.apply();
    auto nDetections = pImpl->mPeakFinder.getNumberOfPeaks();
    if (nDetections < 1){return;} // Swing and a miss
    // Get the detection indices
    auto peaksIndices = pImpl->mPeakFinder.getPeakIndicesPointer();
    // Create detections
    pImpl->mDetections.resize(0);
    pImpl->mDetections.reserve(nDetections); 

    auto signalPtr = pImpl->mSignal;
    auto signalLength = pImpl->mSignalLength;
    for (int i=0; i<nDetections; ++i)
    {
        // Get the detection index in the signal and the corresponding template
        int peakIndex = peaksIndices[i];
        int it = id[peakIndex]; 
        // Set the template ID

        // Get the detected chunk of signal 
        int templateLength = pImpl->mTemplates[it].getSignalLength();
        bool computeAmplitude = true;
        if (peakIndex + templateLength > signalLength)
        {
            fprintf(stderr, "Cannot compute amplitude for detection %d\n", i+1);
            templateLength = signalLength - peakIndex;
            computeAmplitude = false;
        }
        const T *detectedSignalPtr = signalPtr + peakIndex;
        pImpl->mDetections[i].setDetectedSignal(templateLength,
                                                detectedSignalPtr);
        // Compute the onset time
        auto dt = 1/pImpl->mTemplates[it].getSamplingRate(); // Required
        auto detectionTime = peakIndex*dt;
        pImpl->mDetections[i].setDetectionTime(detectionTime);
        // Compute the interpolated onset time
/*
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
        if (pImpl->mTemplates[it].havePhaseOnsetTime())
        {
            auto pickTime = detectionTime
                          + pImpl->mTemplates[it].getPhaseOnsetTime();
            auto intPickTime = intDetTime
                             + pImpl->mTemplates[it].getPhaseOnsetTime();
            pImpl->mDetections[i].setPhaseOnsetTime(pickTime);
            pImpl->mDetections[i].setInterpolatedPhaseOnsetTime(intPickTime);
        }
        // Compute the relative amplitudes / magnitudes
        if (computeAmplitude)
        {
            pImpl->mMagnitudes[it].setDetectedSignal(templateLength,
                                                     detectedSignalPtr);
            auto magType = MFLib::RelativeMagnitudeType::GIBBONS_RINGDAL_2006;
            T alpha
                = pImpl->mMagnitudes[it].computeAmplitudeScalingFactor(magType);
            pImpl->mDetections[it].setAmplitudeScalingFactor(alpha, magType);

            magType = MFLib::RelativeMagnitudeType::SCHAFF_RICHARDS_2014;
            alpha
                = pImpl->mMagnitudes[it].computeAmplitudeScalingFactor(magType);
            pImpl->mDetections[it].setAmplitudeScalingFactor(alpha, magType);
        }
*/
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
