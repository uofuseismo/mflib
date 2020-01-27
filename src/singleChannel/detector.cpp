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
        if (mMaxTemplate){MKL_free(mMaxTemplate);}
        mMaxMF = nullptr;
        mMaxTemplate = nullptr; 
        mDetections.clear();
        mMagnitudes.clear();
        mTemplates.clear();
        mPeakFinder.clear(); 
        mParameters.clear();
        mMaxSizeXC = 0;
        mUseAbsoluteValue = false;
        mInitialized = false;
    }

    T *__attribute__((aligned(64))) mMaxMF = nullptr;
    T *__attribute__((aligned(64))) mMaxTemplate = nullptr;
    std::vector<MFLib::SingleChannel::Detection<T>> mDetections;
    std::vector<MFLib::SingleChannel::RelativeMagnitude<T>> mMagnitudes;
    std::vector<MFLib::WaveformTemplate> mTemplates;
    MFLib::PeakFinder<T> mPeakFinder;
    MFLib::SingleChannel::DetectorParameters mParameters;
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
template<class T>
void Detector<T>::initialize(
    const MFLib::SingleChannel::DetectorParameters &parameters,
    const MFLib::SingleChannel::MatchedFilter<T> &mf)
{
    clear();
    if (!mf.isInitialized())
    {
        throw std::invalid_argument("Matched filter class not initialized\n");
    }
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
    // Class is ready to go
    pImpl->mInitialized = true;
}


/// Compute the detections
template<class T>
void Detector<T>::detect(const MFLib::SingleChannel::MatchedFilter<T> &mf)
{
    pImpl->mDetections.clear();
    if (isInitialized())
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
    // Compute the detections
    auto detectionLength = mf.getFilteredSignalLength();
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

    auto signalPtr = mf.getSignalPointer();
    auto signalLength = mf.getSignalLength();
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
    }
    // Cleanup
    MKL_free(det);
    MKL_free(id);
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
