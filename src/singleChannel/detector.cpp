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
#include "mflib/networkStationPhase.hpp"
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

#pragma omp declare simd
template<class T> T sinc(const T x)
{
    const T tol = std::numeric_limits<T>::epsilon()*100;
    T result = 1;
    if (std::abs(x) > tol)
    {
        T pix = M_PI*x;
        result = std::sin(pix)/pix;
    }
    return result;
} 

#pragma omp declare simd
template<class T> T lanczos(const T x, const int a)
{
    T result = 0;
    if (static_cast<int> (std::abs(x)) <= a){result = sinc(x)*sinc(x/a);}
    return result;
}

template<class T>
T lanczosRefinement(const int n,
                    const T xc[],
                    const int optIndex,
                    const double dt, 
                    const bool lAbs,
                    const int alpha = 19,
                    const int nrefine = 51)
{
    T dtWork = 1/static_cast<double> (nrefine);
    T optShift = 0;
    T yIntMax =-std::numeric_limits<T>::max();
    //for (int idx=1; idx<=2*nrefine-1; idx++)
    for (int idx=1; idx<nrefine; idx++)
    {
        //T shift =-1 + dtWork*idx;
        T shift =-0.5 + dtWork*idx;
        T x = optIndex + shift;
        T yInt = 0;
        #pragma omp simd reduction(+:yInt)
        for (int a=-alpha; a<=alpha; ++a)
        {
            int i = static_cast<int> (std::floor(x)) - a + 1;
            auto xi = x - static_cast<T> (i);
            T xcv = 0;
            if (i >= 0 && i < n){xcv = xc[i];}
            yInt = yInt + xcv*lanczos(xi, alpha);
        }
        if (lAbs){yInt = std::abs(yInt);}
        if (yInt > yIntMax)
        { 
            optShift = dt*shift;
            yIntMax = yInt;
        }
//printf("%d, %lf, %lf, (%lf,%lf,%lf), %lf\n", idx, x, shift, xc[optIndex-1], xc[optIndex], xc[optIndex+1], yInt);
    }
//printf("%lf\n", optShift);
    return optShift;
}

std::vector<int> getUniqueTemplateIDs(const int nDetections,
                                      const int detectionIndex[],
                                      const int ids[])
{
    // Extract the IDs at each detection index
    std::vector<int> uniqueIDs(nDetections);
    #pragma omp simd
    for (int i=0; i<nDetections; ++i)
    {
        int index = detectionIndex[i];
        uniqueIDs[i] = ids[index];
    }
    // Sort (in ascending order) for stl's std::unique
#ifdef USE_PSTL
    std::sort(std::execution::unseq, uniqueIDs.begin(), uniqueIDs.end());
#else
    std::sort(uniqueIDs.begin(), uniqueIDs.end());
#endif
    // Remove duplicates and resize vector
    auto last = std::unique(uniqueIDs.begin(), uniqueIDs.end());
    uniqueIDs.erase(last, uniqueIDs.end()); 
    return uniqueIDs;
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

template<class T>
void cleanReducedTemplate(const int detectionLength, T *det)
{
    const T one = 1;
    #pragma omp simd aligned(det: 64)
    for (int i=0; i<detectionLength; ++i)
    {
        det[i] = std::max(-one, std::min(one, det[i]));
    }
}

template<class T>
void reduceAllCorrelations(
    const bool mUseAbsoluteValue,
    const MFLib::SingleChannel::MatchedFilter<T> &mf,
    T *det, int *id)
{
    auto nt = mf.getNumberOfTemplates();
    auto detectionLength = mf.getFilteredSignalLength();
    auto mfPtr = mf.getMatchedFilterSignalPointer(0);
    if (mUseAbsoluteValue)
    {
        #pragma omp simd aligned(mfPtr, det: 64)
        for (int i=0; i<detectionLength; ++i){det[i] = std::abs(mfPtr[i]);}
    }
    else
    {
#ifdef USE_PSTL
        std::copy(std::execution::unseq, mfPtr, mfPtr+detectionLength, det);
        std::fill(std::execution::unseq, id, id+detectionLength, 0);
#else
        std::copy(mfPtr, mfPtr+detectionLength, det);
        std::fill(id, id+detectionLength, 0);
#endif
    }
    // Now reduce
    for (int it=1; it<nt; ++it)
    {
        mfPtr = mf.getMatchedFilterSignalPointer(it);
        if (mUseAbsoluteValue)
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
    // And clean
    cleanReducedTemplate(detectionLength, det);
}

template<class T>
void setCorrelations(
    const int it,
    const bool mUseAbsoluteValue,
    const MFLib::SingleChannel::MatchedFilter<T> &mf,
    T *det, int *id)
{
    auto detectionLength = mf.getFilteredSignalLength();
    auto mfPtr = mf.getMatchedFilterSignalPointer(it);
    if (mUseAbsoluteValue)
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
    std::fill(id, id+detectionLength, it);
    // And clean
    cleanReducedTemplate(detectionLength, det);
}

}

///--------------------------------------------------------------------------///
///                               Implementation                             ///
///--------------------------------------------------------------------------/// 
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
        mUseAbsoluteValue = false;
        mHaveMatchedFilteredSignals = false;
        mInitialized = false;
    }

    std::vector<MFLib::SingleChannel::Detection<T>> mDetections;
    MFLib::PeakFinder<T> mPeakFinder;
    MFLib::SingleChannel::DetectorParameters mParameters;
 
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

/// Operator
template<class T>
const MFLib::SingleChannel::Detection<T>&
Detector<T>::operator[](const size_t i) const
{
    auto nDetections = static_cast<size_t> (getNumberOfDetections());
    if (i >= nDetections)
    {
        std::string errmsg;
        if (i == 0)
        {
            errmsg = "No detections\n";
        }
        else
        {
            errmsg = "i = " + std::to_string(i) + " cannot exceed "
                  + std::to_string(nDetections) + "\n";
        }
        throw std::invalid_argument(errmsg);
    }
    return pImpl->mDetections[i];
}

template<class T>
MFLib::SingleChannel::Detection<T>
Detector<T>::getDetection(const int id) const
{
    auto nDetections = getNumberOfDetections();
    if (id >= nDetections)
    {
        std::string errmsg;
        if (id == 0)
        {
            errmsg = "No detections\n";
        }
        else
        {
            errmsg = "i = " + std::to_string(id) + " cannot exceed "
                  + std::to_string(nDetections) + "\n";
        }
        throw std::invalid_argument(errmsg);
    }
    return pImpl->mDetections[id];
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
    // Class is ready to go
    pImpl->mInitialized = true;
}

/// Gets the number of detections
template<class T>
int Detector<T>::getNumberOfDetections() const noexcept
{
    return static_cast<int> (pImpl->mDetections.size());
}


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
    auto detectionLength = mf.getFilteredSignalLength();
    if (detectionLength < 1)
    {
        throw std::invalid_argument("Matched filtered signal length is 0\n");
    }
    // Determine if I am making detections for each correlogram or a 
    // a reduced corrrelogram
    bool channelBased = false;
    if (pImpl->mParameters.getMatchedFilteredSignalDetectorPolicy() ==
        MFLib::MatchedFilteredSignalDetectorPolicy::SINGLE)
    {
        channelBased = true;
    }
    // Set workspace
    auto nbytes = sizeof(T)*static_cast<size_t> (detectionLength);
    auto det = static_cast<T *> (MKL_calloc(nbytes, 1, 64));
    nbytes = sizeof(int)*static_cast<size_t> (detectionLength);
    auto id = static_cast<int *> (MKL_calloc(nbytes, 1, 64));
    int nDetections = 0;
    std::vector<MFLib::PeakFinder<T>> peakFinders;
    if (!channelBased)
    {
        peakFinders.resize(1);
        peakFinders[0] = pImpl->mPeakFinder;
        // Initialize
        reduceAllCorrelations(pImpl->mUseAbsoluteValue, mf, det, id);
        // Compute the peaks which are the detections
        peakFinders[0].setSignal(detectionLength, det);
        peakFinders[0].apply();
        nDetections = peakFinders[0].getNumberOfPeaks();
    }
    else
    {
        peakFinders.resize(nt);
        for (int it=0; it<nt; ++it)
        {
            peakFinders[it] = pImpl->mPeakFinder;
            setCorrelations(it, pImpl->mUseAbsoluteValue, mf, det, id);
            // Compute the peaks which are the detections
            peakFinders[it].setSignal(detectionLength, det);
            peakFinders[it].apply();
            nDetections = nDetections + peakFinders[it].getNumberOfPeaks();
        }
    }
    // Swing and a miss
    if (nDetections < 1)
    {
        MKL_free(det);
        MKL_free(id);
        return;
    }
    // Loop through the peaks
    for (int ip=0; ip<static_cast<int> (peakFinders.size()); ++ip)
    {
        // Get the detection indices
        auto peaksIndices = peakFinders[ip].getPeakIndicesPointer();
        // Figure out which templates I'll need and copy them
        auto uniqueTemplateIDs
            = getUniqueTemplateIDs(nDetections, peaksIndices, id);
        std::vector<MFLib::WaveformTemplate> templates(uniqueTemplateIDs.size());
        std::vector<T> templateSignal(uniqueTemplateIDs.size());
        for (int i=0; i<static_cast<int> (uniqueTemplateIDs.size()); ++i)
        {
            auto it = uniqueTemplateIDs[i];
            templates[i] = mf.getWaveformTemplate(it);
            // Ensure we can pair this template to something
            if (!templates[i].haveIdentifier())
            {
                NetworkStationPhase nsp;
                std::pair<NetworkStationPhase, uint64_t> idWork(nsp, it);
                templates[i].setIdentifier(idWork);
            }
        }
        // Initialize the relative magnitudes?
        std::vector<MFLib::SingleChannel::RelativeMagnitude<T>> magnitudes;
        bool lWantMagnitudes = pImpl->mParameters.wantAmplitudeScalingFactor();
        if (lWantMagnitudes)
        {
            magnitudes.resize(templates.size());
            for (int i=0; i<static_cast<int> (templates.size()); ++i)
            {
                magnitudes[i].initialize(templates[i]);
            }
        }
        // Create detections (i.e., make the products)
        auto lGetWaveforms = pImpl->mParameters.wantDetectedWaveform();
        pImpl->mDetections.resize(nDetections); 
        auto signalLength = mf.getSignalLength();
        auto signalPtr = mf.getSignalPointer();
        // Extract the detections
        for (int i=0; i<nDetections; ++i)
        {
            // For simplicity we'll create the detection then copy it later on
            MFLib::SingleChannel::Detection<T> detection; 
            // Get the detection index in the signal and the corresponding
            // template
            int peakIndex = peaksIndices[i];
            int it = id[peakIndex];
            // Get the template ID
            auto jt = getTemplateIndex(it, uniqueTemplateIDs);
            if (jt < 0)
            {
                fprintf(stderr, "Algorithmic failure - negative ID\n");
                continue;
            }
            // Set the basics
            detection.setTemplateIdentifier(templates[jt].getIdentifier());
            detection.setCorrelationCoefficient(det[peakIndex]);
            // Get the detected chunk of signal
            int templateLength = templates[jt].getSignalLength();
            auto lComputeAmplitude = lWantMagnitudes;
            if (peakIndex + templateLength > signalLength)
            {
                fprintf(stderr, "Cannot compute amplitude for detection %d\n",
                        i+1);
                templateLength = signalLength - peakIndex;
                lComputeAmplitude = false;
            }
            const T *detectedSignalPtr = signalPtr + peakIndex;
            if (lGetWaveforms)
            {
                detection.setDetectedSignal(templateLength, detectedSignalPtr);
            }
            // Compute the onset time
            auto dt = 1/templates[jt].getSamplingRate(); // Required
            auto detectionTime = peakIndex*dt;
            detection.setDetectionTime(detectionTime);
            // Compute the interpolated onset time
            auto mfPtr = mf.getMatchedFilterSignalPointer(it);
            /*
            auto shift = quadraticRefinement(detectionLength,
                                             mfPtr,
                                             peakIndex,
                                             dt,
                                             pImpl->mUseAbsoluteValue);
            */
            // printf("%lf, %d, %lf\n", detectionTime, it, shift);
            // TODO 51 and 100 should be parameters
            auto shift = lanczosRefinement(detectionLength, mfPtr,
                                           peakIndex, dt,
                                           pImpl->mUseAbsoluteValue, 51, 100);
            auto intDetTime = detectionTime + shift;
            detection.setDetectionTime(detectionTime);
            detection.setInterpolatedDetectionTime(intDetTime);
            // Try getting the phase onset time
            if (templates[jt].havePhaseOnsetTime())
            {
                auto pickTime = detectionTime
                              + templates[jt].getPhaseOnsetTime();
                auto intPickTime = intDetTime
                                 + templates[jt].getPhaseOnsetTime();
                detection.setPhaseOnsetTime(pickTime);
                detection.setInterpolatedPhaseOnsetTime(intPickTime);
            }
            // Try getting the polarity by investigating the sign of the 
            // detection.  For example, if the Pearson correlation is +0.8,
            // then this is an exact match and the enums will retain their
            // original value.  If the Pearson correlation is, say, -.91, then
            // the detection more than 90 degrees out of phase, so the enums
            // flip sign.  If the original enum was UNKNOWN (0) then multiplying 
            // by the sign of hte detection will leave the enum as UNKNOWN 
            // (e.g., 1*0=0)
            int detectionSign = 1;
            if (mfPtr[peakIndex] < 0){detectionSign =-1;}
            if (mfPtr[peakIndex] == 0){detectionSign = 0;}
            auto polarity = templates[jt].getPolarity();
            auto iPolarity = static_cast<int> (polarity); // Can be: {-1, 0, 1}
            polarity = static_cast<MFLib::Polarity> (detectionSign*iPolarity);
            detection.setPolarity(polarity);
            // Compute the relative amplitudes / magnitudes
            if (lComputeAmplitude)
            {
                MFLib::RelativeMagnitudeType magType;
                magnitudes[jt].setDetectedSignal(templateLength,
                                                 detectedSignalPtr);
                magType = MFLib::RelativeMagnitudeType::GIBBONS_RINGDAL_2006;
                T alpha = magnitudes[jt].computeAmplitudeScalingFactor(magType);
                detection.setAmplitudeScalingFactor(alpha, magType);

                magType = MFLib::RelativeMagnitudeType::SCHAFF_RICHARDS_2014;
                alpha = magnitudes[jt].computeAmplitudeScalingFactor(magType);
                detection.setAmplitudeScalingFactor(alpha, magType);
            }
            // Save this detection
            pImpl->mDetections[i] = detection;
        }
    }
    // Clean up
    MKL_free(det);
    MKL_free(id);
    // Sort the events by onset time if channel-by-channel.  Otherwise, these
    // are already sorted in ascending order.
    if (channelBased)
    {
        std::sort(pImpl->mDetections.begin(), pImpl->mDetections.end(),
                  [](const MFLib::SingleChannel::Detection<T> &a,
                     const MFLib::SingleChannel::Detection<T> &b)
                  {
                     return a.getDetectionTime() < b.getDetectionTime();
                  });
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
