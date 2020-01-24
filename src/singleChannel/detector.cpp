#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <mkl.h>
#include "mflib/singleChannel/detector.hpp"
#include "mflib/singleChannel/detection.hpp"
#include "mflib/singleChannel/matchedFilter.hpp"
#include "private/peakFinder.hpp"

using namespace MFLib::SingleChannel;

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
    pImpl->mDetections.resize(nDetections);
    // Extract the detections

    // Cleanup
    MKL_free(det);
    MKL_free(id);
}

/// Template instantation
template class MFLib::SingleChannel::Detector<double>;
template class MFLib::SingleChannel::Detector<float>;
