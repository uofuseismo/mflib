#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include "mflib/singleChannel/detection.hpp"

using namespace MFLib::SingleChannel;

template<class T>
class Detection<T>::DetectionImpl
{
public:
    std::vector<T> mDetectedSignal;
    double mDetectionTime = 0;
    double mInterpolatedDetectionTime = 0;
    double mPhaseOnsetTime = 0;
    double mInterpolatedPhaseOnsetTime = 0;
    uint64_t mTemplateID = 0;
    bool mHaveDetectedSignal = false;
    bool mHaveDetectionTime = false;
    bool mHaveInterpolatedDetectionTime = false;
    bool mHavePhaseOnsetTime = false;
    bool mHaveInterpolatedPhaseOnsetTime = false;
    bool mHaveTemplateID = false;
    bool mHaveDetection = false;
};

/// Constructor
template<class T>
Detection<T>::Detection() :
    pImpl(std::make_unique<DetectionImpl> ())
{
}

/// Copy c'tor
template<class T>
Detection<T>::Detection(const Detection &detection)
{
    *this = detection;
}

/// Move c'tor
template<class T>
Detection<T>::Detection(Detection &&detection) noexcept
{
    *this = std::move(detection);
}
 
/// Copy assignment operator
template<class T>
Detection<T>& Detection<T>::operator=(const Detection &detection)
{
    if (&detection == this){return *this;}
    pImpl = std::make_unique<DetectionImpl> (*detection.pImpl);
    return *this;
}

/// Move assignment operator
template<class T>
Detection<T>& Detection<T>::operator=(Detection &&detection) noexcept
{
    if (&detection == this){return *this;}
    pImpl = std::move(detection.pImpl);
    return *this;
}

/// Destructor
template<class T>
Detection<T>::~Detection() = default;

/// Clears the class
template<class T>
void Detection<T>::clear() noexcept
{
    pImpl->mDetectedSignal.clear();
    pImpl->mDetectionTime = 0;
    pImpl->mInterpolatedDetectionTime = 0;
    pImpl->mPhaseOnsetTime = 0;
    pImpl->mInterpolatedPhaseOnsetTime = 0;
    pImpl->mTemplateID = 0;
    pImpl->mHaveDetectionTime = false;
    pImpl->mHaveInterpolatedDetectionTime = false;
    pImpl->mHavePhaseOnsetTime = false;
    pImpl->mHaveInterpolatedPhaseOnsetTime = false;
    pImpl->mHaveTemplateID = false;
    pImpl->mHaveDetectedSignal = false;
    pImpl->mHaveDetection = false;
}

/// The template identifier
template<class T>
void Detection<T>::setTemplateIdentifier(const uint64_t id) noexcept
{
    pImpl->mHaveTemplateID = true;
    pImpl->mTemplateID = id;
}

template<class T>
uint64_t Detection<T>::getTemplateIdentifier() const
{
    if (!haveTemplateIdentifier())
    {
        throw std::runtime_error("Template ID not set\n");
    }
    return pImpl->mTemplateID;
}

template<class T>
bool Detection<T>::haveTemplateIdentifier() const noexcept
{
    return pImpl->mHaveTemplateID;
}

/// The detected signal
template<class T>
void Detection<T>::setDetectedSignal(const int n, const T x[]) 
{
    pImpl->mHaveDetectedSignal = false;
    if (n < 1 || x == nullptr)
    {
        if (n < 1){throw std::invalid_argument("n must be positive\n");}
        throw std::invalid_argument("x is NULL\n");
    }
    pImpl->mDetectedSignal.resize(n);
#ifdef USE_PSTL
    std::copy(std::execution::unseq, x, x+n, pImpl->mDetectedSignal.data());
#else
    std::copy(x, x+n, pImpl->mDetectedSignal.data());
#endif
    pImpl->mHaveDetectedSignal = true;
}

template<class T>
void Detection<T>::getDetectedSignal(const int n, T *xIn[]) const
{
    auto npts = getDetectedSignalLength(); // Throws
    auto x = *xIn;
    if (n < npts || x == nullptr)
    {
        if (n < npts)
        {
            throw std::invalid_argument("Signal length = " + std::to_string(n)
                                      + " must be at least = "
                                      + std::to_string(npts) + "\n");
        }
        throw std::invalid_argument("x is NULL\n");
    }
    auto xPtr = getDetectedSignalPointer(); // Throws
#ifdef USE_PSTL
    std::copy(std::execution::unseq, xPtr, xPtr+n, x);
#else
    std::copy(xPtr, xPtr+n, x);
#endif
}

template<class T>
const T* Detection<T>::getDetectedSignalPointer() const
{
    if (!haveDetectedSignal())
    {
        throw std::runtime_error("Detected signal not yet set\n");
    }
    return pImpl->mDetectedSignal.data(); 
}

template<class T>
int Detection<T>::getDetectedSignalLength() const
{
    if (!haveDetectedSignal())
    {
        throw std::runtime_error("Detected signal not yet set\n");
    }
    return static_cast<int> (pImpl->mDetectedSignal.size());
}

template<class T>
bool Detection<T>::haveDetectedSignal() const noexcept
{
    return pImpl->mHaveDetectedSignal;
}

/// The detection time
template<class T>
void Detection<T>::setDetectionTime(const double time)
{
    pImpl->mHaveDetectionTime = false;
    if (time < 0)
    {
        throw std::invalid_argument("Detection onset time must be positive\n");
    }
    pImpl->mDetectionTime = time;
    pImpl->mHaveDetectionTime = true;
}

template<class T>
double Detection<T>::getDetectionTime() const
{
    if (!haveDetectionTime())
    {
        throw std::invalid_argument("Detection time not yet set\n");
    }
    return pImpl->mDetectionTime;
}

template<class T>
bool Detection<T>::haveDetectionTime() const noexcept
{
    return pImpl->mHaveDetectionTime;
}

/// The interpolated detection time
template<class T>
void Detection<T>::setInterpolatedDetectionTime(const double time)
{
    pImpl->mHaveInterpolatedDetectionTime = false;
    if (time < 0)
    {
        throw std::invalid_argument("Detection onset time must be positive\n");
    }
    pImpl->mInterpolatedDetectionTime = time;
    pImpl->mHaveInterpolatedDetectionTime = true;
}

template<class T>
double Detection<T>::getInterpolatedDetectionTime() const
{
    if (!haveInterpolatedDetectionTime())
    {
        throw std::invalid_argument("Detection time not yet set\n");
    }
    return pImpl->mInterpolatedDetectionTime;
}

template<class T>
bool Detection<T>::haveInterpolatedDetectionTime() const noexcept
{
    return pImpl->mHaveInterpolatedDetectionTime;
}

/// The phase onset time
template<class T>
void Detection<T>::setPhaseOnsetTime(const double time)
{
    pImpl->mHavePhaseOnsetTime = false;
    if (time < 0)
    {
        throw std::invalid_argument("Phase onset time must be positive\n");
    }
    pImpl->mPhaseOnsetTime = time;
    pImpl->mHavePhaseOnsetTime = true;
}

template<class T>
double Detection<T>::getPhaseOnsetTime() const
{
    if (!havePhaseOnsetTime())
    {   
        throw std::invalid_argument("Phase onset time not yet set\n");
    }   
    return pImpl->mPhaseOnsetTime;
}

template<class T>
bool Detection<T>::havePhaseOnsetTime() const noexcept
{
    return pImpl->mHavePhaseOnsetTime;
}

/// The interpolated phase onset time
template<class T>
void Detection<T>::setInterpolatedPhaseOnsetTime(const double time)
{
    pImpl->mHaveInterpolatedPhaseOnsetTime = false;
    if (time < 0)
    {
        throw std::invalid_argument("Phase onset time must be positive\n");
    }
    pImpl->mInterpolatedPhaseOnsetTime = time;
    pImpl->mHaveInterpolatedPhaseOnsetTime = true;
}

template<class T>
double Detection<T>::getInterpolatedPhaseOnsetTime() const
{
    if (!haveInterpolatedPhaseOnsetTime())
    {
        throw std::invalid_argument("Phase onset time not yet set\n");
    }
    return pImpl->mInterpolatedPhaseOnsetTime;
}

template<class T>
bool Detection<T>::haveInterpolatedPhaseOnsetTime() const noexcept
{
    return pImpl->mHaveInterpolatedPhaseOnsetTime;
}
/// Determines if I have a detection
/*
template<class T>
bool Detection<T>::haveDetection() const noexcept
{
    return pImpl->mHaveDetection;
}
*/

template class MFLib::SingleChannel::Detection<double>;
template class MFLib::SingleChannel::Detection<float>;
