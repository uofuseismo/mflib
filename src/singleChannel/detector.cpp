#include <cstdio>
#include <cstdlib>
#include <string>
#include "mflib/singleChannel/detector.hpp"
#include "mflib/singleChannel/matchedFilter.hpp"
#include "private/peakFinder.hpp"

using namespace MFLib::SingleChannel;

namespace
{
void detect(const int n, const double x[])
{
    // Compute the zeros
    
}
}

template<class T>
class Detector<T>::DetectorImpl
{
public:
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

/// Compute the detections
