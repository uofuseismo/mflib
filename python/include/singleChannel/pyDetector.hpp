#ifndef PYMFLIB_SINGLECHANNEL_DETECTOR_HPP
#define PYMFLIB_SINGLECHANNEL_DETECTOR_HPP
#include <memory>
#include <pybind11/pybind11.h>
#include "singleChannel/pyDetection.hpp"

namespace MFLib
{
namespace SingleChannel
{
template<class T> class Detector;
}
}

namespace PBMFLib
{
namespace SingleChannel
{

class DetectorParameters;
template<class T> class Detection;
template<class T> class MatchedFilter;

template<class T>
class Detector
{
public:
    Detector();
    ~Detector();
    void initialize(const DetectorParameters &parms);
    void detect(const MatchedFilter<T> &mf);
    int getNumberOfDetections() const noexcept;
    Detection<T> getDetection(int it) const;
private:
    std::unique_ptr<MFLib::SingleChannel::Detector<T>> mDetector; 
};
void initializeDetector(pybind11::module &m);
}
}
#endif
