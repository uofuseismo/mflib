#include <cstdio>
#include <cstdlib>
#include "mflib/singleChannel/detectionTimeInterpolationParameters.hpp"

using namespace MFLib::SingleChannel;

class DetectionTimeInterpolationParameters::InterpParametersImpl
{
public:
    int mLanczosAlpha = 51; 
    int mLanczosInterpPoints = 100;
    DetectionTimeInterpolationType mType
        = DetectionTimeInterpolationType::QUADRATIC;
};

/// Constructor
DetectionTimeInterpolationParameters::DetectionTimeInterpolationParameters() :
    pImpl(std::make_unique<InterpParametersImpl> ())
{
}

/// Copy c'tor
DetectionTimeInterpolationParameters::DetectionTimeInterpolationParameters(
    const DetectionTimeInterpolationParameters &parameters)
{
    *this = parameters;
}

/// Move c'tor
DetectionTimeInterpolationParameters::DetectionTimeInterpolationParameters(
    DetectionTimeInterpolationParameters &&parameters) noexcept
{
    *this = std::move(parameters);
}

/// Copy assignment
DetectionTimeInterpolationParameters&
DetectionTimeInterpolationParameters::operator=(
    const DetectionTimeInterpolationParameters &parameters)
{
    if (&parameters == this){return *this;}
    pImpl = std::make_unique<InterpParametersImpl> (*parameters.pImpl);
    return *this;
}

/// Move assignment
DetectionTimeInterpolationParameters& 
DetectionTimeInterpolationParameters::operator=(
    DetectionTimeInterpolationParameters &&parameters) noexcept
{
    if (&parameters == this){return *this;}
     pImpl = std::move(parameters.pImpl);
     return *this;
}

/// Destructor
DetectionTimeInterpolationParameters::~DetectionTimeInterpolationParameters()
    = default;

/// Resets the class
void DetectionTimeInterpolationParameters::clear() noexcept
{
    pImpl->mLanczosAlpha = 51;
    pImpl->mLanczosInterpPoints = 100;
    pImpl->mType = DetectionTimeInterpolationType::QUADRATIC;
}

/// Lanczos 
void DetectionTimeInterpolationParameters::enableLanczos() noexcept
{
    pImpl->mType = DetectionTimeInterpolationType::LANCZOS;
}

void DetectionTimeInterpolationParameters::setLanczosNumberOfInterpolationPoints(
    const int npts)
{
    if (npts < 1)
    {
        throw std::invalid_argument(
            "Number of interpolation points must be positive\n"); 
    }
    pImpl->mLanczosInterpPoints  = npts;
}

int DetectionTimeInterpolationParameters::getLanczosNumberOfInterpolationPoints() const noexcept
{
    return pImpl->mLanczosInterpPoints;
}

void DetectionTimeInterpolationParameters::setLanczosAlpha(const int alpha)
{
    if (alpha < 1)
    {
        throw std::invalid_argument("alpha must be positive\n");
    }
    pImpl->mLanczosAlpha = alpha;
}

int DetectionTimeInterpolationParameters::getLanczosAlpha() const noexcept
{
    return pImpl->mLanczosAlpha;
}

/// Quadratic
void DetectionTimeInterpolationParameters::enableQuadratic() noexcept
{
    pImpl->mType = DetectionTimeInterpolationType::QUADRATIC;
}

void DetectionTimeInterpolationParameters::disable() noexcept
{
    pImpl->mType = DetectionTimeInterpolationType::NONE;
}

/// Gets the interpolation type
MFLib::DetectionTimeInterpolationType  
DetectionTimeInterpolationParameters::getType() const noexcept
{
    return pImpl->mType;
}
