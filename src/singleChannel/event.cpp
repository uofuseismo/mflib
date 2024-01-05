#include <cstdio>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include "mflib/singleChannel/event.hpp"
#include "private/private.hpp"

using namespace MFLib::SingleChannel;

class Event::EventImpl
{
public:
    double mMagnitude = 0;
    double mLatitude = 0;
    double mLongitude = 0;
    double mDepth = 0;
    double mRMS = 0;
    double mHorizontalError = 0;
    double mVerticalError = 0;
    uint64_t mIdentifier = 0;

    bool mHaveMagnitude = false;
    bool mHaveLatitude = false;
    bool mHaveLongitude = false;
    bool mHaveDepth = false;
    bool mHaveIdentifier = false;
};

/// C'tor
Event::Event() :
    pImpl(std::make_unique<EventImpl> ())
{
}

/// Copy c'tor
Event::Event(const Event &event)
{
    *this = event;
}

/// Move c'tor
Event::Event(Event &&event) noexcept
{
    *this = std::move(event);
}

/// Copy assignment
Event& Event::operator=(const Event &event)
{
    if (&event == this){return *this;}
    pImpl = std::make_unique<EventImpl> (*event.pImpl);
    return *this;
}

/// Move assignment operator
Event& Event::operator=(Event &&event) noexcept
{
    if (&event == this){return *this;}
    pImpl = std::move(event.pImpl);
    return *this;
}

/// Destructor
Event::~Event() = default;

/// Clear the class
void Event::clear() noexcept
{
    pImpl->mMagnitude = 0;
    pImpl->mLatitude = 0;
    pImpl->mLongitude = 0;
    pImpl->mDepth = 0;
    pImpl->mRMS = 0;
    pImpl->mHorizontalError = 0;
    pImpl->mVerticalError = 0;
    pImpl->mIdentifier = 0;
    pImpl->mHaveMagnitude = false;
    pImpl->mHaveLatitude = false;
    pImpl->mHaveLongitude = false;
    pImpl->mHaveDepth = false;
    pImpl->mHaveIdentifier = false;
}

/// Event identifier
void Event::setIdentifier(const uint64_t evid) noexcept
{
    pImpl->mIdentifier = evid;
    pImpl->mHaveIdentifier = true;
}

uint64_t Event::getIdentifier() const
{
    if (!haveIdentifier())
    {
        throw std::runtime_error("Event identifier not yet set\n");
    }
    return pImpl->mIdentifier;
}

bool Event::haveIdentifier() const noexcept
{
    return pImpl->mHaveIdentifier;
}

/// Magnitude
void Event::setMagnitude(const double magnitude) noexcept
{
    pImpl->mMagnitude = magnitude;
    pImpl->mHaveMagnitude = false;
}

double Event::getMagnitude() const
{
    if (!haveMagnitude())
    {
        throw std::runtime_error("Event magnitude not yet set\n");
    }
    return pImpl->mMagnitude;
}

bool Event::haveMagnitude() const noexcept
{
    return pImpl->mHaveMagnitude;
}

/// Latitude 
void Event::setLatitude(const double latitude)
{
    if (latitude <-90 || latitude >+90)
    {
        std::string errmsg = "latitude = " + std::to_string(latitude)
                           + " must be in range [-90,90]\n"; 
        throw std::invalid_argument(errmsg);
    }
    pImpl->mLatitude = latitude;
    pImpl->mHaveLatitude = true;
}

double Event::getLatitude() const
{
    if (!haveLatitude()){throw std::runtime_error("Latitude never set\n");}
    return pImpl->mLatitude;
}

bool Event::haveLatitude() const noexcept
{
    return pImpl->mHaveLatitude;
}

/// Longitude
void Event::setLongitude(const double longitude)
{
    if (longitude <-540 || longitude >+ 540)
    {
        std::string errmsg = "longitude = " + std::to_string(longitude)
                           + " must be in range (-540,540)\n";
        throw std::invalid_argument(errmsg);
    }   
    pImpl->mLongitude = longitude;
    // Force this into bounds [0,360)
    for (int k=0; k<3; ++k)
    {
        if (pImpl->mLongitude >= 0){break;}
        pImpl->mLongitude = pImpl->mLongitude + 360;
    }
    for (int k=0; k<3; ++k)
    {
        if (pImpl->mLongitude < 360){break;}
        pImpl->mLongitude = pImpl->mLongitude - 360;
    }
    pImpl->mHaveLongitude = true;
}

double Event::getLongitude() const
{
    if (!haveLongitude()){throw std::runtime_error("Longitude never set\n");}
    return pImpl->mLongitude;
}

bool Event::haveLongitude() const noexcept
{
    return pImpl->mHaveLongitude;
}

/// Depth
void Event::setDepth(const double depth)
{
    pImpl->mDepth = depth;
    pImpl->mHaveDepth = true;
}

double Event::getDepth() const
{
    if (!haveDepth()){throw std::runtime_error("Depth never set\n");}
    return pImpl->mDepth;
}

bool Event::haveDepth() const noexcept
{
    return pImpl->mHaveDepth;
}

/// Location quality
void Event::setRootMeanSquaredError(const double rms)
{
    if (rms < 0)
    {
        throw std::invalid_argument("RMS = " + std::to_string(rms)
                                  + " must be non-negative\n");
    }
    pImpl->mRMS = rms;
}

double Event::getRootMeanSquaredError() const noexcept
{
    return pImpl->mRMS;
}

void Event::setVerticalError(const double erz)
{
    if (erz < 0)
    {
        throw std::invalid_argument("Vertical error = " + std::to_string(erz)
                                  + " must be non-negative\n");
    }
    pImpl->mVerticalError = erz;
}

double Event::getVerticalError() const noexcept
{
    return pImpl->mVerticalError;
}

void Event::setHorizontalError(const double erh)
{
    if (erh < 0)
    {
        throw std::invalid_argument("Horizontal error = " + std::to_string(erh)
                                  + " must be non-negative\n");
    }
    pImpl->mHorizontalError = erh;
}

double Event::getHorizontalError() const noexcept
{
    return pImpl->mHorizontalError;
}

