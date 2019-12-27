#include <string>
#include "mflib/version.hpp"

using namespace MFLib;

int Version::getMajor() noexcept
{
    return MFLIB_MAJOR;
}

int Version::getMinor() noexcept
{
    return MFLIB_MINOR;
}

int Version::getPatch() noexcept
{
    return MFLIB_PATCH;
}

bool Version::isAtLeast(const int major, const int minor,
                        const int patch) noexcept
{
    if (MFLIB_MAJOR < major){return false;}
    if (MFLIB_MAJOR > major){return true;}
    if (MFLIB_MINOR < minor){return false;}
    if (MFLIB_MINOR > minor){return true;}
    if (MFLIB_PATCH < patch){return false;}
    return true;
}

std::string Version::getVersion() noexcept
{
    std::string version(std::to_string(getMajor()) + "."
                      + std::to_string(getMinor()) + "."
                      + std::to_string(getPatch()));
    return version;
}
