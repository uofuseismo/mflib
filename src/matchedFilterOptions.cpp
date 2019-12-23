#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include <limits>
#include "mflib/waveformTemplate.hpp"
#include "mflib/matchedFilterOptions.hpp"
#include "mflib/enums.hpp"

using namespace MFLib;

namespace
{

int nextPowerOfTwo(const int n)
{
    if (n < 0)
    {   
        throw std::invalid_argument("n = " + std::to_string(n)
                                  + " must be positive\n"); 
    }   
    // Simple base cases
    if (n == 0){return 1;} 
    if (n == 1){return 1;} 
    // General computation
    double dn = static_cast<double> (n);
    // Get log2 and round up
    int itemp = static_cast<int> (std::ceil(std::log2(dn)));
    // Compute next power of 2 greater than n
    uint64_t n2t = static_cast<uint64_t> (std::pow(2, itemp));
    int n2 = static_cast<int> (n2t);
    // Catch any overflow problems associated with int32_t and uint64_t
    if (n2t != static_cast<uint64_t> (n2))
    {   
        throw std::runtime_error("Overflow error");
    }   
    return n2; 
}

/// @brief Computes sizes.
/// @param[in] nb   The filter length.
/// @param[in] nx   The signal length.
/// @result result.first is the FFT size and result.second is the window length.
std::pair<int, int> computeWindowLength(const int nb, const int nx) 
{
    if (nb > nx)
    {
        throw std::invalid_argument("nb cannot exceed nx\n");
    }
    /// Matlab computed floating operation counts for each transform length
    const std::array<double, 21> fftFlops = {
        18, 59, 138, 303, 660, 1441, 3150, 6875, 14952, 32373, 69762,
        149647, 319644, 680105, 1441974, 3047619, 6422736, 13500637, 28311786,
        59244791, 59244791*2.09
    };
    /// Lengths for fast Fourier transforms (must be powers of 2)
    const std::array<uint64_t, 21> fftLengths = {
        2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048,
        4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288,
        1048576, 2097152
    };
    int validSetStart =-1;
    for (int i=0; i<static_cast<int> (fftLengths.size()); ++i)
    {
        if (fftLengths[i] > static_cast<uint64_t> (nb - 1))
        {
            validSetStart = i;
            break;
        }
    }
    if (validSetStart ==-1)
    {
        // throw std::runtime_error("Could not find appropriate padding length\n");
        auto nfft = nextPowerOfTwo(nx);
        auto L =  nfft - (nb - 1); 
        auto result = std::pair(nfft, L); 
        return result;
    }
    else
    {
        //std::vector<int> Nvec(fftLengths.begin() + validSetStart, fftLengths.end());
        // Minimize the (number of blocks) x (number of flops per fft) 
        //std::vector<int> Lvec(Nvec.size(), 0);
        auto objMin = std::numeric_limits<double>::max();
        int nfft = 0;
        int L = 0;
        for (int i=validSetStart; i<static_cast<int> (fftFlops.size()); ++i)
        {
            auto lTemp = fftLengths[i] - (nb - 1); // This is positive
            auto dL = static_cast<double> (lTemp);
            auto obj = std::ceil(nx/dL)*fftFlops[i];
            if (obj < objMin)
            {
                objMin = obj;
                nfft = fftLengths[i];
                L = lTemp;
            }
        }
        return std::pair(nfft, L); 
    }
}
}


class MatchedFilterOptions::MatchedFilterOptionsImpl
{
public:
    /// Collection of templates
    //std::vector<std::vector<T>> mTemplates;
    std::vector<WaveformTemplate> mTemplates;
    /// Matched filter implementation
    MatchedFilterImplementation mImplementation
        = MatchedFilterImplementation::DIRECT;
    /// Sampling rate of templates
    double mSamplingRate =-1;
    /// Signal size
    int mSignalSize = 0;
    /// FFT length
    int mFFTLength = 0;
    /// The block length.  Given a filter length, nb, and FFT length,
    /// the optimal length is  L = mFFTLength - nb - 1.
    int mBlockLength = 0;
    /// Stack absolute value
    bool mStackAbs = true;
};

/// Constructor
MatchedFilterOptions::MatchedFilterOptions() :
    pImpl(std::make_unique<MatchedFilterOptionsImpl> ())
{
}

/// Copy constructor
MatchedFilterOptions::MatchedFilterOptions(
    const MatchedFilterOptions &options)
{
    *this = options;
}

/// Move constructor
MatchedFilterOptions::MatchedFilterOptions(
    MatchedFilterOptions &&options) noexcept
{
    *this = std::move(options);
}

/// Copy assignment
MatchedFilterOptions& MatchedFilterOptions::operator=(
    const MatchedFilterOptions &options)
{
    if (&options == this){return *this;}
    pImpl = std::make_unique<MatchedFilterOptionsImpl> (*options.pImpl);
    return *this;
}

/// Move assignment
MatchedFilterOptions& MatchedFilterOptions::operator=(
    MatchedFilterOptions &&options) noexcept
{
    if (&options == this){return *this;}
    pImpl = std::move(options.pImpl);
    return *this;
}

/// Destructor
MatchedFilterOptions::~MatchedFilterOptions() = default;

/// Clears the memory and resets the class
void MatchedFilterOptions::clear() noexcept
{
    pImpl->mTemplates.clear();
    pImpl->mImplementation = MatchedFilterImplementation::DIRECT;
    pImpl->mSamplingRate =-1;
    pImpl->mSignalSize = 0;
    pImpl->mFFTLength = 0;
    pImpl->mBlockLength = 0;
    pImpl->mStackAbs = true;
}

/// Get the number of templates
int MatchedFilterOptions::getNumberOfTemplates() const noexcept
{
    return static_cast<int> (pImpl->mTemplates.size());
}

/// Clears all the templates
void MatchedFilterOptions::clearTemplates() noexcept
{
    pImpl->mSamplingRate =-1;
    pImpl->mTemplates.clear();
}

void MatchedFilterOptions::addTemplate(const WaveformTemplate &tplate)
{
    pImpl->mFFTLength = 0;
    pImpl->mBlockLength = 0;
    if (!tplate.haveSignal())
    {
        throw std::invalid_argument("Waveform not yet set on template\n");
    }
    // Verify the sampling rates are consistent otherwise Fourier transform
    // will be meaningless.
    double df = tplate.getSamplingRate(); // Throws
    if (pImpl->mSamplingRate < 0){pImpl->mSamplingRate = df;}
    if (std::abs(df - pImpl->mSamplingRate) > 1.e-4)
    {
        throw std::invalid_argument("Template sampling rate = "
                                 + std::to_string(df) + " must equal "
                                 + std::to_string(pImpl->mSamplingRate) + "\n");
    }
    // Check the sampling rate
    pImpl->mTemplates.push_back(tplate); 
}

/// Adds a template
/*
void MatchedFilterOptions::addTemplate(
    const int nSamples, const T templateSignal[])
{
    pImpl->mFFTLength = 0;
    pImpl->mBlockLength = 0;
    if (nSamples < 1 || templateSignal == nullptr)
    {
        if (nSamples < 1)
        {
            throw std::invalid_argument("No samples in template\n");
        }
        throw std::invalid_argument("templateSignal is NULL\n");
    }
    std::vector<T> work(nSamples);
    //std::copy(templateSignal, templateSignal+nSamples, work.begin());
    //pImpl->mTemplates.push_back(work);
} 
*/

/// Gets the template length
int MatchedFilterOptions::getMaxTemplateLength() const
{
    if (getNumberOfTemplates() < 1)
    {
        throw std::runtime_error("No templates set\n");
    }
    int len = 0;
    for (auto &t : pImpl->mTemplates)
    {
        len = std::max(len, t.getSignalLength());
    }
    return len;
}

/// Gets the it'th template
WaveformTemplate MatchedFilterOptions::getTemplate(const int it) const
{
    auto nt = getNumberOfTemplates();
    if (it < 0 || it >= nt)
    {
        throw std::invalid_argument("Template " + std::to_string(it)
                                  + " must be in range [0,"
                                  + std::to_string(nt) + "]\n");
    }
    return pImpl->mTemplates[it];
}

/*
template<class T>
std::vector<T> MatchedFilterOptions<T>::getTemplate(
    const int it) const
{
    auto nt = getNumberOfTemplates();
    if (it < 0 || it >= nt)
    {
        throw std::invalid_argument("Template " + std::to_string(it)
                                  + " must be in range [0,"
                                  + std::to_string(nt) + "]\n");
    }
    // Pad out the it'the template
    auto len = getTemplateLength();
    std::vector<T> work(len, 0);
    //std::copy(pImpl->mTemplates[it].begin(), pImpl->mTemplates[it].end(),
    //          work.begin());
    return work;
}
*/

/// Sets the FFT length
void MatchedFilterOptions::setFFTLength(const int fftLength)
{
    if (getNumberOfTemplates() < 1)
    {
        throw std::runtime_error("No templates set\n");
    }
    int nb = getMaxTemplateLength();
    if (fftLength < nb)
    {
        throw std::invalid_argument("fftLength = " + std::to_string(fftLength)
                                  + " must be at least " + std::to_string(nb)
                                  + "\n");
    }
    pImpl->mFFTLength = fftLength;
    pImpl->mBlockLength = pImpl->mFFTLength - nb + 1;
}

/// Gets the FFT length
int MatchedFilterOptions::getFFTLength() const
{
    if (getNumberOfTemplates() < 1)
    {
        throw std::runtime_error("No templates set\n");
    }
    if (pImpl->mFFTLength < 1)
    {
        auto nx = getSignalSize(); 
        auto nb = getMaxTemplateLength();
        auto result = computeWindowLength(nb, nx);
        pImpl->mFFTLength = result.first;
        pImpl->mBlockLength = result.second;
    }
    return pImpl->mFFTLength;
    
}

/// Gets the block length
int MatchedFilterOptions::getBlockLength() const
{
    if (getNumberOfTemplates() < 1)
    {
        throw std::runtime_error("No templates set\n");
    }
    if (pImpl->mFFTLength < 1)
    {
        auto nx = getSignalSize();
        auto nb = getMaxTemplateLength();
        auto result = computeWindowLength(nb, nx);
        pImpl->mFFTLength = result.first;
        pImpl->mBlockLength = result.second;
    }
    return pImpl->mBlockLength;
}

/// Sets the matched filter implementation
void MatchedFilterOptions::setMatchedFilterImplementation(
    const MatchedFilterImplementation implementation) noexcept
{
    pImpl->mImplementation = implementation;
}

/// Gets the matched filter implementation
MatchedFilterImplementation 
MatchedFilterOptions::getMatchedFilterImplementation() const noexcept
{
    return pImpl->mImplementation;
}

void MatchedFilterOptions::setSignalSize(const int nx) 
{
    pImpl->mFFTLength = 0;
    pImpl->mBlockLength = 0;
    if (nx < 1)
    {
        throw std::invalid_argument("Signal size must be positive\n");
    }
    pImpl->mSignalSize = nx;
}

int MatchedFilterOptions::getSignalSize() const
{
    if (pImpl->mSignalSize < 1)
    {
        throw std::runtime_error("Signal size not set\n");
    } 
    return pImpl->mSignalSize;
}

/// Toggles using the absolute value of the stack
void MatchedFilterOptions::setStackAbsoluteValues(const bool labs) noexcept
{
    pImpl->mStackAbs = labs;
}

bool MatchedFilterOptions::getStackAbsoluteValues() const noexcept
{
    return pImpl->mStackAbs;
}


/// Checks if this is a valid set of template options
bool MatchedFilterOptions::isValid() const noexcept
{
    if (getNumberOfTemplates() < 1){return false;}
    if (getMaxTemplateLength() < 1){return false;}
    if (getSignalSize() < 1){return false;}
    return true;
}

std::pair<int, int>
MFLib::computeOptimalFFTAndBlockLength(const int nb, const int nx)
{
    if (nb < 1)
    {
        throw std::invalid_argument("Filter length must be positive\n");
    }
    if (nx <= nb)
    {
        throw std::invalid_argument("Signal length = " + std::to_string(nx)
                                  + " must exceed filter length = "
                                  + std::to_string(nb) + "\n");
    }
    return computeWindowLength(nb, nx);
}

/// Template class instantiation
//template class MatchedFilterOptions<double>;
//template class MatchedFilterOptions<float>;
