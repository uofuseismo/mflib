#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include "mflib/matchedFilterOptions.hpp"
#include "mflib/enums.hpp"

using namespace MFLib;

namespace
{
/// @brief Computes sizes.
/// @param[in] nb   The filter length.
/// @param[in] nx   The signal length.
/// @result result.first is the FFT size and result.second is the window length.
std::pair<int, int> computeWindowLength(const int nb, const int nx) 
{
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
        if (fftLengths[i] > nb - 1)
        {
            validSetStart = i;
            break;
        }
    }
    if (validSetStart ==-1)
    {
        throw std::runtime_error("Could not find appropriate padding length\n");
    }
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


template<class T>
class MatchedFilterOptions<T>::MatchedFilterOptionsImpl
{
public:
    /// Collection of templates
    std::vector<std::vector<T>> mTemplates;
    /// Matched filter implementation
    MatchedFilterImplementation mImplementation
        = MatchedFilterImplementation::DIRECT;
    /// Signal size
    int mSignalSize = 0;
    /// FFT length
    int mFFTLength = 0;
    /// The block length.  Given a filter length, nb, and FFT length,
    /// the optimal length is  L = mFFTLength - nb - 1.
    int mBlockLength = 0;
};

/// Constructor
template<class T>
MatchedFilterOptions<T>::MatchedFilterOptions() :
    pImpl(std::make_unique<MatchedFilterOptionsImpl> ())
{
}

/// Copy constructor
template<class T>
MatchedFilterOptions<T>::MatchedFilterOptions(
    const MatchedFilterOptions &options)
{
    *this = options;
}

/// Move constructor
template<class T>
MatchedFilterOptions<T>::MatchedFilterOptions(
    MatchedFilterOptions &&options) noexcept
{
    *this = std::move(options);
}

/// Copy assignment
template<class T>
MatchedFilterOptions<T>& MatchedFilterOptions<T>::operator=(
    const MatchedFilterOptions &options)
{
    if (&options == this){return *this;}
    pImpl = std::make_unique<MatchedFilterOptionsImpl> (*options.pImpl);
    return *this;
}

/// Move assignment
template<class T>
MatchedFilterOptions<T>& MatchedFilterOptions<T>::operator=(
    MatchedFilterOptions &&options) noexcept
{
    if (&options == this){return *this;}
    pImpl = std::move(options.pImpl);
    return *this;
}

/// Destructor
template<class T>
MatchedFilterOptions<T>::~MatchedFilterOptions() = default;

/// Clears the memory and resets the class
template<class T>
void MatchedFilterOptions<T>::clear() noexcept
{
    pImpl->mTemplates.clear();
    pImpl->mImplementation = MatchedFilterImplementation::DIRECT;
    pImpl->mSignalSize = 0;
}

/// Get the number of templates
template<class T>
int MatchedFilterOptions<T>::getNumberOfTemplates() const noexcept
{
    return static_cast<int> (pImpl->mTemplates.size());
}

/// Clears all the templates
template<class T>
void MatchedFilterOptions<T>::clearTemplates() noexcept
{
    pImpl->mTemplates.clear();
}

/// Adds a template
template<class T>
void MatchedFilterOptions<T>::addTemplate(
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
    std::copy(templateSignal, templateSignal+nSamples, work.begin());
    pImpl->mTemplates.push_back(work);
} 

/// Gets the template length
template<class T>
int MatchedFilterOptions<T>::getTemplateLength() const
{
    if (getNumberOfTemplates() < 1)
    {
        throw std::runtime_error("No templates set\n");
    }
    size_t len = 0;
    for (auto &t : pImpl->mTemplates)
    {
        len = std::max(len, t.size());
    }
    return static_cast<int> (len); 
}

/// Gets the it'th template
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
    std::copy(pImpl->mTemplates[it].begin(), pImpl->mTemplates[it].end(),
              work.begin());
    return work;
}

/// Sets the FFT length
template<class T>
void MatchedFilterOptions<T>::setFFTLength(const int fftLength)
{
    if (getNumberOfTemplates() < 1)
    {
        throw std::runtime_error("No templates set\n");
    }
    int nb = getTemplateLength();
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
template<class T>
int MatchedFilterOptions<T>::getFFTLength() const
{
    if (getNumberOfTemplates() < 1)
    {
        throw std::runtime_error("No templates set\n");
    }
    if (pImpl->mFFTLength < 1)
    {
        auto nx = getSignalSize(); 
        auto nb = getTemplateLength();
        auto result = computeWindowLength(nb, nx);
        pImpl->mFFTLength = result.first;
        pImpl->mBlockLength = result.second;
    }
    return pImpl->mFFTLength;
    
}

/// Gets the block length
template<class T>
int MatchedFilterOptions<T>::getBlockLength() const
{
    if (getNumberOfTemplates() < 1)
    {
        throw std::runtime_error("No templates set\n");
    }
    if (pImpl->mFFTLength < 1)
    {
        auto nx = getSignalSize();
        auto nb = getTemplateLength();
        auto result = computeWindowLength(nb, nx);
        pImpl->mFFTLength = result.first;
        pImpl->mBlockLength = result.second;
    }
    return pImpl->mBlockLength;
}

/// Sets the matched filter implementation
template<class T>
void MatchedFilterOptions<T>::setMatchedFilterImplementation(
    const MatchedFilterImplementation implementation) noexcept
{
    pImpl->mImplementation = implementation;
}

/// Gets the matched filter implementation
template<class T>
MatchedFilterImplementation 
MatchedFilterOptions<T>::getMatchedFilterImplementation() const noexcept
{
    return pImpl->mImplementation;
}

template<class T>
void MatchedFilterOptions<T>::setSignalSize(const int nx) 
{
    pImpl->mFFTLength = 0;
    pImpl->mBlockLength = 0;
    if (nx < 1)
    {
        throw std::invalid_argument("Signal size must be positive\n");
    }
    pImpl->mSignalSize = nx;
}

template<class T>
int MatchedFilterOptions<T>::getSignalSize() const
{
    if (pImpl->mSignalSize < 1)
    {
        throw std::runtime_error("Signal size not set\n");
    } 
    return pImpl->mSignalSize;
}

/// Checks if this is a valid set of template options
template<class T>
bool MatchedFilterOptions<T>::isValid() const noexcept
{
    if (getNumberOfTemplates() < 1){return false;}
    if (getTemplateLength() < 1){return false;}
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
template class MatchedFilterOptions<double>;
template class MatchedFilterOptions<float>;
