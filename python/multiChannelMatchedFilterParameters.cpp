#include <cstdio>
#include <cstdlib>
#include "mflib/enums.hpp"
#include "mflib/waveformTemplate.hpp"
#include "mflib/matchedFilterParameters.hpp"
#include "pyWaveformTemplate.hpp"
#include "pyMultiChannelMatchedFilterParameters.hpp"

using namespace PBMFLib;

/// Constructor
MultiChannelMatchedFilterParameters::MultiChannelMatchedFilterParameters() :
    mParameters(std::make_unique<MFLib::MatchedFilterParameters> ())
{
    mParameters->setDetectionMode(
        MFLib::MatchedFilterDetectionMode::MULTI_CHANNEL);
}

MultiChannelMatchedFilterParameters::MultiChannelMatchedFilterParameters(
    const MultiChannelMatchedFilterParameters &parms)
{
    *this = parms;
}

MultiChannelMatchedFilterParameters& 
MultiChannelMatchedFilterParameters::operator=(
    const MultiChannelMatchedFilterParameters &parms)
{
    if (&parms == this){return *this;}
    mParameters
        = std::make_unique<MFLib::MatchedFilterParameters> (*parms.mParameters);
    return *this;
}

/// Destructor
MultiChannelMatchedFilterParameters::~MultiChannelMatchedFilterParameters()
    = default;

/// Adds a template
void MultiChannelMatchedFilterParameters::addTemplate(
    const WaveformTemplate &wt)
{
    auto t = wt.getNativeClass();
    mParameters->addTemplate(t);
}

int MultiChannelMatchedFilterParameters::getNumberOfTemplates() const
{
    return mParameters->getNumberOfTemplates();
}

WaveformTemplate 
MultiChannelMatchedFilterParameters::getTemplate(const int it) const
{
    if (it < 0 || it >= getNumberOfTemplates())
    {
        throw std::invalid_argument("template index = "
                                  + std::to_string(it)
                                  + " must be in range [0,"
                                  + std::to_string(getNumberOfTemplates())
                                  + "]\n");
    }
    auto wt = mParameters->getTemplate(it);
    WaveformTemplate wtOut(wt);
    return wtOut;
}

/// FFT Length 
void MultiChannelMatchedFilterParameters::setFFTLength(const int fftLength)
{
    auto nt = mParameters->getNumberOfTemplates();
    if (nt < 1)
    {
        throw std::runtime_error("Templates must be set prior to calling this\n");
    }
    mParameters->setFFTLength(fftLength);
}

int MultiChannelMatchedFilterParameters::getFFTLength() const
{
    auto nt = getNumberOfTemplates(); //mParameters->getNumberOfTemplates();
    if (nt < 1){return 0;}
    return mParameters->getFFTLength();
}

int MultiChannelMatchedFilterParameters::getBlockLength() const
{
    auto nt = getNumberOfTemplates(); //mParameters->getNumberOfTemplates();
    if (nt < 1){return 0;}
    return mParameters->getBlockLength();
}

/// Sets/gets the signals size
void MultiChannelMatchedFilterParameters::setSignalSize(const int n)
{
    if (n < 1)
    {
        throw std::invalid_argument("Signal size must be positive\n");
    }
    mParameters->setSignalSize(n);     
}

int MultiChannelMatchedFilterParameters::getSignalSize() const
{
    if (getNumberOfTemplates() < 1){return 0;}
    return mParameters->getSignalSize();
}

/// Clears the templates from the class
void MultiChannelMatchedFilterParameters::clearTemplates() noexcept
{
    mParameters->clearTemplates();
}

/// Clears the class
void MultiChannelMatchedFilterParameters::clear() noexcept
{
    mParameters->clear();
}

/// Get the native class
MFLib::MatchedFilterParameters
MultiChannelMatchedFilterParameters::getNativeClass() const
{
    auto parmsOut(*mParameters);
    return parmsOut;
}

/// Initialize the class
void PBMFLib::initializeMultiChannelMatchedFilterParameters(pybind11::module &m)
{
    pybind11::class_<PBMFLib::MultiChannelMatchedFilterParameters>
        mfParameters(m, "MatchedFilterParameters");
    /// The default constructor:
    mfParameters.def(pybind11::init<> ());
    mfParameters.doc() = "Defines the options class for matched filtering.";
    /// Adds a template to the class
    mfParameters.def("add_template",
                     &PBMFLib::MultiChannelMatchedFilterParameters::addTemplate,
                     "Adds a waveform template.  The waveform template must at minimum have a signal.");
    mfParameters.def("get_template",
                     &PBMFLib::MultiChannelMatchedFilterParameters::getTemplate,
                     "Gets the it'th waveform template.");
    mfParameters.def("get_signal_size",
                     &PBMFLib::MultiChannelMatchedFilterParameters::getSignalSize,
                     "Gets the expected length of the signal to which the matched filter will be applied.");
    mfParameters.def("set_signal_size",
                     &PBMFLib::MultiChannelMatchedFilterParameters::setSignalSize,
                     "Defines the size of the signals to which the templates will be applied.  This must be positive and should be called after setting all the templates.");
    mfParameters.def("get_number_of_templates",
                     &PBMFLib::MultiChannelMatchedFilterParameters::getNumberOfTemplates,
                     "This is the number of templates that have been set in the class.");
    /// Tuning parameters
    mfParameters.def("get_fft_length",
                     &PBMFLib::MultiChannelMatchedFilterParameters::getFFTLength,
                     "Gets the length of the DFT in the overlap and add method.");
    mfParameters.def("set_fft_length",
                     &PBMFLib::MultiChannelMatchedFilterParameters::setFFTLength,
                     "This is a tuning parameter.  From the signal size and template lengths a nominal FFT length will be computed.  While this length will be a power of 2 the underlying Cooley-Tukey Discrete Fourier Transform may be more efficient for smaller FFT lengths.  This should be set after both the templates and signal size are set."); 
    mfParameters.def("block_length",
                     &PBMFLib::MultiChannelMatchedFilterParameters::getFFTLength,
                     "This is the block length in the overlap-and-add method.  It is computed from the FFT length and the signal size.");
    /// Clear the templates
    mfParameters.def("clear_templates",
                     &PBMFLib::MultiChannelMatchedFilterParameters::clearTemplates,
                     "Clears the existing templates from class's memory.");
    /// Clears the class
    mfParameters.def("clear",
                     &PBMFLib::MultiChannelMatchedFilterParameters::clear,
                     "Clears the class's memory and resets the class.");
}
