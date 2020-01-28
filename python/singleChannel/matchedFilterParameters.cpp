#include <cstdio>
#include <cstdlib>
#include "mflib/enums.hpp"
#include "mflib/waveformTemplate.hpp"
#include "mflib/matchedFilterParameters.hpp"
#include "pyWaveformTemplate.hpp"
#include "singleChannel/pyMatchedFilterParameters.hpp"

using namespace PBMFLib::SingleChannel;

/// Constructor
MatchedFilterParameters::MatchedFilterParameters() :
    mParameters(std::make_unique<MFLib::SingleChannel::MatchedFilterParameters> ())
{
}

MatchedFilterParameters::MatchedFilterParameters(
    const MatchedFilterParameters &parms)
{
    *this = parms;
}

MatchedFilterParameters& MatchedFilterParameters::operator=(
    const MatchedFilterParameters &parms)
{
    if (&parms == this){return *this;}
    mParameters
        = std::make_unique<MFLib::SingleChannel::MatchedFilterParameters>
          (*parms.mParameters);
    return *this;
}

/// Destructor
MatchedFilterParameters::~MatchedFilterParameters() = default;

/// Adds a template
void MatchedFilterParameters::addTemplate(const WaveformTemplate &wt)
{
    auto t = wt.getNativeClass();
    mParameters->addTemplate(t);
}

int MatchedFilterParameters::getNumberOfTemplates() const
{
    return mParameters->getNumberOfTemplates();
}

PBMFLib::WaveformTemplate
MatchedFilterParameters::getTemplate(const int it) const
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
    PBMFLib::WaveformTemplate wtOut(wt);
    return wtOut;
}

/// FFT Length 
void MatchedFilterParameters::setFFTLength(const int fftLength)
{
    auto nt = mParameters->getNumberOfTemplates();
    if (nt < 1)
    {
        throw std::runtime_error("Templates must be set prior to calling this\n");
    }
    mParameters->setFFTLength(fftLength);
}

int MatchedFilterParameters::getFFTLength() const
{
    auto nt = getNumberOfTemplates();
    if (nt < 1){return 0;}
    return mParameters->getFFTLength();
}

int MatchedFilterParameters::getBlockLength() const
{
    auto nt = getNumberOfTemplates();
    if (nt < 1){return 0;}
    return mParameters->getBlockLength();
}

/// Sets/gets the signals size
void MatchedFilterParameters::setSignalSize(const int n)
{
    if (n < 1)
    {
        throw std::invalid_argument("Signal size must be positive\n");
    }
    mParameters->setSignalSize(n);
}

int MatchedFilterParameters::getSignalSize() const
{
    if (getNumberOfTemplates() < 1){return 0;}
    return mParameters->getSignalSize();
}

/// Clears the templates from the class
void MatchedFilterParameters::clearTemplates() noexcept
{
    mParameters->clearTemplates();
}

/// Clears the class
void MatchedFilterParameters::clear() noexcept
{
    mParameters->clear();
}

/// Get the native class
MFLib::SingleChannel::MatchedFilterParameters
MatchedFilterParameters::getNativeClass() const
{
    auto parmsOut(*mParameters);
    return parmsOut;
}

/// Initialize the class
void PBMFLib::SingleChannel::initializeMatchedFilterParameters(pybind11::module &m)
{
    pybind11::class_<PBMFLib::SingleChannel::MatchedFilterParameters>
        mfParameters(m, "MatchedFilterParameters");
    /// The default constructor:
    mfParameters.def(pybind11::init<> ());
    mfParameters.doc() = "Defines the options class for matched filtering.";
    /// Adds a template to the class
    mfParameters.def("add_template",
                     &PBMFLib::SingleChannel::MatchedFilterParameters::addTemplate,
                     "Adds a waveform template.  The waveform template must at minimum have a signal.");
    mfParameters.def("get_template",
                     &PBMFLib::SingleChannel::MatchedFilterParameters::getTemplate,
                     "Gets the it'th waveform template.");
    mfParameters.def("get_signal_size",
                     &PBMFLib::SingleChannel::MatchedFilterParameters::getSignalSize,
                     "Gets the length of the expected signal to which the template will be applied.");
    mfParameters.def("set_signal_size",
                     &PBMFLib::SingleChannel::MatchedFilterParameters::setSignalSize,
                     "Defines the size of the signals to which the templates will be applied.  This must be positive and should be called after setting all the templates.");
    mfParameters.def("get_number_of_templates",
                     &PBMFLib::SingleChannel::MatchedFilterParameters::getNumberOfTemplates,
                     "This is the number of templates that have been set in the class.");
    /// Tuning parameters
    mfParameters.def("get_fft_length",
                     &PBMFLib::SingleChannel::MatchedFilterParameters::getFFTLength,
                     "Gets the length of the Fourier transform in the overlap and add method.");
    mfParameters.def("set_fft_length",
                     &PBMFLib::SingleChannel::MatchedFilterParameters::setFFTLength,
                     "This is a tuning parameter.  From the signal size and template lengths a nominal FFT length will be computed.  While this length will be a power of 2 the underlying Cooley-Tukey Discrete Fourier Transform may be more efficient for smaller FFT lengths.  This should be set after both the templates and signal size are set."); 
    mfParameters.def("get_block_length",
                     &PBMFLib::SingleChannel::MatchedFilterParameters::getFFTLength,
                     "This is the block length in the overlap-and-add method.  It is computed from the FFT length and the signal size.");
    /// Clear the templates
    mfParameters.def("clear_templates",
                     &PBMFLib::SingleChannel::MatchedFilterParameters::clearTemplates,
                     "Clears the existing templates from class's memory.");
    /// Clears the class
    mfParameters.def("clear",
                     &PBMFLib::SingleChannel::MatchedFilterParameters::clear,
                     "Clears the class's memory and resets the class.");
}
