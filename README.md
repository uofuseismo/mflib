# MFLib - A Matched Filtering Library

Contained is the source code for performing matched filtering detection.  In principle, the algorithm is quite simple in that it computes a [Pearson correlation coefficient](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#For_a_sample) at every sample in a time series corresponding to a template.  However, actual implementation can be tedious.  Some things this library provides are:

   1.  Ease of use.  Using any compiled language API is difficult.  However, we look to simplify it to five steps - initialize, set the data, apply the calculator, get the result, release the resources.
   2.  Low-latency Python bindings that directly access the C++ library.
   3.  Performance - this library takes care of the fun things like cache alignment of memory, vectorization, and threading as well as interfacing with high-performance libraries for the heavy computations.

## Installation

See the [installation notes](https://github.com/uofuseismo/mflib/blob/master/INSTALL.md).

## Other Options

   1.  Super-Efficient Cross-Correlation [SEC-C](https://github.com/Naderss/SEC_C) - if you are a Python or Matlab user looking to template matching on a desktop computer and dislike the headache that is compilation then this may be the way to go.  There's actually a good deal of similarity since this library implements in a compiled language many of the speed-ups that the SEC-C group defined in their paper.
   2.  [Fast-Matched-Filter](https://github.com/beridel/fast_matched_filter) - this is another solid time-domain choice.  The core of this appears to be written in C so you'll likely see some pretty good performance.  Additionally, it has Python and Matlab bindings.

## License

This software is distributed under the MIT license.  The details are [here](https://github.com/uofuseismo/mflib/blob/master/LICENSE).
