# MFLib - A Matched Filtering Library

Contained is the source code for performing matched filtering detection.  In principle, the algorithm is quite simple in that it computes a [Pearson correlation coefficient](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#For_a_sample) at every sample in a time series corresponding to a template.  However, the actual implementation in a compiled language is tedious.

## Installation

See the [installation notes](https://github.com/uofuseismo/mflib/blob/master/INSTALL.md).

## Other Options

   1.  Super-Efficient Cross-Correlation [SEC-C](https://github.com/Naderss/SEC_C) - if you are a Python or Matlab user looking to do template matching on a desktop computer and dislike the headache that is compilation then this may be the way to go.  There is actually a good deal of similarity since this library implements in a compiled language many of the speed-ups that the SEC-C group defined in their paper.
   2.  [Fast-Matched-Filter](https://github.com/beridel/fast_matched_filter) - this is a solid time-domain choice.  The core of this appears to be written in C so you'll likely see some pretty good performance.  Additionally, it has Python and Matlab bindings.

## License

This software is distributed under the MIT license.  The details are [here](https://github.com/uofuseismo/mflib/blob/master/LICENSE).
