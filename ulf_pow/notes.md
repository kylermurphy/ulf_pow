# Notes on ULF power

The minimum values for summed and integrated PSD are approximately
- spsd 1E-6
- ipsd 1E-7

This was determined by looking at the GILL spectra and narrowing in on power spectra that had small values. From these spectra it periods where the magnetometer had no power where identified (e.g. flat lines). These days where used to set determine the threshold of summed and integrated PSD. 

An example of this is 2014-03-13. The magnetometer has no bad data however the spectra is flat. The sum power returned is below 1E-17. 

Another example is 2014-06-23. The magnetometer has a short period of data followed by a flat line. The summed power returned is 5.25E-5.

Taking this into account it seems appropriate to drop summed power below 1E-6 and integrated power below 1E-7.