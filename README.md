# Notes on ULF power

The minimum values for summed and integrated PSD are 
- spsd 1E-6
- ipsd 1E-7

This was acomplised by looking at the GILL spectra and narrowing in on power spectra that with small values to identify periods where the magnetometer where showing no power. 

An example of this is 2014-03-13. The magnetometer has no bad data however the spectra is flat. The sum power returned is below 1E-17. 

Another example is 2014-06-23. The magnetometer has a short period of data followed by a flat line. The summed power returned is 5.25E-5.

Taking this into account it seems sensible to drop summed power below 1E-6 and integrated power below 1E-7.