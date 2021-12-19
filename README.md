# FATIMA_152Eu_autocal
Automatic calibration for LaBr3 energy spectra based on 152Eu source measurements.

## Usage:
`./FATIMA_152EU_autocal.py <histogram file name> (<detnumber>)`

### Input:
#### ASCII format histogram file.
This file has to be in two-column format, first column channel number, second
column frequency (counts)
#### Detector number
Optionally a detector number can be given. This will appear in the first
column of the output file.

### Output:
The output is a textfile containing the input file name (or detector number
if given) and the calibration parameters which were determined in the fit.
Additionally a window will open with a plot of the fitted peaks and them
residual of the final calibration. This plot will also be saved as a .png
file for reference.

## Working principle:
The code works in several stages. A brief description (omitting many subtleties)
is given below:

#### Preparation:
First the input histogram is prepared for analysis. It is checked if
the content is sufficient for calibration. Also the specturm is rebinned
and truncated based on some assumptions on the 1408 keV line from the
152Eu decay spectrum.

#### Peak finding:
The code then searches for peaks and tries to identify them with
strong transitions of the 152Eu decay spectrum.

There are parameters which cut peaks of small prominence. These are set
to reasonable values, but have to be adjusted potentially for each case.
Once the peak list is found another check is done. This tries to identify
the peaks with known strong peaks from the 152Eu spectrum, like 121 keV,
779 keV etc. The programme uses approximate relative intensities to do this
check. These can vary for different detector types and might have to be
adjusted. Currently the values are optimised for use with the FATIMA core
LaBr3 detectors.

If the identification fails, the programme will abort.

#### Fitting procedure:
After the peaks are identified the programme will determine each peak's position
using a simple Gauss+BG fit. There are special procedures for peaks that are
close together (411,444 and 1089,1112). With the determined values it will
first attempt a calibration
fit with a 3rd order polynomial. If this doesn't produce a good result it
will do another fit with a 4th order polynomial.

The result is then written to a file. The input spectrum will then be plotted,
along with markers that indicate the identified peaks and the fits of each peak.
Furthermore, the fit residual will be shown in a different panel.

The output can be used to convince yourself that the identification and
fitting has been done properly.
