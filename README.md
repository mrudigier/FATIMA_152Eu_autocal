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

## Working principle:
The code works in several stages.

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
779 keV etc. If this fails, the programme will abort.

#### Fitting procedure:
After the peaks are identified the programme will determine each peaks position
using a simple Gauss+BG fit. There are special procedures for peaks that are
close together (411,444 adn 1089,1112). With the determined values it Will
first attempt a calibration
fit with a 3rd order polynomial. If this doesn't produce a good result it
will do another fit with a 4th order polynomial.

The result is then written to a file. The input spectrum will then be plotted,
along with markers that indicate the identified peaks and the fits of each peak.
Furthermore, the fit residual will be shown in a different panel.

The output can be used to convince yourself that the identification and
fitting has been done properly.
