NOTE: For these you need the latest version of Boris - the source code is available in this repository, but exceeds the currently available 2.3 binaries version.
You could compile the code in the repository or wait for the next binaries release (2.4 or later).
v2.3 doesn't have the following commands:

dp_fitskyrmion command to fit skyrmion profiles in order to extract the radius.
loadovf2mag to load magnetization data from a OVF 2.0 file
v2.3 also has a problem with mxh relaxation when the iDMI module is enabled - this has been fixed now.

To run a batch of simulations which calculates skyrmion diameters for all BSM files run the calculate_skyrmion_diameters.py script in this directory (need Python 2.7).

To process mumax3 data run the calculate_skyrmion_diameters.py script in the respective directory - you need Boris running (> v2.3).