TGS Analyzer
====================
GUI program for extracting and analyzing spectra from astronomical images taken using a
transmission grating spectrometer (TGS).

Using/Installing TGS Analyzer
====================
TGS Analyzer can be run using the Python (2.7) interpreter. Required python libraries are:
- scipy
- matplotlib
- pyfits
- wxpython

Simply copy the script, config, and help files to the desired location to install it. The location of the python
interpreter must be in the path environment variable.

Standalone Executables
====================
PyInstaller has been used to compile TGS Analyzer into standalone executables for some systems. These
can be run without python installed. 
These are found in exec/. PyInstaller can be obtained from https://github.com/pyinstaller/pyinstaller.
To compile for a different system using pyinstaller, clone this repository and pyinstaller on that
machine, and run the command

      pyinstaller tgsa.py -F -i <icon_file>

where *icon_file* is the appropriate icon file found in icons/, if desired.



MIT License
