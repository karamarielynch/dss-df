# DSS Data Fitter

![alt text](https://img.shields.io/badge/License-MIT-blue.svg 'License')
![alt text](https://img.shields.io/badge/Python-3.5-green.svg 'Python version')
![alt text](https://img.shields.io/badge/Tested_on-Linux/Mac-green.svg 'Supported platform')
![alt text](https://img.shields.io/badge/Not_tested_on-Windows-red.svg 'Unsupported platform')

## Purpose
DSS Data Fitter is a GUI for viewing and fitting the DSS data collected with the MIDAS software and CAEN data acquisition hardware. The DSS is the Decay Spectroscopy Station, installed at the [Collinear Resonance Ionization](http://isolde-cris.web.cern.ch/) (CRIS) experiment, at the ISOLDE facility, CERN.

## Contributors
- Kara Marie Lynch (kara.marie.lynch@cern.ch)
- Wouter Gins (wouter.gins@fys.kuleuven.be)

## Features
- Written in Python 3.5 and PyQt4
- Tested on Mac OS X and Linux 

## Installation
- Install DSS Data Fitter directory in location of choice

## Dependencies
- [NumPy](http://www.numpy.org/)
- [PyQt4](https://pypi.python.org/pypi/PyQt4)
- [PyQtGraph](http://www.pyqtgraph.org/)
- [Matplotlib](http://matplotlib.org/)

## Usage
Create `config.txt` file to identify the CAEN channels used when acquiring the data
### Basic usage:
`$ dss_data_fitter.py config.txt`

## Functions:
Tab 1: ‘Plot Data’
- Load data from file
- Choose detector to view
- Change calibration of channel-to-energy conversion

Tab 2: ‘Fit Data’
- Fit 1 or 2 peaks to the data
- Choose line profiles (Gaussian or CrystalBall so far)
- Fit routine outputs fit parameters and reduced chi2

Tab 3: ‘Data Rate’
- Total count rate for data file
- Gated count rate for particular energy range

Tab 4: ‘Coincidences’
- Choose time window for coincidences
- Choose which detector to plot
- Choose which detector to gate on
- Total number of coincidence events
- Create 2D histogram of detector-detector coincidences in time window
- Create histogram of time-delta for successive events in energy range

