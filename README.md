# polair
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21065056.svg)](https://doi.org/10.5281/zenodo.21065056)

Current documentation: https://laurakoehler.github.io/polair/

This package provides a command line tool for processing atmospheric data from the polar research aircrafts Polar 5 and 6. It can be used to process data downloaded from the [DMS](https://dms.awi.de/exportdisplay/) system, convert the values to physical quantities, calibrate, and homogenise the time stamp. Furthermore, it can be used to process the noseboom and the t-bird to obtain corrected physical variables including the wind components. As a final step, it can be used to create standardised netCDF data sets with metadata and flight information. For this, we recommend to document all information needed to process a specific campaign in a dedicated repository. To use the polair package, you need to provide a configuration file with information of the campaign and several yaml dictionaries with information on calibration, variables of interest, units, flight segments, and so on. The BACSAM II campaign repository can serve as a blueprint.

## Setup

Install polair:
```
pip install polair
```

