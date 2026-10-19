# polair
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21065056.svg)](https://doi.org/10.5281/zenodo.21065056)

Current documentation: https://laurakoehler.github.io/polair/

This package provides a command line tool for processing atmospheric data from the polar research aircrafts Polar 5 and 6. It can be used to process data downloaded from the [DMS](https://dms.awi.de/exportdisplay/) system, convert the values to physical quantities, calibrate, and homogenise the time stamp. Furthermore, it can be used to process the noseboom and the t-bird to obtain corrected physical variables including the wind components. As a final step, it can be used to create standardised netCDF data sets with metadata and flight information. For this, we recommend to document all information needed to process a specific campaign in a dedicated repository. To use the polair package, you need to provide a configuration file with information of the campaign and several yaml dictionaries with information on calibration, variables of interest, units, flight segments, and so on. The BACSAM II campaign repository can serve as a blueprint.

## Setup

Install polair:
```
pip install polair
```

### Folder structure

Most of the folders can be defined in the config file. However, if the following structure is kept, not all paths need to be specified in the config file for each flight.

The recommended structure for the input directory for each flight is:
```text
flight_dir/
├── data_dir/
│   ├── prefix_DPM7885_3.SPS_ATEL.dat
│   └── ...
├── dropsondes/
│   ├── main.py
│   └── dropsonde1/
│       ├── CYYYYMMDD_hhmmss.csv
│       └── ...
├── mcpc/
│   ├── YYMMDDA1.TXT
│   └── ...
├── partector/
│   ├── <sn>_hhmmss.txt
│   └── ...
├── uhsas/
│   ├── UHSASYYYYMMDDhhmmss.xls
│   └── ...
├── licor/
│   ├── prefix_LICOR72.CO2D.dat
│   └── ...
├── instrument/
│   └── ...
└── ...
```
  

The recommended structure for the output directory for the campaign is:
```text
processed/
├── noseboom/
├── tbird/
├── dropsondes
├── <instrumen1>
├── <instrumen2>
└── ...
```

The recommended structure for the campaign directory containing all the campaign specific polair information is:
```text
campaign_processing/
├── config.yaml
├── final_variable_attributes.yaml
├── flight_segments.yaml
├── README.md
├── turbulence/
│   ├── calibration.yaml
│   ├── fiveholeprobe_parameters.yaml
│   ├── met_variables.yaml
├── dms_download_templates/
│   ├── template1.xml
│   └── ...
├── notebooks/
│   ├── Manual_processing_steps.ipynb
│   └── Quicklooks.ipynb
└── .gitignore
```

Of course, all these directories can be in one common folder.


### Configuration file

The configuration file contains all campaign specific information. It is a yaml dictionary and for all polair comands, you have to provide the path to the config file. In the following, we show the basic structure with all options. Some dictionary entries are necessary to run the different polair commands, others are optional.

```yaml
campaign:
  name: <campaign_name> # The campaign name used in the file names
  platform: <platform_name> # e.g. Polar 6
  tbird_reversed_antennas: True # required for tbird command, for some campaigns, the antennas in the T-Bird iNAT system were switched which would lead to wrong signs.
  noseboom_twist: 0 # required for noseboom command, angle in degree, leave 0 by default. A twist angle can be defined if the noseboom is twisted in its mounting
  stp: [mcpc, partector, partector_dms] # required for device command: for instruments listed here, the number concentration is reduced to standard temperature and pressure and added as a new variable. Default instruments are mcpc and partector.
  resampling: [<instrument1>, <instrument2>] # optional for device command. Instruments listed here will be resampled to 1 s temporal resolution
  flow_p_dep: [<instrument1>, <instrument2>] # optional for device command. For instruments listed here, the flow correction will include the pressure and use flow/pressure instead of flow as a measure if the flow was out of range or not.
  dev_cut_sec: # optional for device command. For instruments listed here, the time period will be cut at the beginning and end, e.g. if an instrument needs time to start or there are often contaminations from the emissions
    <instrument1>: # instrument where start and stop are modified
      start: <delay_in_sec> # time to start the data set later in seconds
      stop: <ealier_end_in_sec> # time to end the data set earlier in seconds

paths:
  xml_file: <path_to_dms_template> # optional. Path to DMS download template
  variables: <path_to_met_vars> # required for preprocessing command. Yaml file with information on sensors, instruments DMS names, etc.
  calibration: <path_to_calibration_dic> # required for preprocessing command. Contains the calibration information of the DMS sensors
  fiveholeprobe: <path_to_turbulence_calibration> # required for noseboom and tbird command. Contains the turbulence calibration. Calibration parameters are determined iteratively and manually for each campaign, see notebooks/Manual_processing_steps.ipynb.
  processing_log_file: <path_to_log_file> # required. Path of the log file. If it doesn't exist, it will be created.
  processed_variables: <path_to_vars> # required for all commands except preprocessing. Dictionary containing the final varibales with attributes.
  segments: <path_to_segment_yaml> # optional for all commands except preprocessing. Dictionary containing the flight segments.
  outdirs: # required. List of output directories
    raw: <output_dir_preprocessing> # required for preprocessing command. Directory where output files are written to.
    noseboom: <outout_dir_noseboom> # required for noseboom command. Directory where output files are written to.
    tbird: <output_dir_tbird> # required for tbird command. Directory where output files are written to.
    <instrument1>: <output_dir_instrument1> # required for device -i instrument1 command. Directory where output files are written to.
    <instrument2>: <output_dir_instrument2> # required for device -i instrument2 command. Directory where output files are written to. Add all instruments needed for the campaign
    dropsondes: <output_dir_tbird> # required for dropsonde command. Directory where output file is written to.

metadata: # optional but recommended. General metadata for the campaign
  author: <author>
  institution: <institution>
  contact: <email>
  date: <date>
  version: <version>

instrument_metadata: # optional. Metadata will be added as global attributes to the data sets.
  noseboom: # optional for noseboom command.
    title: <title> # end with "... data of", date and research flight number will be added automatically, e.g. for research flight 01 on 2026-01-21 the title will be "... data of RF01 (2026-01-21)"
    author: <authors>
    institution: <institutions>
    contact: <email>
    version: <version>
    doi: <doi>
    conventions: <conventions> # e.g. "CF-1.7"
    processing: <information_on_processing> # e.g. cite polair here: "processed with polair (https://github.com/LauraKoehler/polair, doi: 10.5281/zenodo."
    comment: <comment> 
    date: <processing_date> # optional. If no date is defined, the date of today will be written.
  tbird: # optional for tbird command
    title: <title> # end with "... data of", date and research flight number will be added automatically, e.g. for research flight 01 on 2026-01-21 the title will be "... data of RF01 (2026-01-21)"
    author: <authors>
    institution: <institutions>
    contact: <email>
    version: <version>
    doi: <doi>
    conventions: <conventions> # e.g. "CF-1.7"
    processing: <information_on_processing> # e.g. cite polair here: "processed with polair (https://github.com/LauraKoehler/polair, doi: 10.5281/zenodo."
    comment: <comment> 
    date: <processing_date> # optional. If no date is defined, the date of today will be written.
  <instrument1>: # optional for devive -i instrument1 command.
    title: <title> # end with "... data of", date and research flight number will be added automatically, e.g. for research flight 01 on 2026-01-21 the title will be "... data of RF01 (2026-01-21)"
    author: <authors>
    institution: <institutions>
    contact: <email>
    version: <version>
    doi: <doi>
    conventions: <conventions> # e.g. "CF-1.7"
    processing: <information_on_processing> # e.g. cite polair here: "processed with polair (https://github.com/LauraKoehler/polair, doi: 10.5281/zenodo."
    comment: <comment> 
    date: <processing_date> # optional. If no date is defined, the date of today will be written.
  dropsondes: # optional for the dropsonde command.
    title: <title> # Title of the data set. Since all dropsondes are combined to one output file, the full title can be entered here
    author: <authors>
    institution: <institutions>
    contact: <email>
    version: <version>
    doi: <doi>
    conventions: <conventions> # e.g. "CF-1.7"
    processing: <information_on_processing> # e.g. cite polair here: "processed with polair (https://github.com/LauraKoehler/polair, doi: 10.5281/zenodo."
    comment: <comment> 

flights: # required.
  <number>: # flight number
    date: <flight_date> # required
    flight_dir: <flight_directory> # required for dropsondes command (if there was a sounding at the flight, the data should be in a subfolder of flight_dir called "dropsondes"), recommended for device command. Directory where all raw data of the research flight is collected.
    data_dir: <dms_data_directory> # required for preprocessing, device -i radiation and device -i kt19. Data directory with meteorological, turbulence, irradiance and kt19 data downloaded from DMS.
    <instrument1>: <directory_to_instrument1_data> # optional for device -i instrument1 command. If no directory is specified, the data should be stored in a subfolder named "instrument1" in flight_dir, see folder structure.
    <instrument2>: <directory_to_instrument1_data> # optional for device -i instrument1 command. If no directory is specified, the data should be stored in a subfolder named "instrument1" in flight_dir, see folder structure.
    prefix: <dms_prefix_of_flight> # required for preprocessing command. Specify the common prefix of the files downloaded from DMS.
    tbird: <True or False> # required for preprocessing command. Specify if T-Bird was used during the flight.
    start: <start_time_of_the_flight> # required for all commands except preprocessing
    stop: <stop_time_of_the_flight> # required for all commands except preprocessing
    start_tbird: <start_time_of_the_tbird> # required for tbird command. Only applied if T-Bird was used
    stop_tbird: <stop_time_of_the_tbird> # required for tbird command. Only applied if T-Bird was used
    calibration: # optional. Specify any turubulence calibration segments here.
      out1: # first calibration out
        start: <out_leg_start> # start time of the out leg
        stop: <out_leg_stop> # stop time of the out leg
      ret1: # first calibration ret
        start: <ret_leg_start> # start time of the ret leg
        stop: <ret_leg_stop> # stop time of the ret leg
      out2: # Add as many calibrations as necessary...
        start: <out_leg_start> # start time of the out leg
        stop: <out_leg_stop> # stop time of the out leg
      ret2: # first calibration ret
        start: <ret_leg_start> # start time of the ret leg
        stop: <ret_leg_stop> # stop time of the ret leg
    time_offsets: # optional for device command. Specify time shifts between DMS and instrument times. This is especially relevant for instruments where the times are adjusted manually at the beginning of each flight. For the time shift, the pressure correlations can be used, see notebooks/Manual_processing_steps.ipynb.
      <instrument1>: <time_shift_in_ms> # Time shift of instrument1 with respect to the turbulence data in miliseconds.
      <instrument2>: <time_shift_in_ms> # Time shift of instrument2 with respect to the turbulence data in miliseconds.
```