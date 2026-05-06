# polair

This package provides a command line tool for processing atmospheric data from the polar research aircrafts Polar 5 and 6. It can be used to process data downloaded from the [DMS](https://dms.awi.de/exportdisplay/) system, convert the values to physical quantities, calibrate, and homogenise the time stamp. Furthermore, it can be used to process the noseboom and the t-bird to obtain corrected physical variables including the wind components. As a final step, it can be used to create standardised netCDF data sets with metadata and flight information. For this, we recommend to document all information needed to process a specific campaign in a dedicated repository. To use the polair package, you need to provide a configuration file with information of the campaign and several yaml dictionaries with information on calibration, variables of interest, units, flight segments, and so on. The BACSAM II campaign repository can serve as a blueprint.

## preprocessing

The preprocessing command converts all data in physical quantities. If necessary and provided, it calibrates the data. Finally, they are interpolated on common (default 100 Hz) time stamps. The command also checks for data gaps and sampling errors which are noted in the log file. The calibrated raw data will be saved as a netCDF file. The location needs to be specified in the configuration file of the campaign
```
polair preprocessing -f <flight number> -c <config file>
```

## noseboom

The noseboom command determines the basic meteorological parameters from the calibrated raw data of the noseboom which can be mounted at the nose of Polar 5 and 6. Adiabatic corrections are used and the wind components are determined.
```
polair noseboom -f <flight number> -c <config file>
```

## tbird

The tbird command determines the meteorological parameters from the T-Bird similar to the noseboom command.
```
polair tbird -f <flight number> -c <config file>
```

# Useful code

## Find corrupted lines in DMS files

When downloading the data from the DMS system, it happens that there errors in the files such as broken lines, special charachters, missing entries,... When this is the case, the polair package fails in processing these data. Thus, the corrupted lines have to be removed. Here, we provide some code examples which can be used to quickly identify the corrupted line. In most of the cases, this works:

```python
fn = <corrupted_file>

with open(fn, "rb") as f:
    for i, line in enumerate(f, start=1):
        try:
            line.decode("utf-8")
        except UnicodeDecodeError as e:
            print(f"Unicode error in line {i}: {e}")
            print(line[:200])  # show beginning of bad line
            break
```
This gives the corrupted line number. However, empty lines are not counted. So, if every second line is empty, you need to double the result. If this is not working, try

```python
fn = <corrupted_file>

df = pd.read_csv(fn, header  = 4, sep = r'\s+', names = ["date", "time", f"var"])
bad_idx = [i for i, x in enumerate(df["time"]) if len(str(x)) != 15]
```
bad_idx gives lines with broken timestamps. Duplicated timestamps could also cause problems:

```python
dup_times = df.loc[df["time"].duplicated(), "time"]
```
After removing the corrupted lines manually, the polair commands should work.

## Convert yaml to html