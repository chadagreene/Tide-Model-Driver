# 🚧 Tide Model Driver 3.0 🚧 
Tide Model Driver for MATLAB, version 3.0

🚧 This repo is under construction. 🚧

# Installation 
To start predicting tides, just add TMD to MATLAB and download your favorite tide model data. Each step is described in the links below: 

1. [Add TMD to MATLAB](doc/installing_tmd.md), & 
2. [Get tide model data](doc/tide_model_data.md). 

# TMD Functions
For most applications, **`tmd_predict` is the only function you will need to call directly**, although in some cases you may want `tmd_interp` to retrieve bathymetry, or `tmd_ellipse` to tell you about tidal currents. Documentation is provided for all of the functions with hyperlinked names below. 

#### Main functions:
* `tmd_predict` predicts tidal elevation, transport, or current velocities for given location(s) and time(s). 
* `tmd_interp` provides water column thickness, land/ocean mask, ice shelf flexure, and tidal constiuent parameters at specified geographic locations, for a given tide model. 
* [`tmd_data`](doc/tmd_data_documentation.md) loads gridded tide model data without interpolation. 
* `tmd_ellipse` gives tidal ellipse parameters at specified location(s).  

#### Under-the-hood functions:

* `tidal_range` calculates the peak-to-peak tidal height range. 
* `tmd_astrol` computes the basic astronomical mean longitudes s, h, p, N.
* `tmd_constit` returns amplitude, phase, frequency, alpha, species for tidal constituents. 
* `tmd_harp` predicts tidal time series using harmonic constants. 
* `tmd_InferMinor` returns correction for 16 minor tidal constiuents. 
* `tmd_nodal` calculates the nodal corrections for tidal constituents.
* `tmd_ll2ps` converts geographic coordinates to polar stereographic kilometers for polar models. 
* `tmd_ps2ll` converts polar stereographic kilometers geographic coordinates for polar models. 

# Accessing Help 
For any of the TMD functions listed above, you can get text help in the Command Window by typing `help` followed by the name of the function. For example: 

```matlab
>> help tmd_predict
```
To access formatted documentation with lots of examples for any of the **Main functions** listed above, type `tmd` followed by the function name. For example:

```matlab
>> tmd tmd_predict
```
If you're not sure what function name you're looking for, just type `tmd` into the Command Window, and it will bring up a complete function list:

```matlab
>> tmd 
```

# Tutorials & More Documentation
Be sure to check out the documentation for each of the **Main functions**. There, you'll find multiple examples of how to use each function. In addition, the following tutorials and explanations may be of interest: 

* [TMD Getting Started](doc/tmd_getting_started.md)🚧
* [TMD Model file format](doc/TMD_model_file_format.md)
* [How to: Analyze ADCP current data](doc/tutorial_currents.md)
* [How to: Animate tidal motion](doc/tide_animation.md)🚧
* [How to: Calculate tidal range](doc/tidal_range_calculation.md)🚧

# What's new in TMD 3.0?
TMD3.0 is a rewrite of the [Tide Model Driver for MATLAB v2.5](https://github.com/EarthAndSpaceResearch/TMD_Matlab_Toolbox_v2.5). The goal of this update was to create a cleaner, more efficient set of functions that are better documented, easier to use, easier to debug. We also wanted to reduce the size, complexity, and overall hassle of dealing with multi-file or binary tide model data files. The biggest changes in this update to TMD3.0 are as follows: 

* Switched to a new, [consolidated NetCDF model data format](doc/TMD_model_file_format.md).  
* All functions rewritten for improved performance (*much faster* and lower memory usage than TMD2.5). 
* Most functions were renamed and inputs have changed to make them more intuitive and easier to use.
* All new documentation.
* Made updates to the CATS2008 model, as described [here](doc/cats2008_updates.md).

## Conversion to TMD3.0 NetCDF format

This version of TMD works with NetCDF data, which we have converted from the original binary OTIS format using scripts found in the `tide-model-conversions` folder. 

# Alternatives to MATLAB
Don't like MATLAB? No worries, just try one of these alternatives: 

* **PYTHON:** Tyler Sutterley's [pyTMD](https://github.com/tsutterley/pyTMD) is based on TMD2.5, and reads OTIS and GOT formatted tidal solutions for calculating ocean and load tides.
* **FORTRAN:** A Fortran version of this package is made available through OSU: [OSU Tidal Prediction software (OTPS)](https://www.tpxo.net/otps).
* **Octave** It is *possible* that TMD3.0 works with [Octave](https://www.gnu.org/software/octave/index) as it is currently written. I haven't checked, but most of the TMD functions are pretty simple. If you use Octave and you want to give TMD3.0 a try, please let me know how it goes! 


# Author Info & Citation Information
If you use TMD3.0, please cite the following: 

[Chad A. Greene](https://github.com/chadagreene), [Svetlana Erofeeva](https://ceoas.oregonstate.edu/svetlana-erofeeva), [Laurie Padman](https://github.com/LPadman), [Susan Howard](https://github.com/slhowardESR), [Tyler Sutterley](https://github.com/tsutterley), and [Gary Egbert](https://ceoas.oregonstate.edu/people/gary-egbert). The Tide Model Driver for MATLAB, version 3.0. 