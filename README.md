# Tide Model Driver 3.0
![TMD logo](doc/markdown_figures/tmd_logo.gif)
Tide Model Driver for MATLAB, version 3.0. This is a fast, user-friendly, well-documented rewrite of [ESR](https://www.esr.org/)'s classic Tide Model Driver by [Svetlana Erofeeva](https://ceoas.oregonstate.edu/svetlana-erofeeva).

# Installation 
To start predicting tides, just add TMD to MATLAB and download your favorite tide model data. Each step is described in the links below: 

1. [Add TMD to MATLAB.](doc/installing_tmd.md)
2. [Get tide model data.](doc/tide_model_data.md)

# TMD Functions
For most applications, **`tmd_predict` is the only function you will need to call directly**, although you may occasionally want `tmd_interp` to retrieve bathymetry, `tmd_ellipse` to get tidal current characteristics, or `tmd_conlist` for a quick list of tidal constituents in a given model. Complete documentation is provided for the functions whose names are hyperlinked below, and simple text documentation can be found in the headers of all functions. 

#### Main functions:
* [`tmd_predict`](doc/tmd_predict_documentation.md) predicts tidal elevation, transport, or current velocities for given location(s) and time(s). 
* [`tmd_interp`](doc/tmd_interp_documentation.md) provides water column thickness, land/ocean mask, ice shelf flexure, and tidal constituent parameters at specified geographic locations, for a given tide model. 
* [`tmd_data`](doc/tmd_data_documentation.md) loads gridded tide model data without interpolation. 
* [`tmd_conlist`](doc/tmd_conlist_documentation.md) returns a list of tidal constituents in a TMD3.0 compatible consolidated NetCDF tide model file. 
* [`tmd_ellipse`](doc/tmd_ellipse_documentation.md) gives tidal ellipse parameters at specified location(s).  

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
Be sure to check out the documentation for each of the **Main functions**. There, you'll find multiple examples of how to use each function. In addition, the following tutorials and other documentation may be of interest: 

* TMD Getting Started🚧
* [TMD Model file format](doc/TMD_model_file_format.md)
* [Tide Model Intercomparison](doc/tide_model_intercomparison.md)
* [How to: Analyze ADCP current data](doc/tutorial_currents.md)
* [How to: Calculate tidal range](doc/tutorial_tidal_range.md)
* [How to: Animate tidal motion](doc/tmd_logo_animation.md)

# What's new in TMD 3.0?
TMD3.0 is a rewrite of the [Tide Model Driver for MATLAB v2.5](https://github.com/EarthAndSpaceResearch/TMD_Matlab_Toolbox_v2.5). The goal of this update was to create a clean, efficient set of functions that are well documented, easy to use, easy to debug. We also wanted to reduce the size, complexity, and overall hassle of dealing with multi-file or binary tide model data files. The biggest changes in this update to TMD3.0 are as follows: 

* Introduced a new, [consolidated NetCDF model data format](doc/TMD_model_file_format.md).  
* All functions rewritten for improved performance (*much faster* and lower memory usage than TMD2.5). 
* Most functions were renamed and inputs have changed to make them more intuitive and easier to use.
* All new documentation.
* Updated the CATS2008 model, as described [here](doc/cats2008_updates.md).

# Alternatives to MATLAB
Don't like MATLAB? No worries, just try one of these alternatives: 

* **Python:** [pyTMD](https://github.com/tsutterley/pyTMD) is a tidal prediction software that can read OTIS, GOT and FES formatted tidal solutions for predicting ocean and load tides.
* **Fortran:** A Fortran version of this package is made available through OSU: [OSU Tidal Prediction software (OTPS)](https://www.tpxo.net/otps).
* **Octave:** It is *possible* that TMD3.0 works with [Octave](https://www.gnu.org/software/octave/index) as it is currently written. I haven't checked, but most of the TMD functions are pretty basic, so Octave's functions might be exact clones. If you use Octave and you want to give TMD3.0 a try, please let me know how it goes! 

# Author Info & Citation Information
We are working on a proper literature citation. For now, if you use TMD3.0, please cite the following: 

[Chad A. Greene](https://github.com/chadagreene), [Svetlana Erofeeva](https://ceoas.oregonstate.edu/svetlana-erofeeva), [Laurie Padman](https://github.com/LPadman), [Susan Howard](https://github.com/slhowardESR), [Tyler Sutterley](https://github.com/tsutterley), and [Gary Egbert](https://ceoas.oregonstate.edu/people/gary-egbert) (2022). The Tide Model Driver for MATLAB, version 3.0. [https://github.com/chadagreene/Tide-Model-Driver](https://github.com/chadagreene/Tide-Model-Driver)
