# 🚧 Tide Model Driver 3.0 🚧 
Tide Model Driver for MATLAB, version 3.0

🚧 This repo is under construction. 🚧

# TMD Contents 
### Help 
To access TMD documentation within MATLAB, simply type 

```matlab
tmd 
```
into the Command Window. To access documentation for a specific function, type `tmd` followed by the function name into the Command Window. For example:

```matlab
tmd tmd_predict
```
### Main Functions
For most applications, `tmd_predict` is the only function you will ever need to call directly, although in some cases you may need `tmd_interp` to retrieve bathymetry, or `tmd_ellipse` to tell you about tidal currents. Documentation is provided for all of the functions with hyperlinked names below. 

* `tmd_predict` predicts tidal elevation, transport, or current velocities for given location(s) and time(s). 
* `tmd_interp` provides water column thickness, land/ocean mask, ice shelf flexure, and tidal constiuent parameters at specified geographic locations, for a given tide model. 
* [`tmd_data`](doc/tmd_data_documentation.md) loads gridded tide model data without interpolation. 
* `tmd_ellipse` 
* `tidal_range` 
* `tmd_astrol` computes the basic astronomical mean longitudes s, h, p, N.
* `tmd_constit` returns amplitude, phase, frequency, alpha, species for tidal constituents. 
* `tmd_harp` predicts tidal time series using harmonic constants. 
* `tmd_InferMinor` returns correction for 16 minor tidal constiuents. 
* `tmd_nodal` calculates the nodal corrections for tidal constituents.
* `tmd_ll2ps` converts geographic coordinates to polar stereographic kilometers for polar models. 
* `tmd_ps2ll` converts polar stereographic kilometers geographic coordinates for polar models. 

### Conversion from OTIS to NetCDF format

This version of TMD works with NetCDF data, which have been converted from the original OTIS format using the `tide-model-conversions` folder. 

# Installation 
1. Install TMD 
2. Get tide model data. 

# Tutorials & Further Documentation
* [TMD Getting Started](doc/tmd_getting_started.md).🚧
* [TMD Model file format](doc/TMD_model_file_format.md).
* [ADCP current data analysis](doc/tutorial_currents.md).
* [Tide animation](doc/tide_animation.md).🚧
* [Tidal range calculation](doc/tidal_range_calculation.md).🚧

# What's new in TMD 3.0?

* Switched to consolidated NetCDF model data format. 
* TMD functions rewritten for improved performance. 
* Improved documentation.  
* CATS2008 updates: 
	* CATS resolution increased from 4 km to 2 km. 
	* CATS bathymetry and coastlines adjusted to match BedMachine v2 ([Morlighem et al., 2020](https://doi.org/10.1038/s41561-019-0510-8)). 
	* Ice shelf flexure model included for tidal deflections in grounding zones. 

![Tidal phase of the m2 constituent for the Filchner-Ronne Ice Shelf](tide-model-conversions/CATS2008_update_2022_comparison.png)

# Author Info & Citation Information
If you use TMD3.0, please cite the following: 

[Chad A. Greene](https://github.com/chadagreene), [Svetlana Erofeeva](https://ceoas.oregonstate.edu/svetlana-erofeeva), [Laurie Padman](https://github.com/LPadman), [Susan Howard](https://github.com/slhowardESR), [Tyler Sutterley](https://github.com/tsutterley), and [Gary Egbert](https://ceoas.oregonstate.edu/people/gary-egbert). The Tide Model Driver for MATLAB, version 3.0. 