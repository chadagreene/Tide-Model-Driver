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
### Main Function
Most users will only need to interact with one function. It is: 

* `tmd_predict` predicts tidal elevation, transport, or current velocities for given location(s) and time(s). 

### Other useful functions

* `tmd_interp` provides water column thickness, land/ocean mask, ice shelf flexure, and tidal constiuent parameters at specified geographic locations, for a given tide model. 
* `tmd_data` loads gridded tide model data without interpolation. 
* `tmd_ellipse` 
* `tidal_range` 

### Under-the-hood functions 
You probably won't need to call any of the following functions directly, but they are called by the functions above. 

* `tmd_astrol` computes the basic astronomical mean longitudes s, h, p, N.
* `tmd_constit` returns amplitude, phase, frequency, alpha, species for tidal constituents. 
* `tmd_harp` predicts tidal time series using harmonic constants. 
* `tmd_InferMinor` returns correction for 16 minor tidal constiuents. 
* `tmd_nodal`
* `tmd_ll2ps` converts geographic coordinates to polar stereographic kilometers for polar models. 
* `tmd_ps2ll` converts polar stereographic kilometers geographic coordinates for polar models. 

### Conversion from OTIS to NetCDF format

This version of TMD works with NetCDF data, which have been converted from the original OTIS format using the `tide-model-conversions` folder. 

# Installation 
1. Install TMD 
2. Get tide model data. 

# Tutorials
* [ADCP current data analysis](doc/tutorial_currents.md).

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