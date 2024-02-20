[&larr; Back to TMD3.0 Main Page](../README.md)

# What's new in TMD3.0?
TMD3.0 is a rewrite of the [Tide Model Driver for MATLAB v2.5](https://github.com/EarthAndSpaceResearch/TMD_Matlab_Toolbox_v2.5). The goal of this update was to create a clean, efficient set of functions that are well documented, easy to use, and easy to debug. We also wanted to reduce the size, complexity, and overall hassle of dealing with multi-file or binary tide model data files. The biggest changes in this update to TMD3.0 are as follows: 

* Introduced a new, [consolidated NetCDF model data format](TMD_model_file_format.md).  
* All functions rewritten for improved performance (*much faster* and lower memory usage than TMD2.5). 
* Most functions were renamed and inputs have changed to make them more intuitive and easier to use.
* All new documentation.
* Updated the CATS2008 model, as described [here](cats2008_updates.md).

# Converting from previous TMD releases
Here's a table of how the new functions compare to previous versions of TMD: 

| TMD3.0 function | Previous name | Major Changes                      |
|:---------------:|:-----------:|:------------------------------------:|
|`tmd_predict`  |`tmd_tide_pred   `| Changed the order of function inputs|
|`tmd_interp`   | `tmd_extract_HC`| Changed function inputs & outputs. Phase now in radians.|
|`tmd_data`     | `h_in`, `grd_in`| n/a |
|`tmd_conlist`  |No direct previous equivalent| n/a |
|`tmd_ellipse`|`tmd_get_ellipse`| Now allows input geo coordinates |
|`tidal_range`  |No direct previous equivalent| n/a |
|`tidal_astrol`  |`astrol`| functionally equivalent  |
|`tmd_constit`  |`constit`| functionally equivalent   |
|`tmd_harp`  |`harp`, `harp1`| changed to `datenum` input format |
|`tmd_InferMinor`  |`InferMinor`| functionally equivalent  |
|`tmd_nodal`  |`nodal`, `nodal1`| now requires both p and N as inputs |
|`tmd_ll2ps`  |`mapll`| functionally equivalent  |
|`tmd_ps2ll`  |`mapxy`| functionally equivalent |

#  Comparing TMD3.0 to previous versions 
If you want to compare results of TMD3.0 to previous versions of TMD, feel free to use or modify [this script](../tide-model-conversions/testing/tmd30_vs_tmd24.m) we wrote for that exact purpose. 

# Author Info
This page was written by [Chad A. Greene](https://www.chadagreene.com), December 2023. 
