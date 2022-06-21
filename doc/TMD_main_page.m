%% Tide Model Driver 
% Tide Model Driver for MATLAB, version 3.0. This is a fast, user-friendly, well-documented rewrite of <https://www.esr.org/ ESR>'s classic Tide Model Driver by Svetlana Erofeeva.
%
%% Installation 
% To start predicting tides, just add TMD to MATLAB and download your favorite tide model data. Each step is described in the links below: 
% 
% # <https://github.com/chadagreene/Tide-Model-Driver/blob/main/doc/installing_tmd.md Add TMD to MATLAB.> 
% # <https://github.com/chadagreene/Tide-Model-Driver/blob/main/doc/tide_model_data.md Get tide model data.>
% 
%% TMD Functions
% For most applications, |tmd_predict| is the only function you will need to call directly, although you may occasionally want |tmd_interp| to retrieve bathymetry, |tmd_ellipse| to get tidal current characteristics, or |tmd_conlist| for a quick list of tidal constituents in a given model. Complete documentation is provided for the functions whose names are hyperlinked below, and simple text documentation can be found in the headers of all functions. 
% 
% * Main functions:* 
% * <tmd_predict_documentation.html *|tmd_predict|*> predicts tidal elevation, transport, or current velocities for given location(s) and time(s). 
% * <tmd_interp_documentation.html *|tmd_interp|*> provides water column thickness, land/ocean mask, ice shelf flexure, and tidal constituent parameters at specified geographic locations, for a given tide model. 
% * <tmd_data_documentation.html *|tmd_data|*> loads gridded tide model data without interpolation. 
% * <tmd_conlist_documentation.html *|tmd_conlist|*> returns a list of tidal constituents in a TMD3.0 compatible consolidated NetCDF tide model file. 
% * <tmd_ellipse_documentation.html *|tmd_ellipse|*> gives tidal ellipse parameters at specified location(s).  
% 
% *Under-the-hood functions:*
% 
% * |*tidal_range*| calculates the peak-to-peak tidal height range. 
% * |*tmd_astrol*| computes the basic astronomical mean longitudes s, h, p, N.
% * |*tmd_constit*| returns amplitude, phase, frequency, alpha, species for tidal constituents. 
% * |*tmd_harp*| predicts tidal time series using harmonic constants. 
% * |*tmd_InferMinor*| returns correction for 16 minor tidal constiuents. 
% * |*tmd_nodal*| calculates the nodal corrections for tidal constituents.
% * |*tmd_ll2ps*| converts geographic coordinates to polar stereographic kilometers for polar models. 
% * |*tmd_ps2ll*| converts polar stereographic kilometers geographic coordinates for polar models. 
% 
%% Accessing Help 
% For any of the TMD functions listed above, you can get text help in the Command Window by typing |help| followed by the name of the function. For example: 
% 
%  help tmd_predict
% 
%%
% To access formatted documentation with lots of examples for any of the *Main functions* listed above, type |tmd| followed by the function name. For example:
% 
%  tmd tmd_predict
% 
%%
% If you're not sure what function name you're looking for, just type |tmd| into the Command Window, and it will bring up a complete function list:
% 
%  tmd 
%
%% Tutorials & More Documentation
% Be sure to check out the documentation for each of the *Main functions*. There, you'll find multiple examples of how to use each function. In addition, the following tutorials and other documentation may be of interest: 
% 
% * TMD Getting Started
% * <TMD_model_file_format.html TMD Model file format> 
% * <tide_model_intercomparison.html Tide Model Intercomparison>
% * <tutorial_currents.html How to: Analyze ADCP current data>
% * <tutorial_tidal_range.html How to: Calculate tidal range>
% * <tmd_logo_animation.html How to: Animate Tide Data>

%% Author Info & Citation Information
% We are working on a proper literature citation. For now, if you use TMD3.0, please cite the following: 
%
% Chad A. Greene, Svetlana Erofeeva, Laurie Padman, Susan Howard, Tyler
% Sutterley, and Gary Egbert (2022). The Tide Model Driver for MATLAB, version 3.0. <https://github.com/chadagreene/Tide-Model-Driver>
