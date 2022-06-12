%% TMD3.0 Model File Format
% This page describes the NetCDF data format for TMD3.0 compatible tide
% model files. 
%
% <TMD_Contents.html Back to Tide Model Driver Contents>.

%% Explore file contents 
% Take a look inside any TMD3.0 compatible NetCDF, and the contents should
% be fairly straightforward. Here's an example: 

ncdisp('TPXO9_atlas_v5.nc')

%% Global models vs regional models 
% All model files are either _global_ or _regional_. Global models, like the
% one above, are presented on a regular grid in geographic coordinates (degrees) 
% and have 1D |lon| and |lat| variables. For global models, longitude is
% treated like the x dimension, and latitude is treated like the y
% dimension.
%
% Regional models (CATS, Arctic models, etc) are presented ongrids that are 
% regular in _projected_ coordinates, meaning they are spaced equally in
% meters or kilometers. The proj4 string in the NetCDF file describes the
% projection parameters for any regional model. Regional models also
% include MxN arrays of latitude and longitude values. 
%
% In both global and regional models, all coordinates correspond to grid
% cell centers. 
% 
% Most TMD functions accept latitude and longitude as inputs, and
% automatically project them to regional model coordinates where necessary.

%% Constituents 
% The tidal constituents are easy to miss in the NetCDF. They appear as a
% simple string, in the |constituent_order| attribute of the |constituents|
% variable. Here's how to access them: 

ncreadatt('TPXO9_atlas_v5.nc','constituents','constituent_order')

%%
% Constituents can also be accessed as cell arrays using the |tmd_data|
% function: 

[wct,lon,lat,cons] = tmd_data('TPXO9_atlas_v5.nc','wct'); 

cons

%% 
% I realize this is a somewhat strange way to package the constituent
% names, but it's the simplest way I could figure out how to do it, given
% the limitations of the NetCDF format. 

%% Real and Imaginary Components 
% For tidal height h, zonal transport U, and meridional transport V, you'll
% find real and imaginary components hRe, hIm, URe, UIm, VRe, and VIm in
% the model file. The complex constituent can then be constucted following
% the format 
% 
%  hc = complex(hRe,hIm); 
% 
% or 
% 
%  hc = hRe + 1i*hIm; 
% 
% and the amplitude is then 
% 
%  hAm = abs(hc); 
% 
% and the phase is given by 
% 
%  hPh = angle(hc); 
% 
% Note: Some conventions use the complex conjugate of |hc|, and calculate
% phase angle as |atan2(-imag(h), real(h))|, but TMD3.0 uses the convention
% that's compatible with MATLAB's built-in |complex| and |angle| functions.

%% U and V variables
% Each tide model file contains real and imaginary components of U and V
% transport variables, whose units are m^2/s. 
%
% Regardless of whether the model is global or regional, U and V always
% correspond to zonal (positive pointing geographic east) and meridional
% (positive pointing geographic north) components of transport. 
% 
% Transport estimates are usually pretty good, and represent column-averaged 
% (barotropic) flow. Column-averaged velocities are obtained by dividing
% transport by water column thickness, so the accuracy of velocity
% estimates is subject to the accuracy of the bathymetry estimate. 
% 
% TMD3.0 model file format does not account for vertical variations in
% tidal transport or current velocity. In other words, this is a strictly
% barotropic model, and does not attempt to represent baroclinic flow. 

%% Units 
% TMD3.0 uses meters, seconds, and combinations of meters and seconds. Previous 
% versions of TMD produced tidal currents in cm/s, but we now package
% transports in units of m^2/s, so they may easily be divided by water
% column thickness (m), to get velocities of m/s. 

%% Masks 
% TMD does something new with masking: The scripts that create each model
% file use |regionfill| to interpolate values of tidal constituents
% across all land areas. 
% 
% Here's what the m2 constituent amplitude looks like in the TPXO9 file:

[ph,lon,lat] = tmd_data('TPXO9_atlas_v5.nc','hAm','constituents','m2'); 

figure
h=imagesc(lon,lat,ph);
axis xy image 
caxis([0 2])
cb = colorbar; 
ylabel(cb,'m2 constituent amplitude (m)') 

%% 
% _Why on Earth would we say there's a finite, nonzero tidal amplitude in
% the middle of China, or any other landmass?_ Well, the reason's quite simple: 
% Sometimes you may be interested in tides that are close to shore, for example, 
% if you have a tide gauge at the end of a dock that may lie between
% modeled ocean and land pixels. In such a case, you know that tides exist
% there, but the tide model would only produce NaNs. You could extrapolate
% from the nearest grid cell that produces a finite tide solution, but that'd 
% be inelegant and almost certainly wrong. The |regionfill| approach
% produces smoother, more physical interpolation close to the coast and is
% a reasonable approximation where the underlying tide model cannot offer a
% solution. 
% 
% The <tmd_interp_documentation.html |tmd_interp|> function sets land areas 
% to NaN by default (but offers an option to "unmask"). For your tidal needs, 
% you may wish to mask out land areas using the |mask| variable, which is
% 1 for all ocean grid cells, and 0 everywhe else. Below I'm using the mask
% to set the transparency of the previous image: 

[mask,lon,lat] = tmd_data('TPXO9_atlas_v5.nc','mask'); 
h.AlphaData = mask; 

%% Amplitude, phase, alpha, and omega 
% In each TMD3.0 compatible model file, you will find variables called
% amplitude, phase, omega, and alpha. You probably won't need to interact
% with these variables directly, but they are from the |tmd_constit|
% function. 

%% |scale_factor| 
% The real and imaginary components of h, U, and V are all scaled by a
% strange number listed as the |scale_factor| of each variable. MATLAB
% automatically parses the |scale_factor| when reading in the variables, so
% you don't have to worry about multiplying anything by anything. 
% 
% _Why isn't the scale factor a simple round number?_ Because it doesn't
% have to be. Again, you don't have to worry about the |scale_factor|. The
% value of the scale factor is chosen to take advantage of the full range
% of the int16 Datatype that it's saved as. Any other scale factor would
% either clip the large values, or digitize the data at lower precision
% than necessary. 

%% Extra columns in global models 
% You may notice that global model files have an extra grid cell on each end,
% corresponding to an extra longitude step before 0 and after 360. This
% rationale is that interpolation near 0 degrees longitude would produce
% NaNs, and you'd end up with a seam of missing data at that lontitude. 
% 
% The trick of repeating a few rows or columns around the 0 degree longitude
% is an old one, and it could be performed after loading the data, rather than 
% saving the data. However, |tmd_interp| only loads pixels around query
% points, so the simplest solution is to repeat the data in the data file,
% rather than writing |tmd_interp| to call |tmd_data| multiple times and
% stitch the pieces together before interpolation. 

%% Flexure 
% For the CATS model, we have added a |flexure| variable. This is a
% first-order estimate of tidal flexure of ice shelves in the grounding
% zone, generated by appling a simple 1D linear elastic model to BedMachine
% v2 ice thickness, assuming an elastic modulus of E=4.8 GPa and Poisson's ratio
% nu = 0.4. Values range from zero to about 100, corresponding to the
% percent of tidal range that the ice shelf should exhibit. Values can
% exceed 100 by a few percent near the hydrostatic line. 

%% Author Info
% This document was written by Chad A. Greene, June 2022. 