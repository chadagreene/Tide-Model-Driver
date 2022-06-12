%% |tmd_data|
% |tmd_data| loads tide model data into the Matlab workspace. 
% 
%% Syntax 
% 
%  [Z,x,y] = tmd_data(filename,variable)
%  [...] = tmd_data(...,'constituents',conList)
%  [...] = tmd_data(...,'bounds',[xi yi])
%  [Z,lon,lat] = tmd_data(...,'geo')
%  [...,cons] = tmd_data(...)
%
%% Description 
% 
% |[Z,x,y] = tmd_data(filename,variable)| loads any of the following
% variables from a given tide model file: 
% 
%  * |'h'|   complex tidal height (m)  
%  * |'hRe'| real part of tidal height
%  * |'hIm'| imaginary part of tidal height 
%  * |'hAm'| amplitude of tidal height
%  * |'hPh'| phase of tidal height
%  * |'u'|   complex zonal velocity (m/s) 
%  * |'uRe'| real part of zonal velocity 
%  * |'uIm'| imaginary part of zonal velocity 
%  * |'uAm'| amplitude of zonal velocity
%  * |'uPh'| phase of zonal velocity
%  * |'U'|   complex zonal transport (m^2/s) 
%  * |'URe'|  real part of zonal transport
%  * |'UIm'| imaginary part of zonal transport
%  * |'UAm'| amplitude of zonal transport
%  * |'UPh'| phase of zonal transport 
%  * |'v'|   complex meridional velocity (m/s) 
%  * |'vRe'| real part of meridional velocity 
%  * |'vIm'| imaginary part of meridional velocity
%  * |'vAm'| amplitude of meridional velocity
%  * |'vPh'| phase of meridional velocity 
%  * |'V'|   complex meridional transport (m^2/s)
%  * |'VRe'| real part of meridional transport 
%  * |'VIm'| imaginary part of meridional transport
%  * |'VAm'| amplitude of meridional transport
%  * |'VPh'| phase of meridional transport
%  * |'wct'| water column thickness (m) 
%  * |'mask'| binary land/ocean mask
%  * |'flexure'| ice shelf flexure coefficient from a linear elastic model applied to BedMachine ice thickness (can slightly exceed 1). Only for CATS model. 
%
% |[...] = tmd_data(...,'constituents',conList)| specifies tidal constituents as a 
% cell array (e.g, {'m2','s2'}. If constituents are not specified, all model
% constituents are returned. 
%
% |[...] = tmd_data(...,'bounds',[xi yi])| enter an Mx2 matrix of coordinates 
% to only load data in a rectangle around the specified |[xi(:) yi(:)]| where
% |xi| and |yi| are projected coordinates for tide models that are in projected
% coordinates. For tide models in geo coordinates (global tide models,
% generally), |xi| refers to longitudes of interest and |yi| is latitudes of interest. 
% The advantage of specifying bounds is to minimize the amount of data that
% are loaded, when only a small area within a larger tide model is of
% interest. 
% 
% |[Z,lon,lat] = tmd_data(...,'geo')| returns grid coordinates as
% geographic coordinates. This is the default behavior for global models,
% whereas regional models return projected coordinates by default. 
%
% |[...,cons] = tmd_data(...)| returns a list of constituents in the model.  
% 
%% Example: Global Model water column thickness

[wct,lon,lat] = tmd_data('TPXO9_atlas_v5.nc','wct');

figure
imagesc(lon,lat,wct)
axis xy image 
caxis([0 7000]) 
cb = colorbar; 
ylabel(cb,'water column thickness (m)') 
xlabel('longitude')
ylabel('latitude') 
cmocean deep % optional colormap from Climate Data Toolbox

%% Example: Exploring tide model data
% Explore the model data for the height variable |h| in the CATS model.
% (Here I'm defining |fn| as the filepath to the version of CATS I'm
% working on at the moment, but you may have a different model filename.)

fn = '/Users/cgreene/Downloads/CATS2008/CATS2008_update_2022-06-11.nc';

[h,x,y,cons] = tmd_data(fn,'h'); 

whos h x y cons

%% 
% Above, we see that |h| is a complex variable whose dimensions correspond
% to |y|, |x|, and the constituents inn the model. 

%% Example: Water column thickness
% Using the same CATS filename |fn| we defined above, plot the water column
% thickness in the 

[wct,x,y] = tmd_data(fn,'wct'); 

figure
imagesc(x,y,wct)
axis xy image 
xlabel 'x (kilometers') 
ylabel 'y (kilometers')
cb = colorbar; 
ylabel(cb,'water column thickness (m)') 
cmocean deep % optional colormap 

%% Example
% You see above that the regional CATS model was generated in a nonstandard
% projection. You may prefer to get the geographic coordinates of the grid
% cell centers, like this: (Below I'm plotting with Antarctic Mapping
% Tools' |pcolorps| function): 

[wct,lon,lat] = tmd_data(fn,'wct','geo'); 

figure
pcolorps(lat,lon,wct) % requires AMT
bedmachine % plots Antarctica's coastline, if tyou have the bedmachine toolbox
cb = colorbar; 
ylabel(cb,'water column thickness (m)') 
cmocean deep % optional colormap 

%% Example 3

[hAm,x,y] = tmd_data(fn,'hAm','constituents','s2');  
hPh = tmd_data(fn,'hPh','constituents','s2'); 

figure
imagesc(x,y,hAm); 
axis xy image
hold on
contour(x,y,rad2deg(hPh),-180:30:180,'k')
xlabel 'easting (km)' 
ylabel 'northing (km)' 
caxis([0 1.2])

%%

ocean = tmd_data(fn,'mask'); 

% Set land values to NaN: 
hAm(~ocean) = nan; 
hPh(~ocean) = nan; 

figure
h=imagesc(x,y,hAm); 
h.AlphaData = ocean; % sets transparency
axis xy image
hold on
contour(x,y,rad2deg(hPh),-180:30:180,'k')
xlabel 'easting (km)' 
ylabel 'northing (km)' 
caxis([0 1.2])

%% Example: A global model
% Working with global models is very similar to working with regional
% models. The only real difference is we let lon=x and lat=y. Here's an
% example: 

[M2,lon,lat] = tmd_data


%%

[u1,x1,y1] = tmd_data(fn,'U','bounds',[-2498 3049]);

%%

