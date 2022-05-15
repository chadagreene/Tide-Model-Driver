%% |tmd_interp| documentation
% tmd_interp interpolates tide model data at specified geographic locations. 
% 
%% Syntax 
% 
%  zi = tmd_interp(filename,variable,lati,loni)
%  zi = tmd_interp(...,'constituents',conList) 
%  zi = tmd_interp(...,'coasts',MaskingMethod)
% 
%% Syntax 
% 
% |zi = tmd_interp(filename,variable,lati,loni)| uses the NetCDF tide model
% specified by |filename| to interpolate a specified variable at the
% geographic coordinates |lati,loni|. The variable can be: 
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
%  * |'flexure'| ice shelf flexure coefficient from a linear elastic model applied to BedMachine ice thickness (can slightly exceed 1). 
% 
% |zi = tmd_interp(...,'constituents',conList)| specifies tidal constituents as a 
% cell array (e.g, |{'m2','s2'}|. If constituents are not specified, all constituents 
% from the model are returned. 
% 
% |zi = tmd_interp(...,'coasts',MaskingMethod)| specifies how coastal regions are masked. 
% Can be NaN, 'flexure', or 'unmask'. By default, |MaskingMethod| is |NaN|, meaning outputs 
% are set to |NaN| wherever a nearest-neighbor interpolation of the ocean indicates land. 
% The |'flexure'| option scales tidal constituents by a predicted coefficient of tidal 
% deflection for ice shelf grounding zones. A third option, |'unmask'|, does not apply 
% any masking, which may be preferred close to coasts, where, for example, a tide gauge 
% may exist between land and ocean grid cells. 

%% Example: Water column thickness 

% Create a grid (using the AMT function psgrid): 
[Lat,Lon] = psgrid('scar inlet',1500,1); 

% Get water column thickness at each grid point: 
%wct = tmd_interp('CATS2008_update_2022-04-22.nc','wct',Lat,Lon); 
wct = tmd_interp('/Users/cgreene/Downloads/TPXO9_atlas_v5/TPXO9_atlas30_update_2022-05-03.nc','wct',Lat,Lon); 

% Plot: 
figure
pcolorps(Lat,Lon,wct)
axis tight off
bedmachine   % plots coastline and gl
cmocean deep % colormap
shadem(10)   % hillshade with gain of 10
cb = colorbar; 
ylabel(cb,'water column thickness (m)')

%%


% 600 km wide grid at 0.2 km resolution, centered on (75S,100.5W).
[Lat,Lon] = psgrid(-75,-100.5,100,.2); 

h_default = tmd_interp('CATS2008_update_2022-04-22.nc','hAm',Lat,Lon,'constituents','m2');
h_flexure = tmd_interp('CATS2008_update_2022-04-22.nc','hAm',Lat,Lon,'constituents','m2','coasts','flexure');
h_unmask = tmd_interp('CATS2008_update_2022-04-22.nc','hAm',Lat,Lon,'constituents','m2','coasts','unmask');

figure
subsubplot(1,3,1) 
pcolorps(Lat,Lon,h_default)
axis tight off
bedmachine
caxis([0 0.06])
title 'NaN (default)'

subsubplot(1,3,2) 
pcolorps(Lat,Lon,h_flexure)
axis tight off
bedmachine
caxis([0 0.06])
title 'flexure'

subsubplot(1,3,3) 
pcolorps(Lat,Lon,h_unmask)
axis tight off
bedmachine
caxis([0 0.06])
title 'unmask'

%%

[Lat,Lon] = psgrid('ronne ice shelf',1500,.5); 

h_default = tmd_interp('CATS2008_update_2022-04-22.nc','hAm',Lat,Lon,'constituents','m2');

figure
pcolorps(Lat,Lon,h_default)
axis tight off
bedmachine

%%

%mask = tmd_interp('CATS2008_update_2022-04-11.nc','mask',Lat,Lon); 

h = tmd_interp('CATS2008_update_2022-04-22.nc','hAm',Lat,Lon,'constituent','k2'); 

hold on
contourps(Lat,Lon,h,'k')

%%

U2 = tmd_interp('CATS2008_update_2022-04-22.nc','uAm',Lat,Lon,'constituent','k2'); 
V2 = tmd_interp('CATS2008_update_2022-04-22.nc','vAm',Lat,Lon,'constituent','k2'); 

U(~isfinite(U)) = 0; 
V(~isfinite(V)) = 0; 

[Vx,Vy] = uv2vxvy(Lat,Lon,U,V); 


[X,Y] = ll2ps(Lat,Lon); 

sc = 0.02; 
Xr = imresize(X,sc); 
Yr = imresize(Y,sc); 
Vxr = imresize(Vx,sc); 
Vyr = imresize(Vy,sc); 


hold on
q = quiver(Xr,Yr,Vxr,Vyr,'r'); 
q.AutoScaleFactor = 3; 

%%

