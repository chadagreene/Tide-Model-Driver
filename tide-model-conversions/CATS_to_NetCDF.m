% This file reads the old binary CATS2008 tide model and converts it to NetCDF. 
% This script also densifies the data to 2 km resolution and adjusts coastlines 
% to match BedMachine v2. 
% 
% This script calls some legacy functions from TMD2.5 to load the old data.
% 
% Written by Chad A. Greene, NASA/JPL, February 2022. 

%% Enter initial file info 

% Output file: 
newfilename = ['/Users/cgreene/Documents/data/tides/CATS2008_update_',datestr(now,'yyyy-mm-dd'),'.nc']; 

% Input files: 
filename_grd = '/Users/cgreene/Documents/data/tides/CATS2008/grid_CATS2008'; 
filename_h = '/Users/cgreene/Documents/data/tides/CATS2008/hf.CATS2008.out'; 
filename_u = '/Users/cgreene/Documents/data/tides/CATS2008/uv.CATS2008.out'; 
addpath(genpath('/Users/cgreene/Downloads/TMD3.00_alpha'))
addpath(genpath('/Users/cgreene/Documents/data/tides/CATS2008'))

old_res = 4; % (km) input grid resolution 
new_res = 2; % (km) output grid resolution 

con_string = 'm2 s2 n2 k2 k1 o1 p1 q1 mf mm'; % constituents in original model order

%% Load data

% Read the grid file: 
[ll_lims,wct,mask,~,~] = grd_in(filename_grd); 

% Reorient and reduce data size:
mask = uint8(flipud(mask')); 
wct = uint16(flipud(wct')); 

% Create spatial arrays: 
x = (ll_lims(1)+old_res/2):old_res:(ll_lims(2)-old_res/2); 
y = (ll_lims(4)-old_res/2):-old_res:(ll_lims(3)+old_res/2); 

% Number of constituents: 
conList = strsplit(con_string,' ');
Ncons = length(conList); 

for k=1:Ncons
   
   % Read each tidal height constituent and reorient: 
   tmp = h_in(filename_h,k); 
   h(:,:,k) = flipud(tmp'); 
   
   % Read each transport constituent and reorient: 
   [tmp,tmp2] = u_in(filename_u,k); 
   U(:,:,k) = flipud(tmp'); 
   V(:,:,k) = flipud(tmp2'); 
end

% Mask out erroneous "lakes" in CATS2008 (identified by Tyler Sutterley)  
bad = rms(abs(h),3)<0.003 & mask;
mask(bad) = 0; 

%% Start a phase comparison 
% (Finishes after grid densification) 

oldaz = angle(h(:,:,1)); 
oldaz(~mask) = nan; 

[X,Y] = meshgrid(x,y); 
[Lat,Lon] = mapxy(X,Y,-71,-70,'S'); 

figure
subsubplot(1,2,1) 
pcolorps(Lat,Lon,oldaz); 
axis tight off
bedmachine
cmocean phase 
title 'CATS2008'
ax(1) = gca; 

%% Fill NaNs

for k = 1:Ncons
   tmp = h(:,:,k); 
   tmp = complex(regionfill(real(tmp),mask==0),regionfill(imag(tmp),mask==0));
   h(:,:,k) = tmp; 
   
   tmp = U(:,:,k); 
   tmp = complex(regionfill(real(tmp),mask==0),regionfill(imag(tmp),mask==0));
   U(:,:,k) = tmp; 
   
   tmp = V(:,:,k); 
   tmp = complex(regionfill(real(tmp),mask==0),regionfill(imag(tmp),mask==0));
   V(:,:,k) = tmp; 
   
   k
end

%% Densify the grid 

% The old 4 km resolution data: 
x4 = x; 
y4 = y; 
Z4 = h; 
U4 = U; 
V4 = V; 
clear h U V

% A new resolution grid: 
x = x4(1):new_res:x4(end);
y = y4(1):-new_res:y4(end);
[X,Y] = meshgrid(x,y);

for k = 1:Ncons
   h(:,:,k) = interp2(x4,y4,Z4(:,:,k),X,Y,'spline');
   U(:,:,k) = interp2(x4-old_res/2,y4,U4(:,:,k),X,Y,'spline');
   V(:,:,k) = interp2(x4,y4-old_res/2,V4(:,:,k),X,Y,'spline');
end

oldmask = interp2(x4,y4,mask,X,Y,'nearest'); 
wct_cats = interp2(x4,y4,double(wct),X,Y); 

%% Get the GL to agree with BedMachine

[Lat,Lon] = mapxy(X,Y,-71,-70,'S'); 
Lat(Lat<-90) = -90; % fixes some rounding error

% Load the coastline for context 
[gllat,gllon] = bedmachine_data('coast','geo');

mask_bm = bedmachine_interp('mask',Lat,Lon); 
ground = imfill(ismember(mask_bm,[1 2 4]),'holes'); 

[masktmp,x_bm,y_bm] = bedmachine_data('mask'); 
thickness = bedmachine_data('thickness'); 
dst2ground = double(bwdist(ismember(masktmp,[1 2 4]))) * diff(x_bm(1:2));

% Starting with E=4.8 GPa and nu = 0.4 from Reeh and others (2003) 
rho_sw = 1027; 
gravity = 9.81; 
nu = 0.4; 

beta = (3*rho_sw * gravity*((1-nu^2)./(4.8e9.*thickness.^3))).^(1/4);
flextmp = 1 - exp(-beta.*dst2ground) .* (cos(beta.*dst2ground)+sin(beta.*dst2ground));

flextmp(masktmp==4) = 0; 
flextmp(masktmp==0) = 1; 
flextmp = filt2(flextmp,diff(x_bm(1:2)),new_res*2000,'lp'); % lowpass filter to 4 km before interpolating to 2 km grid 
[Xtmp,Ytmp] = ll2ps(Lat,Lon); 
flexure = interp2(x_bm,y_bm,flextmp,Xtmp,Ytmp); 
flexure(isnan(flexure))=1; 
flexure(oldmask==0 & Lat>max(gllat)) = 0; 
flexure(ground) = 0; 
flexure = uint8(flexure*100); % converts to percent deflection.

wct_bm = filt2(bedmachine_data('wct'),diff(x_bm(1:2)),new_res*2000,'lp'); 
wct_bm = interp2(x_bm,y_bm,wct_bm,Xtmp,Ytmp); 

decay_length = 100; % pixels
weight_cats = exp(-bwdist(isnan(wct_bm))/decay_length); 
weight_bm = 1-weight_cats; 

wct_bm(isnan(wct_bm)) = 0; % just so we can sum without nans bugging us 
wct = (wct_bm.*weight_bm + wct_cats.*weight_cats); % this creates a gentle blending of the two datasets. 

mask = uint8(~ground); 
mask(oldmask==0 & Lat>max(gllat)) = 0; 

wct(mask==0) = 0; 
wct(wct<1 & mask==1) = 1; % just in case of any apparent grounded but not grounded ice 
wct = uint16(wct); 

%%

tmp = angle(h(:,:,1)); 
tmp(~mask) = nan; 

subsubplot(1,2,2) 
hp = pcolorps(Lat,Lon,tmp); 
axis tight off
bedmachine
cmocean phase 
ax(2) = gca; 
title('CATS2008_update_2022','interpreter','none')

linkaxes(ax,'xy'); 
axis([  -1609006.64    -415618.02      63826.14    1152447.96])

% export_fig CATS2008_update_2022_comparison.png -r600 -p0.01

%% Tidal range
% After writing the whole file, calculate the peak-to-peak tidal range and 
% add it to the NetCDF. 

% Takes ~5 hours! 
%h_range = tidal_range(h,conList,mask==1);
% h_range = permute(ncread('CATS2008_update_2022-06-06.nc','h_range'),[2 1]); 

%%

% Scaling factor for saving to NCSHORT (int16):
scale_h = floor(32767/max(abs([real(h(:));imag(h(:))])));
scale_UV = floor(32767/max(abs([real(U(:));real(V(:));imag(U(:));imag(V(:));])));

[ispec,amp,ph,omega,alpha] = tmd_constit(strsplit(con_string));

proj4 = '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=-70 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +type=crs';

%% Write the netcdf 

% 1. Create netCDF file handle and Attributes
mode = netcdf.getConstant('NETCDF4');
mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
ncid=netcdf.create(newfilename,mode);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Conventions','CF-1.7');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Title','CATS2022');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Description','This is CATS2008, but resampled at 2 km resolution, and adjusted to match BedMachine v2 ice masks and water column thickness.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Author','Susan Howard, Svetlana Erofeeva, Laurie Padman, Chad A. Greene, and Tyler Sutterley.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'creation_date',datestr(now,'yyyy-mm-dd'));
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'NetCDF_conversion','Chad A. Greene');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'tmd_version',3.0);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'model_type','ocean');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'license','MIT License');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Data_citation',['Howard, S. L., Erofeeva, S., & Padman, L. (2019) "CATS2008: Circum-Antarctic Tidal Simulation version 2008" u.S. Antarctic Program (USAP) Data Center. doi: https://doi.org/10.15784/601235.'])

% 2. Define dimensions
% Define mapping variable
mapping_var_id= netcdf.defVar(ncid,'mapping','NC_CHAR',[]);
netcdf.putAtt(ncid,mapping_var_id,'grid_mapping_name', 'polar_stereographic');
netcdf.putAtt(ncid,mapping_var_id,'latitude_of_projection_origin',-90.);
netcdf.putAtt(ncid,mapping_var_id,'standard_parallel',-71.);
netcdf.putAtt(ncid,mapping_var_id,'straight_vertical_longitude_from_pole',-70.);
netcdf.putAtt(ncid,mapping_var_id,'false_easting',0.);
netcdf.putAtt(ncid,mapping_var_id,'false_northing',0.); 
netcdf.putAtt(ncid,mapping_var_id,'spatial_proj4',proj4); 

% Define x: 
x_id     = netcdf.defDim(ncid,'x',length(x));
x_var_id = netcdf.defVar(ncid,'x','NC_FLOAT',x_id);
netcdf.putAtt(ncid,x_var_id,'long_name',    'Cartesian x-coordinate, grid cell center');
netcdf.putAtt(ncid,x_var_id,'standard_name','projection_x_coordinate');
netcdf.putAtt(ncid,x_var_id,'units',        'kilometer');

% Define y: 
y_id     = netcdf.defDim(ncid,'y',length(y));
y_var_id = netcdf.defVar(ncid,'y','NC_FLOAT',y_id);
netcdf.putAtt(ncid,y_var_id,'long_name',    'Cartesian y-coordinate, grid cell center');
netcdf.putAtt(ncid,y_var_id,'standard_name','projection_y_coordinate');
netcdf.putAtt(ncid,y_var_id,'units',        'kilometer');

% Define lat
lat_var_id = netcdf.defVar(ncid,'lat','NC_DOUBLE',[x_id y_id]);
netcdf.putAtt(ncid,lat_var_id,'long_name',    'grid cell center latitude');
netcdf.putAtt(ncid,lat_var_id,'standard_name','latitude');
netcdf.putAtt(ncid,lat_var_id,'units',        'degree');

% Define lon
lon_var_id = netcdf.defVar(ncid,'lon','NC_DOUBLE',[x_id y_id]);
netcdf.putAtt(ncid,lon_var_id,'long_name',    'grid cell center longitude');
netcdf.putAtt(ncid,lon_var_id,'standard_name','longitude');
netcdf.putAtt(ncid,lon_var_id,'units',        'degree');

% Define constituents
cons_id = netcdf.defDim(ncid,'constituents',Ncons); 
cons_var_id = netcdf.defVar(ncid,'constituents','NC_BYTE',cons_id); 
netcdf.putAtt(ncid,cons_var_id,'standard_name', 'tidal_constituents');
netcdf.putAtt(ncid,cons_var_id,'long_name','Tidal constituents listed in order in the constituent_order attribute.'); 
netcdf.putAtt(ncid,cons_var_id,'constituent_order',con_string); 

% Define constituent attributes 
amp_var_id = netcdf.defVar(ncid,'amplitude','NC_DOUBLE',cons_id); 
netcdf.putAtt(ncid,amp_var_id,'standard_name',    'amplitude');
netcdf.putAtt(ncid,amp_var_id,'long_name','amplitude of equilibrium tide in for each tidal constituent.'); 
netcdf.putAtt(ncid,amp_var_id,'units',    'meters');

ph_var_id = netcdf.defVar(ncid,'phase','NC_DOUBLE',cons_id); 
netcdf.putAtt(ncid,ph_var_id,'standard_name',    'phase');
netcdf.putAtt(ncid,ph_var_id,'long_name','Astronomical arguments (relative to t0 = 1 Jan 0:00 1992)'); 
netcdf.putAtt(ncid,ph_var_id,'units',    'radians');

om_var_id = netcdf.defVar(ncid,'omega','NC_FLOAT',cons_id); 
netcdf.putAtt(ncid,om_var_id,'standard_name',    'omega');
netcdf.putAtt(ncid,om_var_id,'long_name','frequency'); 
netcdf.putAtt(ncid,om_var_id,'units','1/s'); 

alp_var_id = netcdf.defVar(ncid,'alpha','NC_FLOAT',cons_id); 
netcdf.putAtt(ncid,alp_var_id,'standard_name',    'alpha');
netcdf.putAtt(ncid,alp_var_id,'long_name','loading love number'); 

% Define hRe
hRe_var_id = netcdf.defVar(ncid,'hRe','NC_SHORT',[x_id y_id cons_id]);
netcdf.putAtt(ncid,hRe_var_id,'long_name',    'real component of height constituent');
netcdf.putAtt(ncid,hRe_var_id,'standard_name','height_constituent');
netcdf.putAtt(ncid,hRe_var_id,'grid_mapping', 'polar_stereographic');
netcdf.putAtt(ncid,hRe_var_id,'units',        'm');
netcdf.putAtt(ncid,hRe_var_id,'scale_factor',  1/scale_h);

% Define hIm
hIm_var_id = netcdf.defVar(ncid,'hIm','NC_SHORT',[x_id y_id cons_id]);
netcdf.putAtt(ncid,hIm_var_id,'long_name',    'imaginary component of height constituent');
netcdf.putAtt(ncid,hIm_var_id,'standard_name','height_constituent');
netcdf.putAtt(ncid,hIm_var_id,'grid_mapping', 'polar_stereographic');
netcdf.putAtt(ncid,hIm_var_id,'units',        'm');
netcdf.putAtt(ncid,hIm_var_id,'scale_factor',  1/scale_h);

% Define uRe
uRe_var_id = netcdf.defVar(ncid,'URe','NC_SHORT',[x_id y_id cons_id]);
netcdf.putAtt(ncid,uRe_var_id,'long_name',    'real component of U transport constituent. This is the zonal flow component in geographic coordinates.');
netcdf.putAtt(ncid,uRe_var_id,'standard_name','height_constituent');
netcdf.putAtt(ncid,uRe_var_id,'grid_mapping', 'polar_stereographic');
netcdf.putAtt(ncid,uRe_var_id,'units',        'm^2/s');
netcdf.putAtt(ncid,uRe_var_id,'scale_factor',  1/scale_UV);

% Define uIm
uIm_var_id = netcdf.defVar(ncid,'UIm','NC_SHORT',[x_id y_id cons_id]);
netcdf.putAtt(ncid,uIm_var_id,'long_name',    'imaginary component of U transport constituent. This is the zonal flow component in geographic coordinates.');
netcdf.putAtt(ncid,uIm_var_id,'standard_name','height_constituent');
netcdf.putAtt(ncid,uIm_var_id,'grid_mapping', 'polar_stereographic');
netcdf.putAtt(ncid,uIm_var_id,'units',        'm^2/s');
netcdf.putAtt(ncid,uIm_var_id,'scale_factor',  1/scale_UV);

% Define vRe
vRe_var_id = netcdf.defVar(ncid,'VRe','NC_SHORT',[x_id y_id cons_id]);
netcdf.putAtt(ncid,vRe_var_id,'long_name',    'real component of V transport constituent. This is the meridional flow component in geographic coordinates.');
netcdf.putAtt(ncid,vRe_var_id,'standard_name','height_constituent');
netcdf.putAtt(ncid,vRe_var_id,'grid_mapping', 'polar_stereographic');
netcdf.putAtt(ncid,vRe_var_id,'units',        'm^2/s');
netcdf.putAtt(ncid,vRe_var_id,'scale_factor',  1/scale_UV);

% Define vIm
vIm_var_id = netcdf.defVar(ncid,'VIm','NC_SHORT',[x_id y_id cons_id]);
netcdf.putAtt(ncid,vIm_var_id,'long_name',    'imaginary component of V transport constituent. This is the meridional flow component in geographic coordinates.');
netcdf.putAtt(ncid,vIm_var_id,'standard_name','height_constituent');
netcdf.putAtt(ncid,vIm_var_id,'grid_mapping', 'polar_stereographic');
netcdf.putAtt(ncid,vIm_var_id,'units',        'm^2/s');
netcdf.putAtt(ncid,vIm_var_id,'scale_factor',  1/scale_UV);

% Define wct: 
wct_var_id = netcdf.defVar(ncid,'wct','NC_SHORT',[x_id y_id]);
netcdf.putAtt(ncid,wct_var_id,'long_name','water column thickness');
netcdf.putAtt(ncid,wct_var_id,'standard_name','wct');
netcdf.putAtt(ncid,wct_var_id,'units',    'meters');
netcdf.putAtt(ncid,wct_var_id,'grid_mapping', 'polar_stereographic');

% Define mask
mask_var_id = netcdf.defVar(ncid,'mask','NC_BYTE',[x_id y_id]);
netcdf.putAtt(ncid,mask_var_id,'long_name',    'ocean mask');
netcdf.putAtt(ncid,mask_var_id,'standard_name','ocean_mask');
netcdf.putAtt(ncid,mask_var_id,'grid_mapping', 'polar_stereographic');
netcdf.putAtt(ncid,mask_var_id,'valid_range',  [0 1]);
netcdf.putAtt(ncid,mask_var_id,'flag_values',  [0 1]);
netcdf.putAtt(ncid,mask_var_id,'flag_meanings','land ocean');

% Define flexure
flexure_var_id = netcdf.defVar(ncid,'flexure','NC_BYTE',[x_id y_id]);
netcdf.putAtt(ncid,flexure_var_id,'long_name',    'Forward-modeled coefficient of tidal flexure assuming linear elastic response applied to BedMachine v2 geometry with rho_sw=1027 kg/m3, poisson=0.4, E=4.8 GPa. Can exceed 100% (by a few percent) near the hydrostatic line.');
netcdf.putAtt(ncid,flexure_var_id,'standard_name','ice_flexure_percent');
netcdf.putAtt(ncid,flexure_var_id,'grid_mapping', 'polar_stereographic');
% 
% % define h_range
% R_var_id = netcdf.defVar(ncid,'h_range','NC_SHORT',[x_id y_id]);
% netcdf.putAtt(ncid,R_var_id,'long_name',    'Peak-to-peak tidal range.');
% netcdf.putAtt(ncid,R_var_id,'standard_name','tidal range');
% netcdf.putAtt(ncid,R_var_id,'grid_mapping', 'polar_stereographic');
% netcdf.putAtt(ncid,R_var_id,'units', 'm');
% netcdf.putAtt(ncid,R_var_id,'scale_factor',1/1000);

% Compress and stop variable definition
netcdf.defVarDeflate(ncid,lat_var_id,true,true,9);
netcdf.defVarDeflate(ncid,lon_var_id,true,true,9);
netcdf.defVarDeflate(ncid,hRe_var_id,true,true,9);
netcdf.defVarDeflate(ncid,hIm_var_id,true,true,9);
netcdf.defVarDeflate(ncid,uRe_var_id,true,true,9);
netcdf.defVarDeflate(ncid,uIm_var_id,true,true,9);
netcdf.defVarDeflate(ncid,vRe_var_id,true,true,9);
netcdf.defVarDeflate(ncid,vIm_var_id,true,true,9);
netcdf.defVarDeflate(ncid,wct_var_id,true,true,9);
netcdf.defVarDeflate(ncid,mask_var_id,true,true,9);
netcdf.defVarDeflate(ncid,flexure_var_id,true,true,9);
% netcdf.defVarDeflate(ncid,R_var_id,true,true,9);
netcdf.endDef(ncid);

%3. Place data
netcdf.putVar(ncid,x_var_id,x);
netcdf.putVar(ncid,y_var_id,y);
netcdf.putVar(ncid,lat_var_id,ipermute(single(Lat),[2 1]));
netcdf.putVar(ncid,lon_var_id,ipermute(single(Lon),[2 1]));
netcdf.putVar(ncid,cons_var_id,1:Ncons);
netcdf.putVar(ncid,amp_var_id,amp);
netcdf.putVar(ncid,ph_var_id,ph);
netcdf.putVar(ncid,om_var_id,omega);
netcdf.putVar(ncid,alp_var_id,alpha);
netcdf.putVar(ncid,hRe_var_id,ipermute(int16(scale_h*real(h)),[2 1 3]));
netcdf.putVar(ncid,hIm_var_id,ipermute(int16(scale_h*imag(h)),[2 1 3]));
netcdf.putVar(ncid,uRe_var_id,ipermute(int16(scale_UV*real(U)),[2 1 3]));
netcdf.putVar(ncid,uIm_var_id,ipermute(int16(scale_UV*imag(U)),[2 1 3]));
netcdf.putVar(ncid,vRe_var_id,ipermute(int16(scale_UV*real(V)),[2 1 3]));
netcdf.putVar(ncid,vIm_var_id,ipermute(int16(scale_UV*imag(V)),[2 1 3]));
netcdf.putVar(ncid,wct_var_id,ipermute(wct,[2 1]));
netcdf.putVar(ncid,mask_var_id,ipermute(mask,[2 1]));
netcdf.putVar(ncid,flexure_var_id,ipermute(flexure,[2 1]));
% netcdf.putVar(ncid,R_var_id,ipermute(h_range*1000,[2 1]));

%4. Close file 
netcdf.close(ncid)

disp 'done'


