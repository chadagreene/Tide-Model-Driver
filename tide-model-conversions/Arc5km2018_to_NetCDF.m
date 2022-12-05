% This file reads the old binary Arc5km2018 tide model and converts it to NetCDF. 
% 
% This script calls some legacy functions from TMD2.5 to load the old data.
% 
% Written by Chad A. Greene, NASA/JPL, May 2022. 

addpath(genpath('/Users/cgreene/Documents/MATLAB/TMD3.00_alpha'))
addpath(genpath('/Users/cgreene/Documents/data/tides/Arc5km2018'))

%% Enter initial file info 

% Output file: 
newfilename = ['/Users/cgreene/Documents/data/tides/Arc5km2018.nc']; 

% Input files: 
filename_grd = '/Users/cgreene/Documents/data/tides/Arc5km2018/grid_Arc5km2018'; 
filename_h = '/Users/cgreene/Documents/data/tides/Arc5km2018/h_Arc5km2018'; 
filename_u = '/Users/cgreene/Documents/data/tides/Arc5km2018/UV_Arc5km2018'; 

res = 5; % (km) input grid resolution 

con_string = 'm2 s2 n2 k2 k1 o1 p1 q1 m4 mn4 ms4 2n2'; % constituents in original model order, obtained by rd_con('h_Arc5km2018') 

%% Load data

% Read the grid file: 
[ll_lims,wct,mask,~,~] = grd_in(filename_grd); 

[~,masku,maskv] = Muv(mask);

% Reorient and reduce data size:
mask = uint8(flipud(mask')); 
masku = uint8(flipud(masku')); 
maskv = uint8(flipud(maskv')); 
wct = flipud(wct'); 

% Create spatial arrays: 
x = (ll_lims(1)+res/2):res:(ll_lims(2)-res/2); 
y = (ll_lims(4)-res/2):-res:(ll_lims(3)+res/2); 

% Number of constituents: 
Ncons = length(strsplit(con_string,' ')); 

for k=1:Ncons
   
   % Read each tidal height constituent and reorient: 
   tmp = h_in(filename_h,k); 
   tmp(abs(tmp)>100) = NaN; % fixes a crazy point 
   h(:,:,k) = flipud(tmp'); 
   
   % Read each transport constituent and reorient: 
   [tmp,tmp2] = u_in(filename_u,k); 
   U(:,:,k) = flipud(tmp'); 
   V(:,:,k) = flipud(tmp2'); 
end

%% Look at phase  

oldaz = angle(h(:,:,1)); 

[X,Y] = meshgrid(x,y); 
[Lon,Lat] = xy_ll_Arc5km2018(X,Y,'B');

if true
   figure
   subplot(1,2,1)
   pcolorpsn(Lat,Lon,oldaz); 
   axis tight off
   bedmachine('gl','color',rgb('gray'),'greenland')
   cmocean phase 
   ax(1) = gca; 
end

%% Fill NaNs

% Fix u and v masks where all data are zero amplitude: 
masku(all(abs(U)==0,3) | all(isnan(U),3)) = 0; 
maskv(all(abs(V)==0,3) | all(isnan(V),3)) = 0; 

for k = 1:Ncons
   tmp = h(:,:,k); 
   tmp = complex(regionfill(real(tmp),mask==0),regionfill(imag(tmp),mask==0));
   h(:,:,k) = tmp; 
   
   tmp = U(:,:,k); 
   tmp = complex(regionfill(real(tmp),masku==0),regionfill(imag(tmp),masku==0));
   U(:,:,k) = tmp; 
   
   tmp = V(:,:,k); 
   tmp = complex(regionfill(real(tmp),maskv==0),regionfill(imag(tmp),maskv==0));
   V(:,:,k) = tmp; 
   
   k
end

if true
   subplot(1,2,2)
   pcolorpsn(Lat,Lon,angle(h(:,:,1))); 
   axis tight off
   bedmachine('gl','color',rgb('gray'),'greenland')
   cmocean phase 
   ax(2) = gca; 
   linkaxes(ax,'xy'); 
end

%% Interpolate to ps grid and Fix half-cell offset in transport variables 

% Create a grid in polar stereographic coordinates: 
xps = -2850:res:2495; 
yps  = 2640:-res:-2855; 
% xps = -2645:res:2855; 
% yps = 2495:-res:-2850; 
[Xps,Yps] = meshgrid(xps,yps); 
[Latps,Lonps] = tmd_ps2ll(Xps,Yps,70,0,'N'); 

% Get the Arc5km grid coordinates of the ps grid: 
[Xi,Yi] = xy_ll_Arc5km2018(Lonps,Latps,'F'); 

% Interpolate to new grid
for k = 1:Ncons
   htmp(:,:,k) = interp2(x,y,h(:,:,k),Xi,Yi,'spline'); 
   Utmp(:,:,k) = interp2(x-res/2,y,U(:,:,k),Xi,Yi,'spline'); % half pixel offset for transport variables  
   Vtmp(:,:,k) = interp2(x,y-res/2,V(:,:,k),Xi,Yi,'spline'); 
   k
end

% Interpolate mask and wct to new grid, filling extrapolated values appropriately. 
mask = interp2(x,y,mask,Xi,Yi,'nearest',2);
wct = uint16(interp2(x,y,wct,Xi,Yi,'linear',32767)); 

% Overwrite old h, U, and V with new, interpolated values: 
h = htmp; 
U = Utmp; 
V = Vtmp; 
clear htmp Utmp Vtmp 

% Overwrite old coordinates: 
x = xps; 
y = yps; 
Lat = Latps; 
Lon = Lonps; 

%% Fill missing (extrapolated) values 

for k=1:Ncons
   tmpr = real(h(:,:,k)); 
   tmpi = imag(h(:,:,k)); 
   tmpr(mask==2) = 32767; 
   tmpi(mask==2) = 32767; 
   h(:,:,k) = complex(tmpr,tmpi); 
   
   tmpr = real(U(:,:,k)); 
   tmpi = imag(U(:,:,k)); 
   tmpr(mask==2) = 32767; 
   tmpi(mask==2) = 32767; 
   U(:,:,k) = complex(tmpr,tmpi); 
   
   tmpr = real(V(:,:,k)); 
   tmpi = imag(V(:,:,k)); 
   tmpr(mask==2) = 32767; 
   tmpi(mask==2) = 32767; 
   V(:,:,k) = complex(tmpr,tmpi); 
end

%%

% Maximum real and imaginary values for each constitutent (breaking it up by constituent reduces memory.) 
for k=1:Ncons
   tmp = abs(real(h(:,:,k))); 
   mxhr(k) = max(tmp(mask==1)); 
   tmp = abs(imag(h(:,:,k))); 
   mxhi(k) = max(tmp(mask==1)); 
   
   tmp = abs(real(U(:,:,k))); 
   mxur(k) = max(tmp(mask==1)); 
   tmp = abs(imag(U(:,:,k))); 
   mxui(k) = max(tmp(mask==1)); 
   
   tmp = abs(real(V(:,:,k))); 
   mxvr(k) = max(tmp(mask==1)); 
   tmp = abs(imag(V(:,:,k))); 
   mxvi(k) = max(tmp(mask==1)); 
   k
   
end

% Scaling factor for saving to NCSHORT (int16):
scale_h = 32765/max([mxhr mxhi])
scale_UV= 32765/max([mxur mxui mxvr mxvi])

%%

[ispec,amp,ph,omega,alpha] = tmd_constit(strsplit(con_string));

proj4 = '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +type=crs';

%% Write the netcdf 

% 1. Create netCDF file handle and Attributes
mode = netcdf.getConstant('NETCDF4');
mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
ncid=netcdf.create(newfilename,mode);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Conventions','CF-1.7');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Title','Arc5km2018');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Description','The Arctic Ocean Tidal Inverse Model developed in 2018 (Arc5km2018) is a barotropic tide model on a 5 km polar stereographic grid. It is an update to the original Arctic Ocean Tidal Inverse Model (AOTIM5) model developed in 2004, described by Padman and Erofeeva (2004) (https://doi.org/10.1029/2003GL019003).');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Author','Susan L. Howard & Laurie Padman');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'creation_date',datestr(now,'yyyy-mm-dd'));
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'NetCDF_conversion','Chad A. Greene');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'tmd_version',3.0);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'model_type','ocean');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'license','Creative Commons Attribution 4.0 International License');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Data_citation','Erofeeva S and Egbert G. 2020. Arc5km2018: Arctic Ocean Inverse Tide Model on a 5 kilometer grid, 2018. DOI information can be found at https://arcticdata.io/.')

% 2. Define dimensions
% Define mapping variable
mapping_var_id= netcdf.defVar(ncid,'mapping','NC_CHAR',[]);
netcdf.putAtt(ncid,mapping_var_id,'grid_mapping_name', 'polar_stereographic');
netcdf.putAtt(ncid,mapping_var_id,'latitude_of_projection_origin',90);
netcdf.putAtt(ncid,mapping_var_id,'standard_parallel',70);
netcdf.putAtt(ncid,mapping_var_id,'straight_vertical_longitude_from_pole',0);
netcdf.putAtt(ncid,mapping_var_id,'false_easting',0);
netcdf.putAtt(ncid,mapping_var_id,'false_northing',0); 
netcdf.putAtt(ncid,mapping_var_id,'spatial_proj4',proj4); 

% Define x: 
x_id     = netcdf.defDim(ncid,'x',length(x));
x_var_id = netcdf.defVar(ncid,'x','NC_FLOAT',x_id);
netcdf.putAtt(ncid,x_var_id,'standard_name','projection_x_coordinate');
netcdf.putAtt(ncid,x_var_id,'long_name',    'Cartesian x-coordinate, grid cell center');
netcdf.putAtt(ncid,x_var_id,'units',        'kilometer');

% Define y: 
y_id     = netcdf.defDim(ncid,'y',length(y));
y_var_id = netcdf.defVar(ncid,'y','NC_FLOAT',y_id);
netcdf.putAtt(ncid,y_var_id,'standard_name','projection_y_coordinate');
netcdf.putAtt(ncid,y_var_id,'long_name',    'Cartesian y-coordinate, grid cell center');
netcdf.putAtt(ncid,y_var_id,'units',        'kilometer');

% Define lat
lat_var_id = netcdf.defVar(ncid,'lat','NC_FLOAT',[x_id y_id]);
netcdf.putAtt(ncid,lat_var_id,'standard_name','latitude');
netcdf.putAtt(ncid,lat_var_id,'long_name',    'grid cell center latitude');
netcdf.putAtt(ncid,lat_var_id,'units',        'degree');

% Define lon
lon_var_id = netcdf.defVar(ncid,'lon','NC_FLOAT',[x_id y_id]);
netcdf.putAtt(ncid,lon_var_id,'standard_name','longitude');
netcdf.putAtt(ncid,lon_var_id,'long_name',    'grid cell center longitude');
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
netcdf.putAtt(ncid,hRe_var_id,'standard_name','height_coefficient');
netcdf.putAtt(ncid,hRe_var_id,'long_name',    'real component of height constituent');
netcdf.putAtt(ncid,hRe_var_id,'grid_mapping', 'polar_stereographic');
netcdf.putAtt(ncid,hRe_var_id,'units',        'm');
netcdf.putAtt(ncid,hRe_var_id,'scale_factor',  1/scale_h);
netcdf.defVarFill(ncid,hRe_var_id,false,uint16(32767))

% Define hIm
hIm_var_id = netcdf.defVar(ncid,'hIm','NC_SHORT',[x_id y_id cons_id]);
netcdf.putAtt(ncid,hIm_var_id,'standard_name','height_coefficient');
netcdf.putAtt(ncid,hIm_var_id,'long_name',    'imaginary component of height constituent');
netcdf.putAtt(ncid,hIm_var_id,'grid_mapping', 'polar_stereographic');
netcdf.putAtt(ncid,hIm_var_id,'units',        'm');
netcdf.putAtt(ncid,hIm_var_id,'scale_factor',  1/scale_h);
netcdf.defVarFill(ncid,hIm_var_id,false,uint16(32767))

% Define uRe
uRe_var_id = netcdf.defVar(ncid,'URe','NC_SHORT',[x_id y_id cons_id]);
netcdf.putAtt(ncid,uRe_var_id,'standard_name','transport_coefficient');
netcdf.putAtt(ncid,uRe_var_id,'long_name',    'real component of U transport constituent. This is the zonal (east-west) flow component in geographic coordinates.');
netcdf.putAtt(ncid,uRe_var_id,'grid_mapping', 'polar_stereographic');
netcdf.putAtt(ncid,uRe_var_id,'units',        'm^2/s');
netcdf.putAtt(ncid,uRe_var_id,'scale_factor',  1/scale_UV);
netcdf.defVarFill(ncid,uRe_var_id,false,uint16(32767))

% Define uIm
uIm_var_id = netcdf.defVar(ncid,'UIm','NC_SHORT',[x_id y_id cons_id]);
netcdf.putAtt(ncid,uIm_var_id,'standard_name','transport_coefficient');
netcdf.putAtt(ncid,uIm_var_id,'long_name',    'imaginary component of U transport constituent. This is the zonal (east-west) flow component in geographic coordinates.');
netcdf.putAtt(ncid,uIm_var_id,'grid_mapping', 'polar_stereographic');
netcdf.putAtt(ncid,uIm_var_id,'units',        'm^2/s');
netcdf.putAtt(ncid,uIm_var_id,'scale_factor',  1/scale_UV);
netcdf.defVarFill(ncid,uIm_var_id,false,uint16(32767))

% Define vRe
vRe_var_id = netcdf.defVar(ncid,'VRe','NC_SHORT',[x_id y_id cons_id]);
netcdf.putAtt(ncid,vRe_var_id,'standard_name','transport_coefficient');
netcdf.putAtt(ncid,vRe_var_id,'long_name',    'real component of V transport constituent. This is the meridional (north-south) flow component in geographic coordinates.');
netcdf.putAtt(ncid,vRe_var_id,'grid_mapping', 'polar_stereographic');
netcdf.putAtt(ncid,vRe_var_id,'units',        'm^2/s');
netcdf.putAtt(ncid,vRe_var_id,'scale_factor',  1/scale_UV);
netcdf.defVarFill(ncid,vRe_var_id,false,uint16(32767))

% Define vIm
vIm_var_id = netcdf.defVar(ncid,'VIm','NC_SHORT',[x_id y_id cons_id]);
netcdf.putAtt(ncid,vIm_var_id,'standard_name','transport_coefficient');
netcdf.putAtt(ncid,vIm_var_id,'long_name',    'imaginary component of V transport constituent. This is the meridional (north-south) flow component in geographic coordinates.');
netcdf.putAtt(ncid,vIm_var_id,'grid_mapping', 'polar_stereographic');
netcdf.putAtt(ncid,vIm_var_id,'units',        'm^2/s');
netcdf.putAtt(ncid,vIm_var_id,'scale_factor',  1/scale_UV);
netcdf.defVarFill(ncid,vIm_var_id,false,uint16(32767))

% Define wct: 
wct_var_id = netcdf.defVar(ncid,'wct','NC_SHORT',[x_id y_id]);
netcdf.putAtt(ncid,wct_var_id,'standard_name','wct');
netcdf.putAtt(ncid,wct_var_id,'long_name','water column thickness');
netcdf.putAtt(ncid,wct_var_id,'units',    'meters');
netcdf.putAtt(ncid,wct_var_id,'grid_mapping', 'polar_stereographic');
netcdf.defVarFill(ncid,wct_var_id,false,uint16(32767))

% Define mask
mask_var_id = netcdf.defVar(ncid,'mask','NC_BYTE',[x_id y_id]);
netcdf.putAtt(ncid,mask_var_id,'long_name',    'ocean mask');
netcdf.putAtt(ncid,mask_var_id,'standard_name','ocean_mask');
netcdf.putAtt(ncid,mask_var_id,'grid_mapping', 'polar_stereographic');
netcdf.putAtt(ncid,mask_var_id,'valid_range',  [0 1 2]);
netcdf.putAtt(ncid,mask_var_id,'flag_values',  [0 1 2]);
netcdf.putAtt(ncid,mask_var_id,'flag_meanings','0=land, 1=ocean, 2=out-of-domain');

% % Compress and stop variable definition
% netcdf.defVarDeflate(ncid,lat_var_id,true,true,9);
% netcdf.defVarDeflate(ncid,lon_var_id,true,true,9);
% netcdf.defVarDeflate(ncid,hRe_var_id,true,true,9);
% netcdf.defVarDeflate(ncid,hIm_var_id,true,true,9);
% netcdf.defVarDeflate(ncid,uRe_var_id,true,true,9);
% netcdf.defVarDeflate(ncid,uIm_var_id,true,true,9);
% netcdf.defVarDeflate(ncid,vRe_var_id,true,true,9);
% netcdf.defVarDeflate(ncid,vIm_var_id,true,true,9);
% netcdf.defVarDeflate(ncid,wct_var_id,true,true,9);
% netcdf.defVarDeflate(ncid,mask_var_id,true,true,9);
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

%4. Close file 
netcdf.close(ncid)

disp done
