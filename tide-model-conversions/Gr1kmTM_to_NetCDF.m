% This file reads the old binary Gr1kmTM tide model and converts it to NetCDF. 
% 
% This script calls some legacy functions from TMD2.5 to load the old data.
% 
% Written by Chad A. Greene, NASA/JPL, May 2022. 

addpath(genpath('/Users/cgreene/Documents/MATLAB/TMD3.00_alpha'))

%% Enter initial file info 

% Output file: 
newfilename = ['/Users/cgreene/Documents/data/tides/Gr1kmTM_update_v1.nc']; 

% Input files: 
filename_grd = '/Users/cgreene/Documents/data/tides/Gr1kmTM/grid_Gr1kmTM_v1'; 
filename_h = '/Users/cgreene/Documents/data/tides/Gr1kmTM/h_Gr1kmTM_v1'; 
filename_u = '/Users/cgreene/Documents/data/tides/Gr1kmTM/UV_Gr1kmTM_v1'; 

res = 1; % (km) input grid resolution 

con_string = 'm2 s2 k1 o1 n2 p1 k2 q1'; % constituents in original model order

%% Load data

% Read the grid file: 
[ll_lims,wct,mask,~,~] = grd_in(filename_grd); 

% Reorient and reduce data size:
mask = uint8(flipud(mask')); 
wct = uint16(flipud(wct')); 

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
[Lat,Lon] = mapxy(X,Y,70,-45,'n');

figure
subplot(1,2,1)
pcolorpsn(Lat,Lon,oldaz); 
axis tight off
bedmachine('gl','color',rgb('gray'),'greenland')
cmocean phase 
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

subplot(1,2,2)
pcolorpsn(Lat,Lon,angle(h(:,:,1))); 
axis tight off
bedmachine('gl','color',rgb('gray'),'greenland')
cmocean phase 
ax(2) = gca; 
linkaxes(ax,'xy'); 

%% Fix half-cell offset in transport variables 

for k = 1:Ncons
   U(:,:,k) = interp2(x-res/2,y,U(:,:,k),X,Y); 
   V(:,:,k) = interp2(x,y-res/2,V(:,:,k),X,Y); 
   k
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
scale_h = 32767/max([mxhr mxhi])
scale_UV= 32767/max([mxur mxui mxvr mxvi])

%%

[ispec,amp,ph,omega,alpha] = tmd_constit(strsplit(con_string));

proj4 = '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +type=crs';

%% Write the netcdf 

% 1. Create netCDF file handle and Attributes
mode = netcdf.getConstant('NETCDF4');
mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
ncid=netcdf.create(newfilename,mode);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Conventions','CF-1.7');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Title','Gr1kmTM');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Description','The Greenland 1 kilometer Tide Model (Gr1kmTM) is a barotropic ocean tide model on a 1 km x 1 km polar stereographic grid, developed using the Regional Ocean Modeling System (ROMS). Gr1kmTM consists of spatial grids of complex amplitude coefficients for sea surface height and depth-integrated currents (“volume transports”) for 8 principal tidal constituents: 4 semidiurnal (M2, S2, K2, N2) and 4 diurnal (K1, O1, P1, Q1).');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Author','Susan L. Howard & Laurie Padman');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'creation_date',datestr(now,'yyyy-mm-dd'));
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'NetCDF_conversion','Chad A. Greene');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'tmd_version',3.0);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'license','MIT License');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Data_citation',['Susan L. Howard and Laurie Padman. 2021. Gr1kmTM: Greenland 1 kilometer Tide Model, 2021. urn:node:ARCTIC. doi:10.18739/A2251FM3S.'])

% 2. Define dimensions
% Define mapping variable
mapping_var_id= netcdf.defVar(ncid,'mapping','NC_CHAR',[]);
netcdf.putAtt(ncid,mapping_var_id,'grid_mapping_name', 'polar_stereographic');
netcdf.putAtt(ncid,mapping_var_id,'latitude_of_projection_origin',90);
netcdf.putAtt(ncid,mapping_var_id,'standard_parallel',70);
netcdf.putAtt(ncid,mapping_var_id,'straight_vertical_longitude_from_pole',-45);
netcdf.putAtt(ncid,mapping_var_id,'false_easting',0);
netcdf.putAtt(ncid,mapping_var_id,'false_northing',0); 
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

