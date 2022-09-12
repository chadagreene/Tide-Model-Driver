
% 
% Written by Chad A. Greene, NASA/JPL, June 2022. 
%
% This script takes about 3 minutes to run on my 2019-era laptop. 

%% Enter initial file info 

% Output file: 
newfilename = '/Users/cgreene/Documents/data/tides/EOT20_load.nc'; 

con_string = '2n2 j1 k1 k2 m2 m4 mf mm n2 o1 p1 q1 s1 s2 sa ssa t2'; % constituents in alphabetical order

%% Load data

% Read the grid file: 
ftmp = '/Users/cgreene/Documents/data/tides/EOT20/load_tides/2N2_load_eot20.nc'; 
lon = ncread(ftmp,'lon')'; 
lat = flipud(ncread(ftmp,'lat')); 

mask = uint8(rot90(isfinite(ncread(ftmp,'real')))); 

% Number of constituents: 
cons_cell = strsplit(con_string,' ');
Ncons = length(cons_cell); 

for k=1:Ncons
   
   ftmp = ['/Users/cgreene/Documents/data/tides/EOT20/load_tides/',upper(cons_cell{k}),'_load_EOT20.nc'];
   
   % Read each tidal height constituent and reorient: 
   tmpr = rot90(ncread(ftmp,'real')); 
   tmpi = rot90(ncread(ftmp,'imag')); 
   
   h(:,:,k) = complex(tmpr,tmpi)/100; % divide by 100 to convert cm to m
   k
end

clear tmpr tmpi

%% Extend grid across the prime meridian  
% This is just to allow seamless interpolation. 

lon = [lon(1)-diff(lon(1:2)), lon, lon(end)+diff(lon(1:2))]; 
h = cat(2,h(:,end,:),h,h(:,1,:));
mask = cat(2,mask(:,end),mask,mask(:,1));

%%

% Maximum real and imaginary values for each constitutent (breaking it up by constituent reduces memory.) 
for k=1:Ncons
   tmp = abs(real(h(:,:,k))); 
   mxhr(k) = max(tmp(mask==1)); 
   tmp = abs(imag(h(:,:,k))); 
   mxhi(k) = max(tmp(mask==1)); 
   k
end

% Scaling factor for saving to NCSHORT (int16):
scale_h = 32767/max([mxhr mxhi]);

[ispec,amp,ph,omega,alpha] = tmd_constit(cons_cell);

proj4 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs';

%% Fill NaNs
% You may be asking: "Why are we filling the land areas with tide values??"
% It's because a tide gauge might find itself between an ocean pixel and a 
% land pixel, and in those cases, this regionfill approach should produce 
% the most reasonable estimate of tides. It's only a short extrapolation,
% but this should be smoother and more physical than something like
% nearest-neighbor extrapolation. 
% 
% This section of code might take 40 minutes or so to run.

for k = 1:Ncons
   tmp = h(:,:,k); 
   tmp = complex(regionfill(real(tmp),mask==0),regionfill(imag(tmp),mask==0));
   h(:,:,k) = tmp; 
   
   k
end

%% Write the netcdf 

% 1. Create netCDF file handle and Attributes
mode = netcdf.getConstant('NETCDF4');
mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
ncid=netcdf.create(newfilename,mode);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Conventions','CF-1.7');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Title','DGFI-TUM global empirical ocean tide model EOT20, load tide component.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Description','EOT20 is the latest in a series of global ocean tide models developed at DGFI-TUM. EOT20 is created using residual tidal analysis of multi-mission altimetry data from the period of 1992 to 2019. Eleven satellite altimetry missions are used and the FES2014b tide model used as the reference model for the residual tidal analysis. The model extends throughout the global ocean, with EOT20 ranging from 66°S to 66°N with FES2014b being used to fill in the higher latitudes.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Author','Hart-Davis Michael G., Piccioni Gaia, Dettmering Denise, Schwatke Christian, Passaro Marcello, Seitz Florian.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'creation_date',datestr(now,'yyyy-mm-dd'));
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'tmd_version',3.0);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'model_type','load');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'NetCDF_conversion','Chad A. Greene');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'license','None');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'data_version','1.0.0');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'contact','michael.hart-davis@tum.de');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Data_citation','Hart-Davis Michael, Piccioni Gaia, Dettmering Denise, Schwatke Christian, Passaro Marcello, Seitz Florian (2021). EOT20 - A global Empirical Ocean Tide model from multi-mission satellite altimetry. SEANOE. https://doi.org/10.17882/79489');


% 2. Define dimensions
% Define mapping variable
mapping_var_id= netcdf.defVar(ncid,'mapping','NC_CHAR',[]);
netcdf.putAtt(ncid,mapping_var_id,'grid_mapping_name', 'latitude_longitude');
netcdf.putAtt(ncid,mapping_var_id,'epsg_code',4326);
netcdf.putAtt(ncid,mapping_var_id,'spatial_proj4',proj4); 

% Define lon: 
lon_id     = netcdf.defDim(ncid,'lon',length(lon));
lon_var_id = netcdf.defVar(ncid,'lon','NC_FLOAT',lon_id);
netcdf.putAtt(ncid,lon_var_id,'standard_name','longitude');
netcdf.putAtt(ncid,lon_var_id,'long_name',    'grid cell center longitude (first and last columns are repeats, to enable seamless interpolation)');
netcdf.putAtt(ncid,lon_var_id,'units',        'degrees');

% Define y: 
lat_id     = netcdf.defDim(ncid,'lat',length(lat));
lat_var_id = netcdf.defVar(ncid,'lat','NC_FLOAT',lat_id);
netcdf.putAtt(ncid,lat_var_id,'standard_name','latitude');
netcdf.putAtt(ncid,lat_var_id,'long_name',    'grid cell center latitude');
netcdf.putAtt(ncid,lat_var_id,'units',        'degrees');

% Define constituents
cons_id = netcdf.defDim(ncid,'constituents',Ncons); 
cons_var_id = netcdf.defVar(ncid,'constituents','NC_BYTE',cons_id); 
netcdf.putAtt(ncid,cons_var_id,'standard_name', 'tidal_constituents');
netcdf.putAtt(ncid,cons_var_id,'long_name','Tidal constituents listed in order in the constituent_order attribute.'); 
netcdf.putAtt(ncid,cons_var_id,'constituent_order',con_string); 

% Define constituent attributes 
amp_var_id = netcdf.defVar(ncid,'amplitude','NC_DOUBLE',cons_id); 
netcdf.putAtt(ncid,amp_var_id,'standard_name',    'amplitude');
netcdf.putAtt(ncid,amp_var_id,'long_name','amplitude of equilibrium tide in m for each tidal constituent.'); 
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
hRe_var_id = netcdf.defVar(ncid,'hRe','NC_SHORT',[lon_id lat_id cons_id]);
netcdf.putAtt(ncid,hRe_var_id,'standard_name','height_coefficient');
netcdf.putAtt(ncid,hRe_var_id,'long_name',    'real component of height constituent');
netcdf.putAtt(ncid,hRe_var_id,'units',        'm');
netcdf.putAtt(ncid,hRe_var_id,'scale_factor',  1/scale_h);

% Define hIm
hIm_var_id = netcdf.defVar(ncid,'hIm','NC_SHORT',[lon_id lat_id cons_id]);
netcdf.putAtt(ncid,hIm_var_id,'standard_name','height_coefficient');
netcdf.putAtt(ncid,hIm_var_id,'long_name',    'imaginary component of height constituent');
netcdf.putAtt(ncid,hIm_var_id,'units',        'm');
netcdf.putAtt(ncid,hIm_var_id,'scale_factor',  1/scale_h);

% Define mask
mask_var_id = netcdf.defVar(ncid,'mask','NC_BYTE',[lon_id lat_id]);
netcdf.putAtt(ncid,mask_var_id,'standard_name','ocean_mask');
netcdf.putAtt(ncid,mask_var_id,'long_name',    'ocean mask');
netcdf.putAtt(ncid,mask_var_id,'valid_range',  [0 1]);
netcdf.putAtt(ncid,mask_var_id,'flag_values',  [0 1]);
netcdf.putAtt(ncid,mask_var_id,'flag_meanings','land ocean');

% Compress and stop variable definition
netcdf.defVarDeflate(ncid,lat_var_id,true,true,9);
netcdf.defVarDeflate(ncid,lon_var_id,true,true,9);
netcdf.defVarDeflate(ncid,hRe_var_id,true,true,9);
netcdf.defVarDeflate(ncid,hIm_var_id,true,true,9);
netcdf.defVarDeflate(ncid,mask_var_id,true,true,9);
netcdf.endDef(ncid);

%3. Place data
netcdf.putVar(ncid,lon_var_id,lon);
netcdf.putVar(ncid,lat_var_id,lat);
netcdf.putVar(ncid,cons_var_id,1:Ncons);
netcdf.putVar(ncid,amp_var_id,amp);
netcdf.putVar(ncid,ph_var_id,ph);
netcdf.putVar(ncid,om_var_id,omega);
netcdf.putVar(ncid,alp_var_id,alpha);
netcdf.putVar(ncid,hRe_var_id,ipermute(int16(scale_h*real(h)),[2 1 3]));
netcdf.putVar(ncid,hIm_var_id,ipermute(int16(scale_h*imag(h)),[2 1 3]));
netcdf.putVar(ncid,mask_var_id,ipermute(mask,[2 1]));

%4. Close file 
netcdf.close(ncid)

disp done

