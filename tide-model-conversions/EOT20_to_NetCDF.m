
% This script calls some legacy functions from TMD2.5 to load the old data.
% 
% Written by Chad A. Greene, NASA/JPL, June 2022. 

%% Enter initial file info 

% Output file: 
newfilename = ['/Users/cgreene/Documents/data/tides/EOT20_update_',datestr(now,'yyyy-mm-dd'),'.nc']; 

con_string = 'm2 s2 k1 o1 n2 p1 k2 q1 2n2 mf mm m4 ms4 mn4 s1'; % constituents in proper order
uv = false; 

%% Load data

% Read the grid file: 
[ll_lims,wct,mask,~,~] = grd_in('/Users/cgreene/Downloads/TPXO9_atlas_v5/grid_tpxo9_atlas_30_v5'); 

[lon,lat] = XY(ll_lims,size(wct,1),size(wct,2));
lat = fliplr(lat); 

% Reorient and reduce data size:
mask = uint8(flipud(mask')); 
wct = uint16(flipud(wct')); 

% Number of constituents: 
cons_cell = strsplit(con_string,' ');
Ncons = length(cons_cell); 

for k=1:Ncons
   
   
   ftmp = ['/Users/cgreene/Downloads/TPXO9_atlas_v5/h_',cons_cell{k},'_tpxo9_atlas_30_v5'];
   
   % Read each tidal height constituent and reorient: 
   tmp = h_in(ftmp,k); 
   h(:,:,k) = flipud(tmp'); 
   
   if uv
      % Read each transport constituent and reorient: 
      [tmp,tmp2] = u_in(futmp,k); 
      U(:,:,k) = flipud(tmp'); 
      V(:,:,k) = flipud(tmp2'); 
   end
end

%% Extend grid across the prime meridian  
% This is just to allow seamless interpolation. 

lon = [lon(end)-360, lon, lon(1)+360]; 
h = cat(2,h(:,end,:),h,h(:,1,:));
wct = cat(2,wct(:,end),wct,wct(:,1));
mask = cat(2,mask(:,end),mask,mask(:,1));

if uv 
   U = cat(2,U(:,end,:),U,U(:,1,:));
   V = cat(2,V(:,end,:),V,V(:,1,:));
end

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
if uv
   error 'rewrite scale_UV like scale_h above'; 
   %scale_UV = 32767/max(abs([real(U(:));real(V(:));imag(U(:));imag(V(:));]));
end

[ispec,amp,ph,omega,alpha,constitNum] = tmd_constit(cons_cell);

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
   
   if uv
      tmp = U(:,:,k); 
      tmp = complex(regionfill(real(tmp),mask==0),regionfill(imag(tmp),mask==0));
      U(:,:,k) = tmp; 

      tmp = V(:,:,k); 
      tmp = complex(regionfill(real(tmp),mask==0),regionfill(imag(tmp),mask==0));
      V(:,:,k) = tmp; 
   end
   
   k
end

for k = 1:15
   
end

%% Write the netcdf 

% 1. Create netCDF file handle and Attributes
mode = netcdf.getConstant('NETCDF4');
mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
ncid=netcdf.create(newfilename,mode);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Conventions','CF-1.7');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Title','TPXO9_atlas_v5');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Description','Global tide model at 1/30 degree resolution.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Author','Chad A. Greene');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'creation_date',datestr(now,'yyyy-mm-dd'));
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'tmd_version',3.0);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'license','ask');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Data_citation','Egbert, Gary D., and Svetlana Y. Erofeeva. "Efficient inverse modeling of barotropic ocean tides." Journal of Atmospheric and Oceanic Technology 19.2 (2002): 183-204.')
   
% 2. Define dimensions
% Define mapping variable
mapping_var_id= netcdf.defVar(ncid,'mapping','NC_CHAR',[]);
netcdf.putAtt(ncid,mapping_var_id,'grid_mapping_name', 'latitude_longitude');
netcdf.putAtt(ncid,mapping_var_id,'epsg_code',4326);
netcdf.putAtt(ncid,mapping_var_id,'spatial_proj4',proj4); 

% Define lon: 
lon_id     = netcdf.defDim(ncid,'lon',length(lon));
lon_var_id = netcdf.defVar(ncid,'lon','NC_FLOAT',lon_id);
netcdf.putAtt(ncid,lon_var_id,'long_name',    'grid cell center longitude (first and last columns are repeats, to enable seamless interpolation)');
netcdf.putAtt(ncid,lon_var_id,'standard_name','longitude');
netcdf.putAtt(ncid,lon_var_id,'units',        'degrees');

% Define y: 
lat_id     = netcdf.defDim(ncid,'lat',length(lat));
lat_var_id = netcdf.defVar(ncid,'lat','NC_FLOAT',lat_id);
netcdf.putAtt(ncid,lat_var_id,'long_name',    'grid cell center latitude');
netcdf.putAtt(ncid,lat_var_id,'standard_name','latitude');
netcdf.putAtt(ncid,lat_var_id,'units',        'degrees');

% Define constituents
cons_id = netcdf.defDim(ncid,'cons',Ncons); 
cons_var_id = netcdf.defVar(ncid,'cons','NC_BYTE',cons_id); 
netcdf.putAtt(ncid,cons_var_id,'standard_name',    'tidal_constituents');
netcdf.putAtt(ncid,cons_var_id,'long_name',con_string); 

% Define constituent attributes 
amp_var_id = netcdf.defVar(ncid,'amplitude','NC_DOUBLE',cons_id); 
netcdf.putAtt(ncid,amp_var_id,'standard_name',    'amplitude');
netcdf.putAtt(ncid,amp_var_id,'long_name','constituent amplitude'); 

ph_var_id = netcdf.defVar(ncid,'phase','NC_DOUBLE',cons_id); 
netcdf.putAtt(ncid,ph_var_id,'standard_name',    'phase');
netcdf.putAtt(ncid,ph_var_id,'long_name','Astronomical arguments (relative to t0 = 1 Jan 0:00 1992)'); 

om_var_id = netcdf.defVar(ncid,'omega','NC_FLOAT',cons_id); 
netcdf.putAtt(ncid,om_var_id,'standard_name',    'omega');
netcdf.putAtt(ncid,om_var_id,'long_name','frequency'); 

alp_var_id = netcdf.defVar(ncid,'alpha','NC_FLOAT',cons_id); 
netcdf.putAtt(ncid,alp_var_id,'standard_name',    'alpha');
netcdf.putAtt(ncid,alp_var_id,'long_name','loading love number'); 

% Define hRe
hRe_var_id = netcdf.defVar(ncid,'hRe','NC_SHORT',[lon_id lat_id cons_id]);
netcdf.putAtt(ncid,hRe_var_id,'long_name',    'real component of height constituent');
netcdf.putAtt(ncid,hRe_var_id,'standard_name','height_constituent');
netcdf.putAtt(ncid,hRe_var_id,'grid_mapping', 'polar_stereographic');
netcdf.putAtt(ncid,hRe_var_id,'units',        'm');
netcdf.putAtt(ncid,hRe_var_id,'scale_factor',  1/scale_h);

% Define hIm
hIm_var_id = netcdf.defVar(ncid,'hIm','NC_SHORT',[lon_id lat_id cons_id]);
netcdf.putAtt(ncid,hIm_var_id,'long_name',    'imaginary component of height constituent');
netcdf.putAtt(ncid,hIm_var_id,'standard_name','height_constituent');
netcdf.putAtt(ncid,hIm_var_id,'grid_mapping', 'polar_stereographic');
netcdf.putAtt(ncid,hIm_var_id,'units',        'm');
netcdf.putAtt(ncid,hIm_var_id,'scale_factor',  1/scale_h);

if uv
   % Define uRe
   uRe_var_id = netcdf.defVar(ncid,'URe','NC_SHORT',[lon_id lat_id cons_id]);
   netcdf.putAtt(ncid,uRe_var_id,'long_name',    'real component of U transport constituent.');
   netcdf.putAtt(ncid,uRe_var_id,'standard_name','height_constituent');
   netcdf.putAtt(ncid,uRe_var_id,'grid_mapping', 'polar_stereographic');
   netcdf.putAtt(ncid,uRe_var_id,'units',        'm^2/s');
   netcdf.putAtt(ncid,uRe_var_id,'scale_factor',  1/scale_UV);

   % Define uIm
   uIm_var_id = netcdf.defVar(ncid,'UIm','NC_SHORT',[lon_id lat_id cons_id]);
   netcdf.putAtt(ncid,uIm_var_id,'long_name',    'imaginary component of U transport constituent.');
   netcdf.putAtt(ncid,uIm_var_id,'standard_name','height_constituent');
   netcdf.putAtt(ncid,uIm_var_id,'grid_mapping', 'polar_stereographic');
   netcdf.putAtt(ncid,uIm_var_id,'units',        'm^2/s');
   netcdf.putAtt(ncid,uIm_var_id,'scale_factor',  1/scale_UV);

   % Define vRe
   vRe_var_id = netcdf.defVar(ncid,'VRe','NC_SHORT',[lon_id lat_id cons_id]);
   netcdf.putAtt(ncid,vRe_var_id,'long_name',    'real component of V transport constituent.');
   netcdf.putAtt(ncid,vRe_var_id,'standard_name','height_constituent');
   netcdf.putAtt(ncid,vRe_var_id,'grid_mapping', 'polar_stereographic');
   netcdf.putAtt(ncid,vRe_var_id,'units',        'm^2/s');
   netcdf.putAtt(ncid,vRe_var_id,'scale_factor',  1/scale_UV);

   % Define vIm
   vIm_var_id = netcdf.defVar(ncid,'VIm','NC_SHORT',[lon_id lat_id cons_id]);
   netcdf.putAtt(ncid,vIm_var_id,'long_name',    'imaginary component of V transport constituent.');
   netcdf.putAtt(ncid,vIm_var_id,'standard_name','height_constituent');
   netcdf.putAtt(ncid,vIm_var_id,'grid_mapping', 'polar_stereographic');
   netcdf.putAtt(ncid,vIm_var_id,'units',        'm^2/s');
   netcdf.putAtt(ncid,vIm_var_id,'scale_factor',  1/scale_UV);
end

% Define wct: 
wct_var_id = netcdf.defVar(ncid,'wct','NC_SHORT',[lon_id lat_id]);
netcdf.putAtt(ncid,wct_var_id,'long_name','water column thickness');
netcdf.putAtt(ncid,wct_var_id,'standard_name','wct');
netcdf.putAtt(ncid,wct_var_id,'units',    'meters');
netcdf.putAtt(ncid,wct_var_id,'grid_mapping', 'polar_stereographic');

% Define mask
mask_var_id = netcdf.defVar(ncid,'mask','NC_BYTE',[lon_id lat_id]);
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
if uv
   netcdf.defVarDeflate(ncid,uRe_var_id,true,true,9);
   netcdf.defVarDeflate(ncid,uIm_var_id,true,true,9);
   netcdf.defVarDeflate(ncid,vRe_var_id,true,true,9);
   netcdf.defVarDeflate(ncid,vRe_var_id,true,true,9);
end
netcdf.defVarDeflate(ncid,wct_var_id,true,true,9);
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
if uv
   netcdf.putVar(ncid,uRe_var_id,ipermute(int16(scale_UV*real(U)),[2 1 3]));
   netcdf.putVar(ncid,uIm_var_id,ipermute(int16(scale_UV*imag(U)),[2 1 3]));
   netcdf.putVar(ncid,vRe_var_id,ipermute(int16(scale_UV*real(V)),[2 1 3]));
   netcdf.putVar(ncid,vIm_var_id,ipermute(int16(scale_UV*imag(V)),[2 1 3]));
end
netcdf.putVar(ncid,wct_var_id,ipermute(wct,[2 1]));
netcdf.putVar(ncid,mask_var_id,ipermute(mask,[2 1]));

%4. Close file 
netcdf.close(ncid)

disp done

toc/60