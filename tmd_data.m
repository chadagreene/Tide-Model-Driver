function [Z,x_or_lon,y_or_lat,conList] = tmd_data(filename,variable,varargin) 
% tmd_data loads tide model data into the Matlab workspace. 
% 
%% Syntax 
% 
%  [Z,x_or_lon,y_or_lat] = tmd_data(filename,variable)
%  [...] = tmd_data(...,'constituents',conList)
%  [...] = tmd_data(...,'bounds',[xi yi])
%  [Z,lon,lat] = tmd_data(...,'geo')
%
%% Description 
% 
% [Z,x_or_lon,y_or_lat] = tmd_data(filename,variable)
% 
%  * 'h'   complex tidal height (m)  
%  * 'hRe' real part of tidal height
%  * 'hIm' imaginary part of tidal height 
%  * 'hAm' amplitude of tidal height
%  * 'hPh' phase of tidal height
%  * 'u'   complex zonal velocity (m/s) 
%  * 'uRe' real part of zonal velocity 
%  * 'uIm' imaginary part of zonal velocity 
%  * 'uAm' amplitude of zonal velocity
%  * 'uPh' phase of zonal velocity
%  * 'U'   complex zonal transport (m^2/s) 
%  * 'URe  real part of zonal transport
%  * 'UIm' imaginary part of zonal transport
%  * 'UAm' amplitude of zonal transport
%  * 'UPh' phase of zonal transport 
%  * 'v'   complex meridional velocity (m/s) 
%  * 'vRe' real part of meridional velocity 
%  * 'vIm' imaginary part of meridional velocity
%  * 'vAm' amplitude of meridional velocity
%  * 'vPh' phase of meridional velocity 
%  * 'V'   complex meridional transport (m^2/s)
%  * 'VRe' real part of meridional transport 
%  * 'VIm' imaginary part of meridional transport
%  * 'VAm' amplitude of meridional transport
%  * 'VPh' phase of meridional transport
%  * 'wct' water column thickness (m) 
%  * 'mask' binary land/ocean mask
%  * 'flexure' ice shelf flexure coefficient from a linear elastic model applied to BedMachine ice thickness (can slightly exceed 1). Only for CATS model. 
%
% [...] = tmd_data(...,'constituents',conList) specifies tidal constituents as a 
% cell array (e.g, {'m2','s2'}. If constituents are not specified, all constituents 
% from the model are returned. 
%
% [...] = tmd_data(...,'bounds',[xi yi]) enter an Mx2 matrix of coordinates 
% to only load data in a rectangle around the specified [xi(:) yi(:)] where
% xi and yi are projected coordinates for tide models that are in projected
% coordinates. For tide models in geo coordinates (global tide models,
% generally), xi refers to longitudes of interest and yi is latitudes of interest. 
% The advantage of specifying bounds is to minimize the amount of data that
% are loaded, when only a small area within a larger tide model is of
% interest. 
% 
%% Examples 
% For examples type 
% 
%   tmd tmd_data
% 
%% Author Info 
% This function was written by Chad A. Greene in 2022. 
% 
% See also tmd_interp and tmd_predict. 

%% Error checks 

assert(contains(filename,'.nc'),'Input filename must end in .nc.')
assert(exist(filename,'file'),['Cannot find ',filename,'. Check the path and try again.'])

% Ensure the model file is TMD3.0 compatible: 
try 
   tmd_version = ncreadatt(filename,'/','tmd_version'); 
   assert(tmd_version>=3.0,[filename,' is not compatible with TMD3.0+.'])
catch 
   error([filename,' is not compatible with TMD3.0+.'])
end

%% Input parsing 

narginchk(2,Inf)

% Set defaults 
spatialSubset = false; 
conSubset = false; 
geo = false; % output geographic coordinates 
stride = Inf; 
NCons = 1; 
bounds = [[-Inf;Inf] [-Inf;Inf]];

conList = strsplit(ncreadatt(filename,'constituents','constituent_order')); 

proj4 = ncreadatt(filename,'mapping','spatial_proj4');
if strcmpi(proj4,'+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
   GlobalModel = true; 
else 
   GlobalModel = false; 
end

if nargin>2 
   
   tmp = strncmpi(varargin,'bounds',3); 
   if any(tmp) 
      spatialSubset = true; 
      bounds = varargin{find(tmp)+1}; 
      assert(size(bounds,2)==2,'Spatial bounds must be Mx2 in the form [xi(:) yi(:)] or [loni(:) lati(:)].')
      
      % Wrap longitudes to 360 for global model: 
      if GlobalModel 
         tmp=bounds(:,1); 
         tmp(tmp<0) = tmp(tmp<0)+360; 
         bounds(:,1)  = tmp; 
      end
   end 
   
   tmp = strncmpi(varargin,'constituents',3); 
   if any(tmp) 
      conSubset = true; 
      assert(all(ismember(varargin{find(tmp)+1},conList)),'Error: Some of the requested constituents do not exist in the model file.')
      conList = varargin{find(tmp)+1}; 
      if ischar(conList) 
         conList = cellstr(conList); 
      end
      
      stride = 1; 
      NCons = length(conList); 
   end
   
   if any(strcmpi(varargin,'geo'))
      geo = true; 
   end
      
end

%% Load data

if GlobalModel
   x_or_lon = double(ncread(filename,'lon')); 
   y_or_lat = double(ncread(filename,'lat')); 
else
   x_or_lon = double(ncread(filename,'x')); 
   y_or_lat = double(ncread(filename,'y')); 
end
cons = strsplit(ncreadatt(filename,'constituents','constituent_order')); 

dx = abs(diff(x_or_lon(1:2))); 
dy = abs(diff(y_or_lat(1:2))); 

if spatialSubset
   
    % Get xlimits (xl) and ylimits (yl) of input coordinates + buffer:
    xl = [min(bounds(:,1))-2*dx max(bounds(:,1))+2*dx];
    yl = [min(bounds(:,2))-2*dy max(bounds(:,2))+2*dy];
    
    % Region of rows and columns of pixels to read: 
    ri=find((y_or_lat>=yl(1))&(y_or_lat<=yl(2)));
    ci=find((x_or_lon>=xl(1))&(x_or_lon<=xl(2)));
    
    x_or_lon = x_or_lon(ci); 
    y_or_lat = y_or_lat(ri); 
else
    ri = 1:length(y_or_lat); 
    ci = 1:length(x_or_lon); 
end
 
if conSubset
   Z = NaN(numel(ri),numel(ci),NCons);
end

for k = 1:NCons 
   
   if conSubset
      conind = find(ismember(cons,conList{k}));
      placInd = k;
   else 
      conind = 1; 
      placInd = 1:length(cons);
   end
   
   switch variable
      case 'mask'
         Z = logical(permute(ncread(filename,variable,[ci(1) ri(1)],[numel(ci) numel(ri)]),[2 1]));
      case 'flexure'
         Z = double(permute(ncread(filename,variable,[ci(1) ri(1)],[numel(ci) numel(ri)]),[2 1]))/100;
      case 'wct'
         Z = double(permute(ncread(filename,variable,[ci(1) ri(1)],[numel(ci) numel(ri)]),[2 1]));
      case {'hRe','hIm','URe','UIm','VRe','VIm'}
         Z(:,:,placInd) = double(permute(ncread(filename,variable,[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3]));
      case 'uRe'
         Z(:,:,placInd) = double(permute(ncread(filename,'URe',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3]));
      case 'uIm'
         Z(:,:,placInd) = double(permute(ncread(filename,'UIm',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3]));
      case 'vRe'
         Z(:,:,placInd) = double(permute(ncread(filename,'VRe',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3]));
      case 'vIm'
         Z(:,:,placInd) = double(permute(ncread(filename,'VIm',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3]));
      case 'h'
         Z(:,:,placInd) = complex(double(permute(ncread(filename,'hRe',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])),...
                                  double(permute(ncread(filename,'hIm',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])));
      case 'hAm' 
         Z(:,:,placInd) = abs(complex(double(permute(ncread(filename,'hRe',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])),...
                                      double(permute(ncread(filename,'hIm',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3]))));
      case 'hPh'
          Z(:,:,placInd) = angle(complex(double(permute(ncread(filename,'hRe',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])),...
                                          double(permute(ncread(filename,'hIm',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3]))));
      case {'u','U'}
         Z(:,:,placInd) = complex(double(permute(ncread(filename,'URe',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])),...
                                  double(permute(ncread(filename,'UIm',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])));
      case {'uAm','UAm'} 
         Z(:,:,placInd) = abs(complex(double(permute(ncread(filename,'URe',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])),...
                                      double(permute(ncread(filename,'UIm',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3]))));
      case {'uPh','UPh'}
          Z(:,:,placInd) = angle(complex(double(permute(ncread(filename,'URe',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])),...
                                          double(permute(ncread(filename,'UIm',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3]))));
      case {'v','V'} 
         Z(:,:,placInd) = complex(double(permute(ncread(filename,'VRe',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])),...
                                  double(permute(ncread(filename,'VIm',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])));
      case {'vAm','VAm'} 
         Z(:,:,placInd) = abs(complex(double(permute(ncread(filename,'VRe',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])),...
                                      double(permute(ncread(filename,'VIm',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3]))));
      case {'vPh','VPh'}
          Z(:,:,placInd) = angle(complex(double(permute(ncread(filename,'VRe',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])),...
                                          double(permute(ncread(filename,'VIm',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3]))));
      otherwise
         error(['The requested variable ',variable,' is not a variable in file ',filename,'.'])
   end
   
end

% Convert transports to velocities if requested: 
if ismember(variable,{'uRe','uIm','u','uAm','vRe','vIm','v','vAm'})
   wct = double(permute(ncread(filename,'wct',[ci(1) ri(1)],[numel(ci) numel(ri)]),[2 1]));
   Z = Z./max(wct,10); % divide transports by at least 10 m wct to prevent insanely high velocities.
end

if geo 
   if GlobalModel 
      % do nothing, x_or_lon is already lon and y_or_lat is already lat. 
   else 
      x_or_lon = double(permute(ncread(filename,'lon',[ci(1) ri(1)],[numel(ci) numel(ri)]),[2 1]));
      y_or_lat = double(permute(ncread(filename,'lat',[ci(1) ri(1)],[numel(ci) numel(ri)]),[2 1]));
   end
end

end