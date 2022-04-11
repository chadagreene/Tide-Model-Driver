function [Z,x_or_lon,y_or_lat] = tmd_data(filename,variable,varargin) 
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
%  * 'flexure' ice shelf flexure coefficient from a linear elastic model applied to BedMachine ice thickness (can slightly exceed 1). 
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
%% Example 
% 
% [hAm,x,y] = tmd_data('CATS2022_02-21.nc','hAm');  
% hPh = tmd_data('CATS2022_02-21.nc','hPh'); 
% 
% h = imagesc(x,y,hPh(:,:,1)); 
% cmocean phase % set colormap
% 
%% Author Info 
% This function was written by Chad A. Greene in 2022. 
% 
% See also tmd_interp and tmd_predict. 

%% Input parsing 

narginchk(2,Inf)

% Set defaults 
spatialSubset = false; 
conSubset = false; 
geo = false; % output geographic coordinates 
stride = Inf; 
NCons = 1; 
bounds = [[-Inf;Inf] [-Inf;Inf]];

assert(contains(filename,'.nc'),'Input filename must end in .nc.')
assert(exist(filename,'file'),['Cannot find ',filename,'. Check the path and try again.'])

if nargin>2 
   
   tmp = strncmpi(varargin,'bounds',3); 
   if any(tmp) 
      spatialSubset = true; 
      bounds = varargin{find(tmp)+1}; 
      assert(size(bounds,2)==2,'Spatial bounds must be Mx2 in the form [xi(:) yi(:)] or [loni(:) lati(:)].')
   end 
   
   tmp = strncmpi(varargin,'constituents',3); 
   if any(tmp) 
      conSubset = true; 
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

x_or_lon = double(ncread(filename,'x')); 
y_or_lat = double(ncread(filename,'y')); 
cons = strsplit(ncreadatt(filename,'cons','long_name')); 

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

   switch lower(variable)
      case 'mask'
         Z = logical(permute(ncread(filename,variable,[ci(1) ri(1)],[numel(ci) numel(ri)]),[2 1]));
         
      case {'flexure','wct'}
         Z = double(permute(ncread(filename,variable,[ci(1) ri(1)],[numel(ci) numel(ri)]),[2 1]))/100;

      case {'hre','him','ure','uim','vre','vim'}
         Z(:,:,placInd) = double(permute(ncread(filename,variable,[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3]));

      case 'h'
         Z(:,:,placInd) = complex(double(permute(ncread(filename,'hRe',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])),...
                                  double(permute(ncread(filename,'hIm',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])));
      case 'ham' 
         Z(:,:,placInd) = abs(complex(double(permute(ncread(filename,'hRe',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])),...
                                      double(permute(ncread(filename,'hIm',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3]))));
      case 'hph'
          Z(:,:,placInd) = -angle(complex(double(permute(ncread(filename,'hRe',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])),...
                                          double(permute(ncread(filename,'hIm',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3]))));
      case 'u'
         Z(:,:,placInd) = complex(double(permute(ncread(filename,'URe',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])),...
                                  double(permute(ncread(filename,'UIm',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])));
      case 'uam' 
         Z(:,:,placInd) = abs(complex(double(permute(ncread(filename,'URe',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])),...
                                      double(permute(ncread(filename,'UIm',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3]))));
      case 'uph'
          Z(:,:,placInd) = -angle(complex(double(permute(ncread(filename,'URe',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])),...
                                          double(permute(ncread(filename,'UIm',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3]))));
      case 'v' 
         Z(:,:,placInd) = complex(double(permute(ncread(filename,'VRe',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])),...
                                  double(permute(ncread(filename,'VIm',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])));
      case 'vam' 
         Z(:,:,placInd) = abs(complex(double(permute(ncread(filename,'VRe',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])),...
                                      double(permute(ncread(filename,'VIm',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3]))));
      case 'vph'
          Z(:,:,placInd) = -angle(complex(double(permute(ncread(filename,'VRe',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])),...
                                          double(permute(ncread(filename,'VIm',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3]))));
      otherwise
         error(['The requested variable ',variable,' is not a variable in file ',filename,'.'])
   end
   
end

% Scale variables to SI meters, m/s 
if ismember(lower(variable),{'h','hre','him','ham'})
   Z = Z/1000; 
end

if ismember(variable,{'uRe','uIm','u','uAm','uPh','vRe','vIm','v','vAm','vPh'})
	wct = double(permute(ncread(filename,'wct',[ri(1) ci(1)],[numel(ci) numel(ri)]),[2 1])); 
   Z = Z./wct; 
end

if geo 
   x_or_lon = double(permute(ncread(filename,'lon',[ci(1) ri(1)],[numel(ci) numel(ri)]),[2 1]));
   y_or_lat = double(permute(ncread(filename,'lat',[ci(1) ri(1)],[numel(ci) numel(ri)]),[2 1]));
end

end