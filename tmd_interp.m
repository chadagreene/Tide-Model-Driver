function zi = tmd_interp(filename,variable,lati,loni,varargin)
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
% zi = tmd_interp(filename,variable,lati,loni) uses the NetCDF tide model
% specified by filename to interpolate a specified variable at the
% geographic coordinates lati,loni. The variable can be: 
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
%  * 'URe' real part of zonal transport
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
% zi = tmd_interp(...,'constituents',conList) specifies tidal constituents as a 
% cell array (e.g, {'m2','s2'}). If constituents are not specified, all constituents 
% from the model are returned. 
% 
% zi = tmd_interp(...,'coasts',MaskingMethod) specifies how coastal regions are masked. 
% Can be NaN, 'flexure', or 'unmask'. By default, MaskingMethod is NaN, meaning outputs 
% are set to NaN wherever a nearest-neighbor interpolation of the ocean indicates land. 
% The 'flexure' option scales tidal constituents by a predicted coefficient of tidal 
% deflection for ice shelf grounding zones. A third option, 'unmask', does not apply 
% any masking, which may be preferred close to coasts, where, for example, a tide gauge 
% may exist between land and ocean grid cells. 
% 
%% Examples 
% For examples type 
% 
%   tmd tmd_interp
% 
%% Author Info 
% This function was written by Chad A. Greene in 2022. 
% 
% See also tmd_data and tmd_predict. 

%% Input checks 

narginchk(4,Inf) 
assert(~isnumeric(variable),'Input variable must be a string.')
assert(islatlon(lati,loni),'Inputs lati,loni must be in the range of possible latitude and longitudes.') 
assert(isequal(size(lati),size(loni)),'Dimensions of input lati,loni coordinates must agree.') 
assert(contains(filename,'.nc'),'Input filename must end in .nc.')
assert(exist(filename,'file'),['Cannot find ',filename,'. Check the path and try again.'])

% Ensure the model file is TMD3.0 compatible: 
try 
   tmd_version = ncreadatt(filename,'/','tmd_version'); 
   assert(tmd_version>=3.0,[filename,' is not compatible with TMD3.0+.'])
catch 
   error([filename,' is not compatible with TMD3.0+.'])
end

%% Parse inputs 

MaskingMethod = 'nan'; % default set land areas to nan by nearest-neighbor interpolation of the mask.  
if nargin>4
   tmp = strncmpi(varargin,'coasts',5); 
   if any(tmp)
      MaskingMethod = lower(varargin{find(tmp)+1}); 
      if isnan(MaskingMethod)
         MaskingMethod = 'nan'; % making it a string for later. 
      end
      assert(ismember(MaskingMethod,{'nan','flex','flexure','unmask'}),'Coastal masking method can only be NaN, ''flexure'', or ''unmask''. Try again.')
      tmp(find(tmp)+1) = true; 
      varargin = varargin(~tmp); 
   end
end

%% Project if necessary 

proj4 = ncreadatt(filename,'mapping','spatial_proj4');

switch proj4 
   
   case '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +type=crs'
      [xi,yi] = tmd_ll2ps(lati,loni,70,-45,'N'); 
      
   case '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=-70 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +type=crs'
      [xi,yi] = tmd_ll2ps(lati,loni,-71,-70,'S'); 
      
   case '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +type=crs'
      [xi,yi] = tmd_ll2ps(lati,loni,70,0,'N'); 
      
   otherwise % assume global
      xi = loni; 
      yi = lati; 
      
      % Wrap to 360: 
      xi(xi<0) = xi(xi<0)+360; 
end

%% Load data

% switch lower(variable)
%    case {'h','hre','him'}
%       [zRe,x,y] = tmd_data(filename,'hRe','bounds',[xi(:) yi(:)],varargin{:}); 
%       zIm = tmd_data(filename,'hIm','bounds',[xi(:) yi(:)],varargin{:}); 
%       z = complex(zRe,zIm); 
%    case {'u','ure','uim'}
%       [zRe,x,y] = tmd_data(filename,'URe','bounds',[xi(:) yi(:)],varargin{:}); 
%       zIm = tmd_data(filename,'UIm','bounds',[xi(:) yi(:)],varargin{:}); 
%       z = complex(zRe,zIm); 
%    case {'v','vre','vim'} 
%       [zRe,x,y] = tmd_data(filename,'VRe','bounds',[xi(:) yi(:)],varargin{:}); 
%       zIm = tmd_data(filename,'VIm','bounds',[xi(:) yi(:)],varargin{:}); 
%       z = complex(zRe,zIm); 
%    case {'wct','flexure','mask','ham','hph','uam','uph','vam','vph'} 
%       [z,x,y] = tmd_data(filename,variable,'bounds',[xi(:) yi(:)],varargin{:}); 
%    otherwise 
%       error(['Unrecognized variable ',variable])
% end

switch lower(variable)
   case {'hph'}
      [z,x,y] = tmd_data(filename,'h','bounds',[xi(:) yi(:)],varargin{:}); 
   case {'uph'}
      [z,x,y] = tmd_data(filename,'U','bounds',[xi(:) yi(:)],varargin{:}); 
   case {'vph'} 
      [z,x,y] = tmd_data(filename,'V','bounds',[xi(:) yi(:)],varargin{:}); 
   case {'wct','flexure','mask'} 
      [z,x,y] = tmd_data(filename,variable,'bounds',[xi(:) yi(:)]); 
   otherwise 
      [z,x,y] = tmd_data(filename,variable,'bounds',[xi(:) yi(:)],varargin{:}); 
end

%% Convert transports to velocities if requested

% Lowercase first letter indicates velocity: 
% if ismember(variable(1),{'u','v'})
%    wct = tmd_data(filename,'wct','bounds',[xi(:) yi(:)],varargin{:}); 
%    z = z./max(wct,10); % divide transports by at least 10 m wct to prevent insanely high velocities.
% end

%% Interpolate 

switch variable 
   case 'mask'
      zi = interp2(x,y,z,xi,yi,'nearest'); 
   otherwise 
      zi = nan(size(xi,1),size(xi,2),size(z,3)); 
      for k = 1:size(zi,3) 
         zi(:,:,k) = interp2(x,y,z(:,:,k),xi,yi); 
      end
end

%% Convert complex numbers if requested 

switch lower(variable) 
%    case {'ham','uam','vam'} 
%       zi = abs(zi); 
   case {'hph','uph','vph'} 
      zi = angle(zi); 
   otherwise 
      % I think all the other cases are already taken care of. 

end

%% Mask out land areas 
 
if ~ismember(variable,{'mask','flexure','wct'})
   switch lower(MaskingMethod)
      case 'nan' % default 
         mask = tmd_data(filename,'mask','bounds',[xi(:) yi(:)]); 
         maski = interp2(x,y,uint8(mask~=1),xi,yi,'nearest')==0;
         for k = 1:size(zi,3)
            tmp = zi(:,:,k); 
            tmp(maski==0) = NaN; 
            zi(:,:,k) = tmp; 
         end
      case {'flex','flexure'} 
         flex = tmd_data(filename,'flexure','bounds',[xi(:) yi(:)]); 
         flexi = interp2(x,y,flex,xi,yi);
         for k = 1:size(zi,3)
            zi(:,:,k) = zi(:,:,k).*flexi; 
         end
      case 'unmask'
         % leave it alone 
      otherwise
         error('I do not know how this error is possible.') 
   end
end

end

%% Subfunctions 


function tf = islatlon(lat,lon)
% islatlon determines whether lat,lon is likely to represent geographical
% coordinates. 
% 
%% Citing Antarctic Mapping Tools
% This function was developed for Antarctic Mapping Tools for Matlab (AMT). If AMT is useful for you,
% please cite our paper: 
% 
% Greene, C. A., Gwyther, D. E., & Blankenship, D. D. Antarctic Mapping Tools for Matlab. 
% Computers & Geosciences. 104 (2017) pp.151-157. 
% http://dx.doi.org/10.1016/j.cageo.2016.08.003
% 
% @article{amt,
%   title={{Antarctic Mapping Tools for \textsc{Matlab}}},
%   author={Greene, Chad A and Gwyther, David E and Blankenship, Donald D},
%   journal={Computers \& Geosciences},
%   year={2017},
%   volume={104},
%   pages={151--157},
%   publisher={Elsevier}, 
%   doi={10.1016/j.cageo.2016.08.003}, 
%   url={http://www.sciencedirect.com/science/article/pii/S0098300416302163}
% }
%   
%% Syntax
% 
% tf = islatlon(lat,lon) returns true if all values in lat are numeric
% between -90 and 90 inclusive, and all values in lon are numeric between 
% -180 and 360 inclusive. 
% 
%% Example 1: A single location
% 
% islatlon(110,30)
%    = 0
% 
% because 110 is outside the bounds of latitude values. 
% 
%% Example 2: A grid
% 
% [lon,lat] = meshgrid(-180:180,90:-1:-90); 
% 
% islatlon(lat,lon)
%    = 1 
% 
% because all values in lat are between -90 and 90, and all values in lon
% are between -180 and 360.  What if it's really, really close? What if
% just one value violates these rules? 
% 
% lon(1) = -180.002; 
% 
% islatlon(lat,lon)
%    = 0
% 
%% Author Info
% This function was written by Chad A. Greene of the University of Texas at
% Austin's Institute for Geophysics (UTIG). http://www.chadagreene.com. 
% March 30, 2015. 
% 
% See also wrapTo180, wrapTo360, projfwd, and projinv.  

% Make sure there are two inputs: 
narginchk(2,2)

% Set default output: 
tf = true; 

%% If *any* inputs don't look like lat,lon, assume none of them are lat,lon. 

if ~isnumeric(lat)
    tf = false; 
    return
end

if ~isnumeric(lon)
    tf = false; 
    return
end
if any(abs(lat(:))>90)
    tf = false; 
    return
end

if any(lon(:)>360)
    tf = false; 
    return
end    

if any(lon(:)<-180)
    tf = false; 
end

end