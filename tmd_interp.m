function zi = tmd_interp(filename,variable,lati,loni,varargin)
% 
% 
%% Syntax 
% 
%  zi = tmd_interp(filename,variable,lati,loni)
%  zi = tmd_interp(...,'constituents',conList) 
% 


%% Input checks 

narginchk(4,Inf) 
assert(islatlon(lati,loni),'Inputs lati,loni must be in the range of possible latitude and longitudes.') 
assert(exist(filename,'file')==2,['Cannot find ',filename])
assert(isequal(size(lati),size(loni)),'Dimensions of input coordinates must agree.') 

%% Project if necessary 

proj4 = ncreadatt(filename,'mapping','spatial_proj4');

switch proj4 
   case '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=-70 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +type=crs'
      [xi,yi] = mapll(lati,loni,-71,-70,'S'); 
   otherwise 
      xi = loni; 
      yi = lati; 
end

%% Load data

switch lower(variable)
   case {'h','hre','him'}
      [zRe,x,y] = tmd_data(filename,'hRe','bounds',[xi(:) yi(:)],varargin{:}); 
      zIm = tmd_data(filename,'hIm','bounds',[xi(:) yi(:)],varargin{:}); 
      z = complex(zRe,zIm); 
   case {'u','ure','uim'}
      [zRe,x,y] = tmd_data(filename,'URe','bounds',[xi(:) yi(:)],varargin{:}); 
      zIm = tmd_data(filename,'UIm','bounds',[xi(:) yi(:)],varargin{:}); 
      z = complex(zRe,zIm); 
   case {'v','vre','vim'} 
      [zRe,x,y] = tmd_data(filename,'VRe','bounds',[xi(:) yi(:)],varargin{:}); 
      zIm = tmd_data(filename,'VIm','bounds',[xi(:) yi(:)],varargin{:}); 
      z = complex(zRe,zIm); 
   case {'wct','flexure','mask','ham','hph','uam','uph','vam','vph'} 
      [z,x,y] = tmd_data(filename,variable,'bounds',[xi(:) yi(:)],varargin{:}); 
   otherwise 
      error(['Unrecognized variable ',variable])
end

%% Convert transports to velocities if requested

% Lowercase first letter indicates velocity: 
if ismember(variable(1),{'u','v'})
   wct = tmd_data(filename,'wct','bounds',[xi(:) yi(:)],varargin{:}); 
   z = z./wct; 
end


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
   case {'ham','uam','vam'} 
      zi = abs(zi); 
   case {'hph','uph','vph'} 
      zi = -angle(zi); 
   otherwise 
      % I think all the other cases are already taken care of. 

end

%% NaN out land areas 
 
if ~ismember(variable,{'mask','flexure'})
   mask = tmd_data(filename,'mask','bounds',[xi(:) yi(:)]); 
   maski = interp2(x,y,mask,xi,yi,'nearest'); 
   for k = 1:size(zi,3)
      tmp = zi(:,:,k); 
      tmp(maski==0) = NaN; 
      zi(:,:,k) = tmp; 
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