function [Z,x_or_lon,y_or_lat] = tmd_data(filename,variable,varargin) 
% 
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
%  h 
%  hRe
%  hIm
%  hAm
%  hPh
%  u 
%  uRe
%  uIm
%  uAm
%  uPh
%  U 
%  URe
%  UIm
%  UAm
%  UPh
%  v 
%  vRe
%  vIm
%  vAm
%  vPh
%  V 
%  VRe
%  VIm
%  VAm
%  VPh
%  wct
%  mask
%  flexure 
%
% [...] = tmd_data(...,'constituents',conList)
%
% [...] = tmd_data(...,'bounds',[xi yi])
% 
%% Example 
% 
% [hAm,x,y] = tmd_data('CATS2022_02-21.nc','hAm');  
% hPh = tmd_data('CATS_2022_02-21.nc','hPh'); 
% 
% h = imagesc(x,y,hPh); 
% cmocean phase 
% 
% 


%% Input parsing 

narginchk(2,Inf)

% Set defaults 
spatialSubset = false; 
conSubset = false; 
geo = false; % output geographic coordinates 
stride = Inf; 
NCons = 1; 
bounds = [[-Inf;Inf] [-Inf;Inf]];

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
         Z = double(permute(ncread(filename,variable,[ci(1) ri(1)],[numel(ci) numel(ri)]),[2 1]));

      case {'hre','him','ure','uim','vre','vim'}
         Z(:,:,placInd) = double(permute(ncread(filename,variable,[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3]));

      case 'h'
         Z(:,:,placInd) = complex(double(permute(ncread(filename,'hRe',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])),...
                                  double(permute(ncread(filename,'hIm',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])))/1000;
      case 'ham' 
         Z(:,:,placInd) = abs(complex(double(permute(ncread(filename,'hRe',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3])),...
                                      double(permute(ncread(filename,'hIm',[ci(1) ri(1) conind],[numel(ci) numel(ri) stride]),[2 1 3]))))/1000;
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

if ismember(variable,{'uRe','uIm','u','uAm','uPh','vRe','vIm','v','vAm','vPh'})
	wct = double(permute(ncread(filename,'wct',[ri(1) ci(1)],[numel(ci) numel(ri)]),[2 1])); 
   Z = Z./wct; 
end

if geo 
   x_or_lon = double(permute(ncread(filename,'lon',[ci(1) ri(1)],[numel(ci) numel(ri)]),[2 1]));
   y_or_lat = double(permute(ncread(filename,'lat',[ci(1) ri(1)],[numel(ci) numel(ri)]),[2 1]));
end

end