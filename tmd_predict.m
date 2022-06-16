function z = tmd_predict(filename,lat,lon,t,ptype,varargin)
% tmd_predict constructs tidal predictions from the complex constituent 
% coefficients in consolidated NetCDF tide model files. 
% 
%% Syntax 
% 
%  z = tmd_predict(filename,lat,lon,time)
%  z = tmd_predict(filename,lat,lon,time,ptype)
%  z = tmd_predict(filename,lat,lon,time,ptype,'constituents',conList)
%  z = tmd_predict(...,,'coasts',MaskingMethod)
%  z = tmd_predict(...,'InferMinor',true_or_false)
% 
%% Description 
% 
% z = tmd_predict(filename,lat,lon,time) predicts tide heights a the
% specified lat,lon and time, using the TMD3.0 compatible consolidated
% NetCDF tide model data file. Location(s) lat,lon are decimal degrees, and
% can be scalars, vectors, or MxN arrays. Input time can be datetime or
% MATLAB's datenum format, and must be a scalar or a 1D vector. 
%
% z = tmd_predict(filename,lat,lon,time,ptype) specifies a solution type.
% If ptype is not specified, 'h' is the assumed ptype. Note the ptype is
% case sensitive! Options for ptype are: 
%
%  * 'h' tidal height (m) (default)
%  * 'u' zonal current velocity, positive pointing east (m/s) 
%  * 'v' meridional current velocity, positive pointing north (m/s) 
%  * 'U' zonal transport (m^2/s) 
%  * 'V' meridional height (m^2/s) 
%  
% z = tmd_predict(filename,lat,lon,time,ptype,'constituents',conList) 
% specifies tidal constituents as a cell array (e.g, {'m2','s2'}). If 
% constituents are not specified, all constituents from the model are returned. 
% 
% z = tmd_predict(...,,'coasts',MaskingMethod) specifies how coastal regions are masked. 
% Can be NaN, 'flexure', or 'unmask'. By default, MaskingMethod is NaN, meaning outputs 
% are set to NaN wherever a nearest-neighbor interpolation of the ocean indicates land. 
% The 'flexure' option scales tidal constituents by a predicted coefficient of tidal 
% deflection for ice shelf grounding zones. A third option, 'unmask', does not apply 
% any masking, which may be preferred close to coasts, where, for example, a tide gauge 
% may exist between land and ocean grid cells. 
% 
% z = tmd_predict(...,'InferMinor',true_or_false) specifies whether minor
% constituents should be inferred (true or false). By default, minor constituents are
% inferred unless constituents are specified. 
%
%% Examples 
% For examples type 
% 
%   tmd tmd_predict
% 
%% Author Info 
% This function was written by Chad A. Greene in 2022. 
% 
% See also tmd_data and tmd_predict. 

%% Input checks 
% 

narginchk(4,Inf) 
assert(exist(filename,'file'),['Cannot find model file ',filename,'. Check the file name and make sure it is in the Matlab path or specify the full path to the .nc file.'] ) 
try 
   tmd_version = ncreadatt(filename,'/','tmd_version'); 
   assert(tmd_version>=3.0,[filename,' is not compatible with TMD3.0+.'])
catch 
   error([filename,' is not compatible with TMD3.0+.'])
end
assert(isequal(size(lat),size(lon)),'Dimensions of lat and lon must agree.') 

% Set defaults: 
InferMinorConstituents = true; 
MapSolution = false; 
MaskingMethod = 'nan'; % default set land areas to nan by nearest-neighbor interpolation of the mask.  
conList = strsplit(ncreadatt(filename,'constituents','constituent_order')); 

if nargin<5
   ptype = 'h'; 
end

% Check ptype and define 'z' as default if user has not defined ptype: 
if nargin>4 
    assert(~isnumeric(ptype),'Input error: Solution type ptype must be a string.') 
    assert(ismember(ptype,{'h','z','u','U','v','V'}),'Input error: The fifth input to tmd_predict must be one of the following ptypes: ''h'',''z'',''u'',''U'',''v'',''V''.')
    if strcmpi(ptype,'z')
       ptype='h'; % allows z or h as input. 
    end

end

if nargin>5
   % Which constituents do we solve for? (All by default) 
   tmp = strncmpi(varargin,'constituents',3); 
   if any(tmp)
      conList = varargin{find(tmp)+1}; 
      InferMinorConstituents = false; % If user requests specific constituents, assume they don't want to infer minor constituents. 
   end
   
   % Override InferMinor default behavior if user has a specific preference: 
   tmp = strncmpi(varargin,'inferminor',6); 
   if any(tmp)
      InferMinorConstituents = varargin{find(tmp)+1}; 
      assert(islogical(InferMinorConstituents),'Preference for InferMinor must be either true or false.') 
   end
   
   tmp = strncmpi(varargin,'coasts',5); 
   if any(tmp)
      MaskingMethod = lower(varargin{find(tmp)+1}); 
      if isnan(MaskingMethod)
         MaskingMethod = 'nan'; % making it a string for later. 
      end
      assert(ismember(MaskingMethod,{'nan','flex','flexure','unmask'}),'Coastal masking method can only be NaN, ''flexure'', or ''unmask''. Try again.')
   end

end

if ~InferMinorConstituents
   d_minor = 0; % minor constituent correction. 
end
   
% Parse solution type: 
if numel(lat)>1 & ~isequal(size(lat),size(lon),size(t))
   MapSolution = true; 
end
   
%%

if isdatetime(t)
   t = datenum(t); 
end

% Columnate input coordinates for consistent behavior: 
InputGridSize = size(lat); 
lat = reshape(lat,[],1); 
lon = reshape(lon,[],1); 

InputTimeSize = size(t); 
t = reshape(t,[],1); 

%% Load data

hc = tmd_interp(filename,ptype,lat,lon,'constituents',conList,'coasts',MaskingMethod);
hc = permute(hc,[3 1 2]); % puts constituents in first column. 

[s,h,p,N] = tmd_astrol(t);
[~,~,ph,omega,~] = tmd_constit(conList);

%% Predict tides

if MapSolution
   
   % Preallocate z: 
   z = nan(numel(lat),numel(t));
   
   isf = all(isfinite(hc),1); % continent or out-of-model-domain values are nan.  
   
   % Solve for each timestep:
   for k = 1:length(t)
      
      % Major constiuents: 
      hhat = tmd_harp(t(k),hc(:,isf),conList,p(k),N(k),ph,omega);

      % Minor constiuents:
      if InferMinorConstituents
         d_minor = tmd_InferMinor(hc(:,isf),conList,t(k),s(k),h(k),p(k),N(k)); 
      end
   
      z(isf,k) = d_minor + hhat; 
   end

else % Single-location time series or drift track  

   % Major constituents: 
   hhat = tmd_harp(t,hc,conList,p,N,ph,omega);
   
   if InferMinorConstituents
      d_minor = tmd_InferMinor(hc,conList,t,s,h,p,N); 
   end
   z = d_minor + hhat; 
end

%% Clean up 

if isequal(InputGridSize,[1 1])
   z = reshape(z,InputTimeSize); 
else
   z = reshape(z,InputGridSize(1),InputGridSize(2),[]); 
end

end