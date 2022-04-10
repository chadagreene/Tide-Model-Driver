function z = tmd_predict(filename,lat,lon,t,ptype,varargin)
% 
% 
%% Syntax 
% 
%  z = tmd_predict(filename,lat,lon,time)
%  z = tmd_predict(filename,lat,lon,time,ptype)
%  z = tmd_predict(filename,lat,lon,time,ptype,'constituents',conList)
%  z = tmd_predict(...,'InferMinor',false)
%  z = tmd_predict(...,'flexure')
% 
%% 



%% Input checks 
% 

narginchk(4,Inf) 
assert(exist(filename,'file'),['Cannot find model file ',filename,'. Check the file name and make sure it is in the Matlab path or specify the full path to the .nc file.'] ) 
assert(isequal(size(lat),size(lon)),'Dimensions of lat and lon must agree.') 

% Set defaults: 
InferMinorConstituents = true; 
MapSolution = false; 
flexure = false; % ice shelf flexure zone 
conList = strsplit(ncreadatt(filename,'cons','long_name')); 
if nargin<5
   ptype = 'h'; 
end

% Check ptype and define 'z' as default if user has not defined ptype: 
if nargin>4 
    assert(~isnumeric(ptype),'Input error: Solution type ptype must be a string.') 
    assert(ismember(ptype,{'h','z','u','U','v','V'}),['Input error: Unrecognized ptype ',ptype])
    if strcmpi(ptype,'z')
       ptype='h'; % allows z or h as input. 
    end

end

if nargin>5
   tmp = strncmpi(varargin,'constituents',3); 
   if any(tmp)
      conList = varargin{find(tmp)+1}; 
      InferMinorConstituents = false; 
      d_minor = 0; % minor constituent correction. 
   end
   
   if any(strncmpi(varargin,'flexure',4)) 
      flexure = true; 
   end
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

hc = tmd_interp(filename,ptype,lat,lon,'constituents',conList);
hc = permute(hc,[3 1 2]); % puts constituents in first column. 

if flexure 
   flex = tmd_interp(filename,'flexure',lat,lon); 
end

%%

[astrol_s,astrol_h,astrol_p,astrol_N] = tmd_astrol(t);

if MapSolution
   
   % Preallocate z: 
   z = nan(numel(lat),numel(t));
   
   isf = all(isfinite(hc),1); % continent or out-of-model-domain values are nan.  
   
   % Solve for each timestep:
   for k = 1:length(t)
      
      % Major constiuents: 
      hhat = tmd_harp(t(k),hc(:,isf),conList,astrol_p(k),astrol_N(k));

      % Minor constiuents:
      if InferMinorConstituents
         d_minor = tmd_InferMinor(hc(:,isf),conList,t(k),astrol_s(k),astrol_h(k),astrol_p(k),astrol_N(k)); 
      end
   
      z(isf,k) = d_minor + hhat; 
   end

else % Single-location time series or drift track  

   % Major constituents: 
   hhat = tmd_harp(t,hc,conList,astrol_p,astrol_N);
   
   if InferMinorConstituents
      d_minor = tmd_InferMinor(hc,conList,t,astrol_s,astrol_h,astrol_p,astrol_N); 
   end
   z = d_minor + hhat; 
end

% Account for ice shelf flexure in the grounding zone: 
if flexure 
   z = z.*flex; 
end

%% Clean up 

if isequal(InputGridSize,[1 1])
   z = reshape(z,InputTimeSize); 
else
   z = reshape(z,InputGridSize(1),InputGridSize(2),[]); 
end

end