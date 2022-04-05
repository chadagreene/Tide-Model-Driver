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
assert(~isnumeric(filename),'Model filename must be a string.') 
assert(isequal(size(lat),size(lon)),'Dimensions of lat and lon must agree.') 


% Check ptype and define 'z' as default if user has not defined ptype: 
if nargin>4 
    assert(~isnumeric(ptype),'Input error: Solution type ptype must be a string.') 
    assert(ismember(ptype,{'h','z','u','U','v','V'}),['Input error: Unrecognized ptype ',ptype])
    if strcmpi(ptype,'z')
       ptype='h'; % allows z or h as input. 
    end
else 
    ptype = 'h'; 
end

% Check constituent lists: 
solveall = true; % solve all constituents by default. 
conList = strsplit(ncreadatt(filename,'cons','long_name')); 
%char(pad(conList,4,'right'))

if nargin>5
   tmp = strncmpi(varargin,'constituents',3); 
   if any(tmp)
      conList = varargin{find(tmp)+1}; 
      solveall = false; 
   end
end
   
% Parse solution type: 
if isscalar(lat)
    solntype = 'time series'; 
else
   if isequal(size(lat),size(lon),size(t)) & isvector(t)
      solntype = 'drift track'; 
   else
      solntype = 'map'; 
   end
end
assert(isequal(size(lat),size(lon)),'Inputs lat and lon must have the same dimensions.') 
   
%%

if isdatetime(t)
   t = datenum(t); 
end

% Columnate input coordinates for consistent behavior: 
InputGridSize = size(lat); 
lat = reshape(lat,[],1); 
lon = reshape(lon,[],1); 

%InputTimeSize = size(t); 
t = reshape(t,[],1); 

% % Initialize TS variable: 
% z=[];

%% Solve


if solveall
   hc = tmd_interp(filename,ptype,lat,lon);
else
   hc = tmd_interp(filename,ptype,lat,lon,'constituents',conList);
end

% Switch dimensions bc we've already forced dim2 to be singleton by colunating lat and lon: 
% switch solntype
%    case 'map'
%       hc = permute(hc,[3 2 1]); % Dimensions are now M latpoints by Ncons
%    otherwise
%       hc = permute(hc,[3 1 2]); % Dimensions are now M latpoints by Ncons
% end
%    
hc = permute(hc,[3 1 2]);


[astrol_s,astrol_h,astrol_p,astrol_N] = tmd_astrol(t);

% If drift track or time series: 
d_minor = tmd_InferMinor(hc,conList,t,astrol_s,astrol_h,astrol_p,astrol_N); 

if solveall % solve all constituents: 
    switch solntype 
        case 'time series' 
            z = tmd_harp1(t,hc,conList) + tmd_InferMinor(hc,conList,t);
  
        case 'drift track'
            z = NaN(1,numel(t));  
            for k=1:length(t)
                z(k) = tmd_harp1(t(k),hc(:,k),conList) + tmd_InferMinor(hc(:,k),conList,t(k));
            end
            
        case 'map' 
            hci = NaN(size(hc,2),size(hc,3),size(hc,1)); 
            for k=1:length(t) 
                hci(:,:,k)=hc(k,:,:);
            end
            z=tmd_harp(t,hci,conList);
            dh=tmd_InferMinor(hc,conList,t);
            if ndims(hci)~=ndims(hc)
                dh=dh';
            end
            if ~isequal(size(z),size(dh)) 
               dh = dh'; 
            end
            z=z+dh;
    end
 
else % This is if Cid is defined by user: 


    switch solntype
        case 'time series'
            z=tmd_harp1(t,hc,conList(Cid,:));

        case 'drift track' 
            z = NaN(1,numel(t));  
            for k=1:length(t)
                z(k)=tmd_harp1(t(k),hc(:,k),conList(Cid,:));
            end

        case 'map'
            if length(Cid)>1
                hci = NaN(size(hc,2),size(hc,3),size(hc,1)); 
                for k=1:length(Cid)
                    hci(:,:,k)=hc(k,:,:);
                end
            else
                hci=hc;
            end
            z=tmd_harp(t,hci,conList(Cid,:));
    end
   
end

%% Clean up 

z = reshape(z,InputGridSize(1),InputGridSize(2),[]); 

end