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
    if isscalar(t)
        solntype = 'map'; 
    else
        solntype = 'drift track'; 
        assert(isvector(t),'Input error: Conditions for time series, drift track, or map solution are not met. Check dimensions of lat, lon, and t.')
        assert(isvector(lat),'Input error: Conditions for time series, drift track, or map solution are not met. Check dimensions of lat, lon, and t.')
    end
end

%%

if isdatetime(t)
   t = datenum(t); 
end

%% Get model parameters and prepare inputs: 

% Initialize TS variable: 
z=[];

% Columnate inputs for consistent behavior: 
switch solntype
    case 'time series'
        t = t(:); 
        
    case 'drift track'
        t = t(:); 
        lat = lat(:); 
        lon = lon(:); 
        
    case 'map' 
        % No columnation necessary. 
end

%% Solve

[astrol_s,astrol_h,astrol_p,astrol_N] = tmd_astrol(t);

if solveall
   hc = squeeze(tmd_interp(filename,ptype,lat,lon));
else
   hc = squeeze(tmd_interp(filename,ptype,lat,lon,'constituents',conList));
end
   
NCons = length(hc); 

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
            for k=1:NCons
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

end