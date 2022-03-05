function z = tmd_predict(filename,lat,lon,time,ptype,varargin)
% 
% 
%% Syntax 
% 
%  z = tmd_predict(filename,lat,lon,time)
%  z = tmd_predict(filename,lat,lon,time,ptype)
%  z = tmd_predict(filename,lat,lon,time,ptype,'constituents',conList)
%  z = tmd_predict(...,'InferMinor',false)
%  z = tmd_predict(...,'IceShelfFlexure')
% 
%% 



%% Input checks 
% 

narginchk(4,Inf) 
assert(~isnnumeric(filename),'Model filename must be a string.') 
assert(isequal(size(lat),size(lon)),'Dimensions of lat and lon must agree.') 


% Check ptype and define 'z' as default if user has not defined ptype: 
if nargin>4 
    assert(~isnumeric(ptype),'Input error: Solution type ptype must be a string.') 
    assert(ismember(ptype,{'z','u','U','v','V'}),['Input error: Unrecognized ptype ',ptype])
else 
    ptype = 'z'; 
end

% 
% Three solution types are possible depending on size of inputs 
% t,lat,lon: 
%
%      "Time series" t(N,1),lon(1,1),lat(1,1)
%      "Drift Track" t(N,1),lon(N,1),lat(N,1)
%      "Map"         t(1,1),lon(N,M),lat(N,M) 
% 
%% Syntax 
% 
%  TS = tmd_tide_pred(Model,t,lat,lon,ptype)
%  TS = tmd_tide_pred(Model,t,lat,lon,ptype,Cid)
%  [TS,conList] = tmd_tide_pred(...)
% 
%% Description 
% 
% TS = tmd_tide_pred(Model,t,lat,lon,ptype) returns a tide prediction of type ptype at time
% t for the geolocation(s) lat,lon using a model specified by Model. Models are in OTIS format 
% as follows: 
% 
%       Model: a text file name for a tidal model (e.g. 'Model_CATS2008'), consisting of lines
%             <elevation file name> (e.g., hf.CATS2008.out)
%             <transport file name> (e.g., uv.CATS2008.out)
%             <grid file name>      (e.g., grid_CATS2008b_km) 
%             <function to convert lat,lon to x,y> (e.g., xy_ll_CATS2008) Only include the fourth
%                line for models in cartesian grid (in km) coordinates such as CATS. 
% 
%           t: vector of times expressed in serial days, Matlab's datenum format. See: help datenum. 
% 
%     lat,lon: geo coordinates in decimal degrees. 
% 
%       ptype: solution type can be: 
%                'z' - elevation (m) (default) 
%                'u','v' - velocities (cm/s)
%                'U','V' - transports (m^2/s)
%  
% TS = tmd_tide_pred(Model,t,lat,lon,ptype,Cid) specifies numerical indices of constituents.
% to solve.  By default, all available constituents are solved and minor constituents are inferred.
% Numerical indices of available constituents depend on the Model you are using.  To get a list of 
% available constutents and their numerical indices, use the rd_con function.  For example, 
% conList = rd_con('h_tpxo7.2_load') gives 13 available constituents for the vertical solution of 
% the TPX 7.2 load tide model.  If Cid is specified, minor constituents are NOT inferred.  Minor
% constituents are only inferred if Cid is not specified.  

%% Parse inputs: 


% Check constituent lists: 
solveall = true; % solve all constituents by default. 
if nargin==6
    assert(isnumeric(Cid)==1,'Input error: Constituent IDs Cid must be numeric. Check the table the header of tmd_tide_pred by typing: help tmd_tide_pred.') 
    if ~isempty(Cid)
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

%% Get model parameters and prepare inputs: 

% Initialize TS variable: 
TS=[];

% Get model name and Grid name, check their validity: 
if strcmp(ptype,'z') 
    [ModName,GridName,~] = rdModFile(Model,1);
else
    [ModName,GridName,~] = rdModFile(Model,2);
end

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

% if solveall
%     [amp,pha,~,conList]=tmd_extract_HC(Model,lat,lon,ptype);
% else
%     [amp,pha,~,conList]=tmd_extract_HC(Model,lat,lon,ptype,Cid);
% end
% 
% cph = -1i*pha*pi/180;
% hc = squeeze(amp.*exp(cph));

hc = squeeze(tmd_interp(filename,'h',lat,lon));


[nc,~]=size(conList); % nc here is # of ALL constituents in the model
if strcmp(conList(1:4),'stop')
    TS=NaN;
    return;
end

time=t-datenum(1992,1,1); % corresponds to 48622mjd

if solveall % solve all constituents: 
    switch solntype 
        case 'time series' 
            TS=harp1(time,hc,conList);
            dh = InferMinor(hc,conList,t);
            TS=TS+dh;
  
        case 'drift track'
            TS = NaN(1,numel(t));  
            for k=1:length(t)
                TS(k)=harp1(time(k),hc(:,k),conList);
                dh = InferMinor(hc(:,k),conList,t(k));
                TS(k)=TS(k)+dh;
            end
            
        case 'map' 
            hci = NaN(size(hc,2),size(hc,3),size(hc,1)); 
            for k=1:nc
                hci(:,:,k)=hc(k,:,:);
            end
            TS=harp(time,hci,conList);
            dh=InferMinor(hc,conList,t);
            if ndims(hci)~=ndims(hc)
                dh=dh';
            end
            if ~isequal(size(TS),size(dh)) 
               dh = dh'; 
            end
            TS=TS+dh;
    end
 
else % This is if Cid is defined by user: 

    % It seems odd to change inputs that the user specified explicitly,
    % but I don't want to break anyone's code by changing these next two lines: 
    Cid(Cid<1)=1;
    Cid(Cid>nc)=nc;

    switch solntype
        case 'time series'
            TS=harp1(time,hc,conList(Cid,:));

        case 'drift track' 
            TS = NaN(1,numel(t));  
            for k=1:length(time)
                TS(k)=harp1(time(k),hc(:,k),conList(Cid,:));
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
            TS=harp(time,hci,conList(Cid,:));
    end
   
end

end