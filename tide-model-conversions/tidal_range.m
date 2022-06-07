function R = tidal_range(filename_or_hc,conList,mask)
% tidal_range calculates the peak-to-peak tidal height range. (May take
% hours.) 
% 
%% Syntax
% 
%  R = tidal_range(filename)
%  R = tidal_range(hRe,hIm,conList,mask)
% 
%% Description
% 
% R = tidal_range(filename) returns the peak-to-peak tidal height range for
% a given TMD3.0 tide model file. Tidal range is calculated as the maximum
% minus minimum predicted values over a 14 day period, calculated at 30
% minute temporal resolution, including major and minor constituents. 
% 
% R = tidal_range(hRe,hIm,conList,mask)
%
%% Example 
% 
% R = tidal_range('/Users/cgreene/Downloads/CATS2008/CATS2008_update_2022-04-22.nc');
% 
%% Author Info 
% Written by Chad A. Greene, June 2022. 

%% 

% Two weeks of 30 minute timesteps: 
t = (datenum(2000,1,1):1/48:datenum(2000,1,29))';

%% Load data

if isnumeric(filename_or_hc)
   assert(nargin==3,'If the first input is numeric, the inputs must be hRe, hIm, conList, mask.')
   hc = filename_or_hc; 
else
   conList = strsplit(ncreadatt(filename_or_hc,'cons','long_name')); 
   mask = tmd_data(filename_or_hc,'mask');
   hc = tmd_data(filename_or_hc,'h');
end

hc = cube2rect(hc,mask); % cube2rect is in Climate Data Toolbox for Matlab

[astrol_s,astrol_h,astrol_p,astrol_N] = tmd_astrol(t);
[ispec,~,ph,omega,~,~] = tmd_constit(conList);

%% Solve

% Start with an assumption of zero tides: 
z_min = zeros(size(hc,2),1);
z_max = zeros(size(hc,2),1);

% Solve for each timestep:
w = waitbar(0,'Calculating tidal range (may take a few hours).');
for k = 1:length(t)

   % Major constiuents: 
   hhat = tmd_harp(t(k),hc,conList,astrol_p(k),astrol_N(k),ispec,ph,omega);

   % Minor constiuents:
   d_minor = tmd_InferMinor(hc,conList,t(k),astrol_s(k),astrol_h(k),astrol_p(k),astrol_N(k)); 

   % Just calculate min and max values (thus far) bc all we need is the total range:  
   z_max = max(z_max, d_minor + hhat); 
   z_min = min(z_min, d_minor + hhat); 
   
   waitbar(k/length(t),w,'Calculating tidal range (may take a few hours).')
end

close(w)

R = rect2cube((z_max-z_min)',mask); % rect2cube is in Climate Data Toolbox for Matlab
R(isnan(R)) = 0;

end