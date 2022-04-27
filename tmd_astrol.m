function [s,h,p,N] = tmd_astrol(t)
% astrol computes the basic astronomical mean longitudes s, h, p, N.
% 
%% Syntax
% 
%  [s,h,p,N] = astrol(t)
% 
%% Description 
% 
% [s,h,p,N] = astrol(t) returns astronomical mean longitudes in degrees
% (wrapped 0 to 360). Input times t are datetime or datenum format. 
%   s: moon
%   h: sun
%   p: lunar perigee
%   N: ascending lunar node. Note N is not N', i.e. N is decreasing with time.
% 
%% Example
%
%        
%% Version History: 
% These formulae are for the period 1990 - 2010, and were derived
% by David Cartwright (personal comm., Nov. 1990).
% R. D. Ray    Dec. 1990
% Non-vectorized version. Re-make for Matlab by Lana Erofeeva, 2003
% TMD release 3.00: August 2018, changes by Chad A. Greene include: 
%   - Small changes for code readability. 
%   - Organized and edited this function header. 
%   - Changed the time variable name to t, to remove conflict with inbuilt time function.
%   - Replaced function return with function end.
%   - Removed statements like "if s<0, s = s + circle; end" because they appeared after mod calls, so could not be possible anyway.

%% Get mean longitudes: 

if isdatetime(t)
   t = datenum(t); 
end

t_mjd = t - datenum(1992,1,1) + 48622;

% Convert time from UTC in decimal MJD to ???
T = t_mjd - 51544.4993;

% Mean longitude of moon:
s = 218.3164 + 13.17639648 * T;

% Mean longitude of sun:
h = 280.4661 +  0.98564736 * T;

% Mean longitude of lunar perigee:
p =  83.3535 +  0.11140353 * T;

% Mean longitude of ascending lunar node:
N = 125.0445D0 -  0.05295377D0 * T;

%% Wrap phases: 

s = mod(s,360);
h = mod(h,360);
p = mod(p,360);
N = mod(N,360);

end