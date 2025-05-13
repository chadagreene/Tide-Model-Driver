function [s,h,p,N] = tmd_astrol(t)
% tmd_astrol computes the basic astronomical mean longitudes s, h, p, N, 
% following Foreman & Henry (1989). 
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
% Some guidance from Laurie Padman (because I asked for an explanation, June 2022):
% "The variable names are exactly as used in the Foreman & Henry (1989, FH89)
% paper. The statement about N not being N' can be explained looking at the 
% first page (col 2) of FH89; there, they use n’ as the negative of the longitude 
% of the moon’s ascending node. This statement in the code is telling us that
% ‘n’ in TMD is *not* the negative’d value."    
% 
%% References 
% 
% Foreman, M. G. G., & Henry, R. F. (1989). The harmonic analysis of tidal 
% model time series. Advances in water resources, 12(3), 109-120.
% https://doi.org/10.1016/0309-1708(89)90017-1
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

% Convert time:
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