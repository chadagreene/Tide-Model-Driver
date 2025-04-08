function [s,h,p,N] = tmd_astrol(t, order=1, tt_ut1=0.0007)
% tmd_astrol computes the basic astronomical mean longitudes s, h, p, N, 
% following Foreman & Henry (1989). 
% 
%% Syntax
% 
%  [s,h,p,N] = tmd_astrol(t)
% 
%% Description 
% 
% [s,h,p,N] = tmd_astrol(t) returns astronomical mean longitudes in degrees
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

% convert to J2000 (dynamic time)
J2000 = datenum(2000, 1, 1, 12);
T = t - J2000 + tt_ut1;

if (order == 1)

% Mean longitude of moon:
s = 218.3164 + 13.17639648 * T;

% Mean longitude of sun:
h = 280.4661 +  0.98564736 * T;

% Mean longitude of lunar perigee:
p =  83.3535 +  0.11140353 * T;

% Mean longitude of ascending lunar node:
N = 125.0445D0 -  0.05295377D0 * T;

else

% Meeus arguments:
% Mean longitude of moon:
lunar_longitude = [218.3164591; 13.17639647754579; ...
   -9.9454632D-13; 3.8086292D-20; -8.6184958D-27];
s = polysum(T, lunar_longitude);

% Mean longitude of sun:
solar_longitude = [280.46645; 0.985647360164271; 2.2727347D-13];
h = polysum(T, solar_longitude);

% Mean longitude of lunar perigee:
lunar_perigee = [83.3532430; 0.11140352391786447; ...
   -7.7385418D-12; -2.5636086D-19; 2.95738836D-26];
p = polysum(T, lunar_perigee);

% Mean longitude of ascending lunar node:
lunar_node = [125.0445550; -0.052953762762491446; ...
   1.55628359D-12; 4.390675353D-20; -9.26940435D-27];
N = polysum(T, lunar_node);

end

%% Wrap phases: 

s = mod(s,360);
h = mod(h,360);
p = mod(p,360);
N = mod(N,360);

end

function S = polysum(T, polys)
   S = zeros(size(T));
   for o = 1:length(polys)
      S = S + polys(o)*power(T, o-1);
   end
end
