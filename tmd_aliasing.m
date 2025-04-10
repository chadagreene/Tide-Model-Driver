function period = tmd_aliasing(constituents, sampling, order=1, tt_ut1=0.0007)
% tmd_aliasing returns the tidal aliasing for a repeat period.
% 
%% Syntax
% 
%  period = tmd_aliasing(constituents, sampling)
% 
%% Description 
% 
% period = tmd_conlist(constituents) takes input constituents 
% as cell array (1xN constituents) and input sampling repeat period in seconds.
% Option order is the minimum polynomial order for calculating the astronomical
% mean longitudes. Option tt_ut1 is the delta time between Terrestrial Time (TT)
% and Universal Time (UT1). Output aliasing periods are in seconds.
% 
%% Author Info
% TMD was written by Lana Erofeeva (serofeev@coas.oregonstate.edu) and is maintained by
% Laurie Padman (padman@esr.org). 

% get the angular frequency for tidal constituents
omega = tmd_frequency(constituents, order=order, tt_ut1=tt_ut1);
% convert to cycles per second
f = omega/(2.0*pi);
% calculate the sampling frequency
fs = 1.0/sampling;
% calculate the aliasing frequency
fa = abs(f - fs.*round(f./fs));
% convert aliasing frequency to period
period = 1./fa;

end
