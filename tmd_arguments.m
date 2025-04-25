function [pu, pf, G] = tmd_arguments(t, constituents, order=1, tt_ut1=0.0007)
% Calculates the nodal arguments for the tidal constituents. 
% 
%% Syntax 
% 
%  [pu,pf,G] = tmd_arguments(t, constituents) 
%
%% Description 
% 
% [pu,pf,G] = tmd_arguments(t, constituents) takes input constituents 
% as cell array (1xN constituents). Option order is the minimum polynomial order
% for calculating the astronomical mean longitudes. Option tt_ut1 is the delta
% time between Terrestrial Time (TT) and Universal Time (UT1). Nodal correction
% outputs pu,pf have units of radians and are MxN dimensions, corresponding to M 
% timesteps and N constituents. Arguments G are similarly sizes as the nodal
% corrections and have units of degrees. 
%
%% References 
%
% Doodson, A. T. (1941). Admiralty manual of tides. His Majesty's Stationery Office.
%
% Schureman, P. (1958). Manual of harmonic analysis and prediction of tides.
% Special Edition No. 98. US Coast and Geodetic Survey, United States Government
% Printing Office, Washington, DC.
%
% Foreman, M. G. G., & Henry, R. F. (1989). The harmonic analysis of tidal 
% model time series. Advances in water resources, 12(3), 109-120.
% https://doi.org/10.1016/0309-1708(89)90017-1
%
%% Author Info
% This function is part of the Tide Model Driver (TMD), which was written by Lana Erofeeva
% and is maintained by Laurie Padman. 
% 
% See also tmd_nodal. 

if isdatetime(t)
   t = datenum(t); 
end

% calculate the mean longitudes of the sun and moon
[s, h, p, N] = tmd_astrol(t, order=order, tt_ut1=tt_ut1);

% mean longitude of solar perigee
if order == 1
    ps = 282.8 + 0.0*t;
else
    % convert to J2000 (dynamic time)
    J2000 = datenum(2000, 1, 1, 12);
    T = t - J2000 + tt_ut1;
    % (Simon et al., 1994)
    ps = 282.94 + 1.7192 * T;
end

% initial time conversions
hour = 24.0*np.mod(t, 1);
% convert from hours solar time into mean lunar time in degrees
tau = 15.0*hour - s + h;
% variable for multiples of 90 degrees (Ray technical note 2017)
% full expansion of Equilibrium Tide includes some negative cosine
% terms and some sine terms (Pugh and Woodworth, 2014)
k = 90.0 + 0.0*t;

% determine equilibrium arguments
astro = transpose([tau; s; h; p; N; ps; k]);
G = astro * tmd_doodson(constituents);

% calculate nodal corrections
[pu, pf] = tmd_nodal(constituents, p, N);

end
