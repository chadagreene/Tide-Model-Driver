function omega = tmd_frequency(constituents, order=1, tt_ut1=0.0007)
% Calculates the angular frequency for tidal constituents
% 
%% Syntax 
% 
%  omega = tmd_frequency(constituents) 
%
%% Description 
% 
% omega = tmd_frequency(constituents) takes input constituents 
% as cell array (1xN constituents). Option order is the minimum polynomial order
% for calculating the astronomical mean longitudes. Option tt_ut1 is the delta
% time between Terrestrial Time (TT) and Universal Time (UT1). Output angular
% frequencies omega are in radians per second.
%
%% References 
%
% Doodson, A. T. (1941). Admiralty manual of tides. His Majesty's Stationery Office.
%
%% Author Info
% This function is part of the Tide Model Driver (TMD), which was written by Lana Erofeeva
% and is maintained by Laurie Padman. 
% 
% See also tmd_arguments.

% J2000 (2000-01-01 12:00:00)
J2000 = datenum(2000, 1, 1, 12);
t = J2000 + [0.0, 0.05];
% time interval in seconds
deltat = 86400.0*(t(2) - t(1));
% calculate the mean longitudes of the sun and moon
[s, h, p, n] = tmd_astrol(t, order=order, tt_ut1=tt_ut1);
% mean longitude of solar perigee
if order == 1
    ps = 282.8 + 0.0*t;
else
    % convert to J2000 (dynamic time)
    T = t - J2000 + tt_ut1;
    % (Simon et al., 1994)
    ps = 282.94 + 1.7192 * T;
end

% initial time conversions
hour = 24.0*mod(t, 1);
% convert from hours solar time into mean lunar time in degrees
tau = 15.0*hour - s + h;
% variable for multiples of 90 degrees (Ray technical note 2017)
k = 90.0 + 0.0*t;

% determine differential in equilibrium arguments
astro = [tau; s; h; p; n; ps; k];
rates = transpose(astro(:,2) - astro(:,1))/deltat;
fd = rates * tmd_doodson(constituents);
% convert to radians per second
omega = 2.0*pi*fd/360.0;

end
