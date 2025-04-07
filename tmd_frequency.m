function omega = tmd_frequency(constituents, order=1)

% J2000 (2000-01-01 12:00:00)
t = [730486.5, 730486.55];
% time interval in seconds
deltat = 86400.0*(t(2) - t(1));
% calculate the mean longitudes of the sun and moon
[s, h, p, n] = tmd_astrol(t, order=order);
% mean longitude of solar perigee
if order == 1
    ps = [282.8, 282.8];
else
    % (Simon et al., 1994)
    T = t - 730486.4993;
    ps = 282.94 + 1.7192 * T
end

% initial time conversions
hour = 24.0*mod(t, 1);
% convert from hours solar time into mean lunar time in degrees
tau = 15.0*hour - s + h;
% variable for multiples of 90 degrees (Ray technical note 2017)
k = [90.0, 90.0];

% determine differential in equilibrium arguments
astro = [tau; s; h; p; n; ps; k];
rates = transpose(astro(:,2) - astro(:,1))/deltat;
fd = rates * tmd_doodson(constituents);
% convert to radians per second
omega = 2.0*pi*fd/360.0;

end
