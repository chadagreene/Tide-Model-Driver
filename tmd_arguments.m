function [pu, pf, G] = tmd_arguments(t, constituents, order=1, tt_ut1=0.0007)

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
