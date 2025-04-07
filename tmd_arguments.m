function [pu, pf, G] = tmd_arguments(t, constituents, order=1)

if isdatetime(t)
   t = datenum(t); 
end

% convert times to MJD
t_mjd = t - datenum(1992, 1, 1) + 48622;
% number of temporal values
nt = length(t);

% calculate the mean longitudes of the sun and moon
[s, h, p, N] = tmd_astrol(t, order=order);

% mean longitude of solar perigee
if order == 1
    ps = [282.8, 282.8];
else
    % (Simon et al., 1994)
    T = t_mjd - 51544.4993;
    ps = 282.94 + 1.7192 * T
end

% initial time conversions
hour = 24.0*np.mod(t_mjd, 1);
% convert from hours solar time into mean lunar time in degrees
tau = 15.0*hour - s + h;
% variable for multiples of 90 degrees (Ray technical note 2017)
% full expansion of Equilibrium Tide includes some negative cosine
% terms and some sine terms (Pugh and Woodworth, 2014)
k = 90.0 + zeros(nt, 1);

% determine equilibrium arguments
astro = transpose([tau; s; h; p; N; ps; k]);
G = astro * tmd_doodson(constituents);

% calculate nodal corrections
[pu, pf] = tmd_nodal(constituents, p, N)

end
