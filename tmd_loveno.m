function [h2, k2, l2] = tmd_loveno(omega)

% Love numbers for different frequency bands
if (omega > 1e-4)
    % tides in the semi-diurnal band
    h2 = 0.609;
    k2 = 0.302;
    l2 = 0.0852;
elseif (omega < 2e-5)
    % tides in the long period band
    h2 = 0.606;
    k2 = 0.299;
    l2 = 0.0840;
else
    % calculate love numbers following J. Wahr (1979)
    % use resonance formula for tides in the diurnal band
    % frequency of the o1 tides (radians/second)
    omega_o1 = tmd_frequency('o1');
    % free core nutation frequencies (cycles per sidereal day)
    lambda_fcn = 1.0023214;
    % convert frequency from cycles per sidereal day
    % frequency of free core nutation (radians/second)
    omega_fcn = lambda_fcn*7292115e-11
    % Love number parameters for PREM earth model
    % from Mathews et al. (1995) table 3
    h0 = [5.994e-1, -2.532e-3];
    k0 = [2.962e-1, -1.271e-3];
    l0 = [8.378e-2, 7.932e-5];
    % Love numbers for frequency using equation 4.18 of Wahr (1981)
    % (simplification to use only the free core nutation term)
    ratio = (omega - omega_o1)/(omega_fcn - omega);
    h2 = h0(1) + h0(2)*ratio;
    k2 = k0(1) + k0(2)*ratio;
    l2 = l0(1) + l0(2)*ratio;
end

end
