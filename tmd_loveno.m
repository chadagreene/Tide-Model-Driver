function [h2, k2, l2] = tmd_loveno(omega)
% Compute the body tide Love/Shida numbers for a given frequency. 
% 
%% Syntax 
% 
%  [h2, k2, l2] = tmd_loveno(omega) 
%
%% Description 
% 
% [h2,k2,l2] = tmd_loveno(omega) takes input angular frequencies in radians per second.
% The output Love numbers h2, k2, l2 are dimensionless.
%
%% References 
%
% Wahr, J. M., & Sasao, T. (1981). A diurnal resonance in the ocean tide and in the
% Earth's load response due to the resonant free `core nutation'. Geophysical Journal
% of the Royal Astronomical Society, 64(3), 747-765.
% https://doi.org/10.1111/j.1365-246x.1981.tb02693.x
%
% Mathews, P. M., Buffett, B. A., & Shapiro, I. I. (1995). Love numbers for diurnal
% tides: Relation to wobble admittances and resonance expansions. Journal of
% Geophysical Research: Solid Earth, 100(B6), 9935-9948.
% https://doi.org/10.1029/95jb00670
%
%% Author Info
% This function is part of the Tide Model Driver (TMD), which was written by Lana Erofeeva
% and is maintained by Laurie Padman. 
% 
% See also tmd_nodal. 

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
