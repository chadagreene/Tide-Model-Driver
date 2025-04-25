function doodson = tmd_doodson(constituents)
% tmd_doodson returns a table of Doodson coefficients for tidal constituents
% 
%% Syntax
% 
%  coef = tmd_doodson(constituents)
% 
%% Description 
% 
% coef = tmd_doodson(constituentst) returns the Doodson coefficients
% (Cartwright numbers) for tidal constituents.
% 
% Inputs:
%    constituents: cell array, e.g., {'m2','k1'} (size 1xC), 
% 
%% References 
%
% Doodson, A. T., & Lamb, H. (1921). The harmonic development of the tide-generating potential.
% Proceedings of the Royal Society of London. Series A, Containing Papers of a Mathematical
% and Physical Character, 100(704), 305-329. https://doi.org/10.1098/rspa.1921.0088
%
% Doodson, A. T. (1941). Admiralty manual of tides. His Majesty's Stationery Office.
%
% Cartwright, D. E., & Tayler, R. J. (1971). New computations of the tide-generating potential.
% Geophysical Journal of the Royal Astronomical Society, 23(1), 45-73.
% https://doi.org/10.1111/j.1365-246X.1971.tb01803.x
%
% Cartwright, D. E., & Edden, A. C. (1973). Corrected tables of tidal harmonics.
% Geophysical Journal International, 33(3), 253-264.
% https://doi.org/10.1111/j.1365-246x.1973.tb03420.x
%
%% Author Info
% TMD was written by Lana Erofeeva (serofeev@coas.oregonstate.edu) and is maintained by
% Laurie Padman (padman@esr.org). 
%
% See also tmd_arguments, and tmd_harpd.

c_data = {'sa'; 'ssa'; 'mm'; 'msf'; 'mf'; 'mt'; 'alpha1'; '2q1'; 'sigma1';
    'q1'; 'rho1'; 'o1'; 'tau1'; 'm1'; 'chi1'; 'pi1'; 'p1'; 's1'; 'k1';
    'psi1'; 'phi1'; 'theta1'; 'j1'; 'oo1'; '2n2'; 'mu2'; 'n2'; 'nu2'; 'm2a';
    'm2'; 'm2b'; 'lambda2'; 'l2'; 't2'; 's2'; 'r2'; 'k2'; 'eta2'; 'mns2';
    '2sm2'; 'm3'; 'mk3'; 's3'; 'mn4'; 'm4'; 'ms4'; 'mk4'; 's4'; 's5'; 'm6';
    's6'; 's7'; 's8'; 'm8'; 'mks2'; 'msqm'; 'mtm'; 'n4'; 'eps2'; 'z0'};

% Get the coefficients for the major constituents
coef = tmd_major_table();

[~,kk] = ismember(constituents,c_data);

% Define outputs: 
if any(kk)
   doodson = coef(:, kk);
else
   doodson = zeros(7, length(constituents));
end

end

function coef = tmd_major_table()
% tmd_major_table returns a table of Doodson coefficients for 60 tidal constituents
%
%% Author Info
% TMD was written by Lana Erofeeva (serofeev@coas.oregonstate.edu) and is maintained by
% Laurie Padman (padman@esr.org). 
%
% See also tmd_doodson.

   % modified Doodson coefficients
   % using 7 index variables: tau, s, h, p, n, pp, k
   % tau: mean lunar time
   % s: mean longitude of moon
   % h: mean longitude of sun
   % p: mean longitude of lunar perigee
   % n: mean longitude of ascending lunar node
   % pp: mean longitude of solar perigee
   % k: 90-degree phase
   coef = zeros(7, 60);
   coef(:, 1) = [0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0]; % Sa
   coef(:, 2) = [0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0]; % Ssa
   coef(:, 3) = [0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0]; % Mm
   coef(:, 4) = [0.0, 2.0, -2.0, 0.0, 0.0, 0.0, 0.0]; % MSf
   coef(:, 5) = [0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0]; % Mf
   coef(:, 6) = [0.0, 3.0, 0.0, -1.0, 0.0, 0.0, 0.0]; % Mt
   coef(:, 7) = [1.0, -4.0, 2.0, 1.0, 0.0, 0.0, -1.0]; % alpha1
   coef(:, 8) = [1.0, -3.0, 0.0, 2.0, 0.0, 0.0, -1.0]; % 2q1
   coef(:, 9) = [1.0, -3.0, 2.0, 0.0, 0.0, 0.0, -1.0]; % sigma1
   coef(:, 10) = [1.0, -2.0, 0.0, 1.0, 0.0, 0.0, -1.0];% q1
   coef(:, 11) = [1.0, -2.0, 2.0, -1.0, 0.0, 0.0, -1.0]; % rho1
   coef(:, 12) = [1.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0]; % o1
   coef(:, 13) = [1.0, -1.0, 2.0, 0.0, 0.0, 0.0, 1.0]; % tau1
   coef(:, 14) = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]; % m1
   coef(:, 15) = [1.0, 0.0, 2.0, -1.0, 0.0, 0.0, 1.0]; % chi1
   coef(:, 16) = [1.0, 1.0, -3.0, 0.0, 0.0, 1.0, -1.0]; % pi1
   coef(:, 17) = [1.0, 1.0, -2.0, 0.0, 0.0, 0.0, -1.0]; % p1
   coef(:, 18) = [1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 1.0]; % s1
   coef(:, 19) = [1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]; % k1
   coef(:, 20) = [1.0, 1.0, 1.0, 0.0, 0.0, -1.0, 1.0]; % psi1
   coef(:, 21) = [1.0, 1.0, 2.0, 0.0, 0.0, 0.0, 1.0]; % phi1
   coef(:, 22) = [1.0, 2.0, -2.0, 1.0, 0.0, 0.0, 1.0]; % theta1
   coef(:, 23) = [1.0, 2.0, 0.0, -1.0, 0.0, 0.0, 1.0]; % j1
   coef(:, 24) = [1.0, 3.0, 0.0, 0.0, 0.0, 0.0, 1.0]; % oo1
   coef(:, 25) = [2.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0]; % 2n2
   coef(:, 26) = [2.0, -2.0, 2.0, 0.0, 0.0, 0.0, 0.0]; % mu2
   coef(:, 27) = [2.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0]; % n2
   coef(:, 28) = [2.0, -1.0, 2.0, -1.0, 0.0, 0.0, 0.0]; % nu2
   coef(:, 29) = [2.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0]; % m2a
   coef(:, 30) = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; % m2
   coef(:, 31) = [2.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0]; % m2b
   coef(:, 32) = [2.0, 1.0, -2.0, 1.0, 0.0, 0.0, 2.0]; % lambda2
   coef(:, 33) = [2.0, 1.0, 0.0, -1.0, 0.0, 0.0, 2.0]; % l2
   coef(:, 34) = [2.0, 2.0, -3.0, 0.0, 0.0, 1.0, 0.0]; % t2
   coef(:, 35) = [2.0, 2.0, -2.0, 0.0, 0.0, 0.0, 0.0]; % s2
   coef(:, 36) = [2.0, 2.0, -1.0, 0.0, 0.0, -1.0, 2.0]; % r2
   coef(:, 37) = [2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0]; % k2
   coef(:, 38) = [2.0, 3.0, 0.0, -1.0, 0.0, 0.0, 0.0]; % eta2
   coef(:, 39) = [2.0, -3.0, 2.0, 1.0, 0.0, 0.0, 0.0]; % mns2
   coef(:, 40) = [2.0, 4.0, -4.0, 0.0, 0.0, 0.0, 0.0]; % 2sm2
   coef(:, 41) = [3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; % m3
   coef(:, 42) = coef(:, 19) + coef(:, 30); % mk3
   coef(:, 43) = [3.0, 3.0, -3.0, 0.0, 0.0, 0.0, 0.0]; % s3
   coef(:, 44) = coef(:, 27) + coef(:, 30); % mn4
   coef(:, 45) = [4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; % m4
   coef(:, 46) = coef(:, 30) + coef(:, 35); % ms4
   coef(:, 47) = coef(:, 30) + coef(:, 37); % mk4
   coef(:, 48) = [4.0, 4.0, -4.0, 0.0, 0.0, 0.0, 0.0]; % s4
   coef(:, 49) = [5.0, 5.0, -5.0, 0.0, 0.0, 0.0, 0.0]; % s5
   coef(:, 50) = [6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; % m6
   coef(:, 51) = [6.0, 6.0, -6.0, 0.0, 0.0, 0.0, 0.0]; % s6
   coef(:, 52) = [7.0, 7.0, -7.0, 0.0, 0.0, 0.0, 0.0]; % s7
   coef(:, 53) = [8.0, 8.0, -8.0, 0.0, 0.0, 0.0, 0.0]; % s8
   coef(:, 54) = [8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; % m8
   coef(:, 55) = coef(:, 30) + coef(:, 37) - coef(:, 35); % mks2
   coef(:, 56) = [0.0, 4.0, -2.0, 0.0, 0.0, 0.0, 0.0]; % msqm
   coef(:, 57) = [0.0, 3.0, 0.0, -1.0, 0.0, 0.0, 0.0]; % mtm
   coef(:, 58) = [4.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0]; % n4
   coef(:, 59) = [2.0, -3.0, 2.0, 1.0, 0.0, 0.0, 0.0]; % eps2
   coef(:, 60) = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; % z0
end

function coef = tmd_minor_table()
% tmd_minor_table returns a table of Doodson coefficients for 18 minor tidal constituents
%
%% Author Info
% TMD was written by Lana Erofeeva (serofeev@coas.oregonstate.edu) and is maintained by
% Laurie Padman (padman@esr.org). 
%
% See also tmd_doodson.

   % modified Doodson coefficients for constituents
   % using 7 index variables: tau, s, h, p, n, pp, k
   % tau: mean lunar time
   % s: mean longitude of moon
   % h: mean longitude of sun
   % p: mean longitude of lunar perigee
   % n: mean longitude of ascending lunar node
   % pp: mean longitude of solar perigee
   % k: 90-degree phase
   coef = zeros(7, 18);
   coef(:,1) = [1.0, -3.0, 0.0, 2.0, 0.0, 0.0, -1.0]; % 2q1
   coef(:,2) = [1.0, -3.0, 2.0, 0.0, 0.0, 0.0, -1.0]; % sigma1
   coef(:,3) = [1.0, -2.0, 2.0, -1.0, 0.0, 0.0, -1.0]; % rho1
   coef(:,4) = [1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0]; % m1b
   coef(:,5) = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0]; % m1a
   coef(:,6) = [1.0, 0.0, 2.0, -1.0, 0.0, 0.0, 1.0]; % chi1
   coef(:,7) = [1.0, 1.0, -3.0, 0.0, 0.0, 1.0, -1.0]; % pi1
   coef(:,8) = [1.0, 1.0, 2.0, 0.0, 0.0, 0.0, 1.0]; % phi1
   coef(:,9) = [1.0, 2.0, -2.0, 1.0, 0.0, 0.0, 1.0]; % theta1
   coef(:,10) = [1.0, 2.0, 0.0, -1.0, 0.0, 0.0, 1.0]; % j1
   coef(:,11) = [1.0, 3.0, 0.0, 0.0, 0.0, 0.0, 1.0]; % oo1
   coef(:,12) = [2.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0]; % 2n2
   coef(:,13) = [2.0, -2.0, 2.0, 0.0, 0.0, 0.0, 0.0]; % mu2
   coef(:,14) = [2.0, -1.0, 2.0, -1.0, 0.0, 0.0, 0.0]; % nu2
   coef(:,15) = [2.0, 1.0, -2.0, 1.0, 0.0, 0.0, 2.0]; % lambda2
   coef(:,16) = [2.0, 1.0, 0.0, -1.0, 0.0, 0.0, 2.0]; % l2a
   coef(:,17) = [2.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0]; % l2b
   coef(:,18) = [2.0, 2.0, -3.0, 0.0, 0.0, 1.0, 0.0]; % t2
end
