function h = tmd_lpet(t, lat, order=1, tt_ut1=0.0007)
% long period constituents
c_data = {'node'; 'sa'; 'ssa'; 'msm'; '065.445'; 'mm';
    '065.465'; 'msf'; '075.355'; 'mf'; 'mf+'; '075.575';
    'mst'; 'mt'; '085.465'};

% compute principal mean longitudes
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

% convert to negative mean longitude of the ascending node (N')
n = mod(360.0 - N, 1.0);
% determine equilibrium arguments
astro = transpose([s; h; p; n; ps]);

% Cartwright and Edden potential amplitudes (centimeters)
% assemble long-period tide potential from 15 CTE terms greater than 1 mm
amajor = zeros(15,1);
% group 0,0
% nodal term is included but not the constant term.
amajor(1) = 2.7929;       % node
amajor(2) = -0.4922;      % sa
amajor(3) = -3.0988;      % ssa
% group 0,1
amajor(4) = -0.6728;      % msm
amajor(5) = 0.231;
amajor(6) = -3.5184;      % mm
amajor(7) = 0.228;
% group 0,2
amajor(8) = -0.5837;      % msf
amajor(9) = -0.288;
amajor(10) = -6.6607;     % mf
amajor(11) = -2.763;      % mf+
amajor(12) = -0.258;
% group 0,3
amajor(13) = -0.2422;     % mst
amajor(14) = -1.2753;     % mt
amajor(15) = -0.528;

% Doodson coefficients for 15 long-period terms
coef = zeros(5, 15);
% group 0,0
coef(:,1) = [0.0, 0.0, 0.0, 1.0, 0.0];   % node
coef(:,2) = [0.0, 1.0, 0.0, 0.0, -1.0];  % sa
coef(:,3) = [0.0, 2.0, 0.0, 0.0, 0.0];   % ssa
% group 0,1
coef(:,4) = [1.0, -2.0, 1.0, 0.0, 0.0];   % msm
coef(:,5) = [1.0, 0.0, -1.0, -1.0, 0.0];
coef(:,6) = [1.0, 0.0, -1.0, 0.0, 0.0];   % mm
coef(:,7) = [1.0, 0.0, -1.0, 1.0, 0.0];
% group 0,2
coef(:,8) = [2.0, -2.0, 0.0, 0.0, 0.0];   % msf
coef(:,9) = [2.0, 0.0, -2.0, 0.0, 0.0];
coef(:,10) = [2.0, 0.0, 0.0, 0.0, 0.0];   % mf
coef(:,11) = [2.0, 0.0, 0.0, 1.0, 0.0];   % mf+
coef(:,12) = [2.0, 0.0, 0.0, 2.0, 0.0];
% group 0,3
coef(:,13) = [3.0, -2.0, 1.0, 0.0, 0.0];  % mst
coef(:,14) = [3.0, 0.0, -1.0, 0.0, 0.0];  % mt
coef(:,15) = [3.0, 0.0, -1.0, 1.0, 0.0];

% determine equilibrium arguments
G = astro * coef;
Z = (cosd(G) * amajor)/100.0; % convert to meters

% Love numbers for long-period tides (Wahr, 1981)
k2 = 0.299;
h2 = 0.606;
% tilt factor: response with respect to the solid earth
gamma_2 = (1.0 + k2 - h2);
% 2nd degree Legendre polynomials
P20 = 0.5*(3.0*sind(lat).^2 - 1.0);
% calculate long-period equilibrium tide for latitude
h = gamma_2*sqrt((4 + 1)/(4*pi))*P20.*transpose(Z);

end
